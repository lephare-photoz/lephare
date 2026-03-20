import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import scipy.stats as sc
from scipy.optimize import minimize
import scipy.interpolate as sip

import warnings

class DeltaMaps:
    """
    Handle LSST maglim maps (FITS) and apply them to catalogs.
    """

    def __init__(self, maglim_fits):
        """
        Parameters
        ----------
        maglim_fits : dict
            {band: fits_filename}
            e.g. {"r": "ECDFS_maglim_r.fits", "i": "ECDFS_maglim_i.fits"}
        """
        self.maps = {}

        for band, fname in maglim_fits.items():
            with fits.open(fname) as hdul:
                data = hdul[0].data
                wcs = WCS(hdul[0].header)

            self.maps[band] = {
                "data": data,
                "wcs": wcs,
                "shape": data.shape,
            }

    def bands(self):
        """Return available bands."""
        return list(self.maps.keys())

    def get_map(self, band):
        """Return (data, wcs) for a given band."""
        if band not in self.maps:
            raise ValueError(f"Band '{band}' not loaded.")
        return self.maps[band]["data"], self.maps[band]["wcs"]

    def get_maglim_at_radec(self, band, ra, dec):
        """
        Return maglim values at given RA/Dec positions.

        Parameters
        ----------
        band : str
        ra, dec : array-like (deg)

        Returns
        -------
        maglim : ndarray
        """
        data, wcs = self.get_map(band)

        x, y = wcs.world_to_pixel_values(ra, dec)
        # print(x,y)
        x = np.round(x).astype(int)
        y = np.round(y).astype(int)
        ny, nx = data.shape
        good = (x >= 0) & (x < nx) & (y >= 0) & (y < ny)

        maglim = np.full(len(ra), np.nan)
        maglim[good] = data[y[good], x[good]]

        return maglim

    def add_delta_to_catalog(self, cat, mag_columns, radec_names=["ra","dec"]):
        """
        Add maglim_band and delta_band columns to a catalog.

        Parameters
        ----------
        cat : pandas.DataFrame
            Must contain 'ra', 'dec' and magnitude columns.
        mag_columns : dict
            {band: mag_column_name}

        Returns
        -------
        cat : pandas.DataFrame
        """
        for band, mag_col in mag_columns.items():
            if band not in self.maps:
                raise ValueError(f"Band '{band}' not loaded.")

            maglim = self.get_maglim_at_radec(
                band,
                cat[f"{radec_names[0]}"].values,
                cat[f"{radec_names[1]}"].values,
            )
            cat[f"maglim_{band}"] = maglim
            cat[f"delta_{band}"] = cat[mag_col] - maglim
        # print(cat[f"maglim_{band}"])
        return cat

    def plot_map(self, band, vmin=None, vmax=None, cmap="viridis", figsize=(8, 5)):
        """
        Plot the maglim map for a given band using imshow.
        """
        data, wcs = self.get_map(band)

        plt.figure(figsize=figsize)

        im = plt.imshow(
            data,
            origin="lower",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
        )

        plt.colorbar(im, label=f"PSF mag_lim (5σ, {band}-band)")
        plt.xlabel("Pixel X")
        plt.ylabel("Pixel Y")
        plt.title(f"Maglim map – {band}-band")

        plt.tight_layout()
        plt.show()


class PhotometricErrorModel: 
    """ 
    Empirical gamma model for photometric errors: p(sigma | mag) 
    """ 
    # 1. fit
    def fit(self, mag, sigma, delta_mag=None, min_per_bin=None, sharpness_factor=None, polyfit_deg=4, fit_ends=True): 
        """ 
        Fit conditional gamma model without binning using smoothing splines. 
        Parameters 
        ---------- 
        mag: array-like
            Conditioning variable. 
        sigma: array-like
            Observed uncertainties (>0). 
        delta_mag: float
            mag bins width to compute gamma fit on sigma within.
        min_per_bin: int
            Minimum number of points per bin.
        sharpness_factor: float
            If specified, sharps the local error distributions 
            following loc = mag.min() * sharpness_factor
        polyfit_deg: int
            deg of the polynome to fit on df, loc, scale
        """

        mag = np.asarray(mag) 
        sigma = np.asarray(sigma) 

        # --- Cleaning ---
        good = np.isfinite(mag) & np.isfinite(sigma) & (sigma > 0) 
        mag = mag[good] 
        sigma = sigma[good] 
        # sort
        idx = np.argsort(mag) 
        mag = mag[idx] 
        sigma = sigma[idx]
        N = len(mag)
        if delta_mag is None:
            delta_mag = (mag.max() - mag.min()) / np.sqrt(N)
        if min_per_bin is None:
            min_per_bin = int(np.sqrt(np.log10(N)) * np.sqrt(N)) #this factor has been determined empirically
        if sharpness_factor is True:
            sharpness_factor = 1.1

        # --- Fitting error distribution dynamically with mag bins ---
        # binning in mags
        mag_bins = np.arange(mag.min(), mag.max() + delta_mag, delta_mag)
        n_bins = len(mag_bins)
        eps=1e-12
        i = 0
        gamma_params = []

        while i < n_bins - 1:
            left_idx = np.searchsorted(mag, mag_bins[i], side="left")

            j = 1
            while i + j < n_bins:
                right_idx = np.searchsorted(mag, mag_bins[i + j], side="left")
                if right_idx - left_idx >= min_per_bin:
                    break
                j += 1

            # if reaching the end, take what is remaining
            if i + j >= n_bins and fit_ends == True:
                right_idx = len(mag)
                warnings.warn("fit_ends is True and npoints < min_per_bin in the last bins. This can lead to unexpected behavior.")

            # bin data
            m_bin = mag[left_idx:right_idx]
            s_bin = sigma[left_idx:right_idx]

            # sharp error dist
            if isinstance(sharpness_factor, float):
                s_bin_min = s_bin.min()
                s_bin = s_bin - s_bin_min*sharpness_factor #shift the data to force 'loc' to start higher
                s_bin = s_bin[s_bin>0] + s_bin_min*sharpness_factor #and shift back, do it like this to prevent error

            if len(m_bin) == 0:
                break
            
            # local fit
            middle_mag = 0.5 * (m_bin.min() + m_bin.max())
            # Fit gamma per bin
            df, loc, scale = sc.gamma.fit(s_bin, floc=s_bin.min() - eps)
            # stock in table
            gamma_params.append((middle_mag, df, loc, scale))
            # print(f"bin {i}-{i+j}: {len(m_bin)} samples, mag [{m_bin.min():.3f}, {m_bin.max():.3f}]")
            i += j
            

        # --- Smooth the gamma params VS mag ---
        gamma_params = np.array(gamma_params, dtype=float)
        x = gamma_params[:, 0]
        df_vals = gamma_params[:, 1]
        loc_vals = gamma_params[:, 2]
        scale_vals = gamma_params[:, 3]

        if polyfit_deg == False:
            self.smooth_df = sip.interp1d(x, df_vals, kind=False, bounds_error=False, fill_value=0.0)
            self.smooth_loc = sip.interp1d(x, loc_vals, kind=False, bounds_error=False, fill_value=0.0)
            self.smooth_scale = sip.interp1d(x, scale_vals, kind=False, bounds_error=False, fill_value=0.0)
        else:
            coef_df = np.polyfit(x, df_vals, deg=polyfit_deg)
            self.smooth_df = lambda x: np.polyval(coef_df, x)

            eps_sm = 1e-12
            coef_loc = np.polyfit(x, np.log(loc_vals + eps_sm), deg=polyfit_deg)
            self.smooth_loc = lambda x: np.exp(np.polyval(coef_loc, x))

            coef_scale = np.polyfit(x, np.log(scale_vals + eps_sm), deg=polyfit_deg)
            self.smooth_scale = lambda x: np.exp(np.polyval(coef_scale, x))


        # x_plot = x #np.concatenate((x, np.array([29,30])))
        # plt.figure()
        # plt.plot(x, df_vals)
        # plt.plot(x_plot, self.smooth_df(x_plot))
        # plt.show()
        # plt.figure()
        # plt.plot(x, loc_vals)
        # plt.plot(x_plot, self.smooth_loc(x_plot))
        # plt.show()
        # plt.figure()
        # plt.plot(x, scale_vals)
        # plt.plot(x_plot, self.smooth_scale(x_plot))
        # plt.show()
        return self

    # 2. Prediction
    def predict(self, x, statistic="mean"):
        x = np.asarray(x)

        df = np.maximum(self.smooth_df(x), 1e-8)
        loc = self.smooth_loc(x)
        scale = np.maximum(self.smooth_scale(x), 1e-12)

        if statistic == "mean":
            return loc + df * scale

        elif statistic == "median":
            return sc.gamma.median(a=df, loc=loc, scale=scale)

        elif statistic == "mode":
            mode = loc + (df - 1) * scale
            mode[df <= 1] = loc[df <= 1]  # mode non défini si df ≤ 1
            return mode

        else:
            raise ValueError("statistic must be 'mean', 'median', or 'mode'")

    # 3. Sampling
    def sample(self, x, size=None, random_state=None):
        rng = np.random.default_rng(random_state)
        x = np.asarray(x)

        df = np.maximum(self.smooth_df(x), 1e-6)
        loc = self.smooth_loc(x)
        scale = np.maximum(self.smooth_scale(x), 1e-12)

        return sc.gamma.rvs(a=df, loc=loc, scale=scale, size=size, random_state=rng)

    # 4. Store models in a collection for re-use
    def model_collection(self, x_lists, sigma_lists, model_names, fit_ends=True,
                delta_mag=0.1, min_per_bin=1000, sharpness_factor=None, polyfit_deg=4):
        """
        Build a collection of photometric gamma error models.
        """
        # ---- manage entries ----
        if hasattr(x_lists, "values"):
            x_lists = x_lists.values.T
        else:
            x_lists = np.asarray(x_lists)

        if hasattr(sigma_lists, "values"):
            sigma_lists = sigma_lists.values.T
        else:
            sigma_lists = np.asarray(sigma_lists)

        x_lists = np.atleast_2d(x_lists)
        sigma_lists = np.atleast_2d(sigma_lists)

        if isinstance(model_names, str):
            model_names = [model_names]

        n_models = len(model_names)

        if n_models != x_lists.shape[0]:
            raise ValueError("model_names must match number of x_lists")

        # ---- helper to broadcast parameters ----
        def broadcast_param(param, name):
            if isinstance(param, (list, tuple, np.ndarray)):
                if len(param) != n_models:
                    raise ValueError(f"{name} must have length {n_models}")
                return param
            else:
                return [param] * n_models

        delta_mag_list = broadcast_param(delta_mag, "delta_mag")
        min_per_bin_list = broadcast_param(min_per_bin, "min_per_bin")
        sharpness_list = broadcast_param(sharpness_factor, "sharpness_factor")
        polyfit_deg_list = broadcast_param(polyfit_deg, "polyfit_deg")
        models = {}

        for i, name in enumerate(model_names):
            x = x_lists[i]
            sigma = sigma_lists[i]
            model = PhotometricErrorModel()
            model.fit(x, sigma, delta_mag=delta_mag_list[i],
                min_per_bin=min_per_bin_list[i],
                sharpness_factor=sharpness_list[i],
                polyfit_deg=polyfit_deg_list[i],
                fit_ends=fit_ends)

            models[name] = model

        return models

    def plot_models(self, models, x_lists, sigma_lists, model_names=None, xlimit=None, ylimit=None,
            show_scatter=True, show_model=True, statistic="mean", subsample=None):
        """
        Plot photometric gamma error models.
        """
        # --- Handle DataFrame / array ---
        if hasattr(x_lists, "values"):
            x_lists = x_lists.values.T
        else:
            x_lists = np.asarray(x_lists)

        if hasattr(sigma_lists, "values"):
            sigma_lists = sigma_lists.values.T
        else:
            sigma_lists = np.asarray(sigma_lists)

        x_lists = np.atleast_2d(x_lists)
        sigma_lists = np.atleast_2d(sigma_lists)

        if model_names is None:
            model_names = list(models.keys())

        for i, name in enumerate(model_names):
            model = models[name]
            x = x_lists[i]
            sigma = sigma_lists[i]

            if subsample is not None and subsample < len(x):
                idx = np.random.choice(len(x), subsample, replace=False)
                x = x[idx]
                sigma = sigma[idx]

            plt.figure(figsize=(6, 4))
            if show_scatter:
                plt.scatter(x, sigma, s=5, alpha=0.3, label="Observed")
            if show_model:
                x_grid = np.linspace(np.nanmin(x), np.nanmax(x), 400)
                sigma_model = model.predict(x_grid, statistic=statistic)
                plt.plot(x_grid, sigma_model, color="black", lw=2, label=f"Gamma model ({statistic})")

            if xlimit is not None:
                plt.xlim(xlimit if isinstance(xlimit, tuple) else xlimit[i])
            if ylimit is not None:
                plt.ylim(ylimit if isinstance(ylimit, tuple) else ylimit[i])
            plt.xlabel("Magnitude")
            plt.ylabel(r"$\sigma_{\rm mag}$")
            plt.title(f"Photometric error model – {name}")
            plt.legend()
            plt.tight_layout()
            plt.show()

    def sample_errors_on_catalog(self, cat, mag_columns, models, model_names, 
                    suffix="_err_sim", random_state=None, format="MMEE"):
        """
        Generate photometric errors for multiple magnitude columns
        using corresponding fitted models.

        Parameters
        ----------
        cat : pandas.DataFrame
            Input catalog.
        mag_columns : list of str
            Columns containing magnitudes.
        models : dict
            Dictionary of PhotometricErrorModel instances.
        model_names : list of str
            Model names corresponding to each magnitude column.
        suffix : str
            Suffix for new error columns.
        random_state : int or None
            Seed for reproducibility.
        format : {"MEME", "MMEE", "end"}
            Column ordering:
            - "MEME"  → (mag, err, mag, err, ...)
            - "MMEE" → (mag, mag, mag, err, err, err)
            - "end"     → errors appended at end (default)

        Returns
        -------
        pandas.DataFrame
        """

        import numpy as np

        if len(mag_columns) != len(model_names):
            raise ValueError("mag_columns and model_names must have same length.")

        if format not in ["MEME", "MMEE", "end"]:
            raise ValueError("format must be 'paired', 'grouped', or 'end'.")

        rng = np.random.default_rng(random_state)
        cat_out = cat.copy()

        error_columns = []

        # --- Generate errors
        for col, model_name in zip(mag_columns, model_names):

            if model_name not in models:
                raise ValueError(f"Model '{model_name}' not found in models.")

            model = models[model_name]
            x = np.asarray(cat_out[col])

            sigma_sim = model.sample(x, random_state=rng)

            err_col = col + suffix
            cat_out[err_col] = sigma_sim
            error_columns.append(err_col)

        # --- Reorder columns if needed
        columns = list(cat_out.columns)
        first_mag_index = columns.index(mag_columns[0])

        cols_before = columns[:first_mag_index]
        remaining_cols = [c for c in columns if c not in cols_before]

        if format == "end":
            return cat_out

        if format == "MMEE":

            ordered_main = mag_columns + error_columns
            other_cols = [c for c in remaining_cols if c not in ordered_main]

            final_cols = cols_before + ordered_main + other_cols
            return cat_out[final_cols]


        if format == "MEME":

            paired = []
            for m, e in zip(mag_columns, error_columns):
                paired.extend([m, e])

            other_cols = [c for c in remaining_cols if c not in paired]

            final_cols = cols_before + paired + other_cols
            return cat_out[final_cols]



### TRASH ###


    # def fit(self, mag, sigma, delta_mag=0.1, min_per_bin=1000, sharpness_factor=None, polyfit_deg=4): 
    #     """ 
    #     Fit conditional gamma model without binning using smoothing splines. 
    #     Parameters 
    #     ---------- 
    #     mag: array-like
    #         Conditioning variable. 
    #     sigma: array-like
    #         Observed uncertainties (>0). 
    #     delta_mag: float
    #         mag bins width to compute gamma fit on sigma within.
    #     min_per_bin: int
    #         Minimum number of points per bin.
    #     sharpness_factor: float
    #         If specified, sharps the local error distributions 
    #         following loc = mag.min() * sharpness_factor
    #     polyfit_deg: int
    #         deg of the polynome to fit on df, loc, scale
    #     """ 

    #     import numpy as np
    #     import scipy.stats as sc

    #     mag = np.asarray(mag) 
    #     sigma = np.asarray(sigma) 

    #     # --- Cleaning ---
    #     good = np.isfinite(mag) & np.isfinite(sigma) & (sigma > 0) 
    #     mag = mag[good] 
    #     sigma = sigma[good] 
    #     # sort
    #     idx = np.argsort(mag) 
    #     mag = mag[idx] 
    #     sigma = sigma[idx]

    #     # --- Fitting error distribution dynamically with mag bins ---
    #     # binning in mags
    #     mag_bins = np.arange(mag.min(), mag.max() + delta_mag, delta_mag)  # inclure max
    #     n_bins = len(mag_bins)
    #     eps=1e-12
    #     i = 0
    #     gamma_params = []

    #     while i < n_bins - 1:
    #         left_idx = np.searchsorted(mag, mag_bins[i], side="left")

    #         j = 1
    #         while i + j < n_bins:
    #             right_idx = np.searchsorted(mag, mag_bins[i + j], side="left")
    #             if right_idx - left_idx >= min_per_bin:
    #                 break
    #             j += 1

    #         # if reaching the end, take what is remaining
    #         if i + j >= n_bins:
    #             right_idx = len(mag)

    #         # bin data
    #         m_bin = mag[left_idx:right_idx]
    #         s_bin = sigma[left_idx:right_idx]

    #         # sharp error dist
    #         if isinstance(sharpness_factor, float):
    #             s_bin_min = s_bin.min()
    #             s_bin = s_bin - s_bin_min*sharpness_factor #shift the data to force 'loc' to start higher
    #             s_bin = s_bin[s_bin>0] + s_bin_min*sharpness_factor #and shift back, do it like this to prevent error

    #         if len(m_bin) == 0:
    #             break
            

    #         middle_mag = 0.5 * (m_bin.min() + m_bin.max())
    #         # Fit gamma per bin
    #         df, loc, scale = sc.gamma.fit(s_bin, floc=s_bin.min() - eps)

    #         # stock in table
    #         gamma_params.append((middle_mag, df, loc, scale))

    #         # print(f"bin {i}-{i+j}: {len(m_bin)} samples, mag [{m_bin.min():.3f}, {m_bin.max():.3f}]")
    #         i += j
            

    #     # --- Smooth the gamma params VS mag ---
    #             # --- Store gamma params ---
    #     gamma_params = np.array(gamma_params, dtype=float)
    #     self.gamma_params = gamma_params

    #     x = gamma_params[:, 0]
    #     df_vals = gamma_params[:, 1]
    #     loc_vals = gamma_params[:, 2]
    #     scale_vals = gamma_params[:, 3]

    #     # 
    #     coef_df = np.polyfit(x, df_vals, deg=polyfit_deg)
    #     poly_df = np.poly1d(coef_df)
    #     x_max = x.max()
    #     df_max = poly_df(x_max)
    #     df_slope = np.polyder(poly_df)(x_max)
    #     def smooth_df(m):
    #         m = np.asarray(m)
    #         y = poly_df(m)
    #         mask = m > x_max
    #         if np.any(mask):
    #             y[mask] = df_max + df_slope * (m[mask] - x_max)
    #         return y
    #     self.smooth_df = smooth_df

    #     # 
    #     eps_sm = 1e-12
    #     coef_loc = np.polyfit(x, np.log(loc_vals + eps_sm), deg=polyfit_deg)
    #     poly_loc = np.poly1d(coef_loc)
    #     loc_max = poly_loc(x_max)
    #     loc_slope = np.polyder(poly_loc)(x_max)
    #     def smooth_loc(m):
    #         m = np.asarray(m)
    #         y = poly_loc(m)
    #         mask = m > x_max
    #         if np.any(mask):
    #             y[mask] = loc_max + loc_slope * (m[mask] - x_max)
    #         return np.exp(y)
    #     self.smooth_loc = smooth_loc

    #     #
    #     coef_scale = np.polyfit(x, np.log(scale_vals + eps_sm), deg=polyfit_deg)
    #     poly_scale = np.poly1d(coef_scale)
    #     scale_max = poly_scale(x_max)
    #     scale_slope = np.polyder(poly_scale)(x_max)
    #     def smooth_scale(m):
    #         m = np.asarray(m)
    #         y = poly_scale(m)
    #         mask = m > x_max
    #         if np.any(mask):
    #             y[mask] = scale_max + scale_slope * (m[mask] - x_max)
    #         return np.exp(y)
    #     self.smooth_scale = smooth_scale

    #     x_plot = np.concatenate((x, np.array([29,30])))
    #     plt.figure()
    #     plt.plot(x, df_vals)
    #     plt.plot(x_plot, self.smooth_df(x_plot))
    #     plt.show()
    #     plt.figure()
    #     plt.plot(x, loc_vals)
    #     plt.plot(x_plot, self.smooth_loc(x_plot))
    #     plt.show()
    #     plt.figure()
    #     plt.plot(x, scale_vals)
    #     plt.plot(x_plot, self.smooth_scale(x_plot))
    #     plt.show()
    #     return self


# from scipy.interpolate import interp1d
# from scipy.interpolate import UnivariateSpline
# class PhotometricErrorModel: 
#     """ 
#     Conditional log-normal model for photometric errors: p(sigma | x) 
#     """ 
#     # 1. fit
#     def fit(self, x, sigma, smooth_mu=None, smooth_tau=None): 
#         """ 
#         Fit conditional log-normal model without binning using smoothing splines. 
#         Parameters 
#         ---------- 
#         x : array-like Conditioning variable. 
#         sigma : array-like Observed uncertainties (>0). 
#         smooth_mu : float or None Smoothing factor for mean spline. 
#         smooth_tau : float or None Smoothing factor for variance spline. 
#         """ 
#         x = np.asarray(x) 
#         sigma = np.asarray(sigma) 
#         good = np.isfinite(x) & np.isfinite(sigma) & (sigma > 0) 
#         x = x[good] 
#         sigma = sigma[good] 
#         # Sort for spline stability 
#         idx = np.argsort(x) 
#         x = x[idx] 
#         sigma = sigma[idx] 
#         ln_sigma = np.log(sigma) 

#         # --- Fit mean of log(sigma) 
#         from scipy.ndimage import median_filter 
#         med_mu = median_filter(ln_sigma, size=51) 
#         self.mu_spline = UnivariateSpline(x, med_mu, s=smooth_mu) 

#         # plt.figure()
#         # plt.scatter(x, sigma)
#         # plt.plot(x, np.exp(self.mu_spline(x)),c='r')
#         # plt.show()

#         # --- Estimate local variance 
#         residuals = ln_sigma - self.mu_spline(x) 


#         # Fit variance of residuals 
#         self.tau_spline = UnivariateSpline(x, residuals**2, s=smooth_tau)

#         return self 

#     # 2. Prediction 
#     def predict(self, x, statistic="median"): 
#         x = np.asarray(x) 
#         mu = self.mu_spline(x) 
#         tau = np.sqrt(np.maximum(self.tau_spline(x), 1e-8)) 

#         if statistic == "median": 
#             return np.exp(mu) 
#         elif statistic == "mean": 
#             return np.exp(mu + 0.5 * tau**2) 
#         else: 
#             raise ValueError("statistic must be 'median' or 'mean'") 

#     # 3. Sampling 
#     def sample(self, x): 
#         x = np.asarray(x) 
#         mu = self.mu_spline(x) 
#         tau = np.sqrt(self.tau_spline(x)) 
#         return np.random.lognormal(mu, tau)
#         # return sc.lognorm.rvs(s=tau, loc=0, scale=np.exp(mu))

    
#     # 4. Store models in a collection for re-use
#     def model_collection(self, x_lists, sigma_lists, model_names, smooth_mu=None, smooth_tau=None):
#         """
#         Build a collection of photometric error models.

#         Parameters
#         ----------
#         x_lists : array-like or DataFrame
#             Conditioning variables.
#         sigma_lists : array-like or DataFrame
#             Corresponding photometric uncertainties.
#         model_names : str or list of str
#             Names of the models (e.g. band names).
#         smooth_mu : float or None
#             Smoothing factor for mean spline.
#         smooth_tau : float or None
#             Smoothing factor for variance spline.

#         Returns
#         -------
#         dict
#             Dictionary of models or fitted parameters.
#         """

#         # ---- manage entries ----
#         if hasattr(x_lists, "values"):
#             x_lists = x_lists.values.T
#         else:
#             x_lists = np.asarray(x_lists)

#         if hasattr(sigma_lists, "values"):
#             sigma_lists = sigma_lists.values.T
#         else:
#             sigma_lists = np.asarray(sigma_lists)

#         x_lists = np.atleast_2d(x_lists)
#         sigma_lists = np.atleast_2d(sigma_lists)

#         if isinstance(model_names, str):
#             model_names = [model_names]

#         if len(model_names) != x_lists.shape[0]:
#             raise ValueError("model_names must match number of x_lists")

#         models = {}

#         for i, name in enumerate(model_names):
#             x = x_lists[i]
#             sigma = sigma_lists[i]
#             model = PhotometricErrorModel()

#             model.fit(x, sigma,
#                 smooth_mu=smooth_mu, smooth_tau=smooth_tau)

#             models[name] = model

#         return models

#     def plot_models(self, models, x_lists, sigma_lists, model_names=None,
#         show_scatter=True, show_model=True, statistic="median", subsample=None):
#         """
#         Plot photometric error models.

#         Parameters
#         ----------
#         models : dict
#             Dictionary of PhotometricErrorModel objects.
#         x_lists : array-like or DataFrame
#             Conditioning variables used to fit the models.
#         sigma_lists : array-like or DataFrame
#             Observed photometric uncertainties.
#         model_names : list of str or None
#             Names of models to plot. If None, plot all.
#         show_scatter : bool, optional
#             Whether to show observed data scatter.
#         show_model : bool, optional
#             Whether to show model prediction.
#         statistic : {"median", "mean"}, optional
#             Statistic to plot for the model.
#         subsample : int or None, optional
#             Random subsampling factor for scatter points.
#         """

#         # --- Handle DataFrame / array ---
#         if hasattr(x_lists, "values"):
#             x_lists = x_lists.values.T
#         else:
#             x_lists = np.asarray(x_lists)

#         if hasattr(sigma_lists, "values"):
#             sigma_lists = sigma_lists.values.T
#         else:
#             sigma_lists = np.asarray(sigma_lists)

#         x_lists = np.atleast_2d(x_lists)
#         sigma_lists = np.atleast_2d(sigma_lists)

#         if model_names is None:
#             model_names = list(models.keys())

#         for i, name in enumerate(model_names):
#             model = models[name]
#             x = x_lists[i]
#             sigma = sigma_lists[i]

#             if subsample is not None and subsample < len(x):
#                 idx = np.random.choice(len(x), subsample, replace=False)
#                 x = x[idx]
#                 sigma = sigma[idx]

#             plt.figure(figsize=(6, 4))

#             if show_scatter:
#                 plt.scatter(x, sigma, s=5, alpha=0.3, label="Observed")

#             if show_model:
#                 x_grid = np.linspace( np.nanmin(x), np.nanmax(x), 300)
#                 sigma_model = model.predict(x_grid, statistic=statistic)

#                 plt.plot(x_grid, sigma_model, color="black", lw=2, label=f"Model ({statistic})")

#             plt.xlabel("x")
#             plt.ylabel(r"$\sigma_{\rm mag}$")
#             plt.title(f"Photometric error model – {name}")
#             plt.legend()
#             plt.tight_layout()
#             plt.show()
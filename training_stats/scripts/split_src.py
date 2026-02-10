"""
Here lies the different classes and functions that are used to run SPLIT.

1 - Packages

2- 

"""

#--- Import packages ---
import re
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from lephare import statsplot as lsp

#--- sed grid class ---
# The class used to organize a sed grid thanks to their physical parameters and their Ids.
class SED_GRID:
    def __init__(self, list_path):
        """Load SED list."""
        sed_list = []
        with open(list_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                filepath = line.split()[0]
                sp_type = filepath.split('/')[-1].replace('.sed', "")
                sed_list.append(sp_type)
        self.sed_list = sed_list
        self.sed_grid = None

    def build(self):
        """Build the sed grid [Id, Teff, logg, FeH] from filenames."""
        pattern = r"Teff([-\d\.]+)_logg([-\d\.]+)_FeH([-\d\.]+)"
        sed_grid = []
        for Id, s in enumerate(self.sed_list):
            match = re.search(pattern, s)
            if match:
                Teff, logg, FeH = match.groups()
                Teff, logg, FeH = float(Teff), float(logg), float(FeH)
                sed_grid.append([Id, Teff, logg, FeH])
        self.sed_grid = np.array(sed_grid)
        return self.sed_grid

    def plot(self, cmap='gnuplot'):
        """Plot 3D SED grid with color depending on FeH."""
        if self.sed_grid is None:
            raise ValueError("Call make_sed_grid before plot.")

        Teff, logg, FeH = self.sed_grid[:, 1], self.sed_grid[:, 2], self.sed_grid[:, 3]

        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(projection='3d')

        # Couleur selon la métallicité
        sc = ax.scatter(
            Teff, logg, FeH,
            c=FeH, cmap=cmap,
            alpha=0.8, s=10, edgecolor="none"
        )

        # Colorbar
        cbar = plt.colorbar(sc, ax=ax, pad=0.1, shrink=0.8)
        cbar.set_label('[Fe/H]', rotation=270, labelpad=15)

        ax.set_xlabel('Teff [K]')
        ax.set_ylabel('log(g)')
        ax.set_zlabel('[Fe/H]')
        ax.set_title("SED Grid Colored by Metallicity")

        # Optionnel : meilleure orientation par défaut
        ax.view_init(elev=20, azim=-45)

        plt.tight_layout()
        plt.show()


#--- star pdf class ---
# The class used to project the fitting parameter (Chi2) of a source to each sed of the 
# precedently computed sed grid thanks to the sed Ids. Then the chi2_grid can be marginalised
# and turn to probability values. 


class STAR_PDF:
    def __init__(self, sed_grid):
        if sed_grid is None:
            raise ValueError("SED_GRID: call sed_grid_obj.build() before using STAR_PDF.")

        self.sed_grid = sed_grid   # base grid
        self.len_sed_grid = len(self.sed_grid)
        self.chi2_grid = None

    def build_chi2_grid(self, chi2_array):
        """Attach chi2 array to SED grid."""
        if len(chi2_array) != self.len_sed_grid:
            raise ValueError("chi2 array length does not match sed_grid length.")

        self.chi2_grid = np.column_stack((self.sed_grid, chi2_array))
        return self.chi2_grid

    def project_min_chi2(self, xaxis="Teff", yaxis=None, to_prob=False, fixed=None):
        """
        Project Chi² parallel to one or two axes.
        Uses self.chi2_grid automatically.
        """
        if self.chi2_grid is None:
            raise ValueError("Call build_chi2_grid before projecting.")

        grid = self.chi2_grid
        columns = ["Id", "Teff", "logg", "FeH", "Chi2"]
        col_idx = {k: i for i, k in enumerate(columns)}

        # 1) Apply fixed filters
        mask = np.ones(len(grid), dtype=bool)
        if fixed:
            for key, val in fixed.items():
                mask &= np.isclose(grid[:, col_idx[key]], val)
            grid = grid[mask]
            if len(grid) == 0:
                raise ValueError("No grid point matches the given fixed conditions.")

        xcol = col_idx[xaxis]

        if yaxis:   # ----------- 2D PROJECTION -------------
            ycol = col_idx[yaxis]

            X_vals = np.unique(grid[:, xcol])
            Y_vals = np.unique(grid[:, ycol])
            Z = np.full((len(Y_vals), len(X_vals)), np.nan)

            for i, xv in enumerate(X_vals):
                for j, yv in enumerate(Y_vals):
                    mask_xy = np.isclose(grid[:, xcol], xv) & np.isclose(grid[:, ycol], yv)
                    if np.any(mask_xy):
                        Z[j, i] = np.nanmin(grid[mask_xy, -1])

            if to_prob:
                Zmin = np.nanmin(Z)
                Z = np.exp(-0.5 * (Z - Zmin))

            return X_vals, Y_vals, Z

        else:       # ----------- 1D PROJECTION -------------
            X_vals = np.unique(grid[:, xcol])
            Z = np.full(len(X_vals), np.nan)

            for i, xv in enumerate(X_vals):
                mask_x = np.isclose(grid[:, xcol], xv)
                if np.any(mask_x):
                    Z[i] = np.nanmin(grid[mask_x, -1])

            if to_prob:
                Zmin = np.nanmin(Z)
                Z = np.exp(-0.5 * (Z - Zmin))

            return X_vals, Z

    def plot_min_chi2(self, xaxis="Teff", yaxis=None, to_prob=False,
                      full=False, cmap='gnuplot', fixed=None):

        if self.chi2_grid is None:
            raise ValueError("Call build_chi2_grid before plotting.")

        grid = self.chi2_grid  # just to avoid rewriting this everywhere

        # ---------- FULL CORNER-PLOT ----------
        if full:
            axes = ["Teff", "logg", "FeH"]
            n = len(axes)

            # Uniform axes
            limits = {}
            for ax in axes:
                limits[ax] = np.unique(grid[:, {"Teff": 1, "logg": 2, "FeH": 3}[ax]])
            
            fig = plt.figure(figsize=(8, 8))
            gs = plt.GridSpec(n, n, wspace=0.0, hspace=0.0)

            for i in range(n):  # line
                for j in range(n):  # col
                    if j > i:
                        continue
                    ax = fig.add_subplot(gs[i, j])

                    if i == j:
                        # --- PDF 1D ---
                        X, Z = self.project_min_chi2(xaxis=axes[j], to_prob=to_prob, fixed=fixed)
                        ax.plot(X, Z, color='black', lw=1)
                        # ax.fill_between(X, Z, color='C0', alpha=0.4)
                        ax.set_xlim(limits[axes[j]].min(), limits[axes[j]].max())
                        ax.set_ylim(Z.min() * 0.95, Z.max() * 1.05)
                        if i < n - 1:
                            ax.set_xticklabels([])
                        if i != 0:
                            ax.yaxis.tick_right()
                        else:
                            ax.set_xlabel(axes[j])

                    else:
                        # --- PDF 2D ---
                        X, Y, Z = self.project_min_chi2(xaxis=axes[j], yaxis=axes[i], fixed=fixed)
                        if to_prob:
                            Z = np.exp(-0.5 * (Z - np.nanmin(Z)))
                        im = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cmap)
                        ax.set_xlim(limits[axes[j]].min(), limits[axes[j]].max())
                        ax.set_ylim(limits[axes[i]].min(), limits[axes[i]].max())
                        if i < n - 1:
                            ax.set_xticklabels([])
                        else:
                            ax.set_xlabel(axes[j])
                        if j > 0:
                            ax.set_yticklabels([])
                        else:
                            ax.set_ylabel(axes[i])

            # Colorbar
            cbar_ax = fig.add_axes([1, 0.15, 0.02, 0.7])
            cbar = fig.colorbar(im, cax=cbar_ax)
            cbar.set_label("Probability" if to_prob else "min(χ²)")
            plt.subplots_adjust(left=0.08, right=0.9, top=0.9, bottom=0.08)
            plt.show()
            return


        # --- Single plot mode
        if yaxis:
            X, Y, Z = self.project_min_chi2(xaxis, yaxis, to_prob=to_prob, fixed=fixed)
            plt.figure(figsize=(6, 4))
            plt.pcolormesh(X, Y, Z, shading='auto')
            plt.colorbar(label="Probability" if to_prob else "min(Chi²)")
            plt.xlabel(xaxis)
            plt.ylabel(yaxis)
            plt.title(f"Minimized Chi² map")
            plt.show()
        else:
            X, Z = self.project_min_chi2(xaxis, to_prob=to_prob, fixed=fixed)
            plt.figure(figsize=(6, 3))
            plt.plot(X, Z, marker='o', lw=1)
            plt.xlabel(xaxis)
            plt.ylabel("Probability" if to_prob else "min(Chi²)")
            plt.title("Minimized Chi² PDF")
            plt.grid(True)
            plt.show()

    def plot_full_chi2(self, cmap='gnuplot'):
        """Plot 3D SED grid with color depending on FeH."""
        grid = self.chi2_grid  # just to avoid rewriting this everywhere
        # print(grid)
        Teff, logg, FeH = grid[:, 1], grid[:, 2], grid[:, 3]
        Chi2 = grid[:, 4]
        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(projection='3d')
        # Coulor as a function of Chi2
        sc = ax.scatter(
            Teff, logg, FeH,
            c=Chi2, cmap=cmap,
            alpha=0.8, s=10, edgecolor="none"
        )

        # Colorbar
        cbar = plt.colorbar(sc, ax=ax, pad=0.1, shrink=0.8)
        cbar.set_label('Chi2', rotation=270, labelpad=15)

        ax.set_xlabel('Teff [K]')
        ax.set_ylabel('log(g)')
        ax.set_zlabel('[Fe/H]')
        ax.set_title("SED Grid Colored by Metallicity")
 
        # Optionnel : meilleure orientation par défaut
        ax.view_init(elev=20, azim=-45)
 
        plt.tight_layout()
        plt.show()

    def best_output(self):
        grid = self.chi2_grid
        best_sed_id = np.argmin(grid[:,4])
        chi2_best = min(grid[:,4])
        Teff_best = grid[best_sed_id, 1]
        logg_best = grid[best_sed_id, 2]
        FeH_best = grid[best_sed_id, 3]
        
        chisec_grid = np.delete(grid, best_sed_id, axis=0)
        sec_best_sed_id = np.argmin(chisec_grid[:,4])
        sec_chi2_best = min(chisec_grid[:,4])
        sec_Teff_best = chisec_grid[sec_best_sed_id, 1]
        sec_logg_best = chisec_grid[sec_best_sed_id, 2]
        sec_FeH_best = chisec_grid[sec_best_sed_id, 3]
        sec_best_sed_id = sec_best_sed_id + 1 if sec_best_sed_id >= best_sed_id else sec_best_sed_id
        return [best_sed_id, Teff_best, logg_best, FeH_best, chi2_best, 
                sec_best_sed_id, sec_Teff_best, sec_logg_best, sec_FeH_best, sec_chi2_best]


#--- star pdf analyzer class---
# The class use to compute the different statistical metrics of a star pdf
# in 1, 2 or 3D. The main tool to determine the source type from star pdfs.

from scipy.signal import find_peaks
from scipy.interpolate import interp1d, griddata, RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from numba import njit

class STAR_PDF_ANALYZER:
    """
    Analyzer for chi2_grid (Nx5 array: [Id, Teff, logg, FeH, Chi2]).
    Option A: interpolation is done on chi2 (not on probabilities).
    - to_prob_3D: if True, compute global P3 (full-grid PDF) once at init
    - to_prob: if True, projectors return normalized probabilities per projection
    """
    def __init__(self, chi2_grid, to_prob=True, to_prob_3D=False,
                oversample=1.0, nmax=800):

        chi2_grid = np.asarray(chi2_grid)
        if chi2_grid.ndim != 2 or chi2_grid.shape[1] < 5:
            raise ValueError("chi2_grid must be shape (N,>=5)")

        self.grid = chi2_grid
        self.to_prob = bool(to_prob)
        self.to_prob_3D = bool(to_prob_3D)

        self.col = {"Id": 0, "Teff": 1, "logg": 2, "FeH": 3, "Chi2": 4}

        # precomputation
        self.axis_vals = {}
        self.axis_idx = {}

        for ax in ["Teff", "logg", "FeH"]:
            vals = np.unique(self.grid[:, self.col[ax]])
            vals.sort()
            self.axis_vals[ax] = vals

            # map grid values -> integer indices
            idx = np.searchsorted(vals, self.grid[:, self.col[ax]])
            self.axis_idx[ax] = idx.astype(np.int32)

        self._proj1d_cache = {}
        self._proj2d_cache = {}

        # global chi2
        chi2 = self.grid[:, self.col["Chi2"]]
        self.chi2_min_global = float(np.nanmin(chi2))

        self.P3 = None
        if self.to_prob_3D:
            self._make_prob3d()

        self._regular_cache = {}
        self._oversample = float(oversample)
        self._nmax = int(nmax)

    # --------------------------
    # Internal helpers: PDFs
    # --------------------------
    def _build_fixed_mask(self, fixed):
        if not fixed:
            return np.ones(len(self.grid), dtype=bool)

        mask = np.ones(len(self.grid), dtype=bool)
        for ax, val in fixed.items():
            i = np.where(np.isclose(self.axis_vals[ax], val))[0]
            if len(i) == 0:
                return np.zeros(len(self.grid), dtype=bool)
            mask &= (self.axis_idx[ax] == i[0])
        return mask

    def _make_prob3d(self):
        """Compute full-grid probability P3 normalized."""
        chi2 = self.grid[:, self.col["Chi2"]].astype(float)
        dchi2 = chi2 - np.nanmin(chi2)
        with np.errstate(over="ignore", invalid="ignore"):
            P = np.exp(-0.5 * dchi2)
        P[np.isnan(P)] = 0.0
        s = P.sum()
        if s <= 0:
            P = np.zeros_like(P)
        else:
            P = P / s
        self.P3 = P
        return P
    # --------------------------
    # Regular grid utilities
    # --------------------------
    def _unique_sorted_axis(self, axis):
        vals = np.unique(self.grid[:, self.col[axis]])
        vals = np.sort(vals)
        return vals

    def _is_regular_axis(self, vals, rtol=1e-8):
        """Return True if axis values are regularly spaced (within rtol)."""
        if vals.size < 3:
            return True
        d = np.diff(vals)
        # consider only positive diffs
        dpos = d[d > 0]
        if dpos.size == 0:
            return True
        return np.allclose(dpos, dpos[0], rtol=rtol, atol=0)

    def _optimal_regular_grid(self, vals):
        """Build regular grid using dx = min positive difference (pas minimum observed)."""
        vals = np.sort(np.unique(vals))
        if vals.size < 2:
            return vals.copy()
        d = np.diff(vals)
        dpos = d[d > 0]
        if dpos.size == 0:
            return vals.copy()
        dx_min = float(np.min(dpos))
        xmin, xmax = float(vals[0]), float(vals[-1])
        # number of steps based on dx_min and oversampling factor
        n = int(np.ceil(self._oversample * (xmax - xmin) / dx_min)) + 1
        n = max(2, min(n, self._nmax))
        grid = np.linspace(xmin, xmax, n)
        return grid

    def _get_regular_axis(self, axis):
        """Return cached regular grid for axis or build it. Also return flag regular."""
        if axis in self._regular_cache:
            return self._regular_cache[axis]["vals"], self._regular_cache[axis]["regular"]

        vals = self._unique_sorted_axis(axis)
        regular = self._is_regular_axis(vals)
        if regular:
            grid = vals.copy()
        else:
            grid = self._optimal_regular_grid(vals)
        self._regular_cache[axis] = {"vals": grid, "regular": regular}
        return grid, regular

    # --------------------------
    # Regularize (interpolate) 1D and 2D on chi2 then convert if needed
    # --------------------------
    def _regularize_1d_chi2(self, X, P, axis):
        good = ~np.isnan(P)
        if not np.any(good):
            return np.array([]), np.array([])

        Xg = X[good]
        Pg = P[good]

        grid, is_reg = self._get_regular_axis(axis)

        f = interp1d(Xg, Pg, bounds_error=False, fill_value=np.nan)
        return grid, f(grid)

    def _regularize_2d_chi2(self, X, Y, P, xaxis, yaxis):
        Xr, _ = self._get_regular_axis(xaxis)
        Yr, _ = self._get_regular_axis(yaxis)

        # interpolation séparable (beaucoup plus rapide que griddata)
        Px = np.array([
            np.interp(Xr, X, row, left=np.nan, right=np.nan)
            for row in P
        ])

        Pr = np.array([
            np.interp(Yr, Y, Px[:, i], left=np.nan, right=np.nan)
            for i in range(len(Xr))
        ]).T

        return Xr, Yr, Pr

    # --------------------------
    # Projection on original grid (min chi2 over other dims)
    # --------------------------
    def _project_min_chi2_1d_raw(self, axis, fixed=None):
        idx = self.axis_idx[axis]
        chi2 = self.grid[:, self.col["Chi2"]]

        mask = self._build_fixed_mask(fixed)
        idx = idx[mask]
        chi2 = chi2[mask]

        n = len(self.axis_vals[axis])
        out = np.full(n, np.inf)
        np.minimum.at(out, idx, chi2)

        out[out == np.inf] = np.nan
        return self.axis_vals[axis], out


    def _project_min_chi2_2d_raw(self, xaxis, yaxis, fixed=None):
        xi = self.axis_idx[xaxis]
        yi = self.axis_idx[yaxis]
        chi2 = self.grid[:, self.col["Chi2"]]

        mask = self._build_fixed_mask(fixed)
        xi = xi[mask]
        yi = yi[mask]
        chi2 = chi2[mask]

        nx = len(self.axis_vals[xaxis])
        ny = len(self.axis_vals[yaxis])

        P = np.full((ny, nx), np.inf)
        np.minimum.at(P, (yi, xi), chi2)

        P[P == np.inf] = np.nan
        return self.axis_vals[xaxis], self.axis_vals[yaxis], P

    # --------------------------
    # Public projectors that return either chi2 or probability on regularized grid
    # --------------------------
    def project_1d(self, axis, to_prob=None, fixed=None):
        """
        Return Xr, P where:
         - Xr: regularized axis grid
         - P: if to_prob True -> normalized probability (sum=1); else -> chi2 (min chi2 per X)
        """
        key = (axis, bool(to_prob), None if fixed is None else tuple(sorted(fixed.items())))
        if key in self._proj1d_cache:
            return self._proj1d_cache[key]
        if to_prob is None:
            to_prob = self.to_prob

        Xraw, Pchi2_raw = self._project_min_chi2_1d_raw(axis, fixed=fixed)
        if Xraw.size == 0:
            return np.array([]), np.array([])

        Xr, Pchi2 = self._regularize_1d_chi2(Xraw, Pchi2_raw, axis)

        if Xr.size == 0:
            return np.array([]), np.array([])

        if to_prob:
            # convert chi2 -> probability with local minimum
            if np.all(np.isnan(Pchi2)):
                return Xr, np.zeros_like(Pchi2)
            dchi = Pchi2 - np.nanmin(Pchi2)
            with np.errstate(over="ignore", invalid="ignore"):
                P = np.exp(-0.5 * dchi)
            P[np.isnan(P)] = 0.0
            s = np.nansum(P)
            if s > 0:
                P = P / s
            else:
                P = np.zeros_like(P)
            self._proj1d_cache[key] = (Xr, P)
            return Xr, P
        else:
            self._proj1d_cache[key] = (Xr, P)
            return Xr, Pchi2

    def project_2d(self, xaxis, yaxis, to_prob=None, fixed=None):
        """
        Return Xr, Yr, P where:
         - Xr, Yr: regularized axis grids
         - P: if to_prob True -> normalized probability (sum ~1); else -> chi2 map (min chi2 per cell)
        """
        key = (xaxis, yaxis, bool(to_prob), None if fixed is None else tuple(sorted(fixed.items())))
        if key in self._proj2d_cache:
            return self._proj2d_cache[key]
        if to_prob is None:
            to_prob = self.to_prob

        Xraw, Yraw, Pchi2_raw = self._project_min_chi2_2d_raw(xaxis, yaxis, fixed=fixed)
        if Xraw.size == 0 or Yraw.size == 0:
            return np.array([]), np.array([]), np.array([[]])

        Xr, Yr, Pchi2 = self._regularize_2d_chi2(Xraw, Yraw, Pchi2_raw, xaxis, yaxis)

        if Xr.size == 0 or Yr.size == 0:
            return np.array([]), np.array([]), np.array([[]])

        if to_prob:
            if np.all(np.isnan(Pchi2)):
                return Xr, Yr, np.zeros_like(Pchi2)
            dchi = Pchi2 - np.nanmin(Pchi2)
            with np.errstate(over="ignore", invalid="ignore"):
                P = np.exp(-0.5 * dchi)
            P[np.isnan(P)] = 0.0
            s = np.nansum(P)
            if s > 0:
                P = P / s
            else:
                P = np.zeros_like(P)
            self._proj2d_cache[key] = (Xr, Yr, P)
            return Xr, Yr, P
        else:
            self._proj2d_cache[key] = (Xr, Yr, P)
            return Xr, Yr, Pchi2

    # --------------------------
    # Metrics (1D)
    # --------------------------
    def npeaks_1d(self, axis, threshold=0.2, fixed=None):
        _, P = self.project_1d(axis, to_prob=self.to_prob, fixed=fixed)
        if P.size == 0:
            return -99

        arr = P if self.to_prob else -P
        m = np.nanmax(arr)
        if not np.isfinite(m) or m <= 0:
            return 0

        peaks, _ = find_peaks(arr, height=threshold * m)
        return int(peaks.size)

    def peakratio_1d(self, axis, fixed=None):
        _, P = self.project_1d(axis, to_prob=self.to_prob, fixed=fixed)
        if P.size == 0:
            return -99

        mx = np.nanmax(P)
        return float(np.nanmean(P) / mx) if mx > 0 else -99

    def std_1d(self, axis, fixed=None):
        X, P = self.project_1d(axis, to_prob=self.to_prob, fixed=fixed)
        if P.size == 0:
            return -99

        if self.to_prob:
            W = P
        else:
            dchi = P - np.nanmin(P)
            W = np.exp(-0.5 * dchi)
            W[~np.isfinite(W)] = 0.0

        s = W.sum()
        if s <= 0:
            return -99
        W /= s
        mu = np.sum(X * W)
        return float(np.sqrt(np.sum(W * (X - mu) ** 2)))

    # --------------------------
    # Metrics (2D)
    # --------------------------
    def cov_2d(self, xaxis, yaxis, fixed=None):
        X, Y, P = self.project_2d(xaxis, yaxis, to_prob=True, fixed=fixed)
        if P.size == 0:
            return (-99, -99, -99, -99)

        # adapt proba to the regularized grid
        dx = np.gradient(X)
        dy = np.gradient(Y)
        DX, DY = np.meshgrid(dx, dy)
        W = P * DX * DY

        tot = W.sum()
        if tot <= 0:
            return (-99, -99, -99, -99)
        # make prob grid
        W /= tot
        Xg, Yg = np.meshgrid(X, Y)

        # Compute covariances
        mx = np.sum(W * Xg)
        my = np.sum(W * Yg)

        cov_xx = np.sum(W * (Xg - mx) ** 2)
        cov_yy = np.sum(W * (Yg - my) ** 2)
        cov_xy = np.sum(W * (Xg - mx) * (Yg - my))

        cov_mat = np.array([[cov_xx, cov_xy],
                            [cov_xy, cov_yy]])

        eig = np.linalg.eigvalsh(cov_mat)
        eig = np.maximum(eig, 1e-300)
        axis_ratio = float(np.sqrt(eig[1] / eig[0]))

        return float(cov_xx), float(cov_yy), float(cov_xy), axis_ratio

    def area_highprob_2d(self, xaxis, yaxis, threshold=0.5, fixed=None):
        X, Y, P = self.project_2d(xaxis, yaxis, to_prob=True, fixed=fixed)
        if P.size == 0:
            return -99.0

        Pmax = np.nanmax(P)
        if Pmax <= 0:
            return 0.0
        # adapt proba to the regularized grid
        dx = np.gradient(X)
        dy = np.gradient(Y)
        DX, DY = np.meshgrid(dx, dy)
        A = DX * DY
        
        valid = np.isfinite(P)
        high = P >= threshold * Pmax

        Atot = np.sum(A[valid])
        if Atot <= 0:
            return -99.0

        return float(np.sum(A[high & valid]) / Atot)


    # --------------------------
    # Global metrics
    # --------------------------
    def entropy(self):
        if self.P3 is None:
            # if P3 not precomputed, compute from chi2 on the fly
            self._make_prob3d()
        P = np.array(self.P3, dtype=float)
        mask = P > 0
        if np.sum(mask) == 0:
            return -99.0
        H = -np.sum(P[mask] * np.log(P[mask]))
        return float(H)

    def concentration(self, delta_chi2=1.0):
        chi2 = self.grid[:, self.col["Chi2"]].astype(float)
        if self.P3 is None:
            self._make_prob3d()
        mask = chi2 <= (self.chi2_min_global + float(delta_chi2))
        if np.sum(mask) == 0:
            return 0.0
        return float(np.sum(self.P3[mask]))

    def delta_chi2(self):
        chi2 = np.array(self.grid[:, self.col["Chi2"]], dtype=float)
        if chi2.size < 2:
            return -99.0
        sorted_chi = np.sort(chi2)
        return float(sorted_chi[1] - sorted_chi[0])

    # --------------------------
    # Main entry: compute all features
    # --------------------------
    def compute_all(self, fixed=None):
        """
        Compute all requested features and return as dict.
        fixed: optional dict to restrict projection (e.g. {"FeH":0.0})
        """
        feats = {}

        for ax in ["Teff", "logg", "FeH"]:
            feats[f"npeaks_{ax}"] = self.npeaks_1d(ax, threshold=0.2, fixed=fixed)
            feats[f"peakratio_{ax}"] = self.peakratio_1d(ax, fixed=fixed)
            feats[f"std_{ax}"] = self.std_1d(ax, fixed=fixed)

        pairs = [("Teff", "logg"), ("Teff", "FeH"), ("logg", "FeH")]
        for x, y in pairs:
            cov_xx, cov_yy, cov_xy, axis_ratio = self.cov_2d(x, y, fixed=fixed)
            feats[f"cov_{x}_{y}_xx"] = cov_xx
            feats[f"cov_{x}_{y}_yy"] = cov_yy
            feats[f"cov_{x}_{y}_xy"] = cov_xy
            feats[f"axisratio_{x}_{y}"] = axis_ratio
            feats[f"area_{x}_{y}_p50"] = self.area_highprob_2d(x, y, threshold=0.5, fixed=fixed)
            feats[f"area_{x}_{y}_p90"] = self.area_highprob_2d(x, y, threshold=0.9, fixed=fixed)

        feats["entropy"] = self.entropy()
        feats["concentration_dchi2_1"] = self.concentration(delta_chi2=1.0)
        feats["delta_chi2"] = self.delta_chi2()
        return feats


#--- sample analyzer class
# Run star_pdf_analyzer metrics on a sample

from joblib import Parallel, delayed

class SAMPLE_ANALYZER:
    def __init__(self, chi2_path, sed_grid_path, nrows=None, n_jobs=1):

        # --- Load raw file as strings ---
        raw = pd.read_csv(chi2_path, sep=r"\s+", header=None,
                        dtype="str", comment="#", nrows=nrows)
        raw = raw.to_numpy()
        self.ids = raw[:, 0].astype(np.int64)
        self.chi2_values = raw[:, 1:].astype(float)

        # SED grid (constant for all sources)
        self.sed_grid_obj = SED_GRID(sed_grid_path)
        self.sed_grid = self.sed_grid_obj.build()

        self.n_jobs = int(n_jobs)

    def best_fit_distribution(self, n_rows=None):

        star_pdf = STAR_PDF(self.sed_grid)
        data = []

        N = n_rows if n_rows is not None else len(self.ids)

        for i in range(N):
            _ = star_pdf.build_chi2_grid(self.chi2_values[i])
            df = star_pdf.best_output()
            data.append([self.ids[i]] + df)

        columns = [
            'source_id',
            'best_sed_id', 'Teff_best', 'logg_best', 'FeH_best', 'chi2_best',
            'sec_best_sed_id', 'sec_Teff_best', 'sec_logg_best', 'sec_FeH_best', 'sec_chi2_best'
        ]

        self.sample_df = pd.DataFrame(data, columns=columns)
        return self.sample_df

    def stat_feature_distribution(self, n_rows=None):

        N = n_rows if n_rows is not None else len(self.ids)
        star_pdf = STAR_PDF(self.sed_grid)

        # ---- worker function 
        def onesource(i):
            chi2_grid = star_pdf.build_chi2_grid(self.chi2_values[i])
            analyzer = STAR_PDF_ANALYZER(chi2_grid)
            feats = analyzer.compute_all()
            feats["source_id"] = self.ids[i]
            return feats

        # ---- run
        if self.n_jobs == 1:
            rows = [onesource(i) for i in range(N)]
        else:
            rows = Parallel(n_jobs=self.n_jobs, backend="loky")(
                delayed(onesource)(i) for i in range(N)
            )

        # ---- build dataframe with fixed column order
        columns = ["source_id"] + [k for k in rows[0].keys() if k != "source_id"]
        self.feature_df = pd.DataFrame(rows, columns=columns)

        return self.feature_df


    def plot_best(self, axis='chi2_best'):
        plt.figure()
        plt.hist(self.sample_df[axis].values, bins=50)
        plt.xlabel(axis)
        plt.ylabel("N")
        plt.show()

#--- 
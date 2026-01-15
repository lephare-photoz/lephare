import os

import numpy as np
from scipy.interpolate import interp1d

import lephare as lp

__all__ = [
    "multiply_on_grids",
    "compute_model_reddening",
]


def multiply_on_grids(x1, f1, x2, f2, x=None):
    """
    Interpolate two functions onto the same grid and multiply them elementwise.

    If no target grid `x` is provided, a grid spanning the union of `x1` and `x2`
    is automatically created with the same number of points as `x2`.

    Parameters
    ----------
    x1 : array_like
        Independent variable for the first function.
    f1 : array_like
        Function values corresponding to `x1`.
    x2 : array_like
        Independent variable for the second function.
    f2 : array_like
        Function values corresponding to `x2`.
    x : array_like, optional
        Grid on which to interpolate both functions. If None, an automatic grid
        spanning the min/max of `x1` and `x2` is used.

    Returns
    -------
    numpy.ndarray
        A 2D array with two rows: the first row is the grid `x`, and the second
        row is the product of the interpolated functions `f1 * f2`.

    Notes
    -----
    - Currently uses linear interpolation with `fill_value=0.0` outside the original ranges.
    - This is a temporary implementation; replacing it with native LePHARE classes is recommended.
    """

    if x is None:
        print("Using user defined x grid")
        xmin = min(x1.min(), x2.min())
        xmax = max(x1.max(), x2.max())
        if xmax < xmin:
            # If no overlap return the input filter
            return np.array([x2, f2])
        x = np.linspace(xmin, xmax, len(x2))
    f1i = interp1d(x1, f1, bounds_error=False, fill_value=0.0)(x)
    f2i = interp1d(x2, f2, bounds_error=False, fill_value=0.0)(x)

    return np.array([x, f1i * f2i])  # fftconvolve(f1i, f2i, mode="same") * dx


def compute_model_reddening(config, verbose=False):
    """
    Compute model-dependent reddening following Galametz et al. (2017).

    For each SED in the model library, the function:
    1. Observes the SED (redshift, extinction, etc.).
    2. Convolves it with a set of filters.
    3. Computes the reddening relative to a baseline filter.

    The resulting values can then be applied per model during the C++ fit process.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing model parameters and file paths
        required to generate the SED library and filters.
    verbose : bool, optional
        If True, print additional information during computation. Default is False.

    Returns
    -------
    numpy.ndarray
        A 2D array of reddening values with shape `(n_models, n_filters)`,
        where each row corresponds to a model and each column to a filter.

    Notes
    -----
    - The function assumes the presence of a valid LePHARE environment and
      paths defined by `LEPHAREWORK` or `LEPHAREDIR`.
    - Supports Cardelli, Calzetti, and SMC Prevot extinction laws.
    - Observed SEDs are modified in place using the `observe` function.
    - The output array `values[i, j]` corresponds to the reddening for model `i`
      through filter `j`.
    """
    keymap = lp.all_types_to_keymap(config)

    photz = lp.PhotoZ(keymap)
    # Why are these required?
    keymap["EXT_CURVE"] = lp.keyword("EXT_CURVE", "NONE")
    keymap["GAL_CURVE"] = lp.keyword("GAL_CURVE", "CARDELLI")
    keymap["t"] = lp.keyword("t", "Star")
    star_mag = lp.StarMag(keymap)
    keymap["t"] = lp.keyword("t", "Gal")
    gal_mag = lp.GalMag(keymap)
    keymap["t"] = lp.keyword("t", "QSO")
    qso_mag = lp.QSOMag(keymap)

    # Get the parameters this code is taken from the lp.FitExt runner
    filters = keymap["FILTER_FILE"].split_string("unknown", 1)[0]
    if not os.path.isabs(filters):
        filters = os.path.join(os.environ["LEPHAREWORK"], "filt", filters + ".dat")
    galec = keymap["GAL_CURVE"].split_string("CARDELLI", 1)[0]

    all_filters = lp.GalMag.read_flt(filters)
    original_filters = list([f.data() for f in all_filters])  # These should remain intact always
    print("Number of filters: ", len(original_filters))
    values = np.empty((len(photz.fullLib), len(original_filters)))

    for i, model in enumerate(photz.fullLib):
        for j, filt in enumerate(original_filters):
            if model.is_star():
                mag = star_mag
            elif model.is_gal():
                mag = gal_mag
            elif model.is_qso():
                mag = qso_mag
            # print(i,j)
            # observed_sed=observe(photz, i, mag)
            model.lamb_flux = photz.fullLib[model.index_z0].lamb_flux
            model.generate_spectra(photz.zLib[i], 1.0, mag.opaAll)
            model_filt_product = multiply_on_grids(
                model.data()[0], model.data()[1], filt[0], filt[1], x=filt[0]
            )

            # make a oneflt instance from the updated
            one_filt = all_filters[j]
            for n, lt in enumerate(one_filt.lamb_trans):
                lt.val = model_filt_product[1][n]
            # atmospheric extinction too?
            # returns all_filters, aint, albdav, albd
            if galec == "CARDELLI":
                # If cardelli (hardcoded) with Rv=3.1 by  default
                # output A(lbd)/Av->A(lbd)/(Rv*E(B-V))->A(lbd)/E(B-V)=Rv*A(lbd)/Av
                albdav = lp.cardelli_ext(one_filt)
                albd = albdav * 3.1
            else:
                if not os.path.isabs(galec):
                    galec = os.path.join(os.environ["LEPHAREDIR"], "ext", galec)
                galactic_ext = lp.ext(galec, 1)
                galactic_ext.read(galec)

                #  Rv=3.1 except for Calzetti law (4.05) and SMC Prevot (2.72)
                rv = 3.1
                if "SMC_prevot" in galec:
                    rv = 2.72
                    if "calzetti" in galec:
                        rv = 4.05
                if verbose:
                    print(f"assuming Rv={rv} for this Extinction law {galec}")
                albd = lp.compute_filter_extinction(one_filt, galactic_ext)
            values[i, j] = albd

    # Replace nans with zeros (probably due to non overlap between filter and SED)
    values = np.nan_to_num(values, nan=0.0)
    return values

import numpy as np

import lephare as lp

__all__ = ["compute_model_reddening", "compute_band_pass_correction"]


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
    albd_lib = np.array([g.milky_way_extinction for g in photz.fullLib])
    return np.nan_to_num(albd_lib, nan=0.0)


def compute_band_pass_correction(config, model_number=15, type="star", verbose=False):
    """
    Compute model-dependent band pass correction following Galametz et al. (2017).

    For each SED in the model library, the function:
    1. Observes the SED (redshift, extinction, etc.).
    2. Convolves it with the Johnson B and V filters.
    3. Compares the E(B-V) value relative to a B5 star that defines the E(B-V) value
    4. Returns E(B-V)_SED / E(B-V)_B5 for each model, which is the band pass correction.

    The resulting values can then be applied per model during the C++ fit process.

    Parameters
    ----------
    config : dict
        Configuration dictionary containing model parameters and file paths
        required to generate the SED library and filters.
    model_number : int
        The model number for the reference model. Defaults to 15 in sed/STARS/STAR_MOD_ALL.list
    type : str
        The type of the reference model. Defaults to 'star'.
    verbose : bool, optional
        If True, print additional information during computation. Default is False.

    Returns
    -------
    numpy.ndarray
        A 1D array of band pass correction values with shape `(n_models)`,
        where each row corresponds to a model.

    Notes
    -----
    - The function assumes the presence of a valid LePHARE environment and
      paths defined by `LEPHAREWORK` or `LEPHAREDIR`.
    - Supports Cardelli, Calzetti, and SMC Prevot extinction laws.
    - Observed SEDs are modified in place using the `observe` function.
    - The output array `values[i]` corresponds to the band pass correction for model `i`.
    """
    keymap = lp.all_types_to_keymap(config)
    # Make sure all the config values will work with just B and V
    keymap["FILTER_LIST"] = lp.keyword("FILTER_LIST", "std/B.pb,std/V.pb")
    keymap["FILTER_FILE"] = lp.keyword("FILTER_FILE", "BV_FILTERS_FOR_BAND_PASS_CORRECTION")
    keymap["FILTER_CALIB"] = lp.keyword("FILTER_CALIB", "0,0")
    lp.prepare(keymap)
    photz = lp.PhotoZ(keymap)
    albd_lib = np.array([g.milky_way_extinction for g in photz.fullLib])
    b5_idx = None
    for n, sed in enumerate(photz.fullLib):
        if getattr(sed, f"is_{type}")() and (sed.nummod == model_number):
            b5_idx = n
    # A_B - A_V for the SED, and for the B5 star. The ratio of these
    # is the band pass correction factor to apply to the model reddening.
    ebmv = albd_lib.T[0] - albd_lib.T[1]
    band_pass_correction = ebmv / (albd_lib[b5_idx][0] - albd_lib[b5_idx][1])
    # zeros = np.isclose(values.T[0], 0.0)
    # zeros |= np.isclose(values.T[0], 0.0)
    # band_pass_correction[zeros] = (
    #     1.0  # If there is no reddening in either filter, set the band pass correction to 1 (no correction)
    # )
    return band_pass_correction

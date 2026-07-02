import os

import numpy as np

import lephare as lp

__all__ = ["compute_model_reddening", "compute_band_pass_correction", "classic_extinction_values"]


def compute_model_reddening(config):
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
    keymap = lp.all_types_to_keymap({**config, "APPLY_MW_EXTINCTION": "GALAMETZ"})
    lp.prepare(keymap)
    photz = lp.PhotoZ(keymap)
    albd_lib = np.array([g.milky_way_extinction for g in photz.fullLib])
    return np.nan_to_num(albd_lib, nan=0.0)


def compute_band_pass_correction(config):
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
    # Set to calculate the values if not set
    keymap = lp.all_types_to_keymap({**config, "APPLY_MW_EXTINCTION": "GALAMETZ"})
    lp.prepare(keymap)
    photz = lp.PhotoZ(keymap)
    band_pass_correction = np.array([g.band_pass_correction for g in photz.fullLib])

    # After prepre it is not
    return band_pass_correction


def classic_extinction_values(config, verbose=False):
    """Calculate the extinction values for a set of filters using the classic method

    Parameters
    ==========

    config : dict
        Configuration dictionary containing model parameters and file paths
        required to generate the SED library and filters.

    verbose : bool
        Increase onscreen verbosity

    Returns
    =======

    all_filters : list of lephare.flt
        The list of filter objects
    aint : np.array
        Atmospheric extinction in each filter (mag/airmass)
    albdav : np.array
        Galactic extinction in each filter A(lbd)/Av
    albd : np.array
        Galactic extinction in each filter A(lbd)/E(B-V)
    """
    keymap = lp.all_types_to_keymap(config)
    # Get the parameters
    filters = keymap["FILTER_FILE"].split_string("unknown", 1)[0]
    if not os.path.isabs(filters):  # pragma no cover
        filters = os.path.join(os.environ["LEPHAREWORK"], "filt", filters + ".dat")
    atmec = keymap["EXT_ATMOS_CURVE"].split_string("NONE", 1)[0]
    galec = keymap["EXT_MW_CURVE"].split_string("CARDELLI", 1)[0]

    all_filters, aint, albdav, albd = lp.calculate_extinction_values(filters, atmec, galec, verbose=verbose)
    return all_filters, aint, albdav, albd

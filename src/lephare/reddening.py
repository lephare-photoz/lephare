import os

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
    # Why are these required?
    if "EXT_CURVE" not in keymap:
        keymap["EXT_CURVE"] = lp.keyword("EXT_CURVE", "NONE")
    if "GAL_CURVE" not in keymap:
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

    for j, _ in enumerate(original_filters):
        for i, model in enumerate(photz.fullLib):
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

            # make a oneflt instance from the updated
            one_filt = all_filters[j]
            # for n, lt in enumerate(one_filt.lamb_trans):
            #     lt.val = model_filt_product[1][n]
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
                # albd = lp.compute_filter_sed_extinction(one_filt, galactic_ext)
                albd = lp.compute_filter_sed_extinction(one_filt, galactic_ext, model)
            values[i, j] = albd

    # Replace nans with zeros (probably due to non overlap between filter and SED)
    values = np.nan_to_num(values, nan=0.0)
    return values


def compute_band_pass_correction(config, b5_model=15, verbose=False):
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
    b5_model : int
        The model number for the B5 star. Defaults to 15 in sed/STARS/STAR_MOD_ALL.list
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

    keymap["FILTER_LIST"] = lp.keyword("FILTER_LIST", "std/B.pb,std/V.pb")
    keymap["FILTER_FILE"] = lp.keyword("FILTER_FILE", "BV_FILTERS_FOR_BAND_PASS_CORRECTION")
    keymap["FILTER_CALIB"] = lp.keyword("FILTER_CALIB", "0,0")
    filter_lib = lp.FilterSvc.from_keymap(keymap)
    # Get location to store filter files
    filter_output = os.path.join(os.environ["LEPHAREWORK"], "filt", keymap["FILTER_FILE"].value)
    # Write filter files
    lp.write_output_filter(filter_output + ".dat", filter_output + ".doc", filter_lib)

    photz = lp.PhotoZ(keymap)
    # Why are these required?
    if "EXT_CURVE" not in keymap:
        keymap["EXT_CURVE"] = lp.keyword("EXT_CURVE", "NONE")
    if "GAL_CURVE" not in keymap:
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
    b5_idx = None
    for j, _ in enumerate(original_filters):
        for i, model in enumerate(photz.fullLib):
            if model.is_star():
                mag = star_mag
                if model.nummod == b5_model:
                    b5_idx = i
            elif model.is_gal():
                mag = gal_mag
            elif model.is_qso():
                mag = qso_mag
            # print(i,j)
            # observed_sed=observe(photz, i, mag)
            model.lamb_flux = photz.fullLib[model.index_z0].lamb_flux
            model.generate_spectra(photz.zLib[i], 1.0, mag.opaAll)

            # make a oneflt instance from the updated
            one_filt = all_filters[j]
            # for n, lt in enumerate(one_filt.lamb_trans):
            #     lt.val = model_filt_product[1][n]
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
                # albd = lp.compute_filter_extinction(one_filt, galactic_ext)
                albd = lp.compute_filter_sed_extinction(one_filt, galactic_ext, model)
            values[i, j] = albd

    # A_B - A_V for the SED, and for the B5 star. The ratio of these
    # is the band pass correction factor to apply to the model reddening.
    bmv = values.T[0] - values.T[1]
    band_pass_correction = bmv / (values[b5_idx][0] - values[b5_idx][1])
    # zeros = np.isclose(values.T[0], 0.0)
    # zeros |= np.isclose(values.T[0], 0.0)
    # band_pass_correction[zeros] = (
    #     1.0  # If there is no reddening in either filter, set the band pass correction to 1 (no correction)
    # )
    return band_pass_correction

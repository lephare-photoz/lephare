import datetime

import numpy as np

import lephare as lp

__all__ = ["process", "table_to_data"]


def process(config, input, col_names=None, standard_names=False, filename=None):
    """Run all required steps to produce photometric redshift estimates

    Parameters
    ==========
    config : dict of lephare.keyword
        The configuration for the run
    input : astropy.table.Table
        The input table which must satisfy
    col_names : list
        Input catalogue column names. We will use ordering to determine meaning
    standard_names : bool
        If true we assume standard names.
    filename : str
        Output file name for the output catalogue.

    Returns
    =======
    output : astropy.table.Table
        The output table.
    pdf : np.array
        Array of pdfs for each object
    """
    id, flux, flux_err, context, zspec, string_data = table_to_data(
        config, input, col_names=col_names, standard_names=standard_names
    )
    ng = len(id)
    n_filters = len(config["FILTER_LIST"].value.split(","))
    print(f"Processing {ng} objects with {n_filters} bands")
    photz = lp.PhotoZ(config)
    # Loop over all ng galaxies!
    srclist = []
    for i in range(ng):
        one_obj = lp.onesource(i, photz.gridz)
        one_obj.readsource(str(i), flux[i], flux_err[i], context[i], zspec[i], string_data[i])
        photz.prep_data(one_obj)
        srclist.append(one_obj)

    # If AUTO_ADAPT set compute offsets
    if config["AUTO_ADAPT"].value == "YES":
        a0, a1 = photz.run_autoadapt(srclist)
        offsets = ",".join(np.array(a0).astype(str))
        offsets = "# Offsets from auto-adapt: " + offsets + "\n"
        print(offsets)
    else:
        a1, a0 = np.full(n_filters, 0), np.full(n_filters, 0)  # Do we need to set values?
        print("AUTO_ADAPT set to NO. Using zero offsets.")

    # create the onesource objects
    photozlist = []
    for i in range(ng):
        one_obj = lp.onesource(i, photz.gridz)
        one_obj.readsource(str(i), flux[i], flux_err[i], 63, zspec[i], " ")
        photz.prep_data(one_obj)
        photozlist.append(one_obj)

    # Perform the main run
    photz.run_photoz(photozlist, a0, a1)
    # Get the pdfs
    pdfs = []
    for i in range(ng):
        pdf = photozlist[i].pdfmap[11]
        pdf, zgrid = np.array(pdf.vPDF), np.array(pdf.xaxis)
        pdfs.append(pdf)

    # Loop over objects to compute photoz
    if filename is None:
        now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
        filename = f"lephare_output_{now}.fits"
    output = photz.build_output_tables(photozlist, para_out=None, filename=filename)
    return output, np.array(pdfs), np.array(zgrid)


def table_to_data(config, input, col_names=None, standard_names=False):
    """Take an astropy table and return the arrays required for a run.

    We assume that either the columns are in the standard LePHARE order or
    that the user provides the colnames in the standard order or they
    use a standard naming convention and set the standard_names=True

    The column names are id, f_filtername, ferr_filtername (for all filters)
    and context, zspec, string_input.

    Parameters
    ==========
    config : dict of lephare.keyword
        The config keymap. We need this to know if we have magnitudes or
        fluxes and to get the filter names.
    input : astropy.table.Table
        The input catalogue.
    col_names : list of str
        The column names in the default order.
    standard_names : bool
        If this is set we assume a standard column naming convention.

    Returns
    =======
    id : np.array
        The object ids
    flux : np.array
        The object fluxes or magnitudes
    flux_err : np.array
        The flux or magnitude errors
    context : np.array
        The context determining flux usage
    zspec : np.array
        The spectroscopic redshifts.
    string_data : np.array
        Additional notes as a string.
    """
    n_filters = len(config["FILTER_LIST"].value.split(","))
    if col_names is not None:
        print("Using user defined column names based on ordering.")
        assert len(input.colnames) == 2 * n_filters + 4

    elif standard_names:
        print(
            "Using standard column names id, f_filtername,ferr_filtername,... context, zspec, string_input."
        )
        col_names = ["id"]
        for f in config["FILTER_LIST"].value.split(","):
            col_names += [f"f_{f}"]
            col_names += [f"ferr_{f}"]
    else:
        print("Using user columns from input table assuming they are in the standard order.")
        assert len(input.colnames) == 2 * n_filters + 4
        col_names = list(np.full(2 * n_filters + 4, ""))
        col_names[0] = input.colnames[0]
        col_names[1 : 2 * n_filters : 2] = input.colnames[1 : 2 * n_filters : 2]
        col_names[2 : 2 * n_filters + 1 : 2] = input.colnames[2 : 2 * n_filters + 1 : 2]
        col_names[-3] = input.colnames[-3]
        col_names[-2] = input.colnames[-2]
        col_names[-1] = input.colnames[-1]
        # print(col_names)

    assert len(col_names) == 2 * n_filters + 4
    id = [str(i) for i in input[col_names[0]]]
    flux = input[col_names[1 : 2 * n_filters : 2]]
    flux = np.array([np.array(flux[c]) for c in flux.colnames])
    flux = flux.T
    flux_err = input[col_names[2 : 2 * n_filters + 1 : 2]]
    flux_err = np.array([np.array(flux_err[c]) for c in flux_err.colnames])
    flux_err = flux_err.T
    context = input[col_names[-3]]
    zspec = input[col_names[-2]]
    string_data = input[col_names[-1]]

    # Perform basic checks on table
    # Replace nans with -99

    mask = np.isnan(flux)
    mask |= np.isnan(flux_err)
    flux[mask] = -99.0
    flux_err[mask] = -99.0
    return id, flux, flux_err, context, zspec, string_data

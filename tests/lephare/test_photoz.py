import os

import lephare as lp
import numpy as np
import pytest
from astropy.table import Table


def test_photoz(test_data_dir: str):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    lp.prepare(config)
    photz = lp._lephare.PhotoZ(config)
    # Check the redshift zero case works
    for i, sed in enumerate(photz.fullLib):
        # Check some spectra data is available
        # assert np.array(sed.get_data_vector(0.0, 1000000.0, True, 0.0)).shape != (2, 0)
        # remove for now to use other tests
        # this should eventually test the python code for getting the spectra
        if not np.isclose(photz.zLib[i], 0.0):
            # Check redshifted sources have same model as previous sed
            assert sed.nummod == photz.fullLib[i - 1].nummod


def test_reddening(test_data_dir: str):
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    config["Z_STEP"] = "1.,0.,2."  # Fake star SED gets redshifted out of B band at low z
    # keymap=lp.all_types_to_keymap(config)
    config["APPLY_MW_EXTINCTION"] = "YES"
    config["EXT_MW_CURVE"] = "LMC_Fitzpatrick.dat"
    config["MW_REFERENCE_MODEL"] = "sed/STAR/PICKLES/o5v.sed.ext"
    mw_ebv_test_file = os.path.join(test_data_dir, "examples/mw_ebv.dat")
    config["MW_EBV_FILE"] = mw_ebv_test_file

    lp.prepare(config)

    # The reddening calculator
    albd_lib = lp.compute_model_reddening(config)
    assert albd_lib.shape == (19, 2)
    # test impact of ebv on a source fit
    # Read the test input catalogue
    input_file = os.path.join(test_data_dir, "examples/COSMOS_first100specz.fits")
    input = Table.read(input_file)
    ebv_test = np.linspace(0.0, 0.3, len(input))
    out = np.column_stack([input[input.colnames[0]], ebv_test])
    np.savetxt(mw_ebv_test_file, out, fmt=["%d", "%.8f"])

    test_string = "te s"  # Test with spaces
    input["string_input"][0] = test_string
    # Make a reduced column set for the minimal test
    reduced_cols = []
    for c in input.colnames:
        if not c.startswith("f"):
            reduced_cols.append(c)
        elif "IB527" in c:
            reduced_cols.append(c)
        elif "IB679" in c:
            reduced_cols.append(c)
    with pytest.warns(UserWarning, match="Reddening sent to process"):
        output, photozlist = lp.process(
            config, input[reduced_cols], write_outputs=False, reddening=albd_lib, ebvmw=[0.1] * len(input)
        )
    assert np.isclose(np.sum(albd_lib), 105.84983288457076)
    # Check it gives a warning if no ebv provided
    with pytest.warns(UserWarning, match="No ebv provided. Reddening not applied."):
        lp.process(config, input[reduced_cols], write_outputs=False, reddening=albd_lib)
    # Now test with SMC Prevot curve
    config["GAL_CURVE"] = "MW_seaton.dat"
    albd_lib = lp.compute_model_reddening(config)
    # Check AUTO_ADAPT can run with ebv
    config["AUTO_ADAPT"] = "YES"
    with pytest.warns(UserWarning, match="Reddening sent to process"):
        output, photozlist = lp.process(
            config, input[reduced_cols], write_outputs=False, reddening=albd_lib, ebvmw=[0.1] * len(input)
        )

    # Test the band pass correction
    bpc = lp.compute_band_pass_correction(config, model_number=1)
    print(bpc)
    assert np.isclose(np.sum(bpc), 17.034121587212294)

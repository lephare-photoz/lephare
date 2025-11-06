import inspect

import lephare as lp
import numpy as np
from astropy.table import Table


def test_all_plots():
    """Very simple test that runs all plots on dummy data to check they run without error."""
    n_rows = 100
    n_filters = 6  # u, g, r, z, J, K
    low, high = 0, 30
    n_z = 100

    # Create dictionary for table columns with fake data
    data = {
        "IDENT": np.arange(n_rows),  # ID as integers
        "Z_BEST": np.random.uniform(low, high, n_rows),
        "Z_BEST68_LOW": np.random.uniform(low, high, n_rows),
        "Z_BEST68_HIGH": np.random.uniform(low, high, n_rows),
        "Z_MED": np.random.uniform(low, high, n_rows),
        "Z_MED68_LOW": np.random.uniform(low, high, n_rows),
        "Z_MED68_HIGH": np.random.uniform(low, high, n_rows),
        "CHI_BEST": np.random.uniform(low, high, n_rows),
        "MOD_BEST": np.random.uniform(low, high, n_rows),
        "EXTLAW_BEST": np.random.uniform(low, high, n_rows),
        "EBV_BEST": np.random.uniform(low, high, n_rows),
        "Z_SEC": np.random.uniform(low, high, n_rows),
        "CHI_SEC": np.random.uniform(low, high, n_rows),
        "MOD_SEC": np.random.uniform(low, high, n_rows),
        "EBV_SEC": np.random.uniform(low, high, n_rows),
        "ZQ_BEST": np.random.uniform(low, high, n_rows),
        "CHI_QSO": np.random.uniform(low, high, n_rows),
        "MOD_QSO": np.random.uniform(low, high, n_rows),
        "MOD_STAR": np.random.uniform(low, high, n_rows),
        "CHI_STAR": np.random.uniform(low, high, n_rows),
        "MAG_OBS()": np.random.uniform(low, high, (n_rows, n_filters)),
        "MAG_ABS()": np.random.uniform(low, high, (n_rows, n_filters)),
        "SCALE_BEST": np.random.uniform(low, high, n_rows),
        "NBAND_USED": np.random.randint(1, n_filters + 1, n_rows),
        "CONTEXT": np.random.uniform(low, high, n_rows),
        "ZSPEC": np.random.uniform(low, high, n_rows),
        "AGE_BEST": np.random.uniform(low, high, n_rows),
        "AGE_INF": np.random.uniform(low, high, n_rows),
        "AGE_MED": np.random.uniform(low, high, n_rows),
        "AGE_SUP": np.random.uniform(low, high, n_rows),
        "LDUST_BEST": np.random.uniform(low, high, n_rows),
        "LDUST_INF": np.random.uniform(low, high, n_rows),
        "LDUST_MED": np.random.uniform(low, high, n_rows),
        "LDUST_SUP": np.random.uniform(low, high, n_rows),
        "LUM_TIR_BEST": np.random.uniform(low, high, n_rows),
        "LUM_TIR_INF": np.random.uniform(low, high, n_rows),
        "LUM_TIR_MED": np.random.uniform(low, high, n_rows),
        "LUM_TIR_SUP": np.random.uniform(low, high, n_rows),
        "MASS_BEST": np.random.uniform(low, high, n_rows),
        "MASS_INF": np.random.uniform(low, high, n_rows),
        "MASS_MED": np.random.uniform(low, high, n_rows),
        "MASS_SUP": np.random.uniform(low, high, n_rows),
        "SFR_BEST": np.random.uniform(low, high, n_rows),
        "SFR_INF": np.random.uniform(low, high, n_rows),
        "SFR_MED": np.random.uniform(low, high, n_rows),
        "SFR_SUP": np.random.uniform(low, high, n_rows),
        "SSFR_BEST": np.random.uniform(low, high, n_rows),
        "SSFR_INF": np.random.uniform(low, high, n_rows),
        "SSFR_MED": np.random.uniform(low, high, n_rows),
        "SSFR_SUP": np.random.uniform(low, high, n_rows),
        "LUM_NUV_BEST": np.random.uniform(low, high, n_rows),
        "LUM_R_BEST": np.random.uniform(low, high, n_rows),
        "LUM_K_BEST": np.random.uniform(low, high, n_rows),
        "PDF_BAY_ZG()": np.random.uniform(low, high, (n_rows, n_z)),
    }

    # Create Astropy Table
    t = Table(data)

    test_utils = lp.PlotUtils(
        t,
        sel_filt=3,
        pos_filt=[0, 1, 2, 4, 5, 5],
        range_z=[0, 0.5, 1, 1.5, 3],
        range_mag=[19, 20.5, 21.5, 22.5, 25],
    )

    # Loop over all methods of the instance
    for name, method in inspect.getmembers(test_utils, predicate=inspect.ismethod):
        # Skip private/internal methods (those starting with '_')
        if not name.startswith("_"):
            print(f"Running method: {name}")
            try:
                # Call the method
                method()
            except TypeError as e:
                # Handle methods that require arguments
                print(f"Skipping {name}, requires arguments: {e}")

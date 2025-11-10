import os

import lephare as lp
import numpy as np
import pytest
import scipy.integrate as sciint


def test_filter_extinc():
    """Test the filter_extinc runner makes a file and has good values"""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    out_file = os.path.join(test_dir, "../tmp", "filter_extinc.dat")
    options = {
        "FILTER_FILE": os.path.join(test_dir, "../data", "filt", "LSST_FILTERS.dat"),
        "EXT_CURVE": "SB_calzetti.dat",
        "GAL_CURVE": "CARDELLI",
        "OUTPUT": out_file,
    }
    runner = lp.FiltExt(config_keymap=lp.all_types_to_keymap(lp.default_cosmos_config), **options)
    all_filters, aint, albdav, albd = runner.run()

    # Check it made the file
    assert os.path.exists(out_file)
    with open(out_file, "r") as f:
        contents = f.read()
    assert float(contents.split()[-1]) == pytest.approx(1.3039147701114893)

    # check computation
    atmoext = lp.ext("atmo", 0)
    atmoext.read(os.path.join(os.environ["LEPHAREDIR"], "ext", "SB_calzetti.dat"))
    x1 = [o.lamb for o in atmoext.lamb_ext]
    y1 = [o.val for o in atmoext.lamb_ext]

    def ext(x):
        return np.interp(x, x1, y1, 0, 0)  # 0, 0 are the extrapolated default values

    aint2 = []
    # errs = []
    for f in all_filters:
        x2 = [o.lamb for o in f.lamb_trans]
        y2 = [o.val for o in f.lamb_trans]

        def filt(x):
            # 0, 0 are the extrapolated default values
            return np.interp(x, x2, y2, 0, 0)  # noqa: B023

        def fe(x):
            return ext(x) * filt(x)

        n = sciint.quad(fe, f.lmin(), f.lmax(), limit=200, epsabs=1.0e-3, epsrel=1.0e-4)
        d = sciint.quad(filt, f.lmin(), f.lmax(), limit=200, epsabs=1.0e-3, epsrel=1.0e-4)
        res = n[0] / d[0]
        aint2.append(res)
        # errs.append(sqrt(res**2 *(n[1]**2/n[0]**2 + d[1]**2/d[0]**2)) )

    np.testing.assert_array_almost_equal(aint, aint2, 5.0e-4)

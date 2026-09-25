import lephare as lp
import numpy as np


def test_prior():
    # Instantiate a source (Id, gridz)
    src = lp.onesource(1, np.arange(0, 6.1, 0.1))

    bp = [4, 5]
    prior_blue_low = []
    prior_blue_high = []
    prior_red_low = []
    prior_red_high = []

    # loop over the mag to test 19, 21, 23, 25,27, 29, 31
    test_ab = [
        9.1201084e-28,
        1.4454398e-28,
        2.2908677e-29,
        3.6307805e-30,
        5.7543994e-31,
        9.1201084e-32,
        1.4454398e-32,
    ]
    for k in range(0, 7):
        # define the source
        vals = [-99, -99, -99, -99, test_ab[k], -99]
        err_vals = [-99, -99, -99, -99, 0.01, -99]
        src.readsource("1", vals, err_vals, 0, -99, "test")
        src.fltUsed(-1, -1)
        # Compute the N(z) prior
        prior_blue_low.append(src.nzprior(10, 10, 0.1, bp))
        prior_blue_high.append(src.nzprior(10, 10, 4, bp))
        prior_red_low.append(src.nzprior(6, 10, 0.1, bp))
        prior_red_high.append(src.nzprior(6, 10, 4, bp))

    # 1. prior_blue_low decreasing in [0,1]
    assert all(
        prior_blue_low[i] >= prior_blue_low[i + 1] for i in range(len(prior_blue_low) - 1)
    ), "prior_blue_low not decreasing"
    assert all(0 <= val <= 1 for val in prior_blue_low), "prior_blue_low n'est pas dans [0,1]"
    # 2. prior_blue_high croissant et dans [0,1]
    assert all(
        prior_blue_high[i] <= prior_blue_high[i + 1] for i in range(len(prior_blue_high) - 1)
    ), "prior_blue_high not increasing"
    assert all(0 <= val <= 1 for val in prior_blue_high), "prior_blue_high n'est pas dans [0,1]"
    # 3. prior_red_low décroissant et dans [0,1]
    assert all(
        prior_red_low[i] >= prior_red_low[i + 1] for i in range(len(prior_red_low) - 1)
    ), "prior_red_low not decreasingt"
    assert all(0 <= val <= 1 for val in prior_red_low), "prior_red_low n'est pas dans [0,1]"
    # 4. prior_red_high lower than 0.0001
    assert all(val < 0.0001 for val in prior_red_high), "prior_red_high not lower then 0.0001"

    # case with the first filter solution not available. Should switch to next one
    prior_sbc = []
    prior_scd = []
    # loop over the mag to test 21, 23, 25
    test_ab = [1.4454398e-28, 2.2908677e-29, 3.6307805e-30]
    for k in range(0, 3):
        # define the source
        vals = [-99, -99, -99, -99, -99, test_ab[k]]
        err_vals = [-99, -99, -99, -99, -99, 0.01]
        src.readsource("1", vals, err_vals, 0, -99, "test")
        src.fltUsed(-1, -1)
        # Compute the N(z) prior
        # test the comparison between mod Sbc and Scd
        prior_sbc.append(src.nzprior(10, 11.6, 2, bp))
        prior_scd.append(src.nzprior(10, 10.8, 2, bp))

    # The prior should be higher for Sdc galaxies
    for sbc, scd in zip(prior_sbc, prior_scd):
        assert sbc < scd, f"prior_Sbc ({sbc}) not lower than prior_Scd ({scd})"

import lephare as lp
import numpy as np
from lephare import blackbody, check_first_char, indexes_in_vec


def test_globals_first_char():
    assert check_first_char("test")
    assert check_first_char("   test")
    assert not check_first_char("")
    assert not check_first_char(" ")
    assert not check_first_char("!test")
    assert not check_first_char("#test")
    assert not check_first_char("\t#")


def test_globals_blackbody():
    assert np.isclose(blackbody(10000, 500), 1.018807e-26)


def test_indexes_in_vec():
    v = [0.5, 0.9, 2.9, 3.5]
    val = 0.8
    assert indexes_in_vec(val, v, 0.1) == [1]
    assert indexes_in_vec(val, v, 0.01) == []
    assert indexes_in_vec(val, v, 0.3) == [0, 1]


def test_fast_interpolate():
    x = np.linspace(1, 5, 1000)
    y = x
    z = np.linspace(6, 7, 10)
    # full extrapolation
    assert np.allclose(lp.fast_interpolate(x, y, z, 0), np.zeros_like(z))
    assert np.allclose(lp.fast_interpolate(x, y, z, 1), np.ones_like(z))
    # partial extrapolation
    z = np.linspace(3.1456, 5.4387, 100)
    print(z)
    z2 = lp.fast_interpolate(x, y, z, 0)
    assert np.allclose(z2[:-19], z[:-19])
    assert np.allclose(z2[-19:], np.zeros_like(z2[-19:]))
    z2 = lp.fast_interpolate(x, y, z, 1)
    assert np.allclose(z2[-19:], np.ones_like(z2[-19:]))
    # no extrapolation
    z = np.linspace(2.1456, 4.4387, 100)
    print(z)
    z2 = lp.fast_interpolate(x, y, z, 0)
    assert np.allclose(z2, z)

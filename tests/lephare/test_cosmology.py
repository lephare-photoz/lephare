import numpy as np
import pytest
from lephare import cosmo, indexz, zgrid


def test_cosmology():
    # test default arguments
    c = cosmo()
    assert c == cosmo(70, 0.3, 0.7)
    # test setting arguments
    d = cosmo(h0=67, om0=0.27, l0=0.73)
    assert d != c


def test_cosmology_zgrid():
    # linear grid
    zmin = 0.000
    zmax = 1
    dz = 0.1
    # test linear grid
    grid = zgrid(dz, zmin, zmax)
    dummy = np.arange(zmin, zmax, dz)
    if zmin > 0:
        dummy = np.insert(dummy, 0, 0)
    if dummy[-1] != zmax:
        pygrid = np.insert(dummy, len(dummy), zmax)
    assert np.testing.assert_almost_equal(grid, pygrid) is None


def test_cosmology_indexz():
    grid = np.arange(1, 5, 1)
    # test value outside of grid on the lower side
    assert indexz(0, grid) == 0

    # test value outside of grid on the upper side
    assert indexz(6, grid) == grid.size - 1

    # test exact values in the grid
    assert indexz(1.0, grid) == 0
    assert indexz(2.0, grid) == 1
    assert indexz(3.0, grid) == 2
    assert indexz(4.0, grid) == 3

    # test value greater and close to a grid value
    assert indexz(2.1, grid) == 1

    # test value half to two grid values :
    # convention is to return the value above.
    assert indexz(2.5, grid) == 2
    # test value smaller and close to a grid value
    assert indexz(2.9, grid) == 2


def test_flux_rescaling():
    c = cosmo()
    z1 = 0.1
    dm1 = c.distMod(z1)
    z2 = 0.2
    dm2 = c.distMod(z2)
    assert c.flux_rescaling(z1, z1) == 1.0
    assert c.flux_rescaling(z1, z2) == pytest.approx(np.power(10, 0.4 * (dm2 - dm1)))

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
    with pytest.raises(ValueError):
        zgrid(dz, zmax, zmin)
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


# lephare's h0*1.0224e-12 s->yr conversion constant differs very slightly from
# astropy's (exact Julian year) definition, hence the small but non-zero rtol.
# The cases below cover every branch of cosmo::time(): flat LCDM was already
# indirectly exercised by test_cosmology, but Einstein-de Sitter, empty,
# open and closed matter-dominated universes were not tested at all before.
# This test previously caught a real bug: the open-universe (0<Om0<1, l0=0)
# branch used log10 instead of the natural logarithm required by the
# arccosh(x) = ln(x + sqrt(x^2-1)) identity, which silently produced ages
# wrong by up to a factor ~2 at z=2.
@pytest.mark.parametrize(
    "name,om0,l0",
    [
        ("einstein_de_sitter", 1.0, 0.0),
        ("empty_milne", 0.0, 0.0),
        ("open_matter_only", 0.3, 0.0),
        ("closed_matter_only", 1.5, 0.0),
        ("flat_lcdm", 0.3, 0.7),
    ],
)
@pytest.mark.parametrize("z", [0.0, 0.5, 2.0])
def test_time_against_astropy(name, om0, l0, z):
    astropy_cosmology = pytest.importorskip("astropy.cosmology")
    u = pytest.importorskip("astropy.units")
    h0 = 70.0
    c = cosmo(h0, om0, l0)
    ac = astropy_cosmology.LambdaCDM(H0=h0, Om0=om0, Ode0=l0)
    t_lephare = c.time(z)
    t_astropy = ac.age(z).to(u.yr).value
    assert t_lephare == pytest.approx(t_astropy, rel=1e-3)


def test_time_unsupported_cosmology_raises():
    # l0 !=0 and not flat: not covered by any closed-form branch
    c = cosmo(70, 0.5, 0.2)
    with pytest.raises(RuntimeError):
        c.time(0.5)


# distMet's l0==0 branches are intentionally excluded from coverage
# (LCOV_EXCL_START/STOP) as they are legacy/rarely used; the general
# (non-flat) LCDM branch used for any Om0<1, l0!=0 cosmology was however
# never numerically checked before.
@pytest.mark.parametrize(
    "om0,l0",
    [
        (0.3, 0.7),  # flat, for reference
        (0.3, 0.5),  # non-flat, open-ish LCDM: exercises the general branch
        (0.2, 0.9),  # non-flat, closed-ish LCDM
    ],
)
@pytest.mark.parametrize("z", [0.1, 1.0, 3.0])
def test_distmet_against_astropy(om0, l0, z):
    astropy_cosmology = pytest.importorskip("astropy.cosmology")
    u = pytest.importorskip("astropy.units")
    h0 = 70.0
    c = cosmo(h0, om0, l0)
    ac = astropy_cosmology.LambdaCDM(H0=h0, Om0=om0, Ode0=l0)
    # cosmo::distMet returns the "metric" distance (dlum = dmet*(1+z)),
    # i.e. the comoving distance for a non-flat cosmology.
    d_lephare = c.distMet(z)
    d_astropy = ac.comoving_distance(z).to(u.Mpc).value
    assert d_lephare == pytest.approx(d_astropy, rel=1e-3)

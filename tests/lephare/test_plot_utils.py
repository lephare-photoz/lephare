import inspect
import warnings

import lephare as lp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy.table import Table
from lephare._plot_utils import integrate_pdfs_to_ztrue

matplotlib.use("Agg")

N_FILTERS = 6  # u, g, r, z, J, K
SEL_FILT = 3
RANGE_Z = [0, 0.5, 1, 1.5, 3]
RANGE_MAG = [19, 20.5, 21.5, 22.5, 25]


def _random_table(n_rows=100, low=0, high=30, n_z=100):
    """Build a table of pure noise, spanning values well outside the plotting cuts.

    This is deliberately unphysical: it checks the plot methods survive garbage
    input, where nearly every source is rejected by the internal selections.
    """
    rng = np.random.default_rng(42)

    def col(*shape):
        return rng.uniform(low, high, shape)

    data = {
        "IDENT": np.arange(n_rows),
        "MAG_OBS()": col(n_rows, N_FILTERS),
        "MAG_ABS()": col(n_rows, N_FILTERS),
        "NBAND_USED": rng.integers(1, N_FILTERS + 1, n_rows),
        "PDF_BAY_ZG()": col(n_rows, n_z),
    }
    for key in (
        "Z_BEST Z_BEST68_LOW Z_BEST68_HIGH Z_MED Z_MED68_LOW Z_MED68_HIGH CHI_BEST MOD_BEST "
        "EXTLAW_BEST EBV_BEST Z_SEC CHI_SEC MOD_SEC EBV_SEC ZQ_BEST CHI_QSO MOD_QSO MOD_STAR "
        "CHI_STAR SCALE_BEST CONTEXT ZSPEC AGE_BEST AGE_INF AGE_MED AGE_SUP LDUST_BEST LDUST_INF "
        "LDUST_MED LDUST_SUP LUM_TIR_BEST LUM_TIR_INF LUM_TIR_MED LUM_TIR_SUP MASS_BEST MASS_INF "
        "MASS_MED MASS_SUP SFR_BEST SFR_INF SFR_MED SFR_SUP SSFR_BEST SSFR_INF SSFR_MED SSFR_SUP "
        "LUM_NUV_BEST LUM_R_BEST LUM_K_BEST"
    ).split():
        data[key] = col(n_rows)
    return Table(data)


def _realistic_table(n_rows=400, n_z=100, seed=8):
    """Build a physically plausible LePHARE output catalogue.

    Values are chosen so that a large fraction of sources actually pass the
    ``cond``/``condgal``/``condspec`` selections applied inside every plotting
    method, and so that all four redshift and magnitude panels are populated.
    Without this, the plot bodies are skipped and nothing is really exercised.
    """
    rng = np.random.default_rng(seed)

    # Spectroscopic redshifts spanning the full range_z, inside the outer edges.
    zs = rng.uniform(0.05, 2.9, n_rows)
    # Photo-z scattered about the spec-z, with a realistic (1+z) dependence and
    # a handful of catastrophic outliers.
    zp = np.clip(zs + rng.normal(0, 0.04, n_rows) * (1 + zs), 0.02, 2.95)
    zml = np.clip(zs + rng.normal(0, 0.05, n_rows) * (1 + zs), 0.02, 2.95)
    outlier = rng.random(n_rows) < 0.05
    zp[outlier] = rng.uniform(0.02, 2.95, outlier.sum())

    # Observed magnitudes: one intrinsic brightness per object plus a colour
    # term per filter. The selection filter spans the whole range_mag.
    mag_sel = rng.uniform(19.1, 24.9, n_rows)
    colours = np.array([0.9, 0.45, 0.0, -0.3, -0.7, -1.1])
    mag_obs = mag_sel[:, None] + colours[None, :] + rng.normal(0, 0.15, (n_rows, N_FILTERS))
    mag_obs[:, SEL_FILT] = mag_sel
    # Absolute magnitudes, roughly the observed mag minus a distance modulus.
    mag_abs = mag_obs - 5 * np.log10(3e3 * (1 + zp))[:, None] - 25

    # chi2: most objects are better fit by a galaxy than by a star (condgal),
    # but ~12% are star-like so condstar-gated panels are populated too.
    chi = rng.uniform(3, 40, n_rows)
    chi_star = chi * rng.uniform(1.1, 5.0, n_rows)
    star_like = rng.random(n_rows) < 0.12
    chi_star[star_like] = chi[star_like] * rng.uniform(0.2, 0.9, star_like.sum())

    log_mass = rng.uniform(8.5, 11.5, n_rows)
    log_sfr = rng.uniform(-1.0, 2.2, n_rows)

    def bracket(centre, width):
        """Return (inf, med, sup) triplets bracketing ``centre``."""
        return centre - width, centre, centre + width

    mass_l, mass_m, mass_s = bracket(log_mass, 0.2)
    sfr_l, sfr_m, sfr_s = bracket(log_sfr, 0.3)
    ssfr = log_sfr - log_mass
    ssfr_l, ssfr_m, ssfr_s = bracket(ssfr, 0.3)
    log_ldust = rng.uniform(9.0, 12.0, n_rows)
    ldust_l, ldust_m, ldust_s = bracket(log_ldust, 0.25)
    ltir_l, ltir_m, ltir_s = bracket(log_ldust + 0.1, 0.25)
    age = rng.uniform(1e8, 1e10, n_rows)
    age_l, age_m, age_s = bracket(age, 2e7)

    # Second photo-z peak: present for most objects, flagged absent for the rest.
    zp2 = rng.uniform(0.05, 2.9, n_rows)
    zp2[rng.random(n_rows) < 0.15] = -99.0

    # PDFs: a normalised Gaussian per object, centred on its photo-z.
    zgrid = np.linspace(0, 6, n_z)
    sigma = 0.05 * (1 + zp)
    pdfs = np.exp(-0.5 * ((zgrid[None, :] - zp[:, None]) / sigma[:, None]) ** 2)
    pdfs /= np.trapezoid(pdfs, zgrid, axis=1)[:, None]

    data = {
        "IDENT": np.arange(n_rows),
        "Z_BEST": zp,
        "Z_BEST68_LOW": zp - 0.05 * (1 + zp),
        "Z_BEST68_HIGH": zp + 0.05 * (1 + zp),
        "Z_MED": zml,
        "Z_MED68_LOW": zml - 0.06 * (1 + zml),
        "Z_MED68_HIGH": zml + 0.06 * (1 + zml),
        "CHI_BEST": chi,
        "MOD_BEST": rng.integers(1, 31, n_rows).astype(float),
        "EXTLAW_BEST": rng.integers(0, 3, n_rows).astype(float),
        "EBV_BEST": rng.uniform(0, 0.5, n_rows),
        "Z_SEC": zp2,
        "CHI_SEC": chi * rng.uniform(1.0, 2.0, n_rows),
        "MOD_SEC": rng.integers(1, 31, n_rows).astype(float),
        "EBV_SEC": rng.uniform(0, 0.5, n_rows),
        "ZQ_BEST": rng.uniform(0.05, 2.9, n_rows),
        "CHI_QSO": chi * rng.uniform(1.0, 4.0, n_rows),
        "MOD_QSO": rng.integers(1, 31, n_rows).astype(float),
        "MOD_STAR": rng.integers(1, 31, n_rows).astype(float),
        "CHI_STAR": chi_star,
        "MAG_OBS()": mag_obs,
        "MAG_ABS()": mag_abs,
        "SCALE_BEST": rng.uniform(1e-5, 1e-3, n_rows),
        "NBAND_USED": rng.integers(3, N_FILTERS + 1, n_rows),
        "CONTEXT": np.full(n_rows, 2**N_FILTERS - 1, dtype=float),
        "ZSPEC": zs,
        "AGE_BEST": age,
        "AGE_INF": age_l,
        "AGE_MED": age_m,
        "AGE_SUP": age_s,
        "LDUST_BEST": log_ldust,
        "LDUST_INF": ldust_l,
        "LDUST_MED": ldust_m,
        "LDUST_SUP": ldust_s,
        "LUM_TIR_BEST": log_ldust + 0.1,
        "LUM_TIR_INF": ltir_l,
        "LUM_TIR_MED": ltir_m,
        "LUM_TIR_SUP": ltir_s,
        "MASS_BEST": log_mass,
        "MASS_INF": mass_l,
        "MASS_MED": mass_m,
        "MASS_SUP": mass_s,
        "SFR_BEST": log_sfr,
        "SFR_INF": sfr_l,
        "SFR_MED": sfr_m,
        "SFR_SUP": sfr_s,
        "SSFR_BEST": ssfr,
        "SSFR_INF": ssfr_l,
        "SSFR_MED": ssfr_m,
        "SSFR_SUP": ssfr_s,
        "LUM_NUV_BEST": rng.uniform(9.0, 11.0, n_rows),
        "LUM_R_BEST": rng.uniform(9.0, 11.0, n_rows),
        "LUM_K_BEST": rng.uniform(9.5, 11.5, n_rows),
        "PDF_BAY_ZG()": pdfs,
    }
    return Table(data)


@pytest.fixture
def realistic_table():
    return _realistic_table()


@pytest.fixture
def plot_utils(realistic_table):
    return lp.PlotUtils(
        realistic_table,
        sel_filt=SEL_FILT,
        pos_filt=[0, 1, 2, 4, 5, 5],
        range_z=RANGE_Z,
        range_mag=RANGE_MAG,
    )


def _run_all_public_methods(test_utils):
    """Call every public method with no arguments, returning the names that ran."""
    ran = []
    for name, method in inspect.getmembers(test_utils, predicate=inspect.ismethod):
        # Skip private/internal methods (those starting with '_')
        if name.startswith("_"):
            continue
        print(f"Running method: {name}")
        try:
            # Ignore the Matplotlib legend warning during this method call
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="No artists with labels found to put in legend")
                warnings.filterwarnings("ignore", message="invalid value encountered in divide")
                method()
                ran.append(name)
        except TypeError as e:
            # Handle methods that require arguments
            print(f"Skipping {name}, requires arguments: {e}")
        finally:
            plt.close("all")
    return ran


def test_all_plots(tmp_path, monkeypatch):
    """Very simple test that runs all plots on dummy data to check they run without error."""
    # save_*_plots_pdf write into the working directory, so keep them out of the repo
    monkeypatch.chdir(tmp_path)
    test_utils = lp.PlotUtils(
        _random_table(),
        sel_filt=SEL_FILT,
        pos_filt=[0, 1, 2, 4, 5, 5],
        range_z=RANGE_Z,
        range_mag=RANGE_MAG,
    )
    assert _run_all_public_methods(test_utils)


def test_all_plots_realistic(plot_utils, tmp_path, monkeypatch):
    """Run all plots on a plausible catalogue, so the plot bodies really execute.

    With random data almost every source falls outside the internal z/mag cuts,
    so the ``if len(...) > 0`` guards short-circuit and the drawing code is never
    reached. Here the selections pass, so the histograms and statistics run.
    """
    # save_*_plots_pdf write into the working directory
    monkeypatch.chdir(tmp_path)
    ran = _run_all_public_methods(plot_utils)

    # Sanity check that the fixture really does select a healthy sample; if this
    # regresses, the plots above stop testing anything.
    assert plot_utils.cond.sum() > 100
    assert (plot_utils.cond & plot_utils.condgal & plot_utils.condspec).sum() > 100
    assert (plot_utils.cond & plot_utils.condstar).sum() > 5

    # Every public method should have been runnable without arguments
    assert "zml_zs" in ran
    assert "pit_qq" in ran
    assert (tmp_path / "all_plots.pdf").exists()
    assert (tmp_path / "all_phys_plots.pdf").exists()


def test_plot_utils_defaults(realistic_table):
    """With no binning given, the panel edges are derived from the data quantiles."""
    test_utils = lp.PlotUtils(realistic_table)

    # range_z/range_mag fall back to quartiles of the selected sample
    assert len(test_utils.range_z) == 5
    assert len(test_utils.range_mag) == 5
    assert test_utils.range_z[0] < test_utils.range_z[-1]
    # sel_filt defaults to 0, so the magnitude bins track the first filter
    mag0 = np.asarray(realistic_table["MAG_OBS()"][:, 0])
    expected_mag = np.quantile(mag0[(mag0 > 10) & (mag0 < 40)], [0, 0.25, 0.5, 0.75, 1])
    np.testing.assert_allclose(test_utils.range_mag, expected_mag)
    # pos_filt defaults to all zeros, so every colour uses the first filter
    assert (
        test_utils.uFilt
        == test_utils.bFilt
        == test_utils.rFilt
        == test_utils.zFilt
        == test_utils.jFilt
        == test_utils.kFilt
        == 0
    )
    # Four bins in each of z and mag means a 2x2 panel grid
    assert (test_utils.nbRowZ, test_utils.nbColZ) == (2, 2)
    assert (test_utils.nbRowM, test_utils.nbColM) == (2, 2)


def test_plot_utils_single_panel(realistic_table):
    """A single z/mag bin collapses the panel grid to one row and column."""
    test_utils = lp.PlotUtils(realistic_table, range_z=[0, 3], range_mag=[19, 25])
    assert (test_utils.nbRowZ, test_utils.nbColZ) == (1, 1)
    assert (test_utils.nbRowM, test_utils.nbColM) == (1, 1)
    assert (test_utils.z_min, test_utils.z_max) == (0, 3)
    assert (test_utils.mag_min, test_utils.mag_max) == (19, 25)


def test_plot_utils_out_of_bounds_filters(realistic_table, capsys):
    """Filter indices beyond the number of columns are reset to 0 with a warning."""
    test_utils = lp.PlotUtils(
        realistic_table,
        sel_filt=N_FILTERS + 3,
        pos_filt=[0, N_FILTERS, -2, 1, 2, 3],
        range_z=RANGE_Z,
        range_mag=RANGE_MAG,
    )
    out = capsys.readouterr().out
    assert "sel_filt out of bounds" in out
    assert "pos_filt[1] out of bounds" in out
    assert "pos_filt[2] out of bounds" in out
    # The offending entries were clamped, the valid ones left alone
    assert test_utils.sel_filt == 0
    assert test_utils.bFilt == 0
    assert test_utils.rFilt == 0
    assert test_utils.zFilt == 1
    # Clamping sel_filt means the magnitude column is the first filter's
    np.testing.assert_allclose(test_utils.mag, realistic_table["MAG_OBS()"][:, 0])


def test_title_page_without_network(plot_utils, capsys, monkeypatch):
    """The title page still renders when the logo cannot be downloaded."""

    def no_network(*args, **kwargs):
        raise OSError("network unreachable")

    monkeypatch.setattr(lp._plot_utils.urllib.request, "urlopen", no_network)
    plot_utils.title_page()

    assert "Could not load LePHARE logo: network unreachable" in capsys.readouterr().out
    # The page itself was still produced
    assert plt.get_fignums()
    plt.close("all")


def test_pit_qq_variants(plot_utils, tmp_path, monkeypatch):
    """The PIT/QQ panel can be drawn as QQ only, PIT only, or both."""
    monkeypatch.chdir(tmp_path)

    # Both panels, and a title
    assert plot_utils.pit_qq(title="both") is None

    # QQ only
    plt.close("all")
    assert plot_utils.pit_qq(show_pit=False) is None

    # PIT only, which uses the single-axis histogram branch
    plt.close("all")
    assert plot_utils.pit_qq(show_qq=False) is None

    # Neither, so only the shared axis setup runs
    plt.close("all")
    assert plot_utils.pit_qq(show_pit=False, show_qq=False) is None

    # bins given as explicit edges rather than a count
    plt.close("all")
    assert plot_utils.pit_qq(bins=np.linspace(0, 1, 21)) is None

    # savefig returns the path it wrote
    plt.close("all")
    fig_filename = plot_utils.pit_qq(savefig=True)
    assert fig_filename == "plot_pit_qq_lephare.png"
    assert (tmp_path / fig_filename).exists()
    plt.close("all")


def test_pit_qq_explicit_arrays(plot_utils):
    """PDFs, grid and true redshifts can all be supplied explicitly."""
    zgrid = np.linspace(0, 3, 60)
    ztrue = np.array([0.5, 1.0, 1.5])
    pdfs = np.exp(-0.5 * ((zgrid[None, :] - ztrue[:, None]) / 0.1) ** 2)
    pdfs /= np.trapezoid(pdfs, zgrid, axis=1)[:, None]

    assert plot_utils.pit_qq(pdfs=pdfs, zgrid=zgrid, ztrue=ztrue, bins=10) is None
    plt.close("all")


def test_integrate_pdfs_to_ztrue_well_calibrated():
    """Integrating a PDF up to its own centre gives ~0.5 of the total mass."""
    zgrid = np.linspace(0, 6, 601)
    ztrue = np.array([1.0, 2.0, 3.0])
    pdfs = np.exp(-0.5 * ((zgrid[None, :] - ztrue[:, None]) / 0.1) ** 2)
    pdfs /= np.trapezoid(pdfs, zgrid, axis=1)[:, None]

    pit = integrate_pdfs_to_ztrue(pdfs, zgrid, ztrue)
    assert pit.shape == (3,)
    np.testing.assert_allclose(pit, 0.5, atol=1e-3)


def test_integrate_pdfs_to_ztrue_edges():
    """A ztrue below the grid integrates to 0, above the grid to the full mass."""
    zgrid = np.linspace(1.0, 5.0, 401)
    # A flat, normalised PDF makes the expected integral easy to state exactly.
    pdfs = np.tile(1.0 / (zgrid[-1] - zgrid[0]), (3, len(zgrid)))
    ztrue = np.array([0.5, 3.0, 9.0])

    pit = integrate_pdfs_to_ztrue(pdfs, zgrid, ztrue)
    # Below the first grid point nothing is integrated
    assert pit[0] == 0.0
    # Halfway through a flat PDF is half the mass
    assert pit[1] == pytest.approx(0.5, abs=1e-3)
    # Beyond the last grid point the whole PDF is integrated
    assert pit[2] == pytest.approx(1.0, abs=1e-6)


def test_save_pdf_skips_figureless_methods(plot_utils, tmp_path, monkeypatch, capsys):
    """Plot methods that draw nothing, or draw an empty figure, are skipped.

    The PDF writer collects whatever figures a method left open. A method that
    opens none, or opens one with no axes, must not end up as a blank page.
    """
    monkeypatch.chdir(tmp_path)

    def draws_nothing(**kwargs):
        return None

    def draws_empty_figure(**kwargs):
        plt.figure()  # a figure with no axes at all

    # Replace two of the methods the PDF writer calls
    plot_utils.title_page = draws_nothing
    plot_utils.zml_zs = draws_empty_figure

    plot_utils.save_photoz_plots_pdf(filename="skipped.pdf")

    assert "No figure created by draws_nothing(), skipping." in capsys.readouterr().out
    # The remaining methods still produced a usable PDF
    assert (tmp_path / "skipped.pdf").exists()
    assert (tmp_path / "skipped.pdf").stat().st_size > 0
    plt.close("all")


def test_save_phys_pdf_skips_figureless_methods(plot_utils, tmp_path, monkeypatch, capsys):
    """The physical-parameter PDF writer applies the same skipping rules."""
    monkeypatch.chdir(tmp_path)

    def draws_nothing(**kwargs):
        return None

    def draws_empty_figure(**kwargs):
        plt.figure()

    plot_utils.title_page = draws_nothing
    plot_utils.dist_mass = draws_empty_figure

    plot_utils.save_phys_plots_pdf(filename="skipped_phys.pdf")

    assert "No figure created by draws_nothing(), skipping." in capsys.readouterr().out
    assert (tmp_path / "skipped_phys.pdf").exists()
    assert (tmp_path / "skipped_phys.pdf").stat().st_size > 0
    plt.close("all")


def test_save_pdf_retries_methods_without_kwargs(plot_utils, tmp_path, monkeypatch):
    """Plot methods that take no keyword arguments are retried bare.

    save_*_plots_pdf forwards its **kwargs to every method, but several of them
    (title_page, bzk, ...) accept none, so the call must fall back to method().
    """
    monkeypatch.chdir(tmp_path)
    calls = []

    def takes_no_kwargs():
        calls.append("called")
        plt.figure().add_subplot(111).plot([0, 1], [0, 1])

    plot_utils.title_page = takes_no_kwargs
    # nstep is accepted by the histogram methods but not by title_page
    plot_utils.save_photoz_plots_pdf(filename="retried.pdf", nstep=5)

    # The bare retry happened exactly once for the kwarg-less method
    assert calls == ["called"]
    assert (tmp_path / "retried.pdf").exists()
    plt.close("all")


def test_save_phys_pdf_retries_methods_without_kwargs(plot_utils, tmp_path, monkeypatch):
    """The physical-parameter writer uses the same bare-call fallback."""
    monkeypatch.chdir(tmp_path)
    calls = []

    def takes_no_kwargs():
        calls.append("called")
        plt.figure().add_subplot(111).plot([0, 1], [0, 1])

    plot_utils.title_page = takes_no_kwargs
    plot_utils.save_phys_plots_pdf(filename="retried_phys.pdf", nstep=5)

    assert calls == ["called"]
    assert (tmp_path / "retried_phys.pdf").exists()
    plt.close("all")

import matplotlib

matplotlib.use("Agg")

import os

import lephare as lp
import matplotlib.pyplot as plt

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_spec_plotspec():
    lp._spec.plot_spectrum(os.path.join(TESTDATADIR, "example.spec"))
    assert True


def _write_spec(path, filter_rows, npdf=11, gal_nline=5):
    """Write a minimal .spec file with the given per-filter rows.

    Each row is (mag, err_mag, lambda_eff, width, mag_model). The remaining
    columns match the layout LePHARE writes but are not used by the plot.
    """
    nfilt = len(filter_rows)
    lines = [
        "# Ident Zspec Zphot",
        "1.0 1.500 1.4500",
        "# Mag emag  Lbd_mean  Lbd_width Mag_gal  Mag_FIR  Mag_BCSTOCH",
        f"FILTERS  {nfilt}",
        "# Zstep  PDF",
        f"PDF  {npdf} {npdf}",
        "# Type Nline Model Library Nband  Zphot Zinf Zsup Chi2  PDF  Extlaw EB-V Lir Age Mass SFR SSFR",
        # Only GAL-1 is a used model; the rest carry model = -1 and are skipped
        f"GAL-1 {gal_nline} 30 1 {nfilt} 1.45 1.45 1.45 91.05 -1 0 0.1 -999 -999 -999 -999 -999",
        "GAL-2 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1.",
        "GAL-FIR 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1.",
        "GAL-STOCH 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1.",
        "QSO 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1.",
        "STAR 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1.",
    ]
    # Observed photometry: mag err lbd width mag_gal mag_fir mag_phys pad mag_mod
    for mag, err, lbd, width, mag_mod in filter_rows:
        lines.append(f"{mag} {err} {lbd} {width} {mag_mod} -1 -1 1 {mag_mod}")
    # PDF(z): a simple triangular peak so max() is non-zero
    for i in range(npdf):
        z = 3.0 * i / (npdf - 1)
        prob = 1.0 - abs(i - npdf // 2) / float(npdf)
        lines.append(f"{z} {prob} {prob}")
    # Best-fit template for GAL-1 only
    for i in range(gal_nline):
        lines.append(f"{3000.0 + 1000.0 * i} {24.0 + 0.1 * i}")

    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_plotspec_axis_limits_from_loose_errors(tmp_path):
    """With no error below 0.5 mag, the axis limits fall back to the <10 mag cut."""
    # em is doubled internally, so err = 1.0 gives em1 = 2.0: above the first
    # threshold of 1 but below the fallback threshold of 10.
    rows = [(23.0 + i, 1.0, 4000.0 + 1000.0 * i, 400.0, 23.2 + i) for i in range(4)]
    lp._spec.plotspec(_write_spec(tmp_path / "loose.spec", rows))
    plt.close("all")


def test_plotspec_axis_limits_default(tmp_path):
    """With every error above the fallback cut too, hardcoded axis limits are used."""
    rows = [(23.0 + i, 6.0, 4000.0 + 1000.0 * i, 400.0, 23.2 + i) for i in range(4)]
    lp._spec.plotspec(_write_spec(tmp_path / "wide.spec", rows))
    plt.close("all")


def test_plotspec_no_usable_filters(tmp_path):
    """When no band passes the reliability cuts, a default wavelength range is used."""
    # mag = 1000 marks a saturated/unusable band, as LePHARE writes for masked data
    rows = [(1000.0, 1000.0, 4000.0 + 1000.0 * i, 400.0, 1000.0) for i in range(3)]
    lp._spec.plotspec(_write_spec(tmp_path / "unusable.spec", rows))
    plt.close("all")


def test_plotspec_draws_upper_limits(tmp_path):
    """A negative error flags an upper limit, drawn as a downward arrow."""
    rows = [
        (23.0, 0.1, 4000.0, 400.0, 23.1),
        (24.0, 0.2, 5000.0, 400.0, 24.1),
        # Negative error mag: an upper limit rather than a detection
        (25.0, -1.0, 6000.0, 400.0, 25.1),
    ]
    lp._spec.plotspec(_write_spec(tmp_path / "uplim.spec", rows))

    # The upper limit is drawn with a quiver arrow on the main panel
    main_panel = plt.gcf().axes[0]
    assert any(type(c).__name__ == "Quiver" for c in main_panel.collections)
    plt.close("all")


def test_plotspec_reports_spec_z(tmp_path):
    """A positive spec-z adds a labelled vertical line to the PDF inset."""
    rows = [(23.0 + i, 0.1, 4000.0 + 1000.0 * i, 400.0, 23.1 + i) for i in range(4)]
    lp._spec.plotspec(_write_spec(tmp_path / "specz.spec", rows))

    # The inset is the second axis; it carries the PDF legend
    inset = plt.gcf().axes[1]
    assert inset.get_title() == "PDF(z)"
    labels = [t.get_text() for t in inset.get_legend().get_texts()]
    assert "spec z" in labels
    assert "photo z" in labels
    plt.close("all")

"""Tests of the redshift quality flag (Z_FLAG), computed in C++."""

import os

import lephare as lp
import numpy as np
import pytest
import scipy.signal
from astropy.table import Table


def gauss(x, mu, sig, amp=1.0):
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)


@pytest.fixture
def zgrid():
    return np.linspace(0.0, 6.0, 601)


def test_flagparams_defaults_and_parsing():
    p = lp.FlagParams()
    assert (p.nb_peak, p.height, p.tail, p.peak_ratio, p.error) == (2, 0.43, 0.23, 0.1, 0.1)
    # Empty and NONE give the defaults
    assert lp.FlagParams.from_string("").to_string() == p.to_string()
    assert lp.FlagParams.from_string("NONE").to_string() == p.to_string()
    # Case, spaces and order do not matter, unspecified values keep the default
    q = lp.FlagParams.from_string(" NB_Peak = 3 , error=0.2,, tail=0.5 ")
    assert q.nb_peak == 3
    assert q.error == pytest.approx(0.2)
    assert q.tail == pytest.approx(0.5)
    assert q.height == pytest.approx(0.43)
    # The canonical string can be read back
    r = lp.FlagParams.from_string(q.to_string())
    assert r.to_string() == q.to_string()
    assert "nb_peak=3" in repr(q)


@pytest.mark.parametrize(
    "bad",
    [
        "unknown=1",
        "nb_peak",
        "nb_peak=abc",
        "nb_peak=1.5",
        "nb_peak=0",
        "height=0",
        "height=1.5",
        "tail=2",
        "n_window=0",
        "fwhm_level=1",
        "error=-1",
        "error=nan",
        "error=1,error=2",
    ],
)
def test_flagparams_invalid(bad):
    with pytest.raises(ValueError, match="Z_FLAG_PAR"):
        lp.FlagParams.from_string(bad)


def test_flag_components(zgrid):
    x = zgrid
    # Narrow single peak: perfect
    res = lp.zflag(x, gauss(x, 2.0, 0.03))
    assert res.flag == 0
    assert res.n_peaks == 1
    assert res.error == pytest.approx(0.03, rel=0.05)
    assert res.estimate == pytest.approx(2.0)
    # sigma is the Gaussian equivalent width derived from the FWHM
    assert res.sigma == pytest.approx(0.03, rel=0.05)
    # Larger error only (bit 1)
    assert lp.zflag(x, gauss(x, 2.0, 0.3)).flag & 1
    # Two comparable peaks (bit 8), well separated
    bi = gauss(x, 1.0, 0.03) + gauss(x, 4.0, 0.03)
    res = lp.zflag(x, bi)
    assert res.n_peaks == 2
    assert res.flag & 8
    # A high secondary peak close to the first is merged (distance filter)
    close = gauss(x, 2.0, 0.03) + 0.9 * gauss(x, 2.03, 0.03)
    assert lp.zflag(x, close).n_peaks == 1
    # nb_peak controls bit 8
    assert not (lp.zflag(x, bi, lp.FlagParams.from_string("nb_peak=3")).flag & 8)
    # Heavy tail, bit 4: narrow core plus a broad pedestal holding most of the mass
    tail = gauss(x, 2.0, 0.03, 1.0) + 0.05 * gauss(x, 3.0, 1.0)
    res = lp.zflag(x, tail)
    assert res.tail_mass > 0.23
    assert res.flag & 4
    # Flat PDF: error larger than half the range (bit 16), high mean/max (bit 2)
    res = lp.zflag(x, np.ones_like(x))
    assert res.flag & 16
    assert res.flag & 2
    assert not (res.flag & 1)
    # Thresholds can be changed
    assert lp.zflag(x, gauss(x, 2.0, 0.3), lp.FlagParams.from_string("error=0.5")).flag & 1 == 0


def test_flag_invalid_pdfs(zgrid):
    x = zgrid
    assert lp.zflag(x, np.zeros_like(x)).flag == lp.ZFLAG_INVALID
    bad = gauss(x, 2.0, 0.1)
    bad[10] = np.nan
    assert lp.zflag(x, bad).flag == lp.ZFLAG_INVALID
    bad[10] = np.inf
    assert lp.zflag(x, bad).flag == lp.ZFLAG_INVALID
    bad[10] = -1.0
    assert lp.zflag(x, bad).flag == lp.ZFLAG_INVALID
    assert lp.zflag(x[:2], np.ones(2)).flag == lp.ZFLAG_INVALID
    with pytest.raises(ValueError):
        lp.zflag(x, np.ones(10))


def test_flag_normalisation_and_estimate(zgrid):
    x = zgrid
    p = gauss(x, 2.0, 0.05)
    # The flag does not depend on the normalisation of the PDF
    assert lp.zflag(x, 1.0e-12 * p).flag == lp.zflag(x, 1.0e3 * p).flag
    # The statistics are made around the estimate given
    centred = lp.zflag(x, p, lp.FlagParams(), 2.0)
    off = lp.zflag(x, p, lp.FlagParams(), 3.0)
    assert centred.estimate == 2.0
    assert off.error > 10 * centred.error
    assert off.flag != centred.flag


def test_batch_matches_single(zgrid):
    x = zgrid
    rng = np.random.default_rng(1)
    rows = []
    for _ in range(200):
        n = rng.integers(1, 4)
        p = np.zeros_like(x)
        for _ in range(n):
            p += gauss(x, rng.uniform(0.2, 5.8), 10 ** rng.uniform(-2, 0), rng.uniform(0.1, 1))
        rows.append(p)
    rows[3][:] = 0.0  # invalid
    v = np.array(rows)
    par = lp.FlagParams.from_string("nb_peak=2,height=0.3")
    batch = lp.zflag_batch(x, v, par)
    assert batch.shape == (200,)
    assert batch.dtype == np.int32
    single = np.array([lp.zflag(x, r, par).flag for r in v])
    assert np.array_equal(batch, single)
    assert batch[3] == lp.ZFLAG_INVALID
    assert len(np.unique(batch)) > 3
    with pytest.raises(ValueError):
        lp.zflag_batch(x, v[:, :10])


def test_peak_counting_matches_scipy(zgrid):
    """Peak finding follows scipy.signal.find_peaks(height, distance)."""
    x = zgrid
    rng = np.random.default_rng(2)
    par = lp.FlagParams.from_string("height=0.3,peak_dz=0.05")
    dist = int(np.ceil(0.05 / (x[1] - x[0]) - 1.0e-6))
    for _ in range(300):
        p = np.zeros_like(x)
        for _ in range(rng.integers(1, 7)):
            p += gauss(x, rng.uniform(0.0, 6.0), 10 ** rng.uniform(-2, -0.5), rng.uniform(0.05, 1))
        if rng.random() < 0.3:  # plateaus
            p = np.round(p, 1)
        if p.max() == 0:
            continue
        q = p / np.trapezoid(p, x)
        expected = len(scipy.signal.find_peaks(q, height=par.height * q.max(), distance=dist)[0])
        assert lp.zflag(x, p, par).n_peaks == expected


def test_pdf_quality_flag_method():
    pdf = lp.PDF(0.0, 6.0, 601)
    pdf.setYvals(gauss(np.array(pdf.xaxis), 1.5, 0.04), is_chi2=False)
    res = pdf.quality_flag()
    assert res.flag == 0
    assert res.estimate == pytest.approx(1.5)


def test_onesource_without_pdf_is_invalid():
    src = lp.onesource()
    assert src.zflag == lp.ZFLAG_NOT_COMPUTED
    src.compute_zflag(lp.FlagParams())
    assert src.zflag == lp.ZFLAG_INVALID


def test_zflag_end_to_end(test_data_dir, tmp_path, monkeypatch):
    """Z_FLAG in PARA_OUT: C++ output column, header and table agree."""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    monkeypatch.setenv("LEPHAREDIR", os.path.join(test_dir, "../data"))
    monkeypatch.setenv("LEPHAREWORK", os.path.join(test_dir, "../tmp"))
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    lp.prepare(config)

    # Output parameters: the standard ones plus Z_FLAG
    para = tmp_path / "output.para"
    with open(os.path.join(test_data_dir, "examples/output.para")) as f:
        para.write_text(f.read() + "\nZ_FLAG\n")
    config["PARA_OUT"] = str(para)
    config["CAT_OUT"] = str(tmp_path / "zphot.out")
    config["PDZ_OUT"] = "NONE"
    config["Z_FLAG_PAR"] = "nb_peak=3,error=0.2"

    input = Table.read(os.path.join(test_data_dir, "examples/COSMOS_first100specz.fits"))
    cols = [c for c in input.colnames if not c.startswith("f") or "IB527" in c or "IB679" in c]
    output, photozlist = lp.process(config, input[cols][:20], write_outputs=True)

    flags = np.array(output["Z_FLAG"])
    assert len(flags) == 20
    par = lp.FlagParams.from_string("nb_peak=3,error=0.2")
    for src, f in zip(photozlist, flags):
        assert src.zflag == f
        assert f == src.pdfmap[11].quality_flag(par).flag
    assert output.meta["Z_FLAG_PAR"] == par.to_string()
    assert output.meta["Z_FLAG_DEF"] == lp.ZFLAG_DEFINITION

    # The ascii catalogue has the Z_FLAG column, indexed in the header
    with open(tmp_path / "zphot.out") as f:
        header = [line for line in f if line.startswith("#")]
    text = "".join(header)
    assert "Z_FLAG_PAR" in text
    assert "nb_peak=3" in text
    fmt = header[[i for i, line in enumerate(header) if line.startswith("# Format:")][0] + 1].split()
    col = int(fmt[fmt.index("Z_FLAG") + 1]) - 1
    data = np.loadtxt(tmp_path / "zphot.out", usecols=(col,), dtype=int)
    assert np.array_equal(data, flags)


def test_zflag_requested_only_in_table(test_data_dir, tmp_path, monkeypatch):
    """A table built with a para_out containing Z_FLAG computes it on demand."""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    monkeypatch.setenv("LEPHAREDIR", os.path.join(test_dir, "../data"))
    monkeypatch.setenv("LEPHAREWORK", os.path.join(test_dir, "../tmp"))
    config = lp.read_config(os.path.join(test_data_dir, "examples/COSMOS.para"))
    lp.prepare(config)
    input = Table.read(os.path.join(test_data_dir, "examples/COSMOS_first100specz.fits"))
    cols = [c for c in input.colnames if not c.startswith("f") or "IB527" in c or "IB679" in c]
    output, photozlist = lp.process(config, input[cols][:5], write_outputs=False)
    assert "Z_FLAG" not in output.colnames
    assert all(s.zflag == lp.ZFLAG_NOT_COMPUTED for s in photozlist)

    para = tmp_path / "flag.para"
    para.write_text("IDENT\nZ_BEST\nZ_FLAG\n")
    photz = lp.PhotoZ(config)
    table = photz.build_output_tables(photozlist, para_out=str(para))
    assert list(table.colnames)[:3] == ["IDENT", "Z_BEST", "Z_FLAG"]
    assert all(s.zflag != lp.ZFLAG_NOT_COMPUTED for s in photozlist)

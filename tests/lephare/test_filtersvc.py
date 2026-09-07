import os

import pytest
import yaml
from lephare import LEPHAREDIR
from lephare.filterSvc import FilterSvc

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


@pytest.fixture
def filter_names():
    return ["cosmos/IB527.lowres", "cosmos/IB679.lowres"]


def _write_yaml(tmp_path, config):
    path = tmp_path / "filters.yaml"
    path.write_text(yaml.safe_dump({"filters": config}))
    return str(path)


def test_from_yaml(tmp_path, filter_names):
    """Filters listed in a yaml file are loaded in order, numbered from one."""
    yaml_file = _write_yaml(
        tmp_path,
        {
            "calib": "0",
            "trans": "1",
            "list": [{"name": os.path.join(TESTDATADIR, "filt", name)} for name in filter_names],
        },
    )
    flt_array = FilterSvc.from_yaml(yaml_file)

    assert len(flt_array) == len(filter_names)
    assert [os.path.basename(f.name) for f in flt_array] == [os.path.basename(n) for n in filter_names]
    assert flt_array[0].lambdaMean() == pytest.approx(5262.2831, 1.0e-4)


def test_from_yaml_without_defaults(tmp_path, filter_names):
    """calib/trans may be omitted entirely, in which case they default to zero."""
    yaml_file = _write_yaml(
        tmp_path,
        {"list": [{"name": os.path.join(TESTDATADIR, "filt", filter_names[0])}]},
    )
    (flt_obj,) = FilterSvc.from_yaml(yaml_file)
    expected = FilterSvc.from_file(os.path.join(TESTDATADIR, "filt", filter_names[0]), 1, 0, 0)
    assert flt_obj.lambdaMean() == pytest.approx(expected.lambdaMean())


def test_from_yaml_per_entry_overrides(tmp_path, filter_names):
    """Per-entry calib/trans override the defaults declared at the top level."""
    yaml_file = _write_yaml(
        tmp_path,
        {
            "calib": "0",
            "trans": "0",
            "list": [
                {"name": os.path.join(TESTDATADIR, "filt", filter_names[0]), "trans": "1", "calib": "1"},
                {"name": os.path.join(TESTDATADIR, "filt", filter_names[1])},
            ],
        },
    )
    overridden, default = FilterSvc.from_yaml(yaml_file)

    # The overridden entry uses trans=1/calib=1, which changes the transmission
    # integral relative to the same filter read with the defaults.
    plain = FilterSvc.from_file(os.path.join(TESTDATADIR, "filt", filter_names[0]), 1, 0, 0)
    assert overridden.lambdaMean() != pytest.approx(plain.lambdaMean(), rel=1e-9)
    assert default.lambdaMean() == pytest.approx(
        FilterSvc.from_file(os.path.join(TESTDATADIR, "filt", filter_names[1]), 2, 0, 0).lambdaMean()
    )


def test_from_yaml_expands_lephare_dir(tmp_path, filter_names):
    """A $LEPHAREDIR prefix in a yaml entry is expanded before reading."""
    yaml_file = _write_yaml(
        tmp_path,
        {
            "calib": "0",
            "trans": "1",
            "list": [{"name": os.path.join("$LEPHAREDIR", "filt", filter_names[0])}],
        },
    )
    (flt_obj,) = FilterSvc.from_yaml(yaml_file)
    assert flt_obj.name == os.path.join(LEPHAREDIR, "filt", filter_names[0])


def test_from_file_expands_lephare_dir(filter_names):
    """from_file also expands $LEPHAREDIR before reading the filter."""
    flt_obj = FilterSvc.from_file(os.path.join("$LEPHAREDIR", "filt", filter_names[0]), 7, 1, 0)
    assert flt_obj.name == os.path.join(LEPHAREDIR, "filt", filter_names[0])
    assert flt_obj.lambdaMean() == pytest.approx(5262.2831, 1.0e-4)


def _config_lines(n_filters, n_trans, n_calib):
    """Build a minimal .para body with the given number of comma-separated entries."""
    names = ",".join(["cosmos/IB527.lowres"] * n_filters)
    return (
        f"FILTER_REP {os.path.join(TESTDATADIR, 'filt')}\n"
        f"FILTER_LIST {names}\n"
        f"TRANS_TYPE {','.join(['1'] * n_trans)}\n"
        f"FILTER_CALIB {','.join(['0'] * n_calib)}\n"
        f"FILTER_FILE filter_test\n"
    )


def test_from_config_broadcasts_single_values(tmp_path):
    """A single TRANS_TYPE/FILTER_CALIB is broadcast across all filters."""
    config_file = tmp_path / "broadcast.para"
    config_file.write_text(_config_lines(3, 1, 1))

    flt_array = FilterSvc.from_config(str(config_file))
    assert len(flt_array) == 3


def test_from_config_mismatched_trans(tmp_path):
    """A TRANS_TYPE list of the wrong length is rejected rather than silently truncated."""
    config_file = tmp_path / "bad_trans.para"
    config_file.write_text(_config_lines(3, 2, 1))

    with pytest.raises(RuntimeError, match="FILTER_LIST and FILTER_TRANS do not have the same size"):
        FilterSvc.from_config(str(config_file))


def test_from_config_mismatched_calib(tmp_path):
    """A FILTER_CALIB list of the wrong length is rejected."""
    config_file = tmp_path / "bad_calib.para"
    config_file.write_text(_config_lines(3, 3, 2))

    with pytest.raises(RuntimeError, match="FILTER_LIST and FILTER_CALIB do not have the same size"):
        FilterSvc.from_config(str(config_file))


def _keymap(n_filters, n_trans, n_calib):
    import lephare as lp

    return lp.all_types_to_keymap(
        {
            "FILTER_REP": os.path.join(TESTDATADIR, "filt"),
            "FILTER_LIST": ",".join(["cosmos/IB527.lowres"] * n_filters),
            "TRANS_TYPE": ",".join(["1"] * n_trans),
            "FILTER_CALIB": ",".join(["0"] * n_calib),
        }
    )


def test_from_keymap_broadcasts_single_values():
    """from_keymap broadcasts a single TRANS_TYPE/FILTER_CALIB like from_config does."""
    flt_array = FilterSvc.from_keymap(_keymap(3, 1, 1))
    assert len(flt_array) == 3
    assert all(f.lambdaMean() == pytest.approx(5262.2831, 1.0e-4) for f in flt_array)


def test_from_keymap_mismatched_trans():
    """A TRANS_TYPE list of the wrong length is rejected."""
    with pytest.raises(RuntimeError, match="FILTER_LIST and FILTER_TRANS do not have the same size"):
        FilterSvc.from_keymap(_keymap(3, 2, 1))


def test_from_keymap_mismatched_calib():
    """A FILTER_CALIB list of the wrong length is rejected."""
    with pytest.raises(RuntimeError, match="FILTER_LIST and FILTER_CALIB do not have the same size"):
        FilterSvc.from_keymap(_keymap(3, 3, 2))


class _FakeResponse:
    def __init__(self, content, date="Mon, 01 Jan 2024 00:00:00 GMT"):
        self.content = content
        self.headers = {"Date": date}


def test_svo_request_server_down(monkeypatch, capsys):
    """An unparseable response (e.g. a truncated error page) is reported, not raised."""
    monkeypatch.setattr(
        "lephare.filterSvc.requests.get", lambda *a, **k: _FakeResponse(b"<html>503 unavailable")
    )
    assert FilterSvc.svo_request(1, "Subaru/Suprime.IB527", "AB") is None
    assert "SVO server down" in capsys.readouterr().out


SVO_XML_TEMPLATE = """<?xml version="1.0"?>
<VOTABLE>
  <RESOURCE>
    <INFO name="QUERY_STATUS" value="{status}"/>
    <PARAM name="DetectorType" value="0"/>
    <TABLE>
      <FIELD name="Wavelength" unit="Angstrom" datatype="double"/>
      <FIELD name="Transmission" unit="" datatype="double"/>
      <DATA>
        <TABLEDATA>
          <TR><TD>5000.0</TD><TD>0.1</TD></TR>
          <TR><TD>5100.0</TD><TD>0.5</TD></TR>
          <TR><TD>5200.0</TD><TD>1.0</TD></TR>
          <TR><TD>5300.0</TD><TD>0.5</TD></TR>
          <TR><TD>5400.0</TD><TD>0.1</TD></TR>
        </TABLEDATA>
      </DATA>
    </TABLE>
  </RESOURCE>
</VOTABLE>
"""


def test_svo_request_bad_query_status(monkeypatch):
    """A QUERY_STATUS other than OK means the filter id was wrong."""
    monkeypatch.setattr(
        "lephare.filterSvc.requests.get",
        lambda *a, **k: _FakeResponse(SVO_XML_TEMPLATE.format(status="ERROR").encode()),
    )
    with pytest.raises(AssertionError, match="QUERY_STATUS did not return OK"):
        FilterSvc.svo_request(1, "Subaru/Suprime.NOPE", "AB")


def test_svo_request_parses_votable(monkeypatch, tmp_path):
    """A well-formed VOTable is turned into a filter, and the scratch file removed."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        "lephare.filterSvc.requests.get",
        lambda *a, **k: _FakeResponse(SVO_XML_TEMPLATE.format(status="OK").encode()),
    )
    flt_obj = FilterSvc.svo_request(4, "Subaru/Suprime.IB527", "AB")

    # flt pads the curve with a zero-transmission point at either end, so the
    # support is the tabulated range plus one grid step.
    assert flt_obj.lmin() == pytest.approx(5000.0, abs=15.0)
    assert flt_obj.lmax() == pytest.approx(5400.0, abs=15.0)
    # The symmetric triangular response peaks at the central wavelength
    assert flt_obj.lambdaMean() == pytest.approx(5200.0, rel=1e-3)
    # The SVO metadata is stashed on the returned object
    assert flt_obj.svo_params["DetectorType"] == "0"
    assert "Date" in flt_obj.svo_params
    # The temporary file written while parsing is cleaned up
    assert not (tmp_path / "Suprime.IB527").exists()


def test_from_svo_delegates_to_svo_request(monkeypatch):
    """from_svo is a thin wrapper that forwards counter, id and system."""
    seen = {}

    def fake_request(cls, counter, filter_id, system):
        seen.update(counter=counter, filter_id=filter_id, system=system)
        return "sentinel"

    monkeypatch.setattr(FilterSvc, "svo_request", classmethod(fake_request))
    assert FilterSvc.from_svo(3, "Subaru/Suprime.IB527", "AB", calib=1) == "sentinel"
    assert seen == {"counter": 3, "filter_id": "Subaru/Suprime.IB527", "system": "AB"}


def test_from_yaml_svo_entry(tmp_path, monkeypatch, filter_names):
    """A name prefixed with "svo:" is fetched from the SVO instead of from disk."""
    seen = []

    def fake_from_svo(cls, counter, filter_id, system, calib):
        seen.append((counter, filter_id, system, calib))
        return FilterSvc.from_file(os.path.join(TESTDATADIR, "filt", filter_names[0]), counter, 1, 0)

    monkeypatch.setattr(FilterSvc, "from_svo", classmethod(fake_from_svo))
    yaml_file = _write_yaml(
        tmp_path,
        {
            "calib": "0",
            "trans": "1",
            "list": [
                {"name": "svo:Subaru/Suprime.IB527"},
                {"name": os.path.join(TESTDATADIR, "filt", filter_names[1])},
            ],
        },
    )
    flt_array = FilterSvc.from_yaml(yaml_file)

    assert len(flt_array) == 2
    # The "svo:" prefix is stripped, the system forced to AB, and calib passed on
    assert seen == [(1, "Subaru/Suprime.IB527", "AB", 0)]

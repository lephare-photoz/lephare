"""Tests for how the sedtolib/mag_gal runners validate the object type."""

import os

import lephare as lp
import pytest

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")

BAD_TYPE_MESSAGE = "-t arg must start with G/g Q/q or S/s"


@pytest.fixture
def cosmos_keymap(set_env_vars):
    return lp.all_types_to_keymap(lp.read_config(os.path.join(TESTDATADIR, "examples/COSMOS.para")))


@pytest.mark.parametrize("runner_class", [lp.Sedtolib, lp.MagGal])
@pytest.mark.parametrize("bad_type", ["X", "1", "extragalactic"])
def test_unrecognised_object_type(runner_class, bad_type, cosmos_keymap):
    """A type whose first letter is not G, Q or S is rejected up front."""
    runner = runner_class(config_keymap=cosmos_keymap)
    with pytest.raises(KeyError, match=BAD_TYPE_MESSAGE):
        runner.run(typ=bad_type)


@pytest.mark.parametrize("runner_class", [lp.Sedtolib, lp.MagGal])
def test_empty_object_type(runner_class, cosmos_keymap):
    """An empty type cannot be indexed, so it fails rather than silently defaulting."""
    runner = runner_class(config_keymap=cosmos_keymap)
    with pytest.raises(IndexError):
        runner.run(typ="")


def test_sedtolib_object_type_is_case_insensitive(cosmos_keymap):
    """run() upper-cases the type before dispatching to the library class."""
    runner = lp.Sedtolib(config_keymap=cosmos_keymap)
    runner.run(typ="gal", gal_sed=os.path.join(TESTDATADIR, "sed/GAL/ONE_SED.list"))

    assert runner.typ == "GAL"
    assert isinstance(runner.SEDLib, lp.GalSEDLib)


def test_maggal_object_type_is_case_insensitive(cosmos_keymap):
    """mag_gal applies the same upper-casing to its type argument."""
    runner = lp.MagGal(config_keymap=cosmos_keymap)
    runner.run(typ="gal", VERBOSE="NO")

    assert runner.typ == "GAL"

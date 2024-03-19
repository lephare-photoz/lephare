import numpy as np
import os
import pytest

from lephare._lephare import get_lephare_env, test_first_char, blackbody


def test_globals_first_char():
    assert test_first_char("test")
    assert test_first_char("   test")
    assert not test_first_char("")
    assert not test_first_char(" ")
    assert not test_first_char("!test")
    assert not test_first_char("#test")
    assert not test_first_char("\t#")


def test_globals_blackbody():
    assert np.isclose(blackbody(10000, 500), 1.018807e-26)

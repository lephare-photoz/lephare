import numpy as np
from lephare._lephare import blackbody, test_first_char


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

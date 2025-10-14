import numpy as np
from lephare._lephare import blackbody, check_first_char, indexes_in_vec


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

import os

from lephare import _lephare


def test_one_ei_lambda_constructor():
    el1 = _lephare.oneElLambda(10, 10, 0)
    # test that 1 + 1 = 2
    assert el1.lamb == 10
    assert el1.val == 10
    assert el1.ori == 0
    el2 = _lephare.oneElLambda(el1)
    assert el1.lamb == el2.lamb
    assert el1.val == el2.val
    assert el1.ori == el2.ori


def test_one_ei_lambda_interp():
    el1 = _lephare.oneElLambda(10, 10, 0)
    el2 = _lephare.oneElLambda(20, 20, 0)
    el3 = _lephare.oneElLambda(15, 0, 0)
    el3.interp(el1, el2)
    assert el3.val == 15


def test_one_ei_lambda_constructor2():
    ext = _lephare.ext(name="test", numext=10)
    assert ext.name == "test"
    assert ext.numext == 10
    assert len(ext.lamb_ext) == 0


def test_one_ei_lambda_read_calzetti(test_data_dir):
    filename = os.path.join(test_data_dir, "calzetti.dat")
    ext = _lephare.ext(name="test", numext=10)
    ext.read(filename)
    assert ext.lmin == 400.0
    assert ext.lmax == 40000.0


def test_one_ei_lambda_read_tau10(test_data_dir):
    filename = os.path.join(test_data_dir, "tau10.out")
    ext = _lephare.ext(name="test", numext=10)
    ext.read(filename)

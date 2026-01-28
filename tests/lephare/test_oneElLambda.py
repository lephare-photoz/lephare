import os

import lephare as lp
import numpy as np
import pytest


def test_one_ei_lambda_constructor():
    el1 = lp.oneElLambda(10, 10)
    # test that 1 + 1 = 2
    assert el1.lamb == 10
    assert el1.val == 10
    el2 = lp.oneElLambda(el1)
    assert el1.lamb == el2.lamb
    assert el1.val == el2.val


def test_one_ei_lambda_constructor2():
    ext = lp.ext(name="test", numext=10)
    assert ext.name == "test"
    assert ext.numext == 10
    assert len(ext.lamb_ext) == 0


def test_one_ei_lambda_read_calzetti(test_data_dir):
    filename = os.path.join(test_data_dir, "calzetti.dat")
    ext = lp.ext(name="test", numext=10)
    ext.read(filename)
    assert ext.lmin == 400.0
    assert ext.lmax == 40000.0


def test_one_ei_lambda_read_tau10(test_data_dir):
    filename = os.path.join(test_data_dir, "tau10.out")
    ext = lp.ext(name="test", numext=10)
    ext.read(filename)


def test_concatenate_and_sort():
    x1 = np.linspace(0, 10, 11)
    y1 = np.ones_like(x1)
    x2 = x1 + 0.5
    y2 = np.zeros_like(x2)
    z1 = []
    z2 = []
    for i in range(10):
        z1.append(lp.oneElLambda(x1[i], y1[i]))
        z2.append(lp.oneElLambda(x2[i], y2[i]))
    z = lp.concatenate_and_sort(z1, z2)
    zl = [zz.lamb for zz in z]
    zv = [zz.val for zz in z]
    assert zl == [
        0.0,
        0.5,
        1.0,
        1.5,
        2.0,
        2.5,
        3.0,
        3.5,
        4.0,
        4.5,
        5.0,
        5.5,
        6.0,
        6.5,
        7.0,
        7.5,
        8.0,
        8.5,
        9.0,
        9.5,
    ]
    assert zv == [
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
        1.0,
        0.0,
    ]


def test_make_regular_grid():
    with pytest.raises(RuntimeError):
        grid = lp.make_regular_grid(0, 1, 0)
    with pytest.raises(RuntimeError):
        grid = lp.make_regular_grid(0, 1, -1)
    grid = lp.make_regular_grid(0, 1, 0.1)
    assert np.allclose(np.array(grid), np.linspace(0, 1, 11))


def test_common_interpolate_combined():
    with pytest.raises(RuntimeError):
        lp.common_interpolate_combined(
            np.linspace(10, 20, 10),
            np.linspace(10, 20, 11),
            np.linspace(10, 20, 12),
            np.linspace(10, 20, 12),
            1,
        )
    with pytest.raises(RuntimeError):
        lp.common_interpolate_combined(
            np.linspace(10, 20, 10),
            np.linspace(10, 20, 10),
            np.linspace(10, 20, 10),
            np.linspace(10, 20, 12),
            1,
        )

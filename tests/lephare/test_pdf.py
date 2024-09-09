import lephare as lp
import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")


def test_quadratic_extremum():
    # check alignment assertion
    def line(x):
        return 2 * x + 1

    x = [0, 1, 2]
    y = [line(t) for t in x]
    with pytest.raises(ValueError):
        xm, ym = lp.quadratic_extremum(*x, *y)

    for a in [+1, -1]:

        def parab(x, a):
            return a * x * x

        x = [-1, 0.5, 1]
        y = [parab(t, a) for t in x]
        xm, ym = lp.quadratic_extremum(*x, *y)
        assert xm == pytest.approx(0)
        assert ym == pytest.approx(0)

    # check for both concave and convex (min and max search)
    for a in [-2, 2]:
        b = 1
        c = 2
        x1 = 0.5

        def parab(x, a, b, c):
            return a * x**2 + b * x + c

        x = [x1 - 0.1, x1, x1 + 0.3]
        y = [parab(t, a, b, c) for t in x]
        xm, ym = lp.quadratic_extremum(*x, *y)
        assert xm == pytest.approx(-b / 2.0 / a)
        assert ym == pytest.approx(c - b**2 / 4.0 / a)

    # random parabola and random triplet
    p = -5 + 10 * np.random.random(3)  # a,b,c between -5 and 5

    def parab(x):
        return p[0] * x**2 + p[1] * b * x + p[2] * c

    xt = -p[1] / 2 / p[0]
    yt = parab(xt)
    x = xt - 0.1 + 0.2 * np.random.random(3)  # 3 pts at -0.1 +0.1 from extremum
    y = parab(x)
    xm, ym = lp.quadratic_extremum(*x, *y)
    assert xm == pytest.approx(xt)
    assert ym == pytest.approx(yt)


def test_pdf():
    """Basic test to ensure we can instantiate a PDF object."""
    test_pdf = lp.PDF(0, 3, 10)
    assert test_pdf.size() == 10
    assert len(test_pdf.chi2) == 10


def test_set_yvals():
    """Test to make sure that chi2 and vPDF are set correctly."""
    test_pdf = lp.PDF(0, 3, 10)
    yvals = np.linspace(0.001, 3, 10)
    test_pdf.setYvals(yvals, is_chi2=True)
    np.testing.assert_array_equal(test_pdf.chi2, yvals)
    np.testing.assert_array_almost_equal(test_pdf.vPDF, np.exp(-0.5 * yvals))

    test_pdf.setYvals(yvals, is_chi2=False)
    np.testing.assert_array_equal(test_pdf.vPDF, yvals)
    np.testing.assert_array_almost_equal(test_pdf.chi2, -2 * np.log(yvals))

    # checking that negative value in vPDF raises an exception
    with pytest.raises(ValueError):
        yvals[5] = -1
        test_pdf.setYvals(yvals, is_chi2=False)


def test_plot():
    """Simple test to exercise the plot method."""
    test_pdf = lp.PDF(0, 3, 10)
    yvals = np.linspace(0.001, 3, 10)
    test_pdf.setYvals(yvals, is_chi2=True)
    test_pdf.plot(chi2=True)
    test_pdf.plot(chi2=False)

    test_pdf.setYvals(yvals, is_chi2=False)
    test_pdf.plot(chi2=True)
    test_pdf.plot(chi2=False)


def test_credible_interval():
    test_pdf = lp.PDF(0, 1, 100)
    test_pdf.setYvals(np.ones_like(test_pdf.xaxis), is_chi2=False)
    # check output if argument val is outside of the xgrid range
    a, b = test_pdf.credible_interval(0.2, -1)
    assert a == -1
    assert b == -1
    a, b = test_pdf.credible_interval(0.2, 2)
    assert a == 2
    assert b == 2
    # case where, even as a percentage, level is too high
    a, b = test_pdf.credible_interval(250, 0.5)
    assert a == pytest.approx(0.0)
    assert b == pytest.approx(1)
    # standard cases in the case of a trivial pdf
    a, b = test_pdf.credible_interval(0.2, 0.5)
    assert a == pytest.approx(0.4)
    assert b == pytest.approx(0.6)
    a, b = test_pdf.credible_interval(0.2, 0.9)
    assert a == pytest.approx(0.8)
    assert b == pytest.approx(1.0)
    a, b = test_pdf.credible_interval(0.2, 0.1)
    assert a == pytest.approx(0)
    assert b == pytest.approx(0.2)
    a, b = test_pdf.credible_interval(0.2, 0.0)
    assert a == pytest.approx(0)
    assert b == pytest.approx(0.2)
    # more reasonable pdf
    test_pdf = lp.PDF(0, 10, 1000)

    def gaus(x):
        return np.exp(-((x - 5.0) ** 2) / 2.0)

    test_pdf.setYvals(gaus(np.array(test_pdf.xaxis)), is_chi2=False)
    a, b = test_pdf.credible_interval(68.26, 5)
    assert a == pytest.approx(4, 0.0001)
    assert b == pytest.approx(6, 0.0001)
    a, b = test_pdf.credible_interval(95.46, 5)
    assert a == pytest.approx(3, 0.001)
    assert b == pytest.approx(7, 0.001)

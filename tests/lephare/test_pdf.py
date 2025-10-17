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
    # testing bad input
    with pytest.raises(ValueError):
        test_pdf = lp.PDF(0, 3, 0)
    with pytest.raises(ValueError):
        test_pdf = lp.PDF(3, 0, 10)


def test_normalization():
    test_pdf = lp.PDF(0, 3.0, 101)
    yvals = np.linspace(0, 3.0, 101)
    test_pdf.setYvals(yvals, is_chi2=False)
    integral = test_pdf.normalization()
    assert integral == pytest.approx(4.5)
    integral2 = test_pdf.normalization()
    assert integral2 == pytest.approx(1)


def test_cumulant():
    pdf = lp.PDF(0.0, 1.0, 10)
    yvals = np.ones_like(pdf.xaxis)
    pdf.setYvals(yvals, is_chi2=False)
    cumulant = pdf.cumulant()
    np.testing.assert_almost_equal(cumulant, np.linspace(0, 1, 10))


def test_index():
    # 0, 0.33, 0.66,..., 2.66, 3 (10 values)
    pdf = lp.PDF(0, 3, 10)
    assert pdf.index(-1) == 0
    assert pdf.index(20) == 9
    assert pdf.index(2) == 6
    assert pdf.index(2.1) == 6
    assert pdf.index(2.2) == 7
    assert pdf.index(0) == 0
    assert pdf.index(3) == 9


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
    a, b = test_pdf.credible_interval(0.0, 2)
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
    assert a == pytest.approx(4, 0.01)
    assert b == pytest.approx(6, 0.01)
    a, b = test_pdf.credible_interval(95.46, 5)
    assert a == pytest.approx(3, 0.01)
    assert b == pytest.approx(7, 0.01)

    # analytical pdf
    x = np.linspace(-1, 1, 10001)
    y = np.zeros_like(x)
    y[x >= 0] = -x[x >= 0] + 1
    y[x < 0] = x[x < 0] + 1
    # triangular symmetric, area=1
    pdf = lp.PDF(-1, 1, 10001)
    pdf.setYvals(y, is_chi2=False)
    # central value
    val = 0.0
    for level in np.linspace(0.0, 1.0, 100):
        a, b = pdf.credible_interval(level, val)
        assert a == pytest.approx(-b)
        assert b == pytest.approx(1 - np.sqrt(1 - level))
    for val in [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]:
        for level in [0.1, 0.3, 0.6, 0.9]:
            a, b = pdf.credible_interval(level, val)
            cumul_left = 0.5 - val**2 / 2.0 + val if val > 0 else 0.5 + val**2 / 2.0 + val
            if cumul_left < level / 2.0:
                assert a == -1.0
                if level < 0.5:
                    assert b == pytest.approx(-1 + np.sqrt(2 * level))
                else:
                    assert b == pytest.approx(1 - np.sqrt(2 - 2 * level))
            if cumul_left > 1.0 - level / 2.0:
                assert b == 1.0
                if level < 0.5:
                    assert a == pytest.approx(1 - np.sqrt(2 * level))
                else:
                    assert a == pytest.approx(-1 + np.sqrt(2 - 2 * level))


def test_credible_interval2():
    pdf = lp.PDF(0, 1, 1000)
    pdf.setYvals(np.linspace(2, 1, len(pdf.xaxis)), is_chi2=False)

    # symetric case
    v = 0.5
    for level in [0.01, 0.1, 0.2, 0.34, 0.5]:
        lo, up = pdf.credible_interval(level, v)
        assert lo == pytest.approx(2 - np.sqrt((2 - v) ** 2 + 1.5 * level), 1.0e-4)
        assert up == pytest.approx(2 - np.sqrt((2 - v) ** 2 - 1.5 * level), 1.0e-4)

    # asymetric right case
    v = 0.1
    for level in [0.01, 0.1, 0.2, 0.34, 0.5]:
        lo, up = pdf.credible_interval(level, v)
        if 2 - np.sqrt((2 - v) ** 2 + 1.5 * level) < 0:
            assert lo == 0.0
            assert up == pytest.approx(2 - np.sqrt(4 - 3 * level), 1.0e-4)

    v = 0.9
    for level in [0.01, 0.1, 0.2, 0.34, 0.5]:
        lo, up = pdf.credible_interval(level, v)
        if 0.5 * (2 - v) ** 2 - 0.5 > 1:
            assert up == 0.0
            assert lo == pytest.approx(2 - np.sqrt(3 * level + 1))


def test_credible_interval3():
    # this is a case that caused problems in the past
    pdf = lp.PDF(3, 12.95, 200)
    yvals = np.zeros_like(pdf.xaxis)
    yvals[[119, 122, 123, 124, 126, 127]] = [
        1.90009e01,
        5.56300e-01,
        3.54100e-01,
        8.18000e-02,
        6.70000e-03,
        3.00000e-04,
    ]

    pdf.setYvals(yvals, is_chi2=False)
    a, b = pdf.credible_interval(0.68, 8.9526)
    assert a == pytest.approx(8.9168, 1.0e-4)
    assert b == pytest.approx(8.9884, 1.0e-4)


def test_levelcumu2x():
    pdf = lp.PDF(0, 1, 1000)
    pdf.setYvals(np.linspace(2, 1, len(pdf.xaxis)), is_chi2=False)
    assert pdf.levelCumu2x(0.5) == pytest.approx(2 - np.sqrt(5 / 2))


def test_improve_extremum():
    pdf = lp.PDF(0, 1, 101)
    xaxis = np.array(pdf.xaxis)
    for arg in [True, False]:
        a, b = pdf.improve_extremum(arg)
        assert a == pdf.xaxis[0]
        if arg:
            assert b == lp.HIGH_CHI2
        else:
            assert b == 0.0

    # tweak chi2 to check that parabolic improvement is
    # not done if one of the points is at HIGH_CHI2
    yval = np.array(pdf.chi2)
    yval[50] = 1
    yval[51] = 2
    pdf.setYvals(yval, is_chi2=True)
    a, b = pdf.improve_extremum(True)
    assert a == pdf.xaxis[50]
    assert b == 1
    # tweak vPDF to the same effect, with one of the points at 0
    yval = np.zeros_like(pdf.xaxis)
    yval[50] = 1
    yval[51] = 2
    pdf.setYvals(yval, is_chi2=False)
    a, b = pdf.improve_extremum(False)
    assert a == pdf.xaxis[51]
    assert b == 2

    # 2 normal situations
    pdf.setYvals((xaxis - 0.45) ** 2 + 0.2, True)
    a, b = pdf.improve_extremum(True)
    assert a == pytest.approx(0.45)
    assert b == pytest.approx(0.2)

    y = -((xaxis - 0.45) ** 2) + 0.2
    y[y < 0] = 0.0
    pdf.setYvals(y, False)
    a, b = pdf.improve_extremum(False)
    assert a == pytest.approx(0.45)
    assert b == pytest.approx(0.2)

    # Problematic situations
    pdf = lp.PDF(0, 1, 101)
    xaxis = np.array(pdf.xaxis)
    # extrememum at upper boundary
    pdf.setYvals((xaxis - 1.1) ** 2 + 0.2, True)
    a, b = pdf.improve_extremum(True)
    assert a == 1.0
    assert b == pytest.approx((a - 1.1) ** 2 + 0.2)
    pdf.setYvals(-((xaxis - 1.1) ** 2) + 2, False)
    a, b = pdf.improve_extremum(False)
    assert a == 1.0
    assert b == pytest.approx(-((a - 1.1) ** 2) + 2)
    # extrememum at lower boundary
    pdf = lp.PDF(0, 1, 101)
    xaxis = np.array(pdf.xaxis)
    # extrememum at upper boundary
    pdf.setYvals((xaxis) ** 2 + 0.2, True)
    a, b = pdf.improve_extremum(True)
    assert a == 0.0
    assert b == pytest.approx((a) ** 2 + 0.2)
    pdf.setYvals(-((xaxis) ** 2) + 2, False)
    a, b = pdf.improve_extremum(False)
    assert a == 0.0
    assert b == pytest.approx(-((a) ** 2) + 2)
    # HIGH_CHI2 issue
    pdf = lp.PDF(0, 1, 101)
    xaxis = np.array(pdf.xaxis)
    for i in range(3):
        # chi2 and a val at HIGH_CHI2
        yvals = (xaxis - 0.5) ** 2 + 0.2
        yvals[50 - 1 + i] = lp.HIGH_CHI2
        pdf.setYvals(yvals, True)
        a, b = pdf.improve_extremum(True)
        expected_min = 0.49 if i == 1 else 0.5
        assert a == pytest.approx(expected_min)
        assert b == pytest.approx((a - 0.5) ** 2 + 0.2)
        # vPDF and a val at 0
        yvals = -((xaxis - 0.5) ** 2) + 0.3
        yvals[50 - 1 + i] = 0
        pdf.setYvals(yvals, False)
        a, b = pdf.improve_extremum(False)
        expected_min = 0.49 if i == 1 else 0.5
        assert a == pytest.approx(expected_min)
        assert b == pytest.approx(-((a - 0.5) ** 2) + 0.3)
    # check the delta_x<0.5 condition
    pdf = lp.PDF(0, 1, 3)
    pdf.setYvals([0, 1, 0], False)
    assert a == 0.5
    assert b == 0.3  # normalized


def test_confidence_interval():
    pdf = lp.PDF(-5, 5, 1001)
    pdf.setYvals(np.array(pdf.xaxis) ** 2, is_chi2=True)
    a, b = pdf.confidence_interval(1)
    assert a == pytest.approx(-1)
    assert b == pytest.approx(1)
    a, b = pdf.confidence_interval(2.71)
    assert a == pytest.approx(-np.sqrt(2.71), 1.0e-5)
    assert b == pytest.approx(np.sqrt(2.71), 1.0e-5)

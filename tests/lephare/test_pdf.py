import lephare as lp
import matplotlib
import numpy as np

matplotlib.use("Agg")


def test_pdf():
    """Basic test to ensure we can instantiate a PDF object."""
    test_pdf = lp.PDF(0, 3, 10)
    assert test_pdf.size() == 10
    assert len(test_pdf.chi2) == 10


def test_set_yvals():
    """Test to make sure that chi2 and vPDF are set correctly."""
    test_pdf = lp.PDF(0, 3, 10)
    yvals = np.linspace(0, 3, 10)
    test_pdf.setYvals(yvals, is_chi2=True)
    np.testing.assert_array_equal(test_pdf.chi2, yvals)
    np.testing.assert_array_almost_equal(test_pdf.vPDF, np.exp(-0.5 * yvals))

    test_pdf.setYvals(yvals, is_chi2=False)
    np.testing.assert_array_equal(test_pdf.vPDF, yvals)
    np.testing.assert_array_almost_equal(test_pdf.chi2, -2 * np.log(yvals))


def test_plot():
    """Simple test to exercise the plot method."""
    test_pdf = lp.PDF(0, 3, 10)
    yvals = np.linspace(0, 3, 10)
    test_pdf.setYvals(yvals, is_chi2=True)
    test_pdf.plot(chi2=True)
    test_pdf.plot(chi2=False)

    test_pdf.setYvals(yvals, is_chi2=False)
    test_pdf.plot(chi2=True)
    test_pdf.plot(chi2=False)

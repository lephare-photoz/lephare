import numpy as np
from matplotlib import pylab as plt

from ._lephare import HIGH_CHI2, PDF
from ._utils import continueClass

__all__ = [
    "PDF",
]


@continueClass
class PDF:  # noqa
    def setYvals(self, yvals, is_chi2=False):  # noqa: N802
        if is_chi2:
            yvals[yvals >= HIGH_CHI2] = HIGH_CHI2
            self.chi2 = yvals
            yvals2 = np.exp(-0.5 * yvals)
            yvals2[yvals >= HIGH_CHI2] = 0.0
            self.vPDF = yvals2
        else:
            if np.any(np.array(yvals) < 0):
                raise ValueError("input array has negative values")
            self.vPDF = yvals
            yvals2 = np.zeros_like(yvals)
            yvals2[yvals == 0.0] = HIGH_CHI2
            yvals2[yvals != 0.0] = -2 * np.log(yvals[yvals != 0])
            self.chi2 = yvals2

    def plot(self, chi2=False, kwargs={}):  # noqa: B006
        if chi2:
            plt.plot(self.xaxis, self.chi2, **kwargs)
        else:
            plt.plot(self.xaxis, self.vPDF, **kwargs)

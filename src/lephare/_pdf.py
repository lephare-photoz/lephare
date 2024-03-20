import numpy as np
from matplotlib import pylab as plt

from ._lephare import PDF
from ._utils import continueClass

__all__ = [
    "PDF",
]


@continueClass
class PDF:  # noqa
    def setYvals(self, yvals, is_chi2=False):  # noqa: N802
        if is_chi2:
            self.chi2 = yvals
            self.vPDF = np.exp(-0.5 * yvals)
        else:
            self.vPDF = yvals
            self.chi2 = -2 * np.log(yvals)

    def plot(self, chi2=False, kwargs={}):  # noqa: B006
        if chi2:
            plt.plot(self.xaxis, self.chi2, **kwargs)
        else:
            plt.plot(self.xaxis, self.vPDF, **kwargs)

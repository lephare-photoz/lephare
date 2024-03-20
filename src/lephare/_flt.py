from matplotlib import pylab as plt

from ._lephare import flt
from ._utils import continueClass

__all__ = [
    "flt",
]


@continueClass
class flt:  # noqa
    def plot_filter_curve(self, normed=False):
        filter_name = self.name
        plt.title(filter_name)
        plt.xlabel("wavelength")
        plt.ylabel("value")
        x = self.data()[0]
        y = self.data()[1]
        if normed:
            y /= y.max()
        plt.plot(x, y, label="curve")
        plt.plot([self.lambdaMean(), self.lambdaMean()], [0.9 * y.min(), 1.1 * y.max()], label="lambda mean")
        plt.plot([self.lambdaEff(), self.lambdaEff()], [0.9 * y.min(), 1.1 * y.max()], label="lambda eff")
        plt.legend()

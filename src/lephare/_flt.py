"""Python-only additions to the C++-bound :class:`lephare.flt` class.

This module does not define a new class: ``@continueClass`` (see
:mod:`lephare._utils`) monkey-patches the methods defined below directly
onto the compiled ``lephare._lephare.flt`` class (exposed as
:class:`lephare.flt`). The full attribute/method reference for that class
--- everything defined on the C++ side --- is in the C++ API documentation
(``flt`` in ``src/lib/flt.h``); this module only adds the plotting
convenience method(s) below, which are pure Python and have no C++
equivalent.
"""

from matplotlib import pylab as plt

from lephare import flt

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

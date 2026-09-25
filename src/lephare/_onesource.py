"""Python-only additions to the C++-bound :class:`lephare.onesource` class.

This module does not define a new class: ``@continueClass`` (see
:mod:`lephare._utils`) monkey-patches the methods defined below directly
onto the compiled ``lephare._lephare.onesource`` class (exposed as
:class:`lephare.onesource`). The full attribute/method reference for that
class --- everything defined on the C++ side, including the fit itself,
``ab``/``sab``/``mabs``, the PDFs in ``pdfmap``, etc. --- is in the C++ API
documentation (``onesource`` in ``src/lib/onesource.h``); this module only
adds the pure-Python convenience method(s) below, which have no C++
equivalent.
"""

from ._lephare import onesource
from ._utils import continueClass

__all__ = [
    "onesource",
]


@continueClass
class onesource:  # noqa
    def compute_quality_flag(self):  # pragma no cover
        """
        Compute and return the quality flag for the bayesian PDF and center.

        This method retrieves the PDF object from `self.pdfmap[11]`, uses the first
        element of `self.zgmode` as the center value, and computes a quality flag
        using the PDF's `compute_pdz_flag` method.

        Returns:
        The computed quality flag.
        """
        pdf = self.pdfmap[11]
        center = self.zgmode[0]
        tup = pdf.compute_quality_flag(center)
        self.flag = tup[0]
        return self.flag

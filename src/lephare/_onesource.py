from ._lephare import onesource
from ._utils import continueClass

__all__ = [
    "onesource",
]


@continueClass
class onesource:  # noqa
    def compute_pdz_flag(self):
        """
        Compute and return the quality flag for the bayesian PDF and center.

        This method retrieves the PDF object from `self.pdfmap[1]`, uses the first
        element of `self.zgmode` as the center value, and computes a quality flag
        using the PDF's `compute_pdz_flag` method.

        Returns:
        The computed quality flag.
        """
        pdf = self.pdfmap[11]
        center = self.zgmode[0]
        tup = pdf.compute_pdz_flag(center)
        self.flag = tup[0]
        return self.flag

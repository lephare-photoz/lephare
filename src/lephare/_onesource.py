from ._lephare import onesource
from ._utils import continueClass

__all__ = [
    "onesource",
]


@continueClass
class onesource:  # noqa
    def compute_quality_flag(self):
        """
        Compute and return the quality flag for the current PDF and center.

        This method retrieves the PDF object from `self.pdfmap[9]`, uses the first
        element of `self.zgmin` as the center value, and computes a quality flag
        using the PDF's `compute_quality_flag` method.

        Returns:
        The computed quality flag.
        """
        pdf = self.pdfmap[9]
        center = self.zgmin[0]
        tup = pdf.compute_quality_flag(center)
        self.flag = tup[0]
        return self.flag

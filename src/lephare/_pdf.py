import numpy as np
from matplotlib import pylab as plt
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

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

    def variance(self, estimate):
        """Compute the pseudo variance of the P(z) around estimate"""
        var = trapezoid((np.array(self.xaxis) - estimate) ** 2 * self.vPDF, self.xaxis)
        return np.sqrt(var)

    def approximate_gaussian(self, estimate):
        def gauss(z, a, mu, sigma):
            return a * np.exp(-((z - mu) ** 2) / (2 * sigma**2))

        error = self.variance(estimate)
        mask = (self.xaxis >= estimate - error) & (self.xaxis <= estimate + error)
        z_local = self.xaxis
        pdz_local = np.where(mask, self.vPDF, 0.0)

        p0 = [np.max(pdz_local), estimate, error]
        try:
            popt, _ = curve_fit(gauss, z_local, pdz_local, p0=p0)
            _, mu_fit, siga_fit = popt
            return abs(siga_fit)
        except RuntimeError:
            # Fallback in case fit fails
            return error

    def number_mod(self, threshold=0.75):
        """Count significant local maxima"""
        peaks, _ = find_peaks(self.vPDF, height=threshold * max(self.vPDF))
        return len(peaks)

    def peak_ratio(self):
        """Ratio of max(pdz) to the mean"""
        return np.mean(self.vPDF) / max(self.vPDF)

    def tail_mass(self, estimate, n_window=2, good_sigma=0.01):
        """
        Compute the total probability mass in the tails,
        outside a window around z_best.
        """
        window = self.approximate_gaussian(estimate)
        bound = n_window * window
        mask = (self.xaxis < estimate - bound) | (self.xaxis > estimate + bound)
        if window <= good_sigma:
            return 0.0
        else:
            return trapezoid(np.where(mask, self.vPDF, 0.0), self.xaxis)

    def compute_quality_flag(
        self,
        estimate,
        nb_peak_thresh=1,
        height_thresh=0.5,
        tail_thresh=0.2,
        peak_ratio_thresh=0.25,
        error_thresh=0.2,
    ):
        """
        Compute a 4-bit quality score (0–15) based on PDF shape characteristics.

        The score encodes different quality aspects of a probability distribution,
        where higher values indicate lower quality. Each bit corresponds to a
        specific diagnostic criterion:

        - Bit 0 (value 1): High error variance (above `error_thresh`)
        - Bit 1 (value 2): High peak ratio (above `peak_ratio_thresh`)
        - Bit 2 (value 4): Large tail mass (above `tail_thresh`)
        - Bit 3 (value 8): Multiple peaks (more than `nb_peak_thresh`)
        - Bit 4 (value 16): Extreme error (above `max(xaxis) / 2.5`)

        Parameters
        ----------
        estimate : float
            Central estimate or mean value used for statistics.
        nb_peak_thresh : int, optional
            Maximum acceptable number of peaks. Defaults to 1.
        height_thresh : float, optional
            Threshold for peak detection height. Defaults to 0.5.
        tail_thresh : float, optional
            Threshold for acceptable tail mass. Defaults to 0.2.
        peak_ratio_thresh : float, optional
            Threshold for acceptable peak ratio. Defaults to 0.25.
        error_thresh :float, optional
            Threshold for acceptable variance. Defaults to 0.2.

        Returns a tuple containing:

        - **int**: Quality flag score (0–15)
        - **float**: Estimate value
        - **float**: Error (variance)
        - **float**: Peak ratio
        - **float**: Tail mass
        - **int**: Number of modes (peaks)
        - **float**: Approximate Gaussian sigma

        Notes:
        - A higher score indicates a worse fit or less reliable PDF.
        - The score is cumulative: multiple conditions may be triggered simultaneously.
        """
        score = 0

        # Compute parameters from PDZSTats
        number_mod = self.number_mod(threshold=height_thresh)
        tail_mass = self.tail_mass(estimate)
        peak_ratio = self.peak_ratio()
        error = self.variance(estimate)
        sigma = self.approximate_gaussian(estimate)

        # Bit 0: sigma
        if error > error_thresh and error < np.max(self.xaxis) / 2:
            score += 1
        # Bit 2: peak_ratio
        if peak_ratio > peak_ratio_thresh:
            score += 2
        # Bit 1: tail mass
        if tail_mass > tail_thresh:
            score += 4
        # Bit 0: number of peaks
        if number_mod > nb_peak_thresh:
            score += 8
        # Bit 3: error
        if error > np.max(self.xaxis) / 2.5:
            score += 16

        return int(score), estimate, error, peak_ratio, tail_mass, number_mod, sigma

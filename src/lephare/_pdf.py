"""Python-only additions to the C++-bound :class:`lephare.PDF` class.

This module does not define a new class: ``@continueClass`` (see
:mod:`lephare._utils`) monkey-patches the methods defined below directly
onto the compiled ``lephare._lephare.PDF`` class (exposed as
:class:`lephare.PDF`). The full attribute/method reference for that class
--- ``vPDF``, ``xaxis``, ``secondMax``, ``cumulant``, etc. --- is in the
C++ API documentation (``PDF`` in ``src/lib/PDF.h``); this module only adds
the pure-Python convenience methods below (plotting, peak-fitting helpers),
which have no C++ equivalent.
"""

import warnings

import numpy as np
from matplotlib import pylab as plt
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from ._lephare import HIGH_CHI2, PDF, FlagParams
from ._utils import continueClass

__all__ = [
    "PDF",
]


@continueClass
class PDF:  # noqa
    def setYvals(self, yvals, is_chi2=False):  # noqa: N802
        yvals = np.array(yvals)
        if is_chi2:
            yvals[yvals >= HIGH_CHI2] = HIGH_CHI2
            self.chi2 = yvals
            yvals2 = np.exp(-0.5 * yvals)
            yvals2[yvals >= HIGH_CHI2] = 0.0
            self.vPDF = yvals2
        else:
            if np.any(yvals < 0):
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

    def approximate_gaussian(self, estimate, error=None):
        def gauss(z, a, mu, sigma):
            return a * np.exp(-((z - mu) ** 2) / (2 * sigma**2))

        if error is None:
            error = self.variance(estimate)
        mask = (np.array(self.xaxis) >= estimate - error) & (np.array(self.xaxis) <= estimate + error)
        z_local = self.xaxis
        pdz_local = np.where(mask, self.vPDF, 0.0)

        p0 = [np.max(pdz_local), estimate, error]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, _ = curve_fit(gauss, z_local, pdz_local, p0=p0)

            _, mu_fit, sigma_fit = popt
            return abs(sigma_fit)
        except (RuntimeError, ValueError):
            # Return error if fitting fails
            return error

    def number_mod(self, threshold=0.43, distance=10):
        """Count significant local maxima"""
        if max(self.vPDF) != 0:
            peaks, _ = find_peaks(self.vPDF, height=threshold * max(self.vPDF), distance=distance)
            return len(peaks)
        else:
            return 0

    def peak_ratio(self):
        """Ratio of max(pdz) to the mean"""
        if max(self.vPDF) != 0:
            return np.mean(self.vPDF) / max(self.vPDF)
        else:
            return 0

    def tail_mass(self, estimate, sigma=None, n_window=2, good_sigma=0.01):
        """
        Compute the total probability mass in the tails,
        outside a window around z_best.
        """
        if sigma is None:
            sigma = self.approximate_gaussian(estimate)

        bound = n_window * sigma
        mask = (self.xaxis < estimate - bound) | (self.xaxis > estimate + bound)

        if sigma <= good_sigma:
            return 0.0

        return trapezoid(np.where(mask, self.vPDF, 0.0), self.xaxis)

    def compute_quality_flag(
        self,
        estimate,
        nb_peak_thresh=2,
        height_thresh=0.43,
        tail_thresh=0.23,
        peak_ratio_thresh=0.1,
        error_thresh=0.1,
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

        Args:
        estimate (float): Central estimate or mean value used for statistics.
        nb_peak_thresh (int, optional): Maximum acceptable number of peaks.
            Defaults to 1.
        height_thresh (float, optional): Threshold for peak detection height.
            Defaults to 0.5.
        tail_thresh (float, optional): Threshold for acceptable tail mass.
            Defaults to 0.2.
        peak_ratio_thresh (float, optional): Threshold for acceptable peak ratio.
            Defaults to 0.25.
        error_thresh (float, optional): Threshold for acceptable variance.
            Defaults to 0.2.

        Returns:
        tuple:
            A tuple containing:
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
        good_sigma = self.xaxis[1] - self.xaxis[0]
        # Compute parameters from PDZSTats
        error = self.variance(estimate)
        sigma = self.approximate_gaussian(estimate, error=error)
        tail_mass = self.tail_mass(estimate, sigma=sigma, good_sigma=good_sigma)
        number_mod = self.number_mod(threshold=height_thresh)
        peak_ratio = self.peak_ratio()

        # Bit 0: error
        if error > error_thresh and error < (np.max(self.xaxis) - np.min(self.xaxis)) / 2:
            score += 1
        # Bit 2: peak_ratio
        if peak_ratio > peak_ratio_thresh:
            score += 2
        # Bit 1: tail mass
        if tail_mass > tail_thresh:
            score += 4
        # Bit 0: number of peaks
        if number_mod >= nb_peak_thresh:
            score += 8
        # Bit 3: catastrophic error
        if error >= (np.max(self.xaxis) - np.min(self.xaxis)) / 2:
            score += 16

        return int(score), estimate, error, peak_ratio, tail_mass, number_mod, sigma

    def plot_flag_criteria(self, params=None, estimate=None, ax=None, title=None):
        """Plot the PDF together with the criteria that build the redshift flag.

        The flag is the one computed by :meth:`lephare.PDF.quality_flag` (the
        one written in the ``Z_FLAG`` output column). The figure shows:

        - the PDF, normalised to its maximum;
        - the Gaussian of same peak height and of width ``sigma`` (the
          width of the main peak at ``fwhm_level`` of its maximum, converted
          to a standard deviation), used to define the core of the PDF;
        - the estimate (by default the PDF maximum);
        - the core, ``estimate`` +/- ``n_window`` * ``sigma``; the probability
          mass of the PDF outside of it (shaded) is the tail mass (bit 4);
        - the peak detection threshold ``height`` (fraction of the maximum)
          and the peaks which are counted (bit 8), after discarding peaks
          closer than ``peak_dz``;
        - a text box, next to the axes, with the error, i.e. the pseudo standard
          deviation of the PDF (bits 1 and 16), the ratio mean/max of the PDF
          (bit 2), the tail mass, the number of peaks, and the resulting flag,
          in which the criteria that are triggered are marked.

        Parameters
        ----------
        params : FlagParams or str, optional
            Parameters of the flag, or a ``"name=value,..."`` string as for the
            ``Z_FLAG_PAR`` keyword. Default: the default :class:`FlagParams`.
        estimate : float, optional
            Central value of the statistics. Default: the PDF maximum.
        ax : matplotlib.axes.Axes, optional
            Axes to draw on. Default: the current axes.
        title : str, optional
            Title of the plot.

        Returns
        -------
        ZFlagResult
            The flag and the quantities it is made of.
        """
        if params is None:
            params = FlagParams()
        elif isinstance(params, str):
            params = FlagParams.from_string(params)
        x = np.asarray(self.xaxis, dtype=float)
        v = np.asarray(self.vPDF, dtype=float)
        res = self.quality_flag(params) if estimate is None else self.quality_flag(params, estimate)
        if ax is None:
            ax = plt.gca()
        vmax = np.max(v) if v.size else 0.0
        if res.flag < 0 or not np.isfinite(vmax) or vmax <= 0:
            ax.text(0.5, 0.5, f"invalid PDF (Z_FLAG={res.flag})", transform=ax.transAxes, ha="center")
            return res
        q = v / vmax
        est, sigma = res.estimate, res.sigma
        dx = x[1] - x[0]

        ax.plot(x, q, color="C0", label="P(z)")
        if sigma > 0:
            gauss = np.exp(-0.5 * ((x - est) / sigma) ** 2)
            ax.plot(x, gauss, "--", color="C1", label=f"Gaussian, $\\sigma$={sigma:.3g}")
        ax.axvline(est, color="purple", ls="--", lw=1, label=f"estimate = {est:.3f}")
        ax.axhline(params.height, color="red", lw=1, label=f"peak threshold ({params.height:g} max)")

        # core and tail mass
        lo, hi = est - params.n_window * sigma, est + params.n_window * sigma
        ax.axvspan(lo, hi, color="gray", alpha=0.12, label=f"core $\\pm${params.n_window:g}$\\sigma$")
        tail = (x < lo) | (x > hi)
        ax.fill_between(x, 0, q, where=tail, color="C3", alpha=0.35, label=f"tail mass = {res.tail_mass:.3f}")

        # peaks counted for bit 8 (same convention as scipy.signal.find_peaks)
        dist = max(1, int(np.ceil(params.peak_dz / dx - 1.0e-6)))
        peaks, _ = find_peaks(q, height=params.height, distance=dist)
        if len(peaks):
            ax.plot(x[peaks], q[peaks], "v", color="k", ms=7, label=f"{len(peaks)} peak(s)")

        # criteria and flag
        rng = x[-1] - x[0]
        crit = [
            (1, f"error = {res.error:.3f} > {params.error:g}", params.error < res.error < rng / 2),
            (2, f"mean/max={res.peak_ratio:.3f} > {params.peak_ratio:g}", res.peak_ratio > params.peak_ratio),
            (4, f"tail = {res.tail_mass:.3f} > {params.tail:g}", res.tail_mass > params.tail),
            (8, f"n_peaks = {res.n_peaks} >= {params.nb_peak}", res.n_peaks >= params.nb_peak),
            (16, f"error = {res.error:.3f} >= range/2", res.error >= rng / 2),
        ]
        lines = [f"Z_FLAG = {res.flag}"]
        lines += [f"[{'x' if on else ' '}] {bit:>2}: {txt}" for bit, txt, on in crit]
        ax.text(
            1.02,
            0.0,
            "\n".join(lines),
            transform=ax.transAxes,
            va="bottom",
            ha="left",
            family="monospace",
            fontsize=8,
            bbox={"boxstyle": "round", "fc": "white", "alpha": 0.85},
        )
        ax.set_xlabel("Redshift")
        ax.set_ylabel("P(z) / max")
        ax.set_ylim(0, 1.08)
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8)
        if title is not None:
            ax.set_title(title)
        return res

import argparse
import os
import warnings
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit
from pathlib import Path

###LePHARE sorting script###
'''
A python cript to sort good from bad zphota computed data using :
 - source indents from OUTPUT_CAT[.out]
 - BAY_ZG.prob from PDZ_OUT keyword when running z_phota
'''

def parse_args():
    parser = argparse.ArgumentParser(description="Compute PDZ quality flags and add them to a LePHARE catalog.")

    parser.add_argument("--pdz", required=True, help="Path to the PDZ file (bay_zg.prob)")

    parser.add_argument("--zgrid", required=True, help="Redshift grid as ZSTEP,ZMIN,ZMAX (dz,zmin,zmax)")

    parser.add_argument("--row", type=int, default=0, help="Maximum number of PDZ rows to process (None)")

    return parser.parse_args()

args = parse_args()
PDZ_PATH = args.pdz
row = args.row
ZGRID = args.zgrid
ZSTEP, ZMIN, ZMAX = np.asarray(ZGRID.split(',')).astype(float)
z_grid = np.arange(ZMIN, ZMAX + ZSTEP, ZSTEP)

args = parse_args()


class PDZStats:
    """Compute statistical metrics on LePHARE P(z) distributions."""

    def __init__(self, zgrid, pdz, z_best=None):
        self.zgrid = np.asarray(zgrid, dtype=float)
        pdz = np.asarray(pdz, dtype=float)

        if self.zgrid.ndim != 1 or pdz.ndim != 1:
            raise ValueError("zgrid and pdz must be 1D arrays")
        if len(self.zgrid) != len(pdz):
            raise ValueError("zgrid and pdz must have the same length")
        if np.any(pdz < 0):
            raise ValueError("input PDF contains negative values")

        area = trapezoid(pdz, self.zgrid)
        if area <= 0 or not np.isfinite(area):
            self.pdz = np.zeros_like(pdz)
        else:
            self.pdz = pdz / area

        if z_best is not None:
            self.z_best = z_best
        elif np.any(self.pdz):
            self.z_best = self.zgrid[np.argmax(self.pdz)]
        else:
            self.z_best = self.zgrid[0]

    def zbest(self):
        return self.z_best

    def variance(self, estimate):
        """Compute the pseudo standard deviation of P(z) around estimate."""
        if not np.any(self.pdz):
            return 0.0
        var = trapezoid((self.zgrid - estimate) ** 2 * self.pdz, self.zgrid)
        return np.sqrt(var)

    def approximate_gaussian(self, estimate, error=None):
        """Estimate the local Gaussian sigma around estimate."""
        def gauss(z, a, mu, sigma):
            return a * np.exp(-((z - mu) ** 2) / (2 * sigma**2))

        if error is None:
            error = self.variance(estimate)
        if error <= 0 or not np.any(self.pdz):
            return error

        mask = ((self.zgrid >= estimate - error) &
                (self.zgrid <= estimate + error))
        pdz_local = np.where(mask, self.pdz, 0.0)
        p0 = [np.max(pdz_local), estimate, error]

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, _ = curve_fit(gauss, self.zgrid, pdz_local, p0=p0)
            _, _, sigma_fit = popt
            return abs(sigma_fit)
        except (RuntimeError, ValueError):
            return error

    def number_mod(self, threshold=0.43, distance=5):
        """Count significant local maxima."""
        if not np.any(self.pdz):
            return 0
        peaks, _ = find_peaks(
            self.pdz,
            height=threshold * np.max(self.pdz),
            distance=distance,
        )
        return len(peaks)

    def peak_ratio(self):
        """Ratio of the mean P(z) to its maximum."""
        if not np.any(self.pdz):
            return 0.0
        return np.mean(self.pdz) / np.max(self.pdz)

    def tail_mass(self, estimate, sigma=None, n_window=2, good_sigma=0.01):
        """Compute probability mass outside n_window*sigma around estimate."""
        if sigma is None:
            sigma = self.approximate_gaussian(estimate)
        if sigma <= good_sigma:
            return 0.0

        bound = n_window * sigma
        mask = ((self.zgrid < estimate - bound) |
                (self.zgrid > estimate + bound))
        return trapezoid(np.where(mask, self.pdz, 0.0), self.zgrid)

def compute_pdz_score(pdz, zgrid, nb_peak_thresh=2, height_thresh=0.43,
                    tail_thresh=0.23, peak_ratio_thresh=0.1, error_thresh=0.1, z_best=None):

    """Compute a 5-bit quality score (0-31), higher = worse."""
    if pdz.sum()==0:
        return -99, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, 
    else:
        pdz_stats = PDZStats(zgrid, pdz, z_best)
        zbest = pdz_stats.zbest()
        error = pdz_stats.variance(zbest)
        sigma = pdz_stats.approximate_gaussian(zbest, error=error)
        good_sigma = zgrid[1] - zgrid[0]
        tail_mass = pdz_stats.tail_mass(zbest, sigma=sigma, good_sigma=good_sigma)
        number_mod = pdz_stats.number_mod(threshold=height_thresh, distance=5)
        peak_ratio = pdz_stats.peak_ratio()

        score = 0
        z_range = np.max(zgrid) - np.min(zgrid)
        if error > error_thresh and error < z_range / 2:
            score += 1
        if peak_ratio > peak_ratio_thresh:
            score += 2
        if tail_mass > tail_thresh:
            score += 4
        if number_mod >= nb_peak_thresh:
            score += 8
        if error >= z_range / 2:
            score += 16

        return int(score), zbest, error, peak_ratio, tail_mass, number_mod, sigma

def plot_single_pdz(pdz_path, zgrid, row=None, nb_peak_thresh=2, height_thresh=0.43,
                    tail_thresh=0.23, peak_ratio_thresh=0.1, error_thresh=0.1):
    row = row if row is not None else 0
    pdz_file = np.loadtxt(pdz_path)
    pdz_row = pdz_file[row][1:] / np.max(pdz_file[row][1:])
    # Compute metrics
    score, zbest, error, peak_ratio, tail_mass, number_mod, sigma = compute_pdz_score(pdz_row, zgrid, nb_peak_thresh, height_thresh, tail_thresh, peak_ratio_thresh, error_thresh, z_best=None)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(zgrid, pdz_row, label='PDZ')
    def gaussian(x, A, mu, s):
        return A * np.exp(-((x - mu)**2) / (2 * s**2))
    gauss_fit = gaussian(zgrid, np.max(pdz_row), zbest, sigma)
    plt.plot(zgrid, gauss_fit, label='local gaussian', ls='--')
    plt.axvline(zbest, color='purple', ls='--', label='Mode')
    plt.axhline(height_thresh*np.max(pdz_row), color='red', ls='-', label='peak_threshold')
    plt.axvspan(zbest - 2*sigma, zbest + 2*sigma, 
                color='gray', alpha=0.1, label='±2σ region')

    plt.title(f'PDZ for IDENT={pdz_file[row][0]}')
    plt.xlabel('Redshift (z)')
    plt.ylabel('P(z)')
    plt.legend()

    # Annotate metrics
    plt.text(0.02, 0.95, f"σ ≈ {sigma:.3f}\nTail mass ≈ {tail_mass:.3f}\nFlag = {score}", 
             transform=plt.gca().transAxes, fontsize=10, va='top')


    plt.tight_layout()
    plt.show()

#plot one pdz for example
plot_single_pdz(PDZ_PATH, z_grid, row = row)


# #display statistics from given catalog
# def average_stats(pdz_path, zgrid):
#     pdz_file = np.loadtxt(pdz_path)
#     score_list, zbest_list, error_list, peak_ratio_list, tail_mass_list, number_mod_list, sigma_list =[],[],[],[],[],[],[]
#     for pdz in pdz_file:
#         score, zbest, error, peak_ratio, tail_mass, number_mod, sigma = compute_pdz_score(pdz[1:], zgrid)
#         score_list.append(score)
#         zbest_list.append(zbest)
#         error_list.append(error)
#         peak_ratio_list.append(peak_ratio)
#         tail_mass_list.append(tail_mass)
#         number_mod_list.append(number_mod)
#         sigma_list.append(sigma)
    
#     print('score', np.mean(score_list),
#         '\nzbest', np.mean(zbest_list),
#         '\nerror', np.mean(error_list),
#         '\npeak_ratio', np.mean(peak_ratio_list),
#         '\ntail_mass', np.mean(tail_mass_list),
#         '\nnumber_mod',np.mean(number_mod_list),
#         '\nsigma',np.mean(sigma_list))

#     plt.figure()
#     plt.hist(np.array(sigma_list), bins=100)
#     plt.show()  

# average_stats(PDZ_PATH, z_grid)

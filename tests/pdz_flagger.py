import argparse
import os
import warnings
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit

###LePHARE sorting script###
'''
A python cript to sort good from bad zphota computed data using :
 - source indents from OUTPUT_CAT[.out]
 - BAY_ZG.prob from PDZ_OUT keyword when running z_phota
'''

def parse_args():
    parser = argparse.ArgumentParser(description="Compute PDZ quality flags and add them to a LePHARE catalog.")

    parser.add_argument("--pdz", required=True, help="Path to the PDZ file (bay_zg.prob)")

    parser.add_argument("--cat", required=True, help="Path to the LePHARE catalog")

    parser.add_argument("--zgrid", required=True, help="Redshift grid as ZSTEP,ZMIN,ZMAX (dz,zmin,zmax)")

    parser.add_argument("--cat_flagged", default="cat_flagged.out", help="Output catalog path (cat_flagged.out)")

    parser.add_argument("--max_rows", type=int, default=None, help="Maximum number of PDZ rows to process (None)")

    return parser.parse_args()

args = parse_args()
PDZ_PATH = args.pdz
CAT_PATH = args.cat
CAT_PATH_FLAGGED = args.cat_flagged
max_rows = args.max_rows
ZGRID = args.zgrid
ZSTEP, ZMIN, ZMAX = np.asarray(ZGRID.split(',')).astype(float)
z_grid = np.arange(ZMIN, ZMAX + ZSTEP, ZSTEP)

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

    def number_mod(self, threshold=0.43, distance=10):
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

def load_and_write(catalog_path, flagged_catalog_path, pdz_path, zgrid, max_rows=None):
    """
    Load PDZ and catalog files and rewrite catalog with Z_FLAG column.
    The header is preserved and extended.
    """

    #load PDZ_file
    ids = np.loadtxt(pdz_path, usecols=0, dtype=np.int64)[:max_rows]
    values = np.loadtxt(pdz_path, dtype=np.float64)[:max_rows,1:]
    # pdz_data
    #Compute PDZ_score and store it into an ident dictionary
    pdz_dict = {int(ident): compute_pdz_score(pdf, zgrid)[0] for ident, pdf in zip(ids, values)}
    
    with open(catalog_path, 'r') as fin, open(flagged_catalog_path, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                #append Z_FLAG in header
                if '# IDENT 1 Z_BEST 2' in line:
                    line = line.rstrip() + '  Z_FLAG 25\n'
                if '# IDENT  Z_BEST' in line:
                    line = line.rstrip() + '  Z_FLAG\n'
                fout.write(line)
            else:
                tokens = line.strip().split()
                ident = int(tokens[0])
                score = pdz_dict.get(ident, -1.0)  # use -1.0 if not found
                tokens.append(f"{score}")
                fout.write('          '.join(tokens) + '\n')

load_and_write(CAT_PATH, CAT_PATH_FLAGGED, PDZ_PATH, z_grid, max_rows=max_rows)

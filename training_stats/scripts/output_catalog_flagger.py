import os
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

PDZ_path = '/home/hallouin/Documents/thall_2025/lephare/docs/zphota_training/TRAINING_PDZ_3_MIN_ZG.prob'
CAT_path = '/home/hallouin/Documents/thall_2025/lephare/docs/zphota_training/TRAINING_CAT_3.out'
flagged_CAT_path = '/home/hallouin/Documents/thall_2025/lephare/docs/zphota_training/TRAINING_CAT_3_flag.out'
ZMIN, ZMAX, ZSTEP = 0, 3, 0.01
z_grid = np.arange(ZMIN, ZMAX + ZSTEP, ZSTEP)

class PDZStats:
    """
    A class to compute statistical metrics on LePHARE P(z) distributions
    to estimate quality flags for redshift reliability.
    """

    def __init__(self, zgrid, pdz, z_best=None):
        """
        Initialize with a redshift grid and a normalized P(z).

        Parameters:
        zgrid : array-like
            Array of redshift values corresponding to the PDF bins.
        pdz : array-like
            Probability distribution over zgrid (not normalized).
        z_best : float, optional
            Best-fit redshift (e.g., Z_BEST from LePHARE output).
        """
        self.zgrid = zgrid
        # area = trapezoid(pdz, zgrid)
        # if area == 0:
        #     self.pdz = np.zeros_like(pdz)
        # else:
        #     self.pdz = pdz / area
        self.pdz = pdz / trapezoid(pdz, zgrid)  #Normalize
        self.z_best = z_best if z_best is not None else self.zgrid[np.argmax(self.pdz)]
        
    
    def zbest(self):
        return self.z_best 

    def error(self):
        """Compute sigma of the P(z)"""
        var = trapezoid((self.zgrid - self.z_best) ** 2 * self.pdz, self.zgrid)
        return np.sqrt(var)

    def sigma(self):
        def gauss(z, A, mu, sigma):
            return A * np.exp(-(z - mu)**2 / (2 * sigma**2))

        mask = (self.zgrid >= self.z_best - self.error()) & (self.zgrid <= self.z_best + self.error())
        z_local = self.zgrid
        pdz_local = np.where(mask, self.pdz, 0.0)
        
        p0 = [np.max(pdz_local), self.z_best, self.error()]
        try:
            popt, _ = curve_fit(gauss, z_local, pdz_local, p0=p0)
            _, mu_fit, siga_fit = popt
            return abs(siga_fit)
        except RuntimeError:
            # Fallback in case fit fails
            return self.error()

    def number_mod(self, threshold=0.75):
        """Count significant local maxima"""
        peaks, _ = find_peaks(self.pdz, height=threshold * max(self.pdz))
        return len(peaks)

    def peak_ratio(self):
        """Ratio of max(pdz) to the mean"""
        return np.mean(self.pdz)/ max(self.pdz) 

    def tail_mass(self, n_window = 2, good_sigma =0.01):
        """
        Compute the total probability mass in the tails, outside a window around z_best.
        """
        window = self.sigma()
        mask = (self.zgrid < self.z_best - n_window * window) | (self.zgrid > self.z_best + n_window * window)
        if window <= good_sigma:
            return 0.
        else:
            return trapezoid(np.where(mask, self.pdz, 0.0), self.zgrid)

def compute_pdz_score(pdz, zgrid, nb_peak_thresh=1, height_thresh=0.5, tail_thresh=0.2, peak_ratio_thresh=0.25, error_thresh=0.2,  z_best=None):  
    """
    Returns a 4-bit integer score (0-15), higher = worse.
    impact order: number_mod, tail_mass, peak_ratio, error
    """
    pdz_stats = PDZStats(zgrid, pdz, z_best)

    score = 0

    #Compute parameters from PDZSTats
    zbest = pdz_stats.zbest()
    number_mod = pdz_stats.number_mod(threshold=height_thresh)
    tail_mass = pdz_stats.tail_mass()
    peak_ratio = pdz_stats.peak_ratio()
    error = pdz_stats.error()
    sigma = pdz_stats.sigma()

    # Bit 0: sigma
    if error > error_thresh and error < np.max(zgrid)/2:
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
    if error > np.max(zgrid)/2.5:
        score += 16



    return int(score), zbest, error, peak_ratio, tail_mass, number_mod, sigma

def load_and_write(catalog_path, flagged_catalog_path, pdz_path, zgrid):
    """
    Load PDZ and catalog files and rewrite catalog with PDZ_FLAG column.
    The header is preserved and extended.
    """

    #load PDZ_file
    pdz_data = np.loadtxt(pdz_path)
    #Compute PDZ_score and store it into an ident dictionary
    pdz_dict = {int(row[0]): compute_pdz_score(row[1:], zgrid)[0] for row in pdz_data}
    
    with open(catalog_path, 'r') as fin, open(flagged_catalog_path, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                #append PDZ_FLAG in header
                if '# IDENT  Z_BEST' in line:
                    line = line.rstrip() + '  PDZ_FLAG\n'
                fout.write(line)
            else:
                tokens = line.strip().split()
                ident = int(float(tokens[0]))
                score = pdz_dict.get(ident, -1.0)  # use -1.0 if not found
                tokens.append(f"{score}")
                fout.write('          '.join(tokens) + '\n')

# load_and_write(CAT_path, flagged_CAT_path, PDZ_path, z_grid)

def plot_single_pdz(pdz_path, zgrid, row=None, nb_peak_thresh=1, height_thresh=0.75, tail_thresh=0.5, peak_ratio_thresh=0.4, error_thresh=0.1):
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
plot_single_pdz(PDZ_path, z_grid, row=50000, nb_peak_thresh=1, height_thresh=0.75, tail_thresh=0.25, peak_ratio_thresh=0.75, error_thresh=0.125)


#display statistics from given catalog
def average_stats(pdz_path, zgrid):
    pdz_file = np.loadtxt(pdz_path)
    score_list, zbest_list, error_list, peak_ratio_list, tail_mass_list, number_mod_list, sigma_list =[],[],[],[],[],[],[]
    for pdz in pdz_file:
        score, zbest, error, peak_ratio, tail_mass, number_mod, sigma = compute_pdz_score(pdz[1:], zgrid)
        score_list.append(score)
        zbest_list.append(zbest)
        error_list.append(error)
        peak_ratio_list.append(peak_ratio)
        tail_mass_list.append(tail_mass)
        number_mod_list.append(number_mod)
        sigma_list.append(sigma)
    
    print('score', np.mean(score_list),
        '\nzbest', np.mean(zbest_list),
        '\nerror', np.mean(error_list),
        '\npeak_ratio', np.mean(peak_ratio_list),
        '\ntail_mass', np.mean(tail_mass_list),
        '\nnumber_mod',np.mean(number_mod_list),
        '\nsigma',np.mean(sigma_list))

    plt.figure()
    plt.hist(np.array(sigma_list), bins=100)
    plt.show()  

# average_stats(PDZ_path, z_grid)
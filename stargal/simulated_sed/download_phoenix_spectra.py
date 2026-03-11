#!/usr/bin/env python3
"""
Download PHOENIX HiRes spectra via FTP using the common WAVE file,
convert them into normalized .sed files (at 8004 Å),
and optionally resample them on a custom wavelength grid.

Usage:
- Edit the parameter grids (teff_list, logg_list, feh_list)
- Optionally define lmin, lmax, and dl to crop and resample the spectra
"""

import os
import sys
import io
import ftplib
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
from tqdm import tqdm

# =====================
# USER PARAMETERS
# =====================
OUTDIR = "phoenix_spectra"
os.makedirs(OUTDIR, exist_ok=True)
# sys.appends()

# Grid to download
teff_list = np.concatenate((np.arange(2400, 7000, 200),np.arange(7000, 12200, 200)))
logg_list = np.arange(-0.5, 6.5, 1)
feh_list  = [0.0]

# Optional resampling parameters
LMIN = 1145.0   # Å, set to None to keep full range
LMAX = 25005.0  # Å
DL   = 5.0      # Å, set to None to keep original sampling

# FTP settings
FTP_HOST = "phoenix.astro.physik.uni-goettingen.de"
FTP_BASEDIR = "HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0"
WAVE_FN = "HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"

# =====================
# FTP UTILS
# =====================
def ftp_connect():
    ftp = ftplib.FTP(FTP_HOST, timeout=60)
    ftp.login()  # anonymous
    return ftp

def ftp_download_to_memory(ftp, remote_path):
    bio = io.BytesIO()
    ftp.retrbinary(f"RETR {remote_path}", bio.write)
    bio.seek(0)
    return bio

def ftp_file_exists(ftp, remote_path):
    dirname, fname = os.path.split(remote_path)
    try:
        entries = ftp.nlst(dirname)
    except ftplib.error_perm:
        try:
            ftp.cwd(dirname)
            entries = ftp.nlst()
            ftp.cwd('/')
        except Exception:
            return False
    return fname in [os.path.basename(e) for e in entries]

# =====================
# LOAD OR DOWNLOAD WAVE FILE
# =====================
def get_wave_local(ftp, local_cache_dir=".", force_download=False):
    """Retrieve the common WAVE file (Å) locally, downloading if needed."""
    local_wave_path = os.path.join(local_cache_dir, os.path.basename(WAVE_FN))
    if (not force_download) and os.path.exists(local_wave_path):
        print(f"-> Using local WAVE file {local_wave_path}")
        with fits.open(local_wave_path) as hd:
            wl = hd[0].data.astype(float)
        return wl

    print(f"-> Downloading WAVE file: {WAVE_FN}")
    if not ftp_file_exists(ftp, WAVE_FN):
        raise FileNotFoundError(f"WAVE file not found on FTP: {WAVE_FN}")

    buf = ftp_download_to_memory(ftp, WAVE_FN)
    os.makedirs(os.path.dirname(local_wave_path), exist_ok=True)
    with open(local_wave_path, "wb") as f:
        f.write(buf.getbuffer())
        print(f"WAVE file saved locally: {local_wave_path}")

    with fits.open(local_wave_path) as hd:
        wl = hd[0].data.astype(float)
    print(f"-> WAVE loaded ({len(wl)} points).")
    return wl

# =====================
# RESAMPLING FUNCTION
# =====================

import numpy as np
from scipy.interpolate import interp1d

def resample_spectrum(wl, flux, lmin=None, lmax=None, dl=None, tol=1e-5):
    """
    Crop and/or resample a spectrum to a uniform wavelength grid.
    Exact wavelength matches (within tol) use the original flux value.
    Otherwise, flux is linearly interpolated.

    Parameters
    ----------
    wl : array-like
        Original wavelength array (Å)
    flux : array-like
        Original flux array
    lmin, lmax : float, optional
        Min/max wavelength limits for cropping
    dl : float, optional
        Step size for the new wavelength grid (Å)
    tol : float, optional
        Tolerance (Å) to consider a wavelength as matching an original one

    Returns
    -------
    wl_new, flux_new : np.ndarray
        Resampled wavelength and flux arrays
    """
    # Sort arrays just in case
    order = np.argsort(wl)
    wl = np.array(wl)[order]
    flux = np.array(flux)[order]

    # Crop to limits
    if lmin is not None:
        mask = wl >= lmin
        wl, flux = wl[mask], flux[mask]
    if lmax is not None:
        mask = wl <= lmax
        wl, flux = wl[mask], flux[mask]

    # If no resampling requested, just return cropped data
    if dl is None:
        return wl, flux

    # Build uniform wavelength grid
    wl_new = np.arange(wl[0], wl[-1] + dl/2, dl)

    # Prepare interpolation function for missing points
    interp_func = interp1d(
        wl, flux,
        kind="linear",
        bounds_error=False,
        fill_value=np.nan
    )

    # Allocate flux array
    flux_new = np.empty_like(wl_new)

    # Efficient matching loop
    j = 0
    n = len(wl)
    for i, lam in enumerate(wl_new):
        # Advance j until wl[j] >= lam - tol
        while j < n - 1 and wl[j] < lam - tol:
            j += 1
        # Check for direct match
        if abs(wl[j] - lam) <= tol:
            flux_new[i] = flux[j]
        else:
            # Interpolate
            flux_new[i] = interp_func(lam)

    return wl_new, flux_new


# =====================
# EXPORT
# =====================
def export_sed_file(outdir, teff, logg, feh, wl, flux, norm_wave=5540.0):
    """Normalize at norm_wave and save as .sed text file."""
    order = np.argsort(wl)
    wl, flux = wl[order], flux[order]

    idx = np.argmin(np.abs(wl - norm_wave))
    if flux[idx] <= 0 or np.isnan(flux[idx]):
        good = np.where((~np.isnan(flux)) & (flux > 0))[0]
        if len(good) == 0:
            raise ValueError("Flux array entirely zero or NaN, cannot normalize.")
        idx = good[np.argmin(np.abs(wl[good] - norm_wave))]
    flux /= flux[idx]

    fname = f"Teff{int(teff):04d}_logg{logg:.1f}_feh{feh:+.1f}.sed"
    outpath = os.path.join(outdir, fname)
    with open(outpath, "w") as f:
        f.write(f"# SED PHOENIX Teff={teff} logg={logg:.1f} [Fe/H]={feh}\n")
        f.write("# wl(Å) flux (normalized at 8004 Å)\n")
        for lam, fl in zip(wl, flux):
            f.write(f"{lam:10.6f} {fl:12.6e}\n")
    return outpath

# =====================
# MAIN LOOP
# =====================
if __name__ == "__main__":
    ftp = ftp_connect()
    print("✅ Connected to PHOENIX FTP server.")

    try:
        wl = get_wave_local(ftp, local_cache_dir=OUTDIR, force_download=False)
    except Exception as e:
        ftp.quit()
        raise

    combos = [(t, g, z) for z in feh_list for t in teff_list for g in logg_list]

    for teff, logg, feh in tqdm(combos, desc="Spectra", unit="spec"):
        fn = f"lte{int(teff):05d}-{logg:04.2f}-{feh:.1f}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
        remote_path = f"{FTP_BASEDIR}/{fn}"
        print(f"-> {remote_path}")

        if not ftp_file_exists(ftp, remote_path):
            print(f"⚠️  Not found (skipped): {fn}")
            continue

        try:
            buf = ftp_download_to_memory(ftp, remote_path)
        except Exception as e:
            print(f"⚠️  FTP error for {fn}: {e}")
            continue

        try:
            with fits.open(buf) as hd:
                flux = hd[0].data.astype(float)
        except Exception as e:
            print(f"⚠️  FITS read error for {fn}: {e}")
            continue
            
        # Apply wavelength cropping and resampling
        wl_proc, flux_proc = resample_spectrum(wl, flux, LMIN, LMAX, DL)

        try:
            outp = export_sed_file(OUTDIR, teff, logg, feh, wl_proc, flux_proc, norm_wave=5540.0)
        except Exception as e:
            print(f"⚠️  Export error for {fn}: {e}")
            continue

    ftp.quit()
    print("🌟 Done.")

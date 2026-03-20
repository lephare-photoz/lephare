# Dowloading theoretical star spectra 
# and convert them to .sed as intended for lephare processing

import requests
import xml.etree.ElementTree as ET
import numpy as np
from scipy.interpolate import interp1d
import os
import sys



def download_bt_spectra(
    # === DATABASE & DL DIRECTORY ===
    MODEL = "bt-nextgen-agss2009",
    BASE_URL = None,
    OUTPUT_DIR = os.path.abspath(os.path.join(os.getcwd(), 'stargal/simulated_sed/bt_spectra')), # replace by you output directory

    # === PHYSICAL VALUES SPECTRUM GRID  ===
    teff_values = [20000],
    logg_values = [4.0],
    metallicity_values = [0.0],

    # Optional resampling parameters
    LMIN = 1145.0,   # Å, set to None to keep full range
    LMAX = 25005.0,  # Å
    DL   = 5.0,     # Å, set to None to keep original sampling
    wl_norm = 10000,
    make_sed_list = True,
    overwrite_seds = False,
    list_name = "bt_star_sed_full",
    ):

    if BASE_URL is None:
        BASE_URL = f"https://svo2.cab.inta-csic.es/theory/newov2/ssap.php?model={MODEL}"

    def makeititerable(x):
        if isinstance(x, (list, tuple, np.ndarray)):
            return list(x)
        else:
            return [x]

    teff_values = makeititerable(teff_values)
    logg_values = makeititerable(logg_values)
    metallicity_values = makeititerable(metallicity_values)

    # === DOWNLOAD AND READ VOTABLE ===
    print("Downloading main VOTABLE...")
    resp = requests.get(BASE_URL, verify=False)
    resp.raise_for_status()

    votable_content = resp.text
    root = ET.fromstring(votable_content)

    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}")[0] + "}"

    rows = root.findall(f".//{ns}TR")
    print(f"{len(rows)} <TR> lines found in VOTABLE.")

    # === RETRIEVE INDIVIDUAL SPECTRUM LINKS ===
    spectra_dict = {}  # key = (teff, logg, feh), value = (alpha, url)

    for row in rows:
        cells = [td.text.strip() if td.text else "" for td in row.findall(f"{ns}TD")]
        if len(cells) < 11:
            continue

        try:
            teff = float(cells[0])
            logg = float(cells[1])
            feh = float(cells[2])
            alpha = float(cells[3]) if cells[3] not in ("", "null", "NaN") else 0.0
            url = cells[10].strip()
        except Exception:
            continue

        # Filtrage selon la grille cible
        if (round(teff) in teff_values
            and any(abs(logg - g) < 1e-3 for g in logg_values)
            and any(abs(feh - z) < 1e-3 for z in metallicity_values)):

            key = (round(teff), round(logg, 2), round(feh, 2))
            # Si ce triplet n’existe pas encore ou si ce nouveau alpha est plus proche de 0 → on remplace
            if key not in spectra_dict or abs(alpha) < abs(spectra_dict[key][0]):
                spectra_dict[key] = (alpha, url)

    # Conversion finale en liste
    spectra_links = [(teff, logg, feh, url) for (teff, logg, feh), (alpha, url) in spectra_dict.items()]

    print(f"{len(spectra_links)} spectra found matching with the chosen grid.")

    # === FUNCTION TO LOWER THE SAMPLING OF A SED FLUX = f(WL) ===

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

    # === DOWNLOAD AND CONVERT TO SED ===
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for i, (teff, logg, feh, url) in enumerate(spectra_links, 1):
        fname = f"Teff{int(teff)}_logg{logg:.1f}_FeH{feh:.1f}.sed"
        if os.path.exists(os.path.join(OUTPUT_DIR, fname)) and overwrite_seds == False:
            print(f"T_eff={teff}, logg={logg}, Fe/H={feh} already exists. Pass.")
            continue
        spec = requests.get(url + "&format=ascii", verify=False)
        if spec.status_code != 200:
            print(f"ERROR: {url} not found")
            continue

        print(f"Downloading [{i}/{len(spectra_links)}]... T_eff={teff}, logg={logg}, Fe/H={feh}")

        lines = [l.strip() for l in spec.text.splitlines() if l.strip() and not l.startswith("#")]
        data = np.loadtxt(lines)

        # Norm to wl_norm (in A)
        wl = data[:, 0]
        flux = data[:, 1]
        flux = flux / np.interp(wl_norm, wl, flux)

        # Resample
        wl_proc, flux_proc = resample_spectrum(wl, flux, LMIN, LMAX, DL)

        # Export
        print(os.path.join(OUTPUT_DIR, fname))
        with open(os.path.join(OUTPUT_DIR, fname), "w") as f:
            f.write(f"#SED {MODEL}\n#wl(AA) flux (normed at {wl_norm}AA)\n")
            for w, fl in zip(wl_proc, flux_proc):
                f.write(f"{w:.3f} {fl:.7f}\n")

    print("All available spectrum have been downloaded and formated.")

    # Make sed.list, mandatory for lephare if requested
    if make_sed_list == True:
        SED_list_path = os.path.join(OUTPUT_DIR, f'{list_name}.list')
        
        bt_sed_list = sorted([
            f for f in os.listdir(OUTPUT_DIR)
            if os.path.isfile(os.path.join(OUTPUT_DIR, f)) and f.endswith(".sed")
        ])


        #save to file.list because we must do it like this
        with open(SED_list_path, "w") as f:
            for sed_file in bt_sed_list:
                f.write(OUTPUT_DIR + "/" + sed_file + "\n")


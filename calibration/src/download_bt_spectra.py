# Dowloading theoretical star spectra 
# and convert them to .sed as intended for lephare processing

import requests
import xml.etree.ElementTree as ET
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import sys
import warnings
from io import StringIO
warnings.filterwarnings("ignore")

def download_bt_spectra(
    # === DATABASE & DL DIRECTORY ===
    MODEL = "bt-nextgen-agss2009",
    BASE_URL = None,
    OUTPUT_DIR = os.path.abspath(os.path.join(os.getcwd(), 'stargal/simulated_sed/bt_spectra')),

    # === PHYSICAL VALUES SPECTRUM GRID  ===
    teff_values = [20000],
    logg_values = [4.0],
    metallicity_values = [0.0],

    # Optional resampling parameters
    LMIN = 1145.0,
    LMAX = 25005.0,
    DL   = 5.0,
    wl_norm = 10000,
    make_sed_list = True,
    overwrite_seds = False,
    list_name = "bt_star_sed_full",
    list_folder = None,

    # Parallelization
    max_workers = 8
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

    teff_set = set(teff_values)
    logg_set = set(round(g, 3) for g in logg_values)
    feh_set  = set(round(z, 3) for z in metallicity_values)

    # === DOWNLOAD AND READ VOTABLE ===
    print("Downloading main VOTABLE...")
    session = requests.Session()
    resp = session.get(BASE_URL, verify=False)
    resp.raise_for_status()

    root = ET.fromstring(resp.text)

    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}")[0] + "}"

    rows = root.findall(f".//{ns}TR")
    print(f"{len(rows)} <TR> lines found in VOTABLE.")

    # === RETRIEVE INDIVIDUAL SPECTRUM LINKS ===
    spectra_dict = {}

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

        if (round(teff) in teff_set
            and round(logg, 3) in logg_set
            and round(feh, 3) in feh_set):

            key = (round(teff), round(logg, 2), round(feh, 2))
            if key not in spectra_dict or abs(alpha) < abs(spectra_dict[key][0]):
                spectra_dict[key] = (alpha, url)

    spectra_links = [(teff, logg, feh, url) for (teff, logg, feh), (alpha, url) in spectra_dict.items()]

    print(f"{len(spectra_links)} spectra found matching with the chosen grid.")

    # === FUNCTION TO LOWER THE SAMPLING OF A SED FLUX = f(WL) ===
    def resample_spectrum(wl, flux, lmin=None, lmax=None, dl=None):
        wl = np.asarray(wl)
        flux = np.asarray(flux)

        if lmin is not None:
            mask = wl >= lmin
            wl, flux = wl[mask], flux[mask]
        if lmax is not None:
            mask = wl <= lmax
            wl, flux = wl[mask], flux[mask]

        if dl is None:
            return wl, flux

        wl_new = np.arange(wl[0], wl[-1] + dl/2, dl)
        flux_new = np.interp(wl_new, wl, flux)

        return wl_new, flux_new

    # === PROCESS ONE SPECTRUM (parallel worker) ===
    def process_spectrum(args):
        teff, logg, feh, url = args

        fname = f"Teff{int(teff)}_logg{logg:.1f}_FeH{feh:.1f}.sed"
        outpath = os.path.join(OUTPUT_DIR, fname)

        if os.path.exists(outpath) and not overwrite_seds:
            return f"Skipping existing: {fname}"

        try:
            r = session.get(url + "&format=ascii", verify=False, timeout=30)
            r.raise_for_status()

            data = np.loadtxt(StringIO(r.text))

            wl = data[:, 0]
            flux = data[:, 1]

            flux = flux / np.interp(wl_norm, wl, flux)

            wl_proc, flux_proc = resample_spectrum(wl, flux, LMIN, LMAX, DL)

            np.savetxt(
                outpath,
                np.column_stack([wl_proc, flux_proc]),
                fmt="%.3f %.7f",
                header=f"#SED {MODEL}\n#wl(AA) flux (normed at {wl_norm}AA)",
                comments=""
            )

            return f"Done: T_eff={teff}, logg={logg}, Fe/H={feh}"

        except Exception as e:
            return f"ERROR: {url} ({e})"

    # === DOWNLOAD AND CONVERT TO SED (PARALLEL) ===
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fname_list = [os.path.join(OUTPUT_DIR,
                    f"Teff{int(teff)}_logg{logg:.1f}_FeH{feh:.1f}.sed")
        for (teff, logg, feh, _) in spectra_links
    ]
    print(f"Starting parallel download with {max_workers} workers...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_spectrum, args) for args in spectra_links]

        for i, future in enumerate(as_completed(futures), 1):
            print(f"[{i}/{len(futures)}] {future.result()}")

    print("All available spectrum have been downloaded and formated.")

    # === MAKE SED LIST ===
    if make_sed_list:
        if list_folder is not None:
            os.makedirs(list_folder, exist_ok=True)
            SED_list_path = os.path.join(list_folder, f'{list_name}.list')
        else:
            SED_list_path = os.path.join(OUTPUT_DIR, f'{list_name}.list')

        bt_sed_list = sorted([
            f for f in os.listdir(OUTPUT_DIR)
            if (os.path.isfile(os.path.join(OUTPUT_DIR, f))
                and f.endswith(".sed")
                and os.path.join(OUTPUT_DIR, f) in fname_list)
        ])

        with open(SED_list_path, "w") as f:
            for sed_file in bt_sed_list:
                f.write(OUTPUT_DIR + "/" + sed_file + "\n")

        print('sed list written to ', SED_list_path)
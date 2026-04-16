import yaml
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from astropy.io import fits
from astropy.table import Table
from astropy.io.votable import parse_single_table

def maglib_to_pandas(filename, filters=['u','g','r','i','z','y']):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Identifier la ligne des noms de colonnes
    for i, line in enumerate(lines):
        if not line.startswith('#'):
            header_line_index = i - 1
            break
    
    # Extraire les noms de colonnes
    col_names = lines[header_line_index].strip().lstrip('#').split()
    
    # Détecter les colonnes multi-valeurs
    multi_cols = [c for c in col_names if '[N_filt]' in c]
    
    # Colonnes fixes
    fixed_cols = [c for c in col_names if '[N_filt]' not in c]
    
    data = []
    for line in lines[header_line_index+1:]:
        if line.strip() == '':
            continue
        row = [float(x) for x in line.split()]
        
        # Séparer les valeurs fixes et les multi
        row_fixed = row[:len(fixed_cols)]
        row_multi_start = len(fixed_cols)
        row_multi_values = row[row_multi_start:]
        
        data.append(row_fixed + row_multi_values)
    
    # Créer le DataFrame initial avec uniquement les colonnes fixes + N_filt
    df = pd.DataFrame([d[:len(fixed_cols)] for d in data], columns=fixed_cols)

    # Déplier chaque colonne multi-valeurs
    for col in multi_cols:
        col_idx = col_names.index(col)
        expanded_values = []
        for i, row in enumerate(data):
            Nf = int(df.loc[i, 'N_filt'])
            # Extraire les Nf valeurs correspondantes à cette colonne
            values = row[col_idx: col_idx+Nf]
            # Compléter avec None si Nf < nombre de filtres
            values += [None]*(len(filters)-Nf)
            expanded_values.append(values)
        
        # Ajouter ces colonnes au DataFrame
        new_col_names = [f"{col.split('[')[0]}_{f}" for f in filters]
        df[new_col_names] = pd.DataFrame(expanded_values)
    
    return df

###
def load_hp_file(path, nside_dg=None, nested=True):
    try:
        file_map = hp.read_map(path)
        nside = hp.get_nside(file_map)
    except AttributeError:
        file_map = healsparse.HealSparseMap.read(path)
        nside = file_map.nside_sparse
        npix = hp.nside2npix(nside)

        if nested == True:
            hp_map_nested = np.full(npix, hp.UNSEEN, dtype=float)
            pix = file_map.valid_pixels
            hp_map_nested[pix] = file_map[pix]

            hp_map = hp.reorder(hp_map_nested, n2r=True)
    
    if nside_dg is not None and nside>nside_dg:
        hp_map = hp.ud_grade(file_map, nside_out=nside_dg,
        order_in="RING", order_out="RING", power=0)
        nside = nside_dg
    
    return hp_map, nside


def table_from_file(path, file_type="fits", hdu=1, output="astropy", **read_kwargs):
    """
    Parameters
    ----------
    path : str
        Path to the input file.
    file_type : {'fits', 'csv', 'txt', 'votable'}
        Type of the input file.
    hdu : int or str, optional
        HDU index or name (only used for FITS, default: 1).
    output : {'astropy', 'pandas'}
        Output format.
    **read_kwargs :
        Additional keyword arguments passed to Table.read()
        (useful for csv/txt: delimiter, format, etc.)
    Returns
    -------
    table : astropy.table.Table or pandas.DataFrame
    """

    file_type = file_type.lower()
    if file_type == "fits":
        with fits.open(path, memmap=True) as hdul:
            hdu_obj = hdul[hdu]
            xtension = hdu_obj.header.get("XTENSION", "").upper()

            if xtension in ("BINTABLE", "TABLE"):
                data = Table(hdu_obj.data)
            else:
                raise ValueError(f"HDU {hdu} does not contain a FITS table (XTENSION={xtension})")

    elif file_type in ("csv", "txt"):
        data = Table.read(path, **read_kwargs)

    elif file_type == "votable":
        votable = parse_single_table(path)
        data = votable.to_table()

    else:
        raise ValueError("file_type must be one of {'fits', 'csv', 'txt', 'votable'}")

    if output == "astropy":
        return data
    elif output == "pandas":
        return data.to_pandas()
    else:
        raise ValueError("output must be 'astropy' or 'pandas'")


def _flatten_dict(d, parent_key="", sep="_", ignore_keys=None):
    items = {}
    ignore_keys = ignore_keys or set()
    for k, v in d.items():
        if k in ignore_keys:
            continue
        if v is None:
            continue
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.update(_flatten_dict(v, parent_key=new_key, sep=sep, ignore_keys=ignore_keys))
        else:
            items[new_key] = v
    return items


def yaml_to_pandas(yaml_path, ignore_columns=None, sep="_", object_column_name="object_name"):
    ignore_columns = set(ignore_columns or [])
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)
    if isinstance(data, dict) and len(data) == 1:
        root = next(iter(data.values()))
        if isinstance(root, dict):
            data = root
    rows = []
    for object_name, entry in data.items():
        row = {f"{object_column_name}": object_name}
        if isinstance(entry, dict):
            row.update(_flatten_dict(entry,
                    sep=sep,ignore_keys=ignore_columns))
        rows.append(row)
    return pd.DataFrame(rows)


def join_tables(base_df, tables, how="left", ignore=None, put_match_first=True, drop_unmatched=True, indicator=False):
    """
    Join multiple dataframes onto a base dataframe using specified column mappings,
    with options to reorder columns and remove unmatched rows.

    Parameters
    ----------
    base_df : pd.DataFrame
        The initial dataframe to join others onto.
    tables : tuple or list of tuples
        Each tuple must be (other_df, base_col, other_col).
    how : str
        Join type ("left", "right", "inner", "outer").
    ignore : list[str]
        List of columns to drop.
    put_match_first : bool
        If True, the match key (base_col) is moved to the first column.
    drop_unmatched : bool
        If True, rows where the join key has no match are dropped
        (forces an inner-join behavior even if how != "inner").
    indicator : bool
        If True, return the base_df wtih a joint column present filled with 0 for sources
        present only in the base_df, and 1 if present in both table, according to the matching column.
    
    Returns
    -------
    pd.DataFrame
    """
    # Normalize single join instruction to list
    if isinstance(tables, tuple) and len(tables) == 3:
        tables = [tables]

    result = base_df.copy()
    match_columns = []

    for other_df, base_col, other_col in tables:
        result = result.merge(
            other_df,
            how=how,
            left_on=base_col,
            right_on=other_col,
            indicator=indicator
        )
        match_columns.append(base_col)
    # Drop unmatched rows (forces inner join behavior)
    if drop_unmatched:
        for col in match_columns:
            result = result[result[col].notna()]
    
    if indicator:
        result = result[(result['_merge'] == 'left_only') | (result['_merge'] == 'both')]
        result["present"] = (result['_merge'] == "both").astype(int)
        result = result.drop(columns=["_merge"])
    # Drop ignored columns
    if ignore:
        result = result.drop(columns=[c for c in ignore if c in result.columns])

    # Move match columns to front
    if put_match_first:
        # Keep order if multiple match columns
        front_cols = [c for c in match_columns if c in result.columns]
        other_cols = [c for c in result.columns if c not in front_cols]
        result = result[front_cols + other_cols]

    return result

class SED_GRID:
    '''
    For a SED set dependent on physical parameters, construct a coherent grid.
    '''
    def __init__(self, list_path):
        """Load SED list."""
        self.list_path = list_path
        self.filepaths = []
        self.sed_list = []
        self.sed_grid = None

        with open(list_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                filepath = line.split()[0]
                self.filepaths.append(str(filepath))

                sp_type = os.path.basename(filepath).replace('.sed', "")
                self.sed_list.append(sp_type)
                    
    def build(self, pd_table=False, sort_values=False):
        """
        Build the SED grid [Id, Teff, logg, FeH].
        Parameters
        ----------
        pd_table : bool, optional
            If True, return a pandas DataFrame instead of a numpy array.
        Returns
        -------
        np.ndarray or pandas.DataFrame
        """

        pattern = r"Teff([-\d\.]+)_logg([-\d\.]+)_FeH([-\d\.]+)"
        model_ids = []
        teff = []
        logg = []
        feh = []
        paths = []

        for Id, (s, path) in enumerate(zip(self.sed_list, self.filepaths)):
            match = re.search(pattern, s)
            if match:
                Teff, Logg, FeH = match.groups()
                model_ids.append(float(Id))
                teff.append(float(Teff))
                logg.append(float(Logg))
                feh.append(float(FeH))
                paths.append(path)
                
        sed_grid = np.column_stack((model_ids, teff, logg, feh))

        # Return pandas DataFrame
        if pd_table == True:
            self.sed_grid =  sed_grid
            sed_df = pd.DataFrame(
                np.concatenate([np.array(sed_grid), np.array(self.filepaths).reshape(-1, 1)], axis=1),
                columns=["model", "teff", "logg", "feh", "path"])
            sed_df = sed_df.astype({"model": float, "teff": float, "logg": float, "feh": float})
            sed_df["model"] = sed_df["model"].astype(int)
            if sort_values == True:
                sed_df = sed_df.sort_values(["teff", "logg", "feh"]).reset_index(drop=True)
            return sed_df
        # Return numpy array
        else:
            self.sed_grid = sed_grid
            return sed_grid

    def plot(self, cmap='gnuplot'):
        """Plot 3D SED grid with color depending on FeH."""
        if self.sed_grid is None:
            raise ValueError("Call make_sed_grid before plot.")

        Teff, logg, FeH = self.sed_grid[:, 1], self.sed_grid[:, 2], self.sed_grid[:, 3]

        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(projection='3d')

        # Couleur selon la métallicité
        sc = ax.scatter(
            Teff, logg, FeH,
            c=FeH, cmap=cmap,
            alpha=0.8, s=10, edgecolor="none"
        )

        # Colorbar
        cbar = plt.colorbar(sc, ax=ax, pad=0.1, shrink=0.8)
        cbar.set_label('[Fe/H]', rotation=270, labelpad=15)

        ax.set_xlabel('Teff [K]')
        ax.set_ylabel('log(g)')
        ax.set_zlabel('[Fe/H]')
        ax.set_title("SED Grid Colored by Metallicity")

        # Optionnel : meilleure orientation par défaut
        ax.view_init(elev=20, azim=-45)

        plt.tight_layout()
        plt.show()


class INTERPOLATE_SEDS:
    """
    Lineary interpolate a list of sed on one axis from a list 
    """
    @staticmethod
    def _parse_sed(sed):
        arr = np.array(sed)

        if arr.ndim == 2 and arr.shape[1] == 2:
            lam = arr[:, 0]
            flux = arr[:, 1]
            return lam, flux

        elif arr.ndim == 2 and arr.shape[0] == 2:
            lam = arr[0]
            flux = arr[1]
            return lam, flux

        else:
            raise ValueError("Invalid SED format. Expected [[λ,f],...] or [[λ,...],[f,...]].")

    @staticmethod
    def interpolate_seds(specA, specB, n):
        lambda_A, flux_A = INTERPOLATE_SEDS._parse_sed(specA)
        lambda_B, flux_B = INTERPOLATE_SEDS._parse_sed(specB)

        if not np.allclose(lambda_A, lambda_B):
            raise ValueError("Wavelength grids must be identical for interpolation.")

        lambdas = lambda_A
        ts = np.linspace(0, 1, n + 2)[1:-1]  # exclude endpoints

        return [np.column_stack((lambdas, (1 - t) * flux_A + t * flux_B))
            for t in ts]

    @staticmethod
    def interpolate_from_grid(df, axis, fixed_params,
                            path_col="path", n_interp=1, write=False, output_dir=None, 
                            new_list_path=None, list_path=None, return_results=False):
        """
        Interpolate SEDs along a chosen axis.

        Parameters
        ----------
        df : pandas.DataFrame
        axis : str
            Axis along which to interpolate
        fixed_params : dict or None
            Example: {'logg': 3.0, 'feh': 0.0}
        path_col : str
            Column containing SED file paths
        n_interp : int
            Number of interpolated SEDs between each pair
        write : bool
            If True, write interpolated SEDs to disk
        output_dir : str
            Directory where files are written

        Returns
        -------
        list of dict or None
        """

        df_work = df.copy()

        # --- Apply filtering on other parameters ---
        if fixed_params is not None:
            for key, value in fixed_params.items():
                df_work = df_work[df_work[key] == value]

        # --- Sort along interpolation axis ---
        df_work = df_work.sort_values(axis).reset_index(drop=True)
        results = []
        out_path_list = []
        # --- Loop over consecutive SED pairs ---
        for i in range(len(df_work) - 1):
            rowA = df_work.iloc[i]
            rowB = df_work.iloc[i + 1]

            specA = np.loadtxt(rowA[path_col])
            specB = np.loadtxt(rowB[path_col])

            lam, fluxA = INTERPOLATE_SEDS._parse_sed(specA)
            _, fluxB = INTERPOLATE_SEDS._parse_sed(specB)

            # --- Interpolation ---
            for k in range(1, n_interp + 1):
                t = k / (n_interp + 1)
                flux_interp = (1 - t) * fluxA + t * fluxB
                sed_interp = np.column_stack((lam, flux_interp))
                # Interpolated physical parameters
                new_params = {}

                for col in ["teff", "logg", "feh"]:
                    valA = rowA[col]
                    valB = rowB[col]

                    if col == axis:
                        new_params[col] = (1 - t) * valA + t * valB
                    else:
                        new_params[col] = valA

                # --- Output handling ---
                if write:
                    if output_dir is None:
                        raise ValueError("output_dir must be provided if write=True")

                    filename = (
                        f"Teff{new_params['teff']:.0f}_"
                        f"logg{new_params['logg']:.2f}_"
                        f"FeH{new_params['feh']:.2f}.sed"
                    )
                    os.makedirs(output_dir, exist_ok=True)
                    out_path = os.path.join(output_dir, filename)
                    out_path_list.append(out_path)
                    np.savetxt(out_path, sed_interp)

                if return_results==True:
                    results.append({
                        "sed": sed_interp,
                        "teff": new_params["teff"],
                        "logg": new_params["logg"],
                        "feh": new_params["feh"]
                    })

        # --- remake sed list ---
        if new_list_path and write:
            with open(new_list_path, "w") as fout:
                # Case 1: append to existing list
                if list_path is not None:
                    with open(list_path, "r") as fin:
                        for line in fin:
                            fout.write(line.rstrip() + "\n")
                # Case 2: always append new interpolated SEDs
                for out in out_path_list:
                    fout.write(out + "\n")
            print("SED list written to", new_list_path)

        if return_results==True:
            return results
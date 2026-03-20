import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml
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


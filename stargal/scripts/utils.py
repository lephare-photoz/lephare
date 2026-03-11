import pandas as pd
from astropy.io import fits
from astropy.table import Table

def lephare_to_pandas(CAT_OUT):
    """
    Convert the output from lephare to a pandas data frame.
    Parameters
    ----------
    CAT_OUT: string.
        Output catalog from lephare.
    
    Returns
    -------
    Pandas data frame.
        Re-organized lephare output.
    """

    with open(CAT_OUT, "r") as f:
        lines = f.readlines()
        header_line = None
        for line in lines:
            if line.startswith("# IDENT  Z_BEST"): #line used for the header, always starts like this
                header_line = line
                break

    #add header to column names
    if header_line:
        column_names = header_line.strip("#").strip().split()
    
    zphota = pd.read_csv(CAT_OUT, sep=r'\s+', comment="#", header=None, names=column_names) #zphota dataframe
    return zphota


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
        for c in ignore: 
            if c in result.columns:
                result = result.drop(columns=c)

    # Move match columns to front
    if put_match_first:
        # Keep order if multiple match columns
        front_cols = [c for c in match_columns if c in result.columns]
        other_cols = [c for c in result.columns if c not in front_cols]
        result = result[front_cols + other_cols]
    result = result.loc[:,~result.columns.duplicated()].copy()
    return result

import os
import numpy as np
import pandas as pd
from astropy.io import fits as pf


def apply_context_function(data_array, n_filters, apply_context):
    """
    Handle the context column according to the 'apply_context' flag.
    
    Parameters
    ###
    data_array : np.ndarray
        Input data array (already structured with mags, errors, context, zspec).
    n_filters : int
        Number of photometric filters (used to identify magnitude columns).
    apply_context : str
        "no"   -> do nothing
        "yes"  -> compute context bitmask based on valid magnitudes (<98)
        "null" -> set all context values to 0
    
    Returns
    -------
    np.ndarray
        Modified array with updated or added context values.
    """
    corr_data = np.copy(data_array)

    if apply_context.lower() == "no":
        return corr_data

    elif apply_context.lower() == "null":
        # Set context column (second to last) to zeros
        corr_data[:, -2] = 0
        return corr_data

    elif apply_context.lower() == "yes":
        # Compute context bitmask per source
        for k, source in enumerate(corr_data):
            mags = source[1:1 + n_filters]  # assume mags start at index 1
            context_val = 0
            for i, mag in enumerate(mags):
                if mag < 98:  # 98 or 99 used for non-detections
                    context_val += 2 ** i
            corr_data[k][-2] = context_val  # second to last column = context
        return corr_data

    else:
        raise ValueError("apply_context must be one of: 'no', 'yes', or 'null'")



def format_to_lephareinput(CAT_IN,
                           CAT_OUT,
                           input_columns,
                           n_filters,
                           CAT_TYPE='short',
                           simple_convert=False,
                           shuffle=False,
                           apply_context='no',
                           max_rows=None,
                           error_value_state='default',
                           facticious_specz=None):
    """
    Format input catalog into a LePhare-compatible .dat file.
    
    Parameters
    ----------
    CAT_IN : str
        Path to input catalog (FITS, CSV, or whitespace-separated text).
    CAT_OUT : str
        Path to output .dat file.
    input_columns : list of str
        Ordered list of columns to extract from the input catalog.
    n_filters : int
        Number of photometric filters.
    CAT_TYPE : str, optional
        'short' = Id {mags} 
        'long'  = Id {mags, errs, context, zspec}
    simple_convert : bool, optional
        If True, simply convert the table into a .dat without applying formatting logic.
    shuffle : bool, optional
        If True, shuffle the rows of the catalog.
    apply_context : str, optional
        "no" (do nothing), "yes" (compute context), or "null" (set context to 0).
    max_rows : int, optional
        Maximum number of rows to process (useful for testing).
    error_value_state : str, optional
        "default" → replace missing values with 99.0
        "delete"  → drop rows with any missing value
    facticious_specz : float, optional
        If CAT_TYPE='long' and this is not None, adds a constant 'specz' column
        filled with the provided value.
    """

    ### 1. Read input file ###
    ext = os.path.splitext(CAT_IN)[1].lower()

    if ext == '.fits':
        with pf.open(CAT_IN) as hdul:
            data = hdul[1].data
            arr = np.array(data)
            arr = arr.byteswap().view(arr.dtype.newbyteorder())
            df = pd.DataFrame(arr)

    elif ext == '.csv':
        df = pd.read_csv(CAT_IN)
    else:
        try:
            df = pd.read_table(CAT_IN, delim_whitespace=True, comment='#')
        except Exception as e:
            raise ValueError(f"Unrecognized input format for {CAT_IN}. Error: {e}")

    ### 2. Limit number of rows ###
    if max_rows is not None:
        if max_rows <= 0:
            raise ValueError("max_rows must be a positive integer.")
        df = df.iloc[:max_rows]
        print(f"[INFO] Limiting to first {max_rows} rows.")

    ### 3. Simple conversion ###
    if simple_convert:
        if shuffle:
            df = df.sample(frac=1).reset_index(drop=True)
            print("[INFO] Simple conversion with shuffle enabled.")
        np.savetxt(CAT_OUT, df.values, fmt="%.5f")
        print(f"[OK] Simple conversion completed → {CAT_OUT}")
        return
    
    ### 3bis. Add missing 'context' column automatically ###
    if CAT_TYPE == 'long' and apply_context.lower() in ['yes', 'null']:
        expected_len = 1 + n_filters * 2 + 2  # Id + (mags+errs) + context + zspec
        if (len(input_columns) == expected_len - 1 or expected_len - 2) and apply_context!='no':
            df['context'] = 0
            input_columns.insert(-1, 'context')
            print("[INFO] Added missing 'context' column automatically.")

    ### 4. Validate input columns ###
    if not all(col in df.columns for col in input_columns):
        missing = [c for c in input_columns if c not in df.columns]
        raise ValueError(f"Missing columns in input catalog: {missing}")

    ### 5. Extract selected columns ###
    selected = df[input_columns].copy()

    ### 6. Handle missing values ###
    if error_value_state not in ['default', 'delete']:
        raise ValueError("error_value_state must be 'default' or 'delete'.")

    if error_value_state == 'default':
        n_missing = selected.isna().sum().sum()
        if n_missing > 0:
            print(f"[INFO] Replacing {n_missing} missing values with 99.0")
            selected = selected.fillna(99.0)
    elif error_value_state == 'delete':
        n_before = len(selected)
        selected = selected.dropna().reset_index(drop=True)
        n_after = len(selected)
        print(f"[INFO] Removed {n_before - n_after} rows with missing values")

    ### 7. Check structure for 'long' catalogs ###
    if CAT_TYPE == 'long':
        expected_len = 1 + n_filters * 2 + 2  # Id + (mags+errs) + context + zspec
        print('number columns:',len(selected.columns))
        if facticious_specz is not None:
            selected['specz'] = facticious_specz
            print(f"[INFO] Added constant spec-z column with value {facticious_specz}")
        print('number columns:',len(selected.columns))
        if len(selected.columns) != expected_len:
            raise ValueError(
                f"Column count mismatch ({len(selected.columns)}) ≠ {expected_len} expected "
                f"for CAT_TYPE='long' and n_filters={n_filters}. "
                f"Check input_columns or parameters (context/specz)."
            )

    ### 8. Shuffle rows ###
    if shuffle:
        selected = selected.sample(frac=1).reset_index(drop=True)

    ### 9. Apply context logic ###
    arr = selected.to_numpy(dtype=object)
    if CAT_TYPE == 'long':
        arr = apply_context_function(arr, n_filters, apply_context)

    ### 10. Save formatted output ###
    if CAT_TYPE == 'short':
        fmt = ['%d'] + ['%.5f'] * (n_filters * 2)
    elif CAT_TYPE == 'long':
        fmt = ['%d'] + ['%.5f'] * (n_filters * 2) + ['%d', '%.5f']

    np.savetxt(CAT_OUT, arr, fmt=fmt)
    print(f"[OK] Catalog saved → {CAT_OUT}")

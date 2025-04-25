#Depecrated function, originally used to get the range of a given output in a specific region of the zz-plot

def third_mask_generator(data, pole1, pole2, column, x_col, y_col, boolean_op='or'):
    '''
    From a rectangle defined in (x_col, y_col) space, create a mask using the min/max of a third column.

    Parameters:
    - data: pandas DataFrame
    - pole1, pole2: tuples (x, y) defining opposite corners of a rectangle
    - column: name of column to use for the actual mask (e.g., 'SCALE_BEST')
    - x_col: column used as x-axis in the 2D region (e.g., 'ZSPEC')
    - y_col: column used as y-axis in the 2D region (e.g., 'Z_BEST')
    - boolean_op: 'or' to exclude values within min/max, 'and' to include them

    Returns:
    - Boolean mask (pd.Series)
    '''

    x1, y1 = pole1
    x2, y2 = pole2

    xmin, xmax = min(x1, x2), max(x1, x2)
    ymin, ymax = min(y1, y2), max(y1, y2)

    # Identify points inside the defined 2D region
    region_mask = (
        (data[x_col] >= xmin) & (data[x_col] <= xmax) &
        (data[y_col] >= ymin) & (data[y_col] <= ymax)
    )

    # Get the min/max of the target column within that region
    sub_values = data.loc[region_mask, column]

    if sub_values.empty:
        raise ValueError("No data points found in the defined 2D region.")

    min_val, max_val = sub_values.min(), sub_values.max()

    # Now build the final mask based on min/max of the `column`
    if boolean_op == 'or':
        # Exclude values within the range
        mask = (data[column] <= min_val) | (data[column] >= max_val)
    elif boolean_op == 'and':
        # Include only values within the range
        mask = (data[column] >= min_val) & (data[column] <= max_val)
    else:
        raise ValueError("boolean_op must be 'or' (exclude) or 'and' (include)")

    return mask
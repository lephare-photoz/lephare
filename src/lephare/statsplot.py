# A set of functions to display default interfaced plot.
# Intended to be used with lephare
import numpy as np
import inspect
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import statsmodels.api as sm
from scipy import stats


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


def save_masked_df(masked_output_df, original_input_path, new_input_path):
    """
    Filters the INPUT.dat file using the IDENT values from a masked LePhare output DataFrame
    and writes the result to a new file.

    Parameters
    ----------
    masked_output_df: pandas DataFrame.
        Masked dataframe to save in lephare format.
    original_input_path: str.
        Path to the original INPUT.dat file.
    new_input_path: str.
        Path where the filtered INPUT.dat will be saved.

    Returns
    -------
    None
    """
    # Load original input file
    input_data = []
    with open(original_input_path, "r") as f:
        for line in f:
            input_data.append(line.strip())

    # Extract IDENTs from masked output
    idents_to_keep = set(masked_output_df["IDENT"].astype(float))
    ident_list=[]
    # Write new input file with only the matching IDENTs
    with open(new_input_path, "w") as f:
        for line in input_data:
            ident = float(line.split()[0])
            if ident in idents_to_keep:
                f.write(line + "\n")
    print(f"Filtered input written to: {new_input_path}")


def handle_inputs_for_plot(x=None, y=None, data=None, x_col=None, y_col=None,
                           xlabel=None, ylabel=None):
    """
    Handle DataFrame inputs for plotting or further applications.
    """
    # --- Handle DataFrame input ---
    if data is not None:
        if isinstance(data, (list, tuple)):  # multiple DataFrames
            x_list, y_list = [], []
            for df in data:
                if not isinstance(df, pd.DataFrame):
                    raise TypeError("All elements in `data` must be pandas DataFrames.")
                if x_col not in df.columns or y_col not in df.columns:
                    raise KeyError(f"Columns '{x_col}' and/or '{y_col}' not found in one DataFrame.")
                x_list.append(df[x_col].to_numpy())
                y_list.append(df[y_col].to_numpy())
            x, y = x_list, y_list
        else:  # single DataFrame
            if not isinstance(data, pd.DataFrame):
                raise TypeError("`data` must be a pandas DataFrame if provided.")
            if x_col not in data.columns or y_col not in data.columns:
                raise KeyError(f"Columns '{x_col}' and/or '{y_col}' not found in DataFrame.")
            x = [data[x_col].to_numpy()]
            y = [data[y_col].to_numpy()]
        if xlabel is None:
            xlabel = x_col
        if ylabel is None:
            ylabel = y_col
    else:
        # --- Infer variable names from caller context ---
        frame = inspect.currentframe().f_back
        caller_locals = frame.f_locals
        if xlabel is None:
            for name, val in caller_locals.items():
                if val is x:
                    xlabel = name
                    break
            if xlabel is None:
                xlabel = "x"
        if ylabel is None:
            for name, val in caller_locals.items():
                if val is y:
                    ylabel = name
                    break
            if ylabel is None:
                ylabel = "y"

        # Handle x, y as possibly lists of arrays or single arrays
        if isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)):
            x = [np.asarray(xx) for xx in x]
            y = [np.asarray(yy) for yy in y]
        else:
            x = [np.asarray(x)]
            y = [np.asarray(y)]

    return x, y, xlabel, ylabel


def make_bins(mask_min, mask_max, nbins, logscale=False):
    if logscale:
        return np.logspace(np.log10(max(mask_min, 1e-6)),
                            np.log10(mask_max), nbins)
    else:
        return np.linspace(mask_min, mask_max, nbins)


def scatter_vs_hist2D(x=None, y=None, data=None, x_col='ZSPEC', y_col='Z_BEST',
                      deltaz=None, xlabel=None, ylabel=None, labels=None,
                      zmax=None, xrange=None, yrange=None, cmaps='viridis', dline=True):
    """
    Plot y vs x as scatter + 2D histogram with marginal histograms.
    Supports multiple datasets (lists of arrays or DataFrames).
    """

    # --- Handle input ---
    x_list, y_list, xlabel, ylabel = handle_inputs_for_plot(x, y, data, x_col, y_col, xlabel, ylabel)
    n_series = len(x_list)

    # --- Handle labels ---
    if labels is None:
        labels = [None for _ in range(n_series)]
    elif len(labels) != n_series:
        raise ValueError("Length of 'labels' must match number of datasets.")

    # --- Layout ---
    width_ratios = [2, 0.35, 2, 0.15]
    height_ratios = [0.35, 2]
    fig = plt.figure(constrained_layout=False,
                     figsize=(2.2 * sum(width_ratios), 2.2 * sum(height_ratios)))
    gs = fig.add_gridspec(2, 4, width_ratios=width_ratios, height_ratios=height_ratios)

    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scatter)
    ax_hist2d = fig.add_subplot(gs[1, 2])
    ax_cbar = fig.add_subplot(gs[1, 3])

    # --- Limits and bins ---
    allx = np.concatenate(x_list)
    ally = np.concatenate(y_list)
    x_min, x_max = (np.nanmin(allx), np.nanmax(allx))
    y_min, y_max = (np.nanmin(ally), np.nanmax(ally))

    if xrange is None:
        xrange = (x_min, x_max)
    if yrange is None:
        yrange = (y_min, y_max)

    bins_x = np.linspace(xrange[0], xrange[1], 100)
    bins_y = np.linspace(yrange[0], yrange[1], 100)

    # --- Scatter + histograms ---
    for i, (xx, yy) in enumerate(zip(x_list, y_list)):
        xx = np.ravel(xx)
        yy = np.ravel(yy)
        ax_scatter.scatter(xx, yy, s=10, alpha=0.25, label=labels[i])
        ax_histx.hist(xx, bins=bins_x, alpha=0.35, edgecolor='black')
        ax_histy.hist(yy, bins=bins_y, alpha=0.35, edgecolor='black', orientation='horizontal')

    # --- Scatter decorations ---
    if dline:
        ax_scatter.plot([xrange[0], xrange[1]], [yrange[0], yrange[1]], 'r--', label='y = x')
    if deltaz is not None:
        ax_scatter.plot([0, xrange[1]], [deltaz, yrange[1] + deltaz * (1 + xrange[1])], 'b--', label=r'+Δz')
        ax_scatter.plot([0, xrange[1]], [-deltaz, yrange[1] - deltaz * (1 + xrange[1])], 'b--', label=r'-Δz')

    ax_scatter.set_xlim(xrange)
    ax_scatter.set_ylim(yrange)
    ax_scatter.set_aspect('auto', adjustable='box')
    ax_scatter.set_xlabel(xlabel)
    ax_scatter.set_ylabel(ylabel)
    ax_scatter.legend(fontsize='small', loc='upper left', frameon=False)

    # --- 2D Histogram ---
    h = ax_hist2d.hist2d(allx, ally, bins=[bins_x, bins_y], cmap=cmaps)
    if dline:
        ax_hist2d.plot([xrange[0], xrange[1]], [yrange[0], yrange[1]], 'r--')
    ax_hist2d.set_xlim(xrange)
    ax_hist2d.set_ylim(yrange)
    ax_hist2d.set_aspect('auto', adjustable='box')
    ax_hist2d.set_xlabel(xlabel)
    ax_hist2d.tick_params(axis="y", left=False)

    # --- Colorbar ---
    plt.colorbar(h[3], cax=ax_cbar, label="Counts")

    # --- Clean up ---
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    plt.setp(ax_histy.get_yticklabels(), visible=False)
    plt.setp(ax_hist2d.get_yticklabels(), visible=False)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.05, hspace=0.02)



# z_spec = np.random.uniform(0, 2, (1000, 3))
# z_photo = z_spec + np.random.normal(0, 0.1, (1000, 3))
# labels = ["Training", "Validation", "Test"]
# scatter_vs_hist2D(z_spec, z_photo, deltaz=0.15, labels=labels)


def pit_qqplot(x=None, y=None, data=None, x_col='ZSPEC', y_col='Z_BEST',
               dist='norm', xlabel=None, ylabel=None, bins=30):
    """
    Display QQ-plot, QQ residuals, histograms of both distributions,
    and the Probability Integral Transform (PIT) in a compact GridSpec layout.

    Can work directly with a pandas DataFrame (by specifying or defaulting to 'ZSPEC' and 'Z_BEST').

    Parameters
    ----------
    data : pandas.DataFrame, optional
        Input table. If provided, x and y are taken from columns `x_col` and `y_col`.
    x, y : array-like, optional
        The two distributions to compare (e.g. true vs predicted).
    x_col, y_col : str, optional
        Column names to use from the DataFrame. Defaults are 'ZSPEC' and 'Z_BEST'.
    dist : str or scipy.stats distribution, optional
        Theoretical distribution to fit for PIT and PDF overlay (default 'norm').
    xlabel, ylabel : str, optional
        Axis labels for QQ-plot. If None, inferred from variable or column names.
    bins : int, optional
        Number of bins for histograms and PIT histogram.
    """

    # --- Handle input ---
    x, y, xlabel, ylabel = handle_inputs_for_plot(x, y, data, x_col, y_col, xlabel, ylabel)

    # --- Clean and align data ---
    x = np.asarray(x)
    y = np.asarray(y)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    n = min(len(x), len(y))
    x = np.sort(x)[:n]
    y = np.sort(y)[:n]

    # --- Fit distribution for PIT and PDF overlay ---
    x_df, x_loc, x_scale = stats.chi2.fit(x)
    pit_values = stats.chi2.cdf(y, x_df, loc=x_loc, scale=x_scale)

        # --- Create layout ---
    fig = plt.figure(figsize=(10, 7))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1.5, 1.0])
    plt.subplots_adjust(wspace=0.3, hspace=0.35)

    # === QQ-Plot ===
    ax_qq = fig.add_subplot(gs[0, 0])
    sm.qqplot_2samples(x, y, line="45", ax=ax_qq)
    ax_qq.set_title("QQ-Plot: Quantile Comparison", fontsize=13, weight='bold')
    ax_qq.set_xlabel(f"{xlabel} quantiles", fontsize=11)
    ax_qq.set_ylabel(f"{ylabel} quantiles", fontsize=11)
    ax_qq.grid(alpha=0.3)

    # === QQ Residuals ===
    ax_resid = fig.add_subplot(gs[1, 0], sharex=ax_qq)
    ax_resid.plot(x, y - x, 'o', alpha=0.4, markersize=3)
    ax_resid.axhline(0, color='red', linestyle='--', lw=1.5)
    ax_resid.set_xlabel(f"{xlabel} quantiles", fontsize=11)
    ax_resid.set_ylabel("Residuals", fontsize=11)
    ax_resid.set_title("QQ Residuals", fontsize=13, weight='bold')
    ax_resid.grid(alpha=0.3)

    # === Histograms ===
    ax_hist = fig.add_subplot(gs[0, 1])
    ax_hist.hist(x, bins=bins, density=True, alpha=0.6,
                 edgecolor='black', label=f"{xlabel}", range=(0,2))
    ax_hist.hist(y, bins=bins, density=True, alpha=0.6,
                 edgecolor='black', label=f"{ylabel}", range=(0,2))

    xx = np.linspace(min(x), max(x), 1000)
    x_pdf = stats.chi2.pdf(xx, x_df, loc=x_loc, scale=x_scale)
    ax_hist.plot(xx, x_pdf, 'r-', lw=1.5, label=f'{xlabel} fit')
    ax_hist.set_xlim(0, 2)
    
    ax_hist.set_title("Distribution Comparison", fontsize=13, weight='bold')
    ax_hist.set_xlabel("Value", fontsize=11)
    ax_hist.set_ylabel("Density (normalized)", fontsize=11)
    ax_hist.legend(frameon=False, fontsize=10)
    ax_hist.grid(alpha=0.3)

    # === PIT histogram ===
    ax_pit = fig.add_subplot(gs[1, 1])
    ax_pit.hist(pit_values, bins=bins, edgecolor='black', alpha=0.6, range=(0,1))
    ax_pit.axhline(len(pit_values) / bins, linestyle='--', color='red', lw=1.5)
    ax_pit.set_xlim(0, 1)
    ax_pit.set_xlabel("PIT value", fontsize=11)
    ax_pit.set_ylabel("Frequency", fontsize=11)
    ax_pit.set_title("Probability Integral Transform", fontsize=13, weight='bold')
    ax_pit.grid(alpha=0.3)

    # --- Final touches ---
    for ax in [ax_qq, ax_resid, ax_hist, ax_pit]:
        ax.tick_params(labelsize=10)

    plt.tight_layout()
    plt.show()
# z_spec = np.random.chisquare(df=4, size=1000)
# z_photo = z_spec + np.random.normal(0, 0.3, 1000)
# pit_qqplot(z_spec, z_photo, dist='chi2', xlabel="z_spec", ylabel="z_photo")


def chi_stats(x=None, y=None, data=None, x_col='CHI_STAR', y_col='CHI_BEST',
              xlabel=None, ylabel=None, labels=None,
              mask_min=None, mask_max=None, bins=100, log=False,
              scatter_xlim=(None, None), scatter_ylim=(None, None)):
    """
    Plot Chi² statistics (e.g. CHI_STAR vs CHI_BEST) for one or several datasets.
    Shows histograms and scatter plots, supports log scaling and multiple inputs.

    Compatible with handle_inputs_for_plot() and scatter_vs_hist2D().

    Parameters
    ----------
    x, y : array-like, list of arrays, or None
        Chi² values for comparison (if DataFrame not provided).
    data : pandas.DataFrame or list of DataFrames, optional
        Source data for x and y columns.
    x_col, y_col : str, optional
        Column names in the DataFrame(s).
    xlabel, ylabel : str, optional
        Axis labels.
    labels : list of str, optional
        Names for each dataset.
    mask_min, mask_max : float, optional
        Clip values below/above these limits (for histograms only).
    bins : int, optional
        Number of histogram bins.
    log : bool, optional
        If True, use log-scale on x-axes for histograms and scatter plot.
    scatter_xlim, scatter_ylim : tuple, optional
        Limits for scatter plot axes.
    """

    import numpy as np
    import matplotlib.pyplot as plt

    # --- Handle input ---
    x_list, y_list, xlabel, ylabel = handle_inputs_for_plot(
        x, y, data, x_col, y_col, xlabel, ylabel
    )
    n_series = len(x_list)

    # --- Labels ---
    if labels is None:
        labels = [None for i in range(n_series)]
    elif len(labels) != n_series:
        raise ValueError("Length of 'labels' must match number of datasets.")

    # --- Figure layout ---
    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)
    ax_star = fig.add_subplot(gs[0, 0])
    ax_best = fig.add_subplot(gs[0, 1])
    ax_delta = fig.add_subplot(gs[1, 0])
    ax_scatter = fig.add_subplot(gs[1, 1])

    # --- Loop over datasets ---
    for i, (xx, yy) in enumerate(zip(x_list, y_list)):
        xx = np.ravel(xx)
        yy = np.ravel(yy)

        # Always scatter valid finite points (no masking)
        mask_scatter = np.isfinite(xx) & np.isfinite(yy)
        xx_scat = xx[mask_scatter]
        yy_scat = yy[mask_scatter]
        ax_scatter.scatter(xx_scat, yy_scat, s=5, alpha=0.3, label=labels[i])

        # For histograms: apply mask_min/max
        mask_x_hist = np.isfinite(xx)
        mask_y_hist = np.isfinite(yy)
        if mask_min is not None:
            mask_x_hist &= xx >= mask_min
            mask_y_hist &= yy >= mask_min
        if mask_max is not None:
            mask_x_hist &= xx <= mask_max
            mask_y_hist &= yy <= mask_max
        
        xx_hist = xx[mask_x_hist]
        yy_hist = yy[mask_y_hist]

        delta_chi = None
        if len(xx) > 0 and len(yy) > 0:
            delta_chi = xx - yy
        
        # --- Try to plot each histogram independently ---
        try:
            if len(xx_hist) > 0:
                ax_star.hist(xx_hist, bins=make_bins(0, mask_max, bins), alpha=0.6,
                             edgecolor='black', label=labels[i])
        except Exception as e:
            print(f"No CHI_STAR values in range for dataset'{labels[i]}': {e}")

        try:
            if len(yy_hist) > 0:
                ax_best.hist(yy_hist, bins=make_bins(0, mask_max, bins), alpha=0.6,
                             edgecolor='black', label=labels[i])
        except Exception as e:
            print(f"No CHI_BEST values in range for dataset '{labels[i]}': {e}")

        try:
            if delta_chi is not None and len(delta_chi) > 0:
                ax_delta.hist(delta_chi, bins=make_bins(mask_min, mask_max, bins), alpha=0.6,
                              edgecolor='black', label=labels[i])
        except Exception as e:
            print(f"No ΔChi values in range for dataset '{labels[i]}': {e}")


    # --- Axes formatting ---
    for ax in [ax_star, ax_best, ax_delta]:
        ax.set_xlabel("Chi² value")
        ax.set_ylabel("Count")
        if log:
            ax.set_xscale("log")

    ax_star.set_title(f"{xlabel} distribution", fontsize=13, weight='bold')
    ax_best.set_title(f"{ylabel} distribution", fontsize=13, weight='bold')
    ax_delta.set_title(f"ΔChi² = {xlabel} - {ylabel}", fontsize=13, weight='bold')
    ax_delta.axvline(0, color='gray', linestyle='--')

    # --- Scatter formatting ---
    ax_scatter.set_title(f"{xlabel} vs {ylabel}", fontsize=13, weight='bold')
    ax_scatter.set_xlabel(xlabel)
    ax_scatter.set_ylabel(ylabel)
    ax_scatter.plot([0, 1e9], [0, 1e9], 'k--', lw=1)
    if scatter_xlim[0] is not None and scatter_xlim[1] is not None:
        ax_scatter.set_xlim(scatter_xlim)
    if scatter_ylim[0] is not None and scatter_ylim[1] is not None:
        ax_scatter.set_ylim(scatter_ylim)
    if log:
        ax_scatter.set_xscale("log")
        ax_scatter.set_yscale("log")

    # --- Legends ---
    for ax in [ax_star, ax_best, ax_delta, ax_scatter]:
        ax.legend(fontsize=9, frameon=False)
    plt.tight_layout()
    plt.show()


# np.random.seed(0)
# chi_star = np.random.chisquare(df=4, size=3000)
# chi_best = np.random.chisquare(df=3, size=3000)
# chi_stats(x=chi_star, y=chi_best, log=False)
# # chi_stats(x=chi_star, y=chi_best, mask_max=20, log=True)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def histograms(data_list, col='Z_BEST', labels=None, bins=100, 
               xrange=(0,2), log=False, density=False, xlabel=None, 
               ylabel="Count", title=None, figsize=(8, 5), alpha=0.6):
    """
    Plot multiple histograms on the same figure.
    Supports lists, numpy arrays, or pandas DataFrames.

    Parameters
    ----------
    data_list : list
        List of arrays, Series, or DataFrames.
    col : str, optional
        Column name to extract if elements of `data_list` are DataFrames.
    labels : list of str, optional
        Labels for each series.
    bins : int, optional
        Number of histogram bins (default 100).
    range : tuple, optional
        (min, max) range for histogram.
    log : bool, optional
        Use log scale on y-axis.
    density : bool, optional
        Normalize histograms to unit area.
    xlabel, ylabel, title : str, optional
        Axis labels and title.
    figsize : tuple, optional
        Figure size.
    alpha : float, optional
        Transparency of histograms.
    """

    # --- Input normalization ---
    if not isinstance(data_list, (list, tuple)):
        data_list = [data_list]
    n_series = len(data_list)

    bins=make_bins(xrange[0], xrange[1], bins)

    # --- Labels ---
    if labels is None:
        labels = [None for i in range(n_series)]
    elif len(labels) != n_series:
        raise ValueError("Length of 'labels' must match number of datasets.")

    # --- Extract arrays ---
    values = []
    for data in data_list:
        if isinstance(data, pd.DataFrame):
            if col is None:
                raise ValueError("If passing DataFrames, specify `col` to select a column.")
            values.append(data[col].to_numpy())
        elif isinstance(data, (pd.Series, np.ndarray, list)):
            values.append(np.asarray(data))
        else:
            raise TypeError("Unsupported type in data_list; must be DataFrame, Series, ndarray, or list.")

    # --- Plot ---
    plt.figure(figsize=figsize)
    for vals, label in zip(values, labels):
        vals = vals[np.isfinite(vals)]  # ignore NaN/inf
        plt.hist(vals, bins=bins, range=xrange, alpha=alpha,
                 edgecolor='black', label=label, density=density)

    # --- Style ---
    plt.xlabel(xlabel or col or "Value", fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    if log:
        plt.yscale("log")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.show()



def more_chi_stats(x=None, y=None, data=None, x_col='CHI_STAR', y_col='CHI_BEST',
              xlabel=None, ylabel=None, mask_min=None, mask_max=None, bins=100, log=False):
    """
    Plot the Chi² statistics of a LePhare run (template fitting on galaxies and stars),
    showing model comparison between galaxy and star fits.

    Supports direct array input (x, y) or a pandas DataFrame with column names.

    Parameters
    ----------
    x, y : array-like, optional
        Chi² values for star and best-fit galaxy models.
    data : pandas.DataFrame, optional
        Input DataFrame containing chi² columns.
    x_col, y_col : str, optional
        Column names to use from the DataFrame (defaults: 'CHI_STAR' and 'CHI_BEST').
    xlabel, ylabel : str, optional
        Axis labels (inferred automatically if None).
    mask_min, mask_max : float, optional
        Minimum/maximum χ² values to keep. If None, no clipping.
    bins : int, optional
        Number of bins for histograms. Default = 100.
    log : bool, optional
        If True, sets log-scale on x-axis for histograms and scatter.

    Returns
    -------
    None
        Displays a 2x2 grid:
            [0,0] ΔCHI histogram
            [0,1] CHI_BEST vs ΔCHI scatter
            [1,0] CHI_STAR vs ΔCHI scatter
            [1,1] 
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # --- Handle input ---
    x, y, xlabel, ylabel = handle_inputs_for_plot(x, y, data, x_col, y_col, xlabel, ylabel)

    # --- Clean / mask data ---
    x = np.asarray(x)
    y = np.asarray(y)

    # Set sensible defaults if None
    if mask_min is None:
        mask_min = 1e-3 if log else np.nanmin([np.nanmin(x), np.nanmin(y)])
    if mask_max is None:
        mask_max = 1e9 if log else np.nanmax([np.nanmax(x), np.nanmax(y)])

    # Apply mask
    mask = np.isfinite(x) & np.isfinite(y) & (x >= mask_min) & (y <= mask_max)
    x, y = x[mask], y[mask]
    delta_chi = x - y


    # --- Determine plotting ranges dynamically ---
    hist_range = (mask_min, mask_max)
    scatter_xlim = (
        np.nanmin(x) if mask_min is None else mask_min,
        np.nanmax(x) if mask_max is None else mask_max
    )
    scatter_ylim = (
        np.nanmin(y) if mask_min is None else mask_min,
        np.nanmax(y) if mask_max is None else mask_max
    )

    # --- Setup figure ---
    fig = plt.figure(constrained_layout=False, figsize=(10, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1], hspace=0.35, wspace=0.3)

    # --- [0,0] Histogram of DELTA CHI ---
    ax_delta = fig.add_subplot(gs[0, 0])
    ax_delta.hist(delta_chi, bins=bins, alpha=0.7, edgecolor='black', color='purple', range=hist_range)
    ax_delta.axvline(0, color='gray', linestyle='--')
    ax_delta.set_title(f"ΔChi² = {xlabel} - {ylabel}")
    ax_delta.set_xlabel("ΔChi²")
    ax_delta.set_ylabel("Count")
    if log:
        ax_delta.set_xscale("log")

    # --- [0,1] CHI_BEST vs DELTA_CHI ---
    ax_scatter = fig.add_subplot(gs[0, 1])
    ax_scatter.scatter(y, x-y, s=5, alpha=0.5, c='orange')
    ax_scatter.set_xlim(*scatter_xlim)
    ax_scatter.set_ylim(*scatter_ylim)
    ax_scatter.set_title(f"{ylabel} vs ΔChi²")
    ax_scatter.set_xlabel(ylabel)
    ax_scatter.set_ylabel("ΔChi²")
    ax_scatter.legend(loc='upper left', fontsize=8, frameon=False)
    if log:
        ax_scatter.set_xscale("log")
        ax_scatter.set_yscale("log")

    # --- [1,0] CHI_STAR vs DELTA_CHI ---
    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_scatter.scatter(x, x-y, s=5, alpha=0.5, c='blue')
    ax_scatter.set_xlim(*scatter_xlim)
    ax_scatter.set_ylim(*scatter_ylim)
    ax_scatter.set_title(f"{xlabel} vs ΔChi²")
    ax_scatter.set_xlabel(xlabel)
    ax_scatter.set_ylabel("ΔChi²")
    ax_scatter.legend(loc='upper left', fontsize=8, frameon=False)
    if log:
        ax_scatter.set_xscale("log")
        ax_scatter.set_yscale("log")

    # --- [1,1] Scatter CHI_STAR vs CHI_BEST ---
    ax_scatter = fig.add_subplot(gs[1, 1])
    ax_scatter.scatter(x, y, s=5, alpha=0.5, c='teal')
    ax_scatter.set_xlim(-scatter_xlim[1]/50, scatter_xlim[1])
    ax_scatter.set_ylim(-scatter_ylim[1]/50, scatter_ylim[1])
    ax_scatter.set_title(f"{xlabel} vs {ylabel}")
    ax_scatter.set_xlabel(xlabel)
    ax_scatter.set_ylabel(ylabel)
    ax_scatter.legend(loc='upper left', fontsize=8, frameon=False)
    if log:
        ax_scatter.set_xscale("log")
        ax_scatter.set_yscale("log")

    plt.tight_layout()
    plt.show()
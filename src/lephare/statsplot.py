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
    Parameters
    ----------
    data : pandas.DataFrame, optional
        Input table. If provided, x and y are taken from columns `x_col` and `y_col`.
    x, y : array-like or 2D arrays, optional
        Input data. Ignored if `data` is provided.
    x_col, y_col : str, optional
        Column names to use from the DataFrame. Defaults are 'ZSPEC' and 'Z_BEST'.
    xlabel, ylabel : str, optional
        Axis labels. If not provided, inferred from variable or column names.

    
    Returns
    -------
    x, y : array-like or 2D arrays, optional
        Input data for the scatter plot and 2D histogram.
    xlabel, ylabel : str, optional
        Axis labels.
    """

    # --- Handle DataFrame input ---
    if data is not None:
        if not isinstance(data, pd.DataFrame):
            raise TypeError("`data` must be a pandas DataFrame if provided.")
        if x_col not in data.columns or y_col not in data.columns:
            raise KeyError(f"Columns '{x_col}' and/or '{y_col}' not found in DataFrame.")
        x = data[x_col].to_numpy()
        y = data[y_col].to_numpy()
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


    # --- Prepare data shape ---
    x = np.asarray(x)
    y = np.asarray(y)

    return x, y, xlabel, ylabel


def scatter_vs_hist2D(x=None, y=None, data=None, x_col='ZSPEC', y_col='Z_BEST',
                      deltaz=None, xlabel=None, ylabel=None, labels=None, cmaps='viridis'):
    """
    Plot y vs x as scatter + 2D histogram with marginal histograms.

    Supports both 1D and 2D input arrays, or a pandas DataFrame.
    If x and y are 2D (same shape), each pair (x_i, y_i) is plotted 
    with a different color (shared across scatter and histograms).
    
    Optionally, you can provide custom labels for each data component.

    Parameters
    ----------
    data : pandas.DataFrame, optional
        Input table. If provided, x and y are taken from columns `x_col` and `y_col`.
    x, y : array-like or 2D arrays, optional
        Input data for the scatter plot and 2D histogram. Ignored if `data` is provided.
    x_col, y_col : str, optional
        Column names to use from the DataFrame. Defaults are 'ZSPEC' and 'Z_BEST'.
    deltaz : float, optional
        If given, adds dashed blue lines representing ±Δz tolerance around the y=x line.
    xlabel, ylabel : str, optional
        Axis labels. If not provided, inferred from variable or column names.
    labels : list of str, optional
        Labels for each series (used in legend). Matches number of columns in x/y.

    Returns
    -------
    None
        Displays the combined figure.
    """

    # --- Handle input ---
    x, y, xlabel, ylabel = handle_inputs_for_plot(x, y, data, x_col, y_col, xlabel, ylabel)


    # --- Prepare data shape ---
    x = np.asarray(x)
    y = np.asarray(y)

    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape.")

    # Convert 1D inputs to 2D for consistent looping
    if x.ndim == 1:
        x = x[:, None]
        y = y[:, None]

    n_series = x.shape[1]

    # --- Handle custom labels ---
    if labels is None:
        labels = [f"Set {i+1}" for i in range(n_series)]
    else:
        if len(labels) != n_series:
            raise ValueError("Length of 'labels' must match the number of series (columns in x).")

    # --- Colors ---
    cmap = plt.cm.tab10  # consistent color palette
    colors = cmap(np.linspace(0, 1, n_series))

    # --- Figure layout ---
    width_ratios = [2, 0.35, 2, 0.15]
    height_ratios = [0.35, 2]

    fig = plt.figure(
        constrained_layout=False,
        figsize=(2.2 * sum(width_ratios), 2.2 * sum(height_ratios))
    )
    gs = fig.add_gridspec(
        2, 4, width_ratios=width_ratios, height_ratios=height_ratios
    )

    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0])
    ax_histy = fig.add_subplot(gs[1, 1])
    ax_hist2d = fig.add_subplot(gs[1, 2])
    ax_cbar = fig.add_subplot(gs[1, 3])


    # --- Scatter plots ---
    for i in range(n_series):
        ax_scatter.scatter(
            x[:, i], y[:, i], s=10, alpha=0.15, color=colors[i], label=labels[i]
        )

    ax_scatter.plot([0, np.max(x)], [0, np.max(x)], 'r--', label='y = x')

    if deltaz is not None:
        ax_scatter.plot(
            [0, np.max(x)],
            [deltaz, np.max(x) + deltaz * (1 + np.max(x))],
            'b--', label=r'Upper cut ($+\Delta z$)'
        )
        ax_scatter.plot(
            [0, np.max(x)],
            [-deltaz, np.max(x) - deltaz * (1 + np.max(x))],
            'b--', label=r'Lower cut ($-\Delta z$)'
        )

    ax_scatter.set_xlim(0, np.max(x))
    ax_scatter.set_ylim(0, np.max(x))
    ax_scatter.set_aspect('equal', adjustable='box')
    ax_scatter.set_xlabel(xlabel)
    ax_scatter.set_ylabel(ylabel)
    ax_scatter.legend(fontsize='small', loc='upper left', frameon=False)

    # --- Marginal histograms (same colors + labels) ---
    for i in range(n_series):
        ax_histx.hist(
            x[:, i], bins=40, alpha=0.3, edgecolor='black',
            color=colors[i], label=labels[i]
        )
        ax_histy.hist(
            y[:, i], bins=40, alpha=0.3, edgecolor='black',
            orientation='horizontal', color=colors[i], label=labels[i]
        )

    ax_histx.tick_params(axis="x", bottom=False)
    ax_histx.set_ylabel('Counts')

    ax_histy.tick_params(axis="y", left=False)
    ax_histy.set_xlabel('Counts')

    # --- 2D histogram of all points combined ---
    h = ax_hist2d.hist2d(x.flatten(), y.flatten(), bins=80, cmap=cmaps)
    ax_hist2d.plot([0, np.max(x)], [0, np.max(x)], 'r--')
    ax_hist2d.set_xlim(0.01, np.max(x))
    ax_hist2d.set_ylim(0, np.max(x))
    ax_hist2d.set_aspect('equal', adjustable='box')
    ax_hist2d.set_xlabel(xlabel)
    ax_hist2d.tick_params(axis="y", left=False)

    # --- Colorbar ---
    plt.colorbar(h[3], cax=ax_cbar, label="Counts")

    # --- Final layout tweaks ---
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    plt.setp(ax_histy.get_yticklabels(), visible=False)
    plt.setp(ax_hist2d.get_yticklabels(), visible=False)
    plt.tight_layout()

    # --- Final manual layout tuning ---
    plt.subplots_adjust(wspace=0.05, hspace=0.02)

    for ax in [ax_histx, ax_histy]:
        ax.margins(0)
    for ax in [ax_scatter, ax_histx, ax_histy, ax_hist2d]:
        ax.set_anchor('C')

    pos_scatter = ax_scatter.get_position()
    pos_histy = ax_histy.get_position()
    ax_histy.set_position([
        pos_histy.x0, pos_scatter.y0,
        pos_histy.width, pos_scatter.height
        
    ])
    pos_histx = ax_histx.get_position()
    ax_histx.set_position([
        pos_scatter.x0, pos_histx.y0,
        pos_scatter.width, pos_histx.height-0.01
    ])

    pos_scatter = ax_scatter.get_position()
    pos_hist2d = ax_hist2d.get_position()
    pos_cbar = ax_cbar.get_position()

    # aligne verticalement le 2D hist sur le scatter
    ax_hist2d.set_position([
        pos_histy.x1 + 0.03,    # décalé vers la droite
        pos_scatter.y0,             # même bas que le scatter
        pos_hist2d.width,           # garde la largeur du 2D hist
        pos_scatter.height          # même hauteur que le scatter
    ])

    new_pos_hist2d = ax_hist2d.get_position()
    ax_cbar.set_position([
        new_pos_hist2d.x1 + 0.01,
        new_pos_hist2d.y0,
        pos_cbar.width * 0.7,
        new_pos_hist2d.height
    ])




    plt.show()

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
    ax_resid.set_ylabel("Residuals (y - x)", fontsize=11)
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
            [0,0] CHI_STAR histogram
            [0,1] CHI_BEST histogram
            [1,0] ΔCHI histogram
            [1,1] CHI_STAR vs CHI_BEST scatter
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

    # --- [0,0] Histogram of CHI_STAR ---
    ax_star = fig.add_subplot(gs[0, 0])
    print(x)
    if log:
        logbins = np.logspace(np.log10(np.min(x)),
                            np.log10(np.max(x)), int(np.log10(np.max(x))-np.log10(np.min(x))+1))
        ax_star.set_xscale("log")
    if log:
        print(logbins)
    ax_star.hist(x, bins=logbins if log else bins, alpha=0.7, edgecolor='black', range=(0, mask_max))
    ax_star.set_title(f"{xlabel} distribution")
    ax_star.set_xlabel("Chi² value")
    ax_star.set_ylabel("Count")

    # --- [0,1] Histogram of CHI_BEST ---
    ax_best = fig.add_subplot(gs[0, 1])
    ax_best.hist(y, bins=bins, alpha=0.7, edgecolor='black', color='orange', range=(0, mask_max))
    ax_best.set_title(f"{ylabel} distribution")
    ax_best.set_xlabel("Chi² value")
    ax_best.set_ylabel("Count")
    if log:
        ax_best.set_xscale("log")

    # --- [1,0] Histogram of ΔCHI ---
    ax_delta = fig.add_subplot(gs[1, 0])
    ax_delta.hist(delta_chi, bins=bins, alpha=0.7, edgecolor='black', color='purple', range=hist_range)
    ax_delta.axvline(0, color='gray', linestyle='--')
    ax_delta.set_title(f"ΔChi² = {xlabel} - {ylabel}")
    ax_delta.set_xlabel("ΔChi²")
    ax_delta.set_ylabel("Count")
    if log:
        ax_delta.set_xscale("log")

    # --- [1,1] Scatter CHI_STAR vs CHI_BEST ---
    ax_scatter = fig.add_subplot(gs[1, 1])
    ax_scatter.scatter(x, y, s=5, alpha=0.5, c='teal')
    # ax_scatter.plot([0, max(scatter_xlim[1], scatter_ylim[1])],
    #                 [0, max(scatter_xlim[1], scatter_ylim[1])],
    #                 'r--', label="y=x")
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


# np.random.seed(0)
# chi_star = np.random.chisquare(df=4, size=3000)
# chi_best = np.random.chisquare(df=3, size=3000)
# chi_stats(x=chi_star, y=chi_best, log=False)
# # chi_stats(x=chi_star, y=chi_best, mask_max=20, log=True)


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
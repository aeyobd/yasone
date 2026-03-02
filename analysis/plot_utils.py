import matplotlib.pyplot as plt
import numpy as np
import sys
from ana_utils import fit_err

def plot_mag_comparison(
    table,
    filters=("g", "r", "i"),
    mag_ref_fmt="{filt}_MAG",
    mag_cmp_fmt="{filt}_MAG_2",
    mag_ref_err_fmt=None,
    mag_cmp_err_fmt=None,
    ref_mag_upper = True,
    cmp_mag_upper = True,
    ref_label="reference",
    cmp_label="comparison",
    colors=None,
    xlim=(26, 17),
    ylim_mag=(26, 17),
    ylim_res=(1, -1),
    figsize=(7, 4),
    scatter_kwargs=None,
):
    """
    Plot magnitude comparison and residual panels for multiple filters.

    Parameters
    ----------
    table : astropy Table
        Joined catalogue containing magnitude columns.
    filters : tuple
        Filters to plot (e.g., ("g","r","i")).
    mag_ref_fmt : str
        Format string for reference magnitude column.
    mag_cmp_fmt : str
        Format string for comparison magnitude column.
    ref_label : str
        Y-label for top-left magnitude panel.
    cmp_label_fmt : str
        X-label format string for bottom panels.
    colors : list or None
        List of colors (length must match filters).
    scatter_kwargs : dict or None
        Extra kwargs passed to plt.scatter.
    """

    if scatter_kwargs is None:
        scatter_kwargs = dict(s=1, alpha=0.1)

    if colors is None:
        colors = ["C0"] * len(filters)

    fig, axs = plt.subplots(
        2,
        len(filters),
        sharex=True,
        sharey="row",
        figsize=figsize,
        gridspec_kw={"hspace": 0, "wspace": 0},
    )

    # Ensure axs is always 2D (even if one filter)
    axs = np.atleast_2d(axs)

    for i, filt in enumerate(filters):

        if ref_mag_upper:
            mag_ref_col = mag_ref_fmt.format(filt=filt.upper())
        else:
            mag_ref_col = mag_ref_fmt.format(filt=filt)

        if cmp_mag_upper:
            mag_cmp_col = mag_cmp_fmt.format(filt=filt.upper())
        else:
            mag_cmp_col = mag_cmp_fmt.format(filt=filt)

        mag_ref = table[mag_ref_col]
        mag_cmp = table[mag_cmp_col]


        color = colors[i]

        # ---- Top row: magnitude comparison ----
        ax = axs[0, i]
        ax.scatter(mag_ref, mag_cmp, c=color, **scatter_kwargs)

        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim_mag)
        ax.plot(xlim, xlim, color="black", zorder=-1)

        if i == 0:
            ax.set_ylabel(f"mag ({cmp_label})")

        # ---- Bottom row: residual ----
        ax_res = axs[1, i]
        ax_res.scatter(mag_ref, mag_cmp - mag_ref, c=color, **scatter_kwargs)

        if mag_ref_err_fmt is not None:
            if ref_mag_upper:
                mag_ref_err_col = mag_ref_err_fmt.format(filt=filt.upper())
            else:
                mag_ref_err_col = mag_ref_err_fmt.format(filt=filt)

            if cmp_mag_upper:
                mag_cmp_err_col = mag_cmp_err_fmt.format(filt=filt.upper())
            else:
                mag_cmp_err_col = mag_cmp_err_fmt.format(filt=filt)

            mag_ref_err = table[mag_ref_err_col]
            mag_cmp_err = table[mag_cmp_err_col]
            mag_err = np.sqrt(mag_ref_err**2 + mag_cmp_err**2)
            plot_sliding_window_mag_error(ax_res, mag_ref, mag_err)

        ax_res.set_ylim(*ylim_res)
        ax_res.axhline(0, color="black", zorder=-1)
        ax_res.set_xlabel(f"{filt} mag ({ref_label})")

        if i == 0:
            ax_res.set_ylabel("residual")

    plt.tight_layout()
    return fig, axs



def plot_sliding_window_mag_error(ax, mags, mag_errors, window_size=100):
    """
    Plot sliding-window mean magnitude error vs magnitude.

    window_size = number of objects per window.
    """

    mags = np.asarray(mags)
    mag_errors = np.asarray(mag_errors)

    # Sort by magnitude
    sort_idx = np.argsort(mags)
    mags_sorted = mags[sort_idx]
    errors_sorted = mag_errors[sort_idx]

    centers = []
    means = []

    for i in range(len(mags_sorted) - window_size):
        window_mags = mags_sorted[i:i+window_size]
        window_errs = errors_sorted[i:i+window_size]

        centers.append(np.nanmedian(window_mags))
        means.append(np.nanmedian(window_errs))

    ax.plot(centers, means, '-k')
    ax.plot(centers, [-m for m in means], '-k')

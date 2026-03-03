import matplotlib.pyplot as plt
import numpy as np
import sys
from ana_utils import fit_err

from mpl_toolkits.axes_grid1 import make_axes_locatable



def plot_mag_2d_residual(
        table,
        filters=("g", "r", "i"),
        mag_ref_fmt="{filt}_MAG",
        mag_cmp_fmt="{filt}_MAG_2",
        max_diff=1,
        figsize=(4, 5)
    ):

    fig, axs = plt.subplots(
        1, len(filters),
        sharey=True,
        sharex=True,
        figsize=figsize,
        gridspec_kw={"hspace": 0, "wspace": 0},
    )

    im = None

    for i, filt in enumerate(filters):
        plt.sca(axs[i])
        FILT = filt.upper()
        mag_ref_col = mag_ref_fmt.format(FILT=FILT, filt=filt)
        mag_cmp_col = mag_cmp_fmt.format(FILT=FILT, filt=filt)

        mag_ref = table[mag_ref_col]
        mag_cmp = table[mag_cmp_col]

        resid = mag_ref - mag_cmp
        filt_good = np.isfinite(resid)
        filt_good &= ~resid.mask

        im = plt.scatter(table["xi"][filt_good], table["eta"][filt_good],
                         c=resid[filt_good],
                    cmap="RdBu", vmin=-max_diff, vmax=max_diff, s=1, lw=0)
        axs[i].set_aspect(1)
        plt.xlabel("xi")
        plt.title(filt)
        if i == 0:
            plt.ylabel("eta")


    divider = make_axes_locatable(axs[-1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
       
    plt.colorbar(im, cax=cax, label="residual (mag)")


    plt.tight_layout()


    


def plot_mag_residual_hist(
        table,
        filters=("g", "r", "i"),
        mag_ref_fmt="{filt}_MAG",
        mag_cmp_fmt="{filt}_MAG_2",
        mag_ref_err_fmt=None,
        mag_cmp_err_fmt=None,
        xlim=(-5, 5),
        sigma_scale=1,
        bins = 50,
    ):
    bins = np.linspace(*xlim, bins)


    for i, filt in enumerate(filters):
        FILT = filt.upper()
        mag_ref_col = mag_ref_fmt.format(FILT=FILT, filt=filt)
        mag_cmp_col = mag_cmp_fmt.format(FILT=FILT, filt=filt)

        mag_ref = table[mag_ref_col]
        mag_cmp = table[mag_cmp_col]

        mag_ref_err_col = mag_ref_err_fmt.format(filt=filt, FILT=FILT)
        mag_cmp_err_col = mag_cmp_err_fmt.format(filt=filt, FILT=FILT)

        mag_ref_err = table[mag_ref_err_col]
        mag_cmp_err = table[mag_cmp_err_col]
        mag_err = np.sqrt(mag_ref_err**2 + mag_cmp_err**2)

        z = (mag_cmp - mag_ref) / mag_err

        filt_good = np.isfinite(z)
        filt_good &= ~z.mask
        plt.hist(z[filt_good], bins, label=filt, histtype="step", density=True)

    x = np.linspace(*xlim, 1000)
    y = 1/np.sqrt(2*sigma_scale*np.pi) * np.exp(-x**2 / 2 / sigma_scale**2)
    plt.plot(x, y, color="k")

    plt.xlabel("z-score")
    plt.legend()


def plot_mag_comparison(
        table,
        filters=("g", "r", "i"),
        mag_ref_fmt="{filt}_MAG",
        mag_cmp_fmt="{filt}_MAG_2",
        mag_ref_err_fmt=None,
        mag_cmp_err_fmt=None,
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
        3,
        len(filters),
        sharex=True,
        sharey="row",
        figsize=figsize,
        gridspec_kw={"hspace": 0, "wspace": 0},
    )

    # Ensure axs is always 2D (even if one filter)
    axs = np.atleast_2d(axs)

    for i, filt in enumerate(filters):
        FILT = filt.upper()

        mag_ref_col = mag_ref_fmt.format(FILT=FILT, filt=filt)
        mag_cmp_col = mag_cmp_fmt.format(FILT=FILT, filt=filt)
        mag_cmp_col = mag_cmp_fmt.format(FILT=FILT, filt=filt)

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

        # ---- Middle row: residual ----
        ax_res = axs[1, i]
        ax_res.scatter(mag_ref, mag_cmp - mag_ref, c=color, **scatter_kwargs)

        mag_err = None
        if mag_ref_err_fmt is not None:
            mag_ref_err_col = mag_ref_err_fmt.format(filt=filt, FILT=FILT)
            mag_cmp_err_col = mag_cmp_err_fmt.format(filt=filt, FILT=FILT)

            mag_ref_err = table[mag_ref_err_col]
            mag_cmp_err = table[mag_cmp_err_col]
            mag_err = np.sqrt(mag_ref_err**2 + mag_cmp_err**2)
            plot_sliding_window_mag_error(ax_res, mag_ref, mag_err)

        ax_res.set_ylim(*ylim_res)
        ax_res.axhline(0, color="black", zorder=-1)


        ax_z = axs[2, i]
        if mag_ref_err_fmt is not None:
            # ---- Bottom row: residual ----
            ax_z.scatter(mag_ref, (mag_cmp - mag_ref)/mag_err, c=color, **scatter_kwargs)

            ax_z.set_ylim(-5, 5)
            ax_z.axhline(0, color="black", zorder=-1)

            ax_z.set_xlabel(f"{filt} mag ({ref_label})")
        else:
            ax_z.remove()
            ax_res.set_xlabel(f"{filt} mag ({ref_label})")

        if i == 0:
            ax_res.set_ylabel("residual")
            ax_z.set_ylabel("z-score")

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

        centers.append(np.ma.median(window_mags))
        means.append(np.ma.median(window_errs))

    ax.plot(centers, means, '-k')
    ax.plot(centers, [-m for m in means], '-k')

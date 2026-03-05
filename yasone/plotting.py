import matplotlib.pyplot as plt
import numpy as np
from .analysis import fit_err

from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import visualization as aviz
from astropy.nddata.blocks import block_reduce
from astropy.nddata.utils import Cutout2D
from astropy.nddata import CCDData

from matplotlib import pyplot as plt
import numpy as np
import seaborn



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
        axs[i].invert_xaxis()
        plt.title(filt)
        if i == 0:
            plt.ylabel("eta")


    plt.tight_layout()

    plt.colorbar(im, ax=(axs), label="residual (mag)",fraction=0.02, pad=0.04)

    


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

            ax_z.set_ylim(5, -5)
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



cmap = seaborn.color_palette("mako", as_cmap=True)

def show_images(images, titles=None, figscale=3, **kwargs):
    """
    Shows a series of images as a row
    """

    num_images = len(images)

    fig, axs = plt.subplots(1, num_images, 
                            figsize=(num_images * figscale, figscale))

    for i in range(num_images):
        show_image(images[i], fig=fig, ax=axs[i], **kwargs)
        if titles is not None:
            axs[i].set_title(titles[i])

    plt.tight_layout()

    

def show_image_residual(image, image2, perc=99, clim=None, **kwargs):
    diff = image - image2
    if clim is None:
        val_low = np.abs(np.percentile(diff, 100 - perc))
        val_high = np.abs(np.percentile(diff, perc))
        val = np.maximum(val_low, val_high)
        clim = (-val, val)

    return show_image(diff, clim=clim, cmap="RdBu")


def show_image(image,
               percl=99, percu=None, is_mask=False,
               figsize=(10, 10),
               dpi = None,
               cmap=cmap, log=False, clip=True,
               clabel=None,
               clim=None,
               show_colorbar=True, show_ticks=True,
               fig=None, ax=None, input_ratio=None):
    """
    Show an image in matplotlib with some basic astronomically-appropriate stretching.

    Parameters
    ----------
    image
        The image to show
    percl : number
        The percentile for the lower edge of the stretch (or both edges if ``percu`` is None)
    percu : number or None
        The percentile for the upper edge of the stretch (or None to use ``percl`` for both)
    figsize : 2-tuple
        The size of the matplotlib figure in inches
    """
    if percu is None:
        percu = percl
        percl = 100 - percl

    if dpi is None:
        dpi = image.shape[0] / figsize[0]

    if (fig is None and ax is not None) or (fig is not None and ax is None):
        raise ValueError('Must provide both "fig" and "ax" '
                         'if you provide one of them')
    elif fig is None and ax is None:
        if figsize is not None:
            # Rescale the fig size to match the image dimensions, roughly
            image_aspect_ratio = image.shape[0] / image.shape[1]
            figsize = (max(figsize) * image_aspect_ratio, max(figsize))

        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)


    # To preserve details we should *really* downsample correctly and
    # not rely on matplotlib to do it correctly for us (it won't).

    # So, calculate the size of the figure in pixels, block_reduce to
    # roughly that,and display the block reduced image.

    # Thanks, https://stackoverflow.com/questions/29702424/how-to-get-matplotlib-figure-size

    fig_size_pix = fig.get_size_inches() * fig.dpi

    ratio = (image.shape // fig_size_pix).max()

    if ratio < 1:
        ratio = 1

    ratio = input_ratio or ratio

    reduced_data = block_reduce(image, ratio)

    if not is_mask:
        # Divide by the square of the ratio to keep the flux the same in the
        # reduced image. We do *not* want to do this for images which are
        # masks, since their values should be zero or one.
         reduced_data = reduced_data / ratio**2

    # Of course, now that we have downsampled, the axis limits are changed to
    # match the smaller image size. Setting the extent will do the trick to
    # change the axis display back to showing the actual extent of the image.
    extent = [0, image.shape[1], 0, image.shape[0]]

    if log:
        stretch = aviz.LogStretch()
    else:
        stretch = aviz.LinearStretch()

    if clim is None:
        clim = aviz.AsymmetricPercentileInterval(percl, percu)
    else:
        clim = aviz.ManualInterval(*clim)

    norm = aviz.ImageNormalize(reduced_data,
                               interval=clim,
                               stretch=stretch, clip=clip)

    if is_mask:
        # The image is a mask in which pixels should be zero or one.
        # block_reduce may have changed some of the values, so reset here.
        reduced_data = reduced_data > 0
        # Set the image scale limits appropriately.
        scale_args = dict(vmin=0, vmax=1)
    else:
        scale_args = dict(norm=norm)

    im = ax.imshow(reduced_data, origin='lower',
                   cmap=cmap, extent=extent, aspect='equal', **scale_args)

    if show_colorbar:
        # I haven't a clue why the fraction and pad arguments below work to make
        # the colorbar the same height as the image, but they do....unless the image
        # is wider than it is tall. Sticking with this for now anyway...
        # Thanks: https://stackoverflow.com/a/26720422/3486425
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=clabel)
        # In case someone in the future wants to improve this:
        # https://joseph-long.com/writing/colorbars/
        # https://stackoverflow.com/a/33505522/3486425
        # https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes

    if not show_ticks:
        ax.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False)


def image_snippet(image, center, width=50, axis=None, fig=None,
                  is_mask=False, pad_black=False, **kwargs):
    """
    Display a subsection of an image about a center.

    Parameters
    ----------

    image : numpy array
        The full image from which a section is to be taken.

    center : list-like
        The location of the center of the cutout.

    width : int, optional
        Width of the cutout, in pixels.

    axis : matplotlib.Axes instance, optional
        Axis on which the image should be displayed.

    fig : matplotlib.Figure, optional
        Figure on which the image should be displayed.

    is_mask : bool, optional
        Set to ``True`` if the image is a mask, i.e. all values are
        either zero or one.

    pad_black : bool, optional
        If ``True``, pad edges of the image with zeros to fill out width
        if the slice is near the edge.
    """
    if pad_black:
        sub_image = Cutout2D(image, center, width, mode='partial', fill_value=0)
    else:
        # Return a smaller subimage if extent goes out side image
        sub_image = Cutout2D(image, center, width, mode='trim')
    show_image(sub_image.data, cmap='gray', ax=axis, fig=fig,
               show_colorbar=False, show_ticks=False, is_mask=is_mask,
               **kwargs)


def _mid(sl):
    return (sl.start + sl.stop) // 2


from .plotting import show_image

from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.table import Table, join, vstack

plate_scale = 0.254 # arcmin / pixel

"""
The airmass to extinction coefficient wth uncertainty from the OSIRIS/GTC
Broad Band Imaging calibration documentation.
"""
EXTINCTION_DICT = {
    'Sloan_u': [0.45, 0.02],
    'Sloan_g': [0.15, 0.02],
    'Sloan_r': [0.07, 0.01],
    'Sloan_i': [0.04, 0.01],
    'Sloan_z': [0.03, 0.01]
}


def plot_sources(objects, img, scale=1):
    """
    Plots the sources from a sep.extract run as red ellipses
    on the given image, scaled up by the provided scale.
    """
    show_image(img)

    ax = plt.gca()
    
    # plot an ellipse for each object
    for i in range(len(objects)):
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=scale*objects['a'][i],
                    height=scale*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi, linewidth=0.5)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)


def swap_byteorder(data):
    """
    Swaps the byteorder of the data. Necessary to work with sep.
    """
    data = data.astype(data.dtype.newbyteorder('='))
    return data


def to_mag(flux, flux_err, zero_point=0, atm_extinction=0):
    mag =  zero_point - 2.5*np.log10(flux)  - atm_extinction
    magerr = 2.5 / np.log(10) * flux_err / flux
    return mag, magerr


def get_atm_extinction(airmass, filt):
    A, A_err = EXTINCTION_DICT[filt]
    return  A * airmass, A_err * airmass


def xmatch(ra1, dec1, ra2, dec2, rmax=1*u.arcsec):
    """
    Crossmatch two RA/Dec catalogs and keep only unique matches
    in catalog2, retaining the closest match.
    Assumes RA and Dec are in degrees.

    Returns
    -------
    idx1 : array
        Indices in catalog1 of surviving matches
    idx2 : array
        Corresponding unique indices in catalog2
    sep  : Quantity
        Angular separations of surviving matches
    """

    coord1 = coordinates.SkyCoord(ra1, dec1, unit=u.deg)
    coord2 = coordinates.SkyCoord(ra2, dec2, unit=u.deg)

    xmatch_indicies, xmatch_sep, _ = coord1.match_to_catalog_sky(coord2)

    # Initial radius cut
    has_xmatch = xmatch_sep < rmax

    # Work only on valid matches
    idx1_valid = np.where(has_xmatch)[0]
    idx2_valid = xmatch_indicies[has_xmatch]
    sep_valid = xmatch_sep[has_xmatch]

    # Sort valid matches by separation (closest first)
    order = np.argsort(sep_valid)
    idx1_sorted = idx1_valid[order]
    idx2_sorted = idx2_valid[order]
    sep_sorted = sep_valid[order]

    # Keep only first occurrence of each catalog2 index
    _, unique_keep = np.unique(idx2_sorted, return_index=True)

    idx1_keep = idx1_sorted[unique_keep]

    # Reset duplicates to "no match"
    has_xmatch_unique = np.zeros_like(has_xmatch)
    has_xmatch_unique[idx1_keep] = True

    # For remved duplicates, set separation to large value
    xmatch_sep_unique = xmatch_sep.copy()
    xmatch_sep_unique[~has_xmatch_unique] = np.inf * xmatch_sep.unit

    Ndup = np.sum(has_xmatch) - len(has_xmatch_unique)
    if Ndup > 0:
        print(f"removing {Ndup} duplicate xmatches")

    return xmatch_indicies, xmatch_sep_unique, has_xmatch_unique



def outer_join_xmatch(cat1, cat2, ra1="ra", dec1="dec",
                      ra2="ra", dec2="dec",
                      xmatch_radius = 1 * u.arcsec,
                      lsuffix="_1",
                      rsuffix="_2",
                      **kwargs,
                      ):
    xmatch_idx, xmatch_sep, has_xmatch = xmatch(
            cat1[ra1], cat1[dec1], cat2[ra2], cat2[dec2],
            xmatch_radius)

    cat1 = cat1.copy()
    cat1.rename_columns(cat1.colnames, [name + lsuffix for name in cat1.colnames])
    cat2 = cat2.copy()
    cat2.rename_columns(cat2.colnames, [name + rsuffix for name in cat2.colnames])

    # create indexes to use for xmatch
    cat1["_idx1"] = np.arange(len(cat1))
    cat2["_idx2"] = np.arange(len(cat2))
    matched_cat1 = cat1[has_xmatch]
    matched_cat2 = cat2[xmatch_idx[has_xmatch]]

    # rename matched indicies
    matched_cat2["_idx1"] = matched_cat1["_idx1"]

    matched_join = join(
        matched_cat1, 
        matched_cat2,
        keys = "_idx1",
        join_type="inner",
        **kwargs
    )

    assert len(matched_join) == np.sum(has_xmatch)

    unmatched_cat1 = cat1[~has_xmatch]
    matched_idx2 = xmatch_idx[has_xmatch]
    unmatched_mask2 = np.ones(len(cat2), dtype=bool)
    unmatched_mask2[matched_idx2] = False
    unmatched_cat2 = cat2[unmatched_mask2]


    print("matched count", len(matched_join))
    print("unmatched left", len(unmatched_cat1))
    print("total left", len(cat1))
    print("unmatched right", len(unmatched_cat2))
    print("total right", len(cat2))

    assert len(matched_join) + len(unmatched_cat1) == len(cat1)
    assert len(matched_join) + len(unmatched_cat2) == len(cat2)

    # Remove temporary keys
    for t in (matched_join, unmatched_cat1, unmatched_cat2):
        for col in ['_idx1', '_idx2']:
            if col in t.colnames:
                t.remove_column(col)

    # Stack everything (outer behavior)
    result = vstack([matched_join, unmatched_cat1, unmatched_cat2],
                    join_type='outer')

    for col in result.colnames:
        if col.endswith(lsuffix):
            col_new = col[:-len(lsuffix)]
            if col_new +rsuffix not in result.colnames:
                result.rename_column(col, col_new)
        if col.endswith(rsuffix):
            col_new = col[:-len(rsuffix)]
            if col_new + lsuffix not in result.colnames:
                result.rename_column(col, col_new)
            
    return result


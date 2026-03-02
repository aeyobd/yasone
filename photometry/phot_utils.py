import sys
sys.path.append("../imaging")
from convenience_functions import show_image
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.table import Table, join, vstack

plate_scale = 0.254 # arcmin / pixel

# From documentation (OSIRIS/GTC Broad Band Imaging calibration)
extinction_dict = {
    'Sloan_u': [0.45, 0.02],
    'Sloan_g': [0.15, 0.02],
    'Sloan_r': [0.07, 0.01],
    'Sloan_i': [0.04, 0.01],
    'Sloan_z': [0.03, 0.01]
}

# There are also colour corrections
#u’ – u’0 = -25.807(±0.053) - 0.071 (±0.023) (u’0 – g’0)
#g’ – g’0 = -28.823 (±0.040) - 0.078 (±0.013) (g’0 – r’0)
#r’ – r’0 = -29.291 (±0.017) - 0.114 (±0.028) (r’0 – i’0)
#i’ – i’0 = -28.857 (±0.015) - 0.079 (±0.041) (i’0 – z’0)
#z’ – z’0 = -28.231(±0.031) - 0.072 (±0.052) (i’0 – z’0)


def plot_sources(objects, img, scale=6):
    show_image(img,)

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
    data = data.astype(data.dtype.newbyteorder('='))
    return data



def calc_zero_point(flux, mag, atm_extinction):
    zero_point = mag_std + 2.5*np.log10(counts_std / exp_time_std)


def to_mag(flux, flux_err, zero_point, atm_extinction):
    mag =  zero_point - 2.5*np.log10(flux)  - atm_extinction
    magerr = 2.5 / np.log(10) * flux_err / flux
    return mag, magerr

def get_atm_extinction(airmass, filt):
    A, A_err = extinction_dict[filt]
    return  A * airmass, A_err * airmass


def xmatch(ra1, dec1, ra2, dec2, rmax=1*u.arcsec):
    coord1 = coordinates.SkyCoord(ra1, dec1, unit=u.deg)
    coord2 = coordinates.SkyCoord(ra2, dec2, unit=u.deg)

    xmatch_indicies, xmatch_sep, _ = coord1.match_to_catalog_sky(coord2)

    has_xmatch = xmatch_sep < rmax
    return xmatch_indicies, xmatch_sep, has_xmatch



def outer_join_xmatch(cat1, cat2, ra1="ALPHA_J2000", dec1="DELTA_J2000",
                      ra2="ALPHA_J2000", dec2="DELTA_J2000",
                      lprefix="1_", rprefix="2_",
                      xmatch_radius = 1 * u.arcsec,
                      **kwargs
                      ):
    xmatch_idx, xmatch_sep, has_xmatch = xmatch(
            cat1[ra1], cat1[dec1], cat2[ra2], cat2[dec2],
            xmatch_radius)

    cat1 = cat1.copy()
    cat2 = cat2.copy()

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
    #assert len(matched_join) + len(unmatched_cat2) == len(cat2)

        # Remove temporary keys
    for t in (matched_join, unmatched_cat1, unmatched_cat2):
        for col in ['_idx1', '_idx2']:
            if col in t.colnames:
                t.remove_column(col)

    # Stack everything (outer behavior)
    result = vstack([matched_join, unmatched_cat1, unmatched_cat2],
                    join_type='outer')

    #assert len(result) == len(matched_join) + len(unmatched_cat1) + len(unmatched_cat2)
    return result


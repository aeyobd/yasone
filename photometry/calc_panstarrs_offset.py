from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

import tomli_w

import arya
from astropy.table import Table
from astropy import units as u

import astropy
import tomllib

import sys
sys.path.append(".")
sys.path.append("../imaging")

from phot_utils import xmatch, get_atm_extinction
from astropy.nddata import CCDData



def get_panstarrs(objname):
    panstarrs = Table.read(f"../survey_data/{objname}_panstarrs.fits")
    for col in panstarrs.colnames:
        panstarrs[col] = panstarrs[col].reshape(-1)

    panstarrs_filt = (panstarrs["qualityFlag"] & 4) > 0

    for filt in ["G", "R", "I"]:
        panstarrs_filt &= panstarrs[f"{filt.lower()}Flags"] & (1 + 2 + 4) == 0
        panstarrs_filt &= panstarrs[f"{filt.lower()}MeanPSFMag"] > -10
        # panstarrs_filt &= panstarrs[f"{filt}MeanPSFMagNpt"].reshape(-1) > 2

    # apply corrections from Tonry et al. 2012, Apj.759.99, table 6,
    # the SDSS - PS1 vs g-r PS1 quadratics

    gr = panstarrs["gMeanPSFMag"] - panstarrs["rMeanPSFMag"]

    panstarrs["G_MAG"] = panstarrs["gMeanPSFMag"] + 0.013 + 0.145*gr + 0.019*gr**2
    panstarrs["R_MAG"] = panstarrs["rMeanPSFMag"] - 0.001 + 0.004*gr + 0.007*gr**2
    panstarrs["I_MAG"] = panstarrs["iMeanPSFMag"] - 0.005 + 0.011*gr + 0.010*gr**2

    # technically should also add ri uncertainty here, but this doesnt matter
    # too much
    panstarrs["G_MAG_ERR"] = panstarrs["gMeanPSFMagErr"] + 0.008
    panstarrs["R_MAG_ERR"] = panstarrs["rMeanPSFMagErr"] + 0.004
    panstarrs["I_MAG_ERR"] = panstarrs["iMeanPSFMagErr"] + 0.004



    print("median correction shifts")
    print("g", np.ma.median(panstarrs["G_MAG"] - panstarrs["gMeanPSFMag"]))
    print("r", np.ma.median(panstarrs["R_MAG"] - panstarrs["rMeanPSFMag"]))
    print("i", np.ma.median(panstarrs["I_MAG"] - panstarrs["iMeanPSFMag"]))
    return panstarrs[panstarrs_filt]


def determine_zeropoints(cat):
    panstarrs_zeropoints = {}

    for filtname in ["G", "R", "I"]:
        mag0 = cat[f"{filtname}_MAG"]
        mag_ps = cat[f"{filtname}_MAG_PS1"]
        filt = cat["R_FLAGS"] == 0
        filt &= ~mag0.mask
        filt &= mag_ps > -10

        w = 1/(cat[f"{filtname}_MAG_PS1_ERR"]**2)
        residuals = mag0[filt] - mag_ps[filt]

        m_weighted = weighted_median(residuals, w[filt])

        m = np.median(np.array(residuals)) # convert to arra
        # first to ignore mask warning
        m_sc, med_sc, err_sc = astropy.stats.sigma_clipped_stats(residuals)
        print("filter: ", filtname)
        print("weighted median, median", m_weighted, m)
        print("sigmaclipped mean, median: ", m_sc, ",", med_sc)
        panstarrs_zeropoints[filtname.lower()] = med_sc

    return panstarrs_zeropoints




def match_catalogues(cat_gtc, panstarrs):
    xmatch_idx, xmatch_dist, xmatch_filt_panstarrs = xmatch(
        cat_gtc["ra"], cat_gtc["dec"],
        panstarrs["ra"].reshape(-1), panstarrs["dec"].reshape(-1), 
        1*u.arcsec
    )

    panstarrs_xmatch = panstarrs[xmatch_idx]
    osiris_x_panstarrs = cat_gtc[xmatch_filt_panstarrs]

    for filt in ["G", "R", "I"]:
        cat_gtc[f"{filt}_MAG_PS1"] = panstarrs_xmatch[f"{filt}_MAG"].reshape(-1)
        cat_gtc[f"{filt}_MAG_PS1_ERR"] = panstarrs_xmatch[f"{filt}_MAG_ERR"].reshape(-1)

    return cat_gtc[xmatch_filt_panstarrs]


def main():
    if len(sys.argv) != 3:
        print("usage: scriptname.py objname catname")
        sys.exit(1)

    objname, catname = sys.argv[1:]

    cat_gtc = Table.read(f"{objname}/{catname}.cat")
    cat_panstarrs = get_panstarrs(objname)
    cat_matched = match_catalogues(cat_gtc, cat_panstarrs)

    cat_matched.write(f"{objname}/{catname}_x_PS.cat", format="fits", overwrite=True)
    zeropoints = determine_zeropoints(cat_matched)

    with open(Path(objname) / f"{catname}_panstarrs_shift.toml", "wb") as f:
        tomli_w.dump(zeropoints, f)

    plot_residuals(cat_matched, zeropoints)
    plt.savefig(f"{objname}/figures/{catname}_panstarrs_residuals.pdf")



def plot_residuals(cat, panstarrs_zeropoints):
    fig, axs = plt.subplots(2, 3, sharex=True, sharey="row", figsize=(7, 4), gridspec_kw={"hspace": 0, "wspace": 0})
    c = None
    for i in range(3):
        plt.sca(axs[0][i])
        cs = arya.COLORS[0] #np.log10(cat["FWHM_WORLD"].to("arcsec") / u.arcsec)

        filt = ["g", "r", "i"][i]
        mag0 = cat[f"{filt.upper()}_MAG"] - panstarrs_zeropoints[filt]
        mag_ps = cat[f"{filt.upper()}_MAG_PS1"]
        mag_ps_err = cat[f"{filt.upper()}_MAG_PS1_ERR"]
        plt.scatter(mag0, mag_ps, s=1,  c=cs, alpha=0.1)
        plt.xlim(24, 17)
        plt.ylim(24, 17)
        plt.plot([26, 10], [26, 10], color="black", zorder=-1)
        if i == 0:
            plt.ylabel(f"mag (PS)")
        
        plt.sca(axs[1][i])
        
        yerr = np.maximum(mag_ps_err, 0)
        # plt.errorbar(mag0[1:-1:20], (mag_ps-mag0)[1:-1:20], yerr=yerr[1:-1:20], fmt=".", capsize=0, )
        c = plt.scatter(mag0, mag_ps - mag0, s=1, c=cs, alpha=0.1)

        plt.ylim(3, -3)
        plt.axhline(0, color="black", zorder=-1)
        
        plt.xlabel(f"{filt} (osiris)")
        if i == 0:
            plt.ylabel("residual")


    plt.tight_layout()



def weighted_median(values, weights):
    """
    Compute the weighted median of an array of values.
    
    This implementation sorts values and computes the cumulative
    sum of the weights. The weighted median is the smallest value for
    which the cumulative sum is greater than or equal to half of the
    total sum of weights.
    Parameters
    ----------
    values : array-like
        List or array of values on which to calculate the weighted median.
    weights : array-like
        List or array of weights corresponding to the values.
    Returns
    -------
    float
        The weighted median of the input values.
    """
    # Convert input values and weights to numpy arrays
    values = np.array(values)
    weights = np.array(weights)
    
    # Get the indices that would sort the array
    sort_indices = np.argsort(values)
    
    # Sort values and weights according to the sorted indices
    values_sorted = values[sort_indices]
    weights_sorted = weights[sort_indices]  

    # Compute the cumulative sum of the sorted weights
    cumsum = weights_sorted.cumsum()
    
    # Calculate the cutoff as half of the total weight sum
    cutoff = weights_sorted.sum() / 2.
    
    # Return the smallest value for which the cumulative sum is greater
    # than or equal to the cutoff
    return values_sorted[cumsum >= cutoff][0]



if __name__ == "__main__":
    main()

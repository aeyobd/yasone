import os

from pathlib import Path

from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt
import arya
import numpy as np


import tomllib


from astropy.stats import sigma_clipped_stats


from astropy.coordinates import SkyCoord
import sys
sys.path.append(".")
sys.path.append("../imaging")

from phot_utils import outer_join_xmatch, get_atm_extinction, to_mag


airmass = 1.2 # calibrated out anyways



def get_zeropoint(filt, exposure=190, gain=1.9):
    d = None
    with open(f"std1/img_{filt}_11/flat_fielded-astrom-zeropoint.toml", "rb") as f:
        d = tomllib.load(f)

    return d["zeropoint"] + d["ap_corr"] + 2.5 * np.log10(exposure) - get_atm_extinction(airmass, f"Sloan_{filt}")[0] - 2.5 * np.log10(gain)



def get_airmass(file):
    path = Path(file).parent
    img = CCDData.read(path / "flat_fielded.fits")
    return img.header["AIRMASS"]


def calibrate_mag(cat, filt):
    idx_aper = 3

    filt_bad = (cat["MAG_APER_3"] > 50) | (cat["FLAGS"] & 64 + 128 + 32 + 16 > 0) 
    cat["MAG"] = cat["MAG_APER_3"] + get_zeropoint(filt) 
    cat["MAG_ERR"] = cat["MAGERR_APER_3"]

    cat["MAG"][filt_bad] = np.nan
    cat["MAG_ERR"][filt_bad] = np.nan


def flatten_arrays(cat):
    cat_out = Table()

    for col in cat.colnames:
        if len(cat[col].shape) > 1:
            for i in range(cat[col].shape[1]):
                cat_out[col + f"_{i}"] = cat[col][:, i]
        else:
            cat_out[col] = cat[col]
    return cat_out


def main():
    assert len(sys.argv) in [5, 6]
    if len(sys.argv) == 5:
        objname, catname, outname, foldername = sys.argv[1:]
        suffix=""
    else:
        objname, catname, outname, foldername, suffix = sys.argv[1:]




    filename_g = f"{objname}/{foldername}_g{suffix}/{catname}"
    filename_r = f"{objname}/{foldername}_r{suffix}/{catname}"
    filename_i = f"{objname}/{foldername}_i{suffix}/{catname}"


    cat_g = Table.read(filename_g, hdu=2)
    cat_r = Table.read(filename_r, hdu=2)
    cat_i = Table.read(filename_i, hdu=2)
    cat_g = flatten_arrays(cat_g)
    cat_r = flatten_arrays(cat_r)
    cat_i = flatten_arrays(cat_i)


    calibrate_mag(cat_g, "g")
    calibrate_mag(cat_r, "r")
    calibrate_mag(cat_i, "i")

    cat_g.rename_columns(cat_g.colnames, ["G_" + name for name in
                                          cat_g.colnames])
    cat_r.rename_columns(cat_r.colnames, ["R_" + name for name in
                                          cat_r.colnames])

    cat_i.rename_columns(cat_i.colnames, ["I_" + name for name in
                                          cat_i.colnames])

    combined = outer_join_xmatch(cat_r, cat_g, 
                                 ra1 = "R_ALPHA_J2000", dec1="R_DELTA_J2000",
                                 ra2 = "G_ALPHA_J2000", dec2="G_DELTA_J2000",
                                 )


    combined["ra"] = combined["R_ALPHA_J2000"].copy()
    combined["dec"] = combined["R_DELTA_J2000"].copy()

    filt = combined["ra"].mask | combined["dec"].mask
    combined["ra"][filt] = combined["G_ALPHA_J2000"][filt]
    combined["dec"][filt] = combined["G_DELTA_J2000"][filt]


    combined = outer_join_xmatch(combined, cat_i, 
                                 ra1 = "ra", dec1="dec",
                                 ra2 = "I_ALPHA_J2000", dec2="I_DELTA_J2000",
                                 )


    filt = combined["ra"].mask | combined["dec"].mask
    combined["ra"][filt] = combined["I_ALPHA_J2000"][filt]
    combined["dec"][filt] = combined["I_DELTA_J2000"][filt]

    combined.write(f"{objname}/{outname}{suffix}.cat", format="fits", overwrite=True)


if __name__ == "__main__":
    main()

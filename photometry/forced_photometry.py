import sys
from dataclasses import dataclass, field
from pathlib import Path

from typing import List
import numpy as np

from photutils import aperture as phot_aperture
from photutils import background
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.table import Table
import astropy.units as u

APERTURE_RADII = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0]

@dataclass
class aperture_params:
    aperture_radii:list = field(default_factory=lambda: APERTURE_RADII)
    back_filter_size: int = 5
    back_size: int = 128
    xy_bounds: float = 4.0
    psf_size: int = 25



def create_aperture_phot(img, positions, params):
    apertures = [phot_aperture.CircularAperture(positions, r=r) for r in
                 params.aperture_radii]
    phot_table = phot_aperture.aperture_photometry(img, apertures)

    # rename header info for FITS standard
    phot_table.meta = {k.replace("aperture_photometry", "ap").replace("aperture", "ap"): v 
                       for k, v in phot_table.meta.items()}
    return phot_table



def load_image(objname, filtname, imgname):
    img = CCDData.read(f"./{objname}/img_{filtname}_{imgname}/flat_fielded-astrom.fits")
    img.wcs = CCDData.read(f"./{objname}/img_{filtname}_{imgname}/nobkg.fits").wcs
    img_err = CCDData.read(f"./{objname}/img_{filtname}_{imgname}/flat_fielded.weight.fits")
    img_flag = CCDData.read(f"./{objname}/img_{filtname}_{imgname}/flag.fits",
                            unit="adu")

    img.uncertainty = StdDevUncertainty(img_err.data)
    img.mask = img_flag.data > 4
    return img


def load_ref_cat(obj, filt):
    cat = Table.read(f"./{obj}/coadd_median_{filt}/detection.cat", hdu=2)
    img_coadd = CCDData.read(f"./{obj}/coadd_median_{filt}/coadd.fits", unit="adu")
    cat_coords = img_coadd.wcs.pixel_to_world(cat["X_IMAGE"] , cat["Y_IMAGE"])
    return cat, cat_coords


def to_image_coords(img, cat, coords, pad=0):
    x_cat_all, y_cat_all = img.wcs.world_to_pixel(coords)

    filt_inbounds = x_cat_all >= 0 + pad
    filt_inbounds &= x_cat_all <= img.shape[1] - pad - 1
    filt_inbounds &= y_cat_all >= 0 + pad
    filt_inbounds &= y_cat_all <= img.shape[0] - pad - 1
    filt_inbounds &= cat["FLAGS"] <= 64
    filt_inbounds &= cat["IMAFLAGS_ISO"] <= 4

    index = np.arange(len(coords))

    x_cat = x_cat_all[filt_inbounds]
    y_cat = y_cat_all[filt_inbounds]
    index_cat = index[filt_inbounds]

    return [a for a in zip(x_cat, y_cat)], index_cat


def subtract_background(img, params):
    bkg_estimator = background.SExtractorBackground()
    bkg = background.Background2D(img, (params.back_size, params.back_size), 
        filter_size=(params.back_filter_size, params.back_filter_size),
        bkg_estimator=bkg_estimator)

    img_nobkg = CCDData(img-bkg.background, unit="adu")
    tot_err = np.sqrt(bkg.background_rms**2 + img.uncertainty.array**2 * u.adu**2)

    print("median rms:", bkg.background_rms_median)

    img_nobkg.uncertainty = StdDevUncertainty(tot_err)
    return img_nobkg


params_dict = {
    "": aperture_params(),
    "_small_bkg": aperture_params(back_size=64),
    "_large_bkg": aperture_params(back_size=256),
    "_small_bkg_filt": aperture_params(back_filter_size=3),
    "_large_bkg_filt": aperture_params(back_filter_size=11),
}


def main(objname, filtname, params, suffix=""):
    cat, coords = load_ref_cat(objname, filtname)

    for imgpath in Path(f"{objname}").glob(f"img_{filtname}_*"):
        imgname = imgpath.stem[-2:]
        img = load_image(objname, filtname, imgname)
        img_coords, index = to_image_coords(img, cat, coords)

        img_nobkg = subtract_background(img, params)

        phot_result = create_aperture_phot(img_nobkg, img_coords, params)
        phot_result["index"] = index

        phot_result.write(imgpath / f"forced_ap_phot{suffix}.fits",
                          overwrite=True)


    
if __name__ == "__main__":
    objname, filtname = sys.argv[1:]

    for name, params in params_dict.items():
        main(objname, filtname, params, suffix=name)

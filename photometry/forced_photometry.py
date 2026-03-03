import sys
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

from photutils import aperture as phot_aperture
from photutils import background

import astropy
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.table import Table
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy import units as u

APERTURE_RADII = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0]


@dataclass
class PhotometryParams:
    """
    A class containing the parameters for use with aperture and psf
    photometry.

    Attributes:
        aperture_radii (list): The radii of apertures to calculate photometry
        for
        back_size (int): The size of the SourceExtractor-like background mesh
        back_filter_size (odd int): The size of the additional median filter on
        the background
        xy_bounds (float): the maximum that coordinates are allowed to vary in x and y
        (pixels)
        psf_size (int): The size to use for creating the emperical psf and
        fitting psf
        psf_oversample (int): How much to oversample the psf (for the emperical
        psf estimation)
        psf_aperture_radius (float): The aperture to use for the initial psf
        flux guess. 
    """

    aperture_radii:list = field(default_factory=lambda: APERTURE_RADII)
    back_filter_size: int = 5
    back_size: int = 154 # correct for slightly smaller back_filter_size
    xy_bounds: float = 4.0
    local_background: bool = False
    local_background_r_in: float = 15
    local_background_r_out: float = 20
    psf_size: int = 25
    psf_oversample: int = 3
    psf_aperture_radius: float = 5
    outdir:str = "."

def get_median_background(img, positions, params):
    annulus_aperture = phot_aperture.CircularAnnulus(positions,
                                       r_in=params.local_background_r_in, 
                                       r_out=params.local_background_r_out)

    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    bkg_stats = phot_aperture.ApertureStats(img, annulus_aperture,
                                            sigma_clip=sigclip)

    return bkg_stats.median


def create_aperture_phot(img, positions, params):
    """
    Given an image, positions (as a list of tuples), and PhotometryParams,
    returns a Table containing the aperture photometry results from phot_utils.
    """
    apertures = [phot_aperture.CircularAperture(positions, r=r) for r in
                 params.aperture_radii]
    phot_table = phot_aperture.aperture_photometry(img, apertures)

    # local background
    phot_table["median_local_background"] = get_median_background(img,
                                                                  positions,
                                                                  params)


    for i, ap in enumerate(params.aperture_radii):
        phot_table[f"aperture_sum_lb_{i}"] = (
                phot_table[f"aperture_sum_{i}"] - np.pi *
                params.aperture_radii[i]**2 *
                phot_table["median_local_background"]
                )
    # rename header info for FITS standard
    phot_table.meta = {k.replace("aperture_photometry", "ap").replace("aperture", "ap"): v 
                       for k, v in phot_table.meta.items()}
    return phot_table



def load_image(objname, filtname, imgname):
    """
    Given the object name (yasone1-3), the filter name (g, r, i) and the image
    name (usually 00 through 15), returns a CCDData containing the image (with
    background), but updating the image mask and uncertainty.
    """

    img_path = Path(f"./{objname}/img_{filtname}_{imgname}")
    img = CCDData.read(img_path / "flat_fielded-astrom.fits")

    head = fits.Header.fromtextfile(img_path / "nobkg.head")
    img.wcs = astropy.wcs.WCS(head)

    img_err = CCDData.read(img_path / "flat_fielded.weight.fits")
    img_flag = CCDData.read(img_path / "flag.fits", unit="adu")

    img.uncertainty = StdDevUncertainty(img_err.data)
    img.mask = img_flag.data > 4
    return img


def load_ref_cat_coords(objname, filtname):
    """
    Based on the object name (obj) and the filter name (filt),
    returns the reference (stacked) catalogue and the 
    """
    ref_path = Path(f"./{objname}/coadd_median_{filtname}")
    cat = Table.read(ref_path / "detection.cat", hdu=2)
    img_coadd = CCDData.read(ref_path / "coadd.fits", unit="adu")

    cat_coords = img_coadd.wcs.pixel_to_world(cat["X_IMAGE"]-1,
                                              cat["Y_IMAGE"]-1)
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
    "": PhotometryParams(),
    #"_small_bkg": PhotometryParams(back_size=64),
    #"_large_bkg": PhotometryParams(back_size=256),
    #"_small_bkg_filt": PhotometryParams(back_filter_size=3),
    #"_large_bkg_filt": PhotometryParams(back_filter_size=11),
}


def main(objname, filtname, params, suffix=""):
    cat, coords = load_ref_cat_coords(objname, filtname)

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

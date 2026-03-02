import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.segmentation import detect_sources, SourceCatalog
from photutils.psf import PSFPhotometry, SourceGrouper, MoffatPSF
from photutils.background import LocalBackground

# ------------------------------------------------------------
# Main pipeline
# ------------------------------------------------------------

def run_pipeline(imgname, weightname, flagsname,
                 threshold_sigma=5,
                 npixels=20,
                 fit_shape=(401, 401),
                 ):

    # --------------------------------------------------------
    # Read data
    # --------------------------------------------------------

    image = fits.getdata(imgname)
    weight = fits.getdata(weightname)
    flags = fits.getdata(flagsname)

    error = 1.0 / np.sqrt(weight)

    # Example flag mask (same logic as your bitmask)
    flag_sat_extended = 8
    mask = (flags & flag_sat_extended) > 0

    image = np.array(image, dtype=float)
    image[mask] = np.nan

    # --------------------------------------------------------
    # Background statistics
    # --------------------------------------------------------

    mean, median, std = sigma_clipped_stats(image, mask=mask)
    threshold = median + threshold_sigma * std

    # --------------------------------------------------------
    # Source detection (segmentation)
    # --------------------------------------------------------

    segm = detect_sources(image, threshold, npixels=npixels)
    catalog = SourceCatalog(image, segm)

    positions = np.transpose((catalog.xcentroid,
                              catalog.ycentroid))

    init_table = Table()
    init_table['x_0'] = positions[:, 0]
    init_table['y_0'] = positions[:, 1]

    # --------------------------------------------------------
    # PSF model definition
    # --------------------------------------------------------

    psf_model = MoffatPSF(
    )

    # Optional parameter bounds
    psf_model.alpha.bounds = (0.01, 30) # size
    psf_model.beta.bounds = (0.5, 2) # powerlaw

    fitter = LevMarLSQFitter()

    # Optional: grouping for overlapping sources
    grouper = SourceGrouper(min_separation=20)

    # Optional: local background estimation
    local_bkg = LocalBackground(inner_radius=10, outer_radius=20)

    phot = PSFPhotometry(
        psf_model=psf_model,
        fit_shape=fit_shape,
        #fitter=fitter,
        grouper=grouper,
        #localbkg_estimator=local_bkg,
        aperture_radius=40,
        xy_bounds = 10,
    )

    # --------------------------------------------------------
    # Perform PSF photometry
    # --------------------------------------------------------

    print(len(init_table), "saturated groups")
    result_tab = phot(
        image,
        init_params=init_table,
        error=error,
        mask=mask
    )
    print("results completed")
    print(result_tab)

    # --------------------------------------------------------
    # Model and residual images
    # --------------------------------------------------------

    model_image = phot.make_model_image(image.shape, psf_shape=fit_shape)
    residual_image = phot.make_residual_image(image, psf_shape=fit_shape)

    return model_image, residual_image, result_tab


# ------------------------------------------------------------
# Script-style execution
# ------------------------------------------------------------

if __name__ == "__main__":

    imgname = "nobkg.fits"
    weightname = "nobkg.weight.fits"
    flagsname = "flag.fits"

    model_img, resid_img, catalog = run_pipeline(
        imgname,
        weightname,
        flagsname
    )

    fits.writeto("bright_model.fits", model_img, overwrite=True)
    fits.writeto("bright_residual.fits", resid_img, overwrite=True)
    catalog.write("bright_catalogue.fits", overwrite=True)

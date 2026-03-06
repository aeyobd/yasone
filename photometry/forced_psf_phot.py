from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from photutils import background
from photutils import psf

import astropy
from astropy.nddata import CCDData
from astropy.table import Table
from astropy import units as u
from astropy.visualization import simple_norm

import sys
sys.path.append(".")
from forced_photometry import PhotometryParams, load_image, load_ref_cat_coords, subtract_background, to_image_coords
from yasone.plotting import show_image



def main(objname, filtname, params, suffix=""):
    for imgpath in Path(f"{objname}").glob(f"img_{filtname}_*"):
        params.outdir = imgpath
        imgname = imgpath.stem[-2:]
        phot_result = derive_psf_phot(objname, filtname, imgname, params)
        break


################ Analysis 
def get_inbounds_filter(cat, img, params):
    psf_size = params.psf_size
    
    cat_coords = img_coadd.wcs.pixel_to_world(cat["X_IMAGE"]-1,
                                              cat["Y_IMAGE"]-1)
    x_cat_all, y_cat_all = img.wcs.world_to_pixel(cat_coords)

    filt_inbounds = x_cat_all > 0 + psf_size/2
    filt_inbounds &= x_cat_all < img.shape[1] - psf_size/2 - 1
    filt_inbounds &= y_cat_all > 0 + psf_size/2
    filt_inbounds &= y_cat_all < img.shape[0] - psf_size/2 - 1
    filt_inbounds &= cat["FLAGS"] < 64
    filt_inbounds &= cat["IMAFLAGS_ISO"] < 4

    x_cat = x_cat_all[filt_inbounds]
    y_cat = y_cat_all[filt_inbounds]
    return x_cat, y_cat


def get_bright_stars(cat_frame, img,  params):
    psf_size = params.psf_size
    hsize = (psf_size - 1) / 2
    x = cat_frame['X_IMAGE'] - 1
    y = cat_frame['Y_IMAGE'] - 1
    mask_edge = ((x > hsize) & (x < (img.shape[1] - 1 - hsize)) &
            (y > hsize) & (y < (img.shape[0] - 1 - hsize)))


    filt_good = cat_frame["FLAGS"] == 0
    filt_good &= cat_frame["FLAGS_WEIGHT"] == 0
    filt_good &= cat_frame["FWHM_IMAGE"] < 6
    filt_good &= cat_frame["ELLIPTICITY"] < 0.3

    filt_good &= cat_frame["MAGERR_APER"][:, 3] < 0.01
    filt_good &=  mask_edge

    good_stars_tbl = Table()
    good_stars_tbl["x"] = x[filt_good]
    good_stars_tbl["y"] = y[filt_good]
    return good_stars_tbl


def fit_ana_psf(epsf, params):

    psf_model = psf.MoffatPSF()
    psf_model.alpha.bounds = [0, 30] # radius
    psf_model.beta.bounds = [0.5, 6] # power
    psf_model.alpha.fixed = False
    psf_model.beta.fixed = False


    shape = epsf.data.shape
    psfphot = psf.PSFPhotometry(psf_model, shape, aperture_radius=5,  xy_bounds=5)
    init_params = Table()
    init_params["x"] = [shape[0]/2]
    init_params["y"] = [shape[1]/2]

    result = psfphot(epsf.data, init_params=init_params) 


    psf_model_fit = psf.MoffatPSF(alpha =
                                  result["alpha_fit"][0]/params.psf_oversample, beta=result["beta_fit"][0])
    return psf_model_fit
    

def build_psf_model(img_nobkg, cat_frame, params):
    good_stars_tbl = get_bright_stars(cat_frame, img_nobkg, params)
    stars = psf.extract_stars(img_nobkg, good_stars_tbl, size=params.psf_size)

    epsf_builder = psf.EPSFBuilder(oversampling=params.psf_oversample, maxiters=5,
                               progress_bar=False)
    epsf, fitted_stars = epsf_builder(stars)


    psf_model_ana = fit_ana_psf(epsf, params)
    psf_model_ana.fixed["alpha"] = True
    psf_model_ana.fixed["beta"] = True

    return epsf, stars



def derive_psf_phot(obj, filt, imgname, params):
    img = load_image(obj, filt, imgname)

    cat, coords = load_ref_cat_coords(objname, filtname)
    img_coords, index = to_image_coords(img, cat, coords, pad=params.psf_size)

    img_nobkg = subtract_background(img, params)

    # use the SE on the image to build quality sample
    cat_frame = Table.read(f"./{obj}/img_{filt}_{imgname}/detection.cat", hdu=2)
    epsf, fitted_stars = build_psf_model(img_nobkg, cat_frame, params)

    psf_ana = fit_ana_psf(epsf, params)

    init_params = Table()
    init_params["x"] = [c[0] for c in img_coords]
    init_params["y"] = [c[1] for c in img_coords]

    localbkg_estimator  = background.LocalBackground(params.local_background_r_in, 
                                       params.local_background_r_out)

    psfphot = psf.PSFPhotometry(epsf, 
        (params.psf_size, params.psf_size),
        aperture_radius=params.psf_aperture_radius,  xy_bounds=params.xy_bounds,
        progress_bar=True,
        localbkg_estimator=localbkg_estimator
    )
    phot_result = psfphot(img_nobkg, init_params = init_params,
                   error=img_nobkg.uncertainty) 


    psfphot_ana = psf.PSFPhotometry(psf_ana,
        (params.psf_size, params.psf_size),
        aperture_radius=params.psf_aperture_radius,  xy_bounds=params.xy_bounds,
        progress_bar=True,
        localbkg_estimator=localbkg_estimator
    )

    phot_result_ana = psfphot_ana(img_nobkg, init_params = init_params,
                                error=img_nobkg.uncertainty,)
                                
    

    for col in phot_result_ana.colnames:
        phot_result[col + "_ana"] = phot_result_ana[col]

    phot_result["index"] = index

    phot_result.write(params.outdir / "forced_psf_phot.fits", overwrite=True)

    all_results = {
        "img": img,
        "img_nobkg": img_nobkg,
        "epsf": epsf,
        "stars": fitted_stars,
        "phot_result": phot_result,
        "phot_result_ana": phot_result_ana,
        "psfphot": psfphot,
        "psfphot_ana": psfphot_ana,
    }

    plot_results(all_results, params)

    return 





################ Plotting 
def plot_cutouts(cutouts, params):
    nrows = 10
    ncols = 5
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20),
                           squeeze=True)
    ax = ax.ravel()
    idxs = np.random.choice(np.arange(len(cutouts)), nrows*ncols)
    for i in range(nrows * ncols):
        norm = simple_norm(cutouts[idxs[i]], 'log', percent=99.0)
        ax[i].imshow(cutouts[idxs[i]], norm=norm, origin='lower', cmap='viridis')

    savefig("cutouts", params)



def plot_results(all_results, params):

    if not (params.outdir / "figures").is_dir():
        (params.outdir / "figures").mkdir()

    epsf = all_results["epsf"]
    psfphot = all_results["psfphot"]
    psfphot_ana = all_results["psfphot_ana"]
    img = all_results["img"]
    img_nobkg = all_results["img_nobkg"]

    plot_cutouts(all_results["stars"], params)

    show_image(epsf.data, log=True, figsize=(5, 5), dpi=100, clim=(0, None))
    savefig("epsf", params)

    #show_image(psfphot.make_model_image(epsf.data.shape) / u.adu, dpi=100, log=True, clim=(None, None))
    #savefig("anapsf_model", params)

    #show_image(psfphot.make_residual_image(epsf.data), dpi=100, cmap="RdBu")
    #savefig("anapsf_model_residual", params)

    show_image(img_nobkg, log=True)
    phot = all_results["phot_result"]
    plt.scatter(phot["x_fit"], phot["y_fit"], lw=0.5, s=10, color="none", ec="C1")
    savefig("bkg_subtracted", params)

    show_image(img.data - img_nobkg.data, log=True)
    savefig("bkg", params)

    img_model = psfphot.make_model_image(img.shape)
    show_image(img_model.data, log=True)
    savefig("epsf_model_image", params)

    img_model = psfphot_ana.make_model_image(img.shape)
    show_image(img_model.data, log=True)
    savefig("anapsf_model_image", params)

    img_resid = psfphot_ana.make_residual_image(img_nobkg)
    show_image(img_resid, log=True, clim=(2000, 0))
    savefig("anapsf_residual_image", params)

    img_resid = psfphot.make_residual_image(img_nobkg)
    show_image(img_resid, log=True, clim=(2000, 0))
    savefig("epsf_model_image", params)

    fig, axs = plt.subplots(1, 2, figsize=(5, 2))
    plt.sca(axs[0])
    plt.hist(phot["x_fit"] - phot["x_init"], 100)
    plt.xlabel("xfit - x0")

    plt.sca(axs[0])
    plt.hist(phot["y_fit"] - phot["y_init"], 100)
    plt.xlabel("yfit - y0")
    savefig("psf_centroid_shifts", params)
 

def savefig(figstem, params):
    plt.savefig(params.outdir / "figures" / f"{figstem}.jpg")


if __name__ == "__main__":
    objname = sys.argv[1]

    params = PhotometryParams()
    for filtname in ["g", "r", "i"]:
        main(objname, filtname, params)

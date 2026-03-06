from dataclasses import dataclass, asdict
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

from shapely.geometry import Point
from shapely import Polygon

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats, mad_std
from astropy.table import Table

import arya
from . import analysis
from . import read_iso



isos_fe_h = ["p0.00", "m0.25", "m0.50", "m0.75", "m1.00", "m1.25", "m1.50", "m1.75", "m2.00", "m2.50", "m3.00"]

isonames = [f"../MIST/MIST_v1.2_vvcrit0.4_SDSSugriz/MIST_v1.2_feh_{iso_fe_h}_afe_p0.0_vvcrit0.4_SDSSugriz.iso.cmd" for iso_fe_h in isos_fe_h]
isochrones = None


BKG_KWARGS = dict(s=3, lw=0, color="k", label="all")
SELECTED_KWARGS = dict(s=9, lw=0, color="C3", label="selected")

def load_isochrones():
    global isochrones 
    isochrones = {fe_h: read_iso.ISOCMD(name) for fe_h, name in zip(isos_fe_h, isonames)}



@dataclass
class filter_params:
    "the name of the object (yasone1, yasone2, yasone3)"
    objname: str = ""

    """The size of circles to use for selecting central sources """
    r_cen: float = 25/60

    """The maximum distance from the centre in xi and eta, used for some
    background/plotting statistics"""
    xymax: float = 3

    "The maximum value of the `R_FLAGS` column"
    flags_max: int = 15
    "The maximum value of the `R_FLAGS_WEIGHT` column"
    flags_weight_max: int = 3
    """The maximum value of the `R_IMAFLAGS_ISO` column. By default, 
    filters out bad pixels (flag 16), and saturated regions (flags 8 and 4)"""
    flags_iso_max: int = 3
    "The maximum value of the `R_FLAGS_WEIGHT` column"
    detection_sigma: float = 1.5

    """The E(B-V) extinction for the isochrone selection. 
    Note that catalogues are usually already extinction corrected"""
    E_BV: float = 0

    "Extra reddining in the i band?"
    A_i_extra: float = 0

    """The distance modulus for the isochrone selection"""
    dm: float = 0

    """The intrinsic width of the isochrone selection"""
    iso_width: float = 0.2

    """The metallicity of the isochrone. Should be one of `isos_fe_h`"""
    iso_fe_h: str = 0

    """The log10 of the isochrone age in years"""
    iso_log_age: float = 0

    """Scale the colour uncertainties for the isochrone selection""" 
    err_scale: float = 1

    """The mode to filter stars on
    - "n": no magnitude based filtering
    - "r": exists in r-band
    - "gr": exists in g & r-band, passes gr-isochrone
    - "ri": exists in r & i-band, passes ri-isochrone
    - "gri": exists all three bands and passes gr and ri-isochrones.
    """ 
    mode: str = "ri"

    """The column of the catalogue to use for g-mag. 
    The column+_ERR must also exist"""
    gcol:str = "G_MAG"

    """The column of the catalogue to use for r-mag. 
    The column+_ERR must also exist"""
    rcol: str = "R_MAG"

    """The column of the catalogue to use for i-mag. 
    The column+_ERR must also exist"""
    icol: str = "I_MAG"


    # star galaxy methods
    "The maximum ellipticity (R_ELLIPTICITY) of selected members"
    ellipticity_max: float = None

    "The maximum half-light radius (R_FLUX_RADIUS) of selected members"
    r_half_max: float = None

    "The maximum FWHM (R_FWHM_IMAGE) of selected members"
    fwhm_max: float = None

    "The maximum MAG_23 parameter (i.e. 5px - 3px magnitudes) in R-band"
    r23_max: float = None

    "Set r23_max to be this sigma above the median for well-measured stars "
    r23_max_sigma: float = None

    
    "the position of a background circle in arcminutes"
    xi_bkg: float = 1

    "the position of a background circle in arcminutes"
    eta_bkg: float = -1

    gr_err = None
    ri_err = None


def calibrate_mag_col(cat, col, mag_min = 19, mag_max=22):
    """
    TODO: use function for filtering by flags / mask
    """
    for filt in ["G", "R", "I"]:
        filt_good = cat[f"{filt}_MAG"] > mag_min
        filt_good &= cat[f"{filt}_MAG"] < mag_max
        filt_good &= ~cat[f"{filt}_MAG"].mask
        filt_good &= cat[f"{filt}_FLAGS"] == 0
        if hasattr(cat[f"{filt}_{col}"], "mask"):
            filt_good &= ~cat[f"{filt}_{col}"].mask

        residual = np.ma.median((cat[f"{filt}_{col}"] - cat[f"{filt}_MAG"])[filt_good])

        print(residual)
        cat[f"{filt}_{col}"] -= residual




def filter_not_bright_point(cat):
    filt_good = cat["R_FLUX_RADIUS"] < 3
    filt_good &= cat["R_MAG"] < 21
    filt_good &= cat["R_ELLIPTICITY"] < 0.3
    return filt_good


def derive_threshhold(flux_param, cat, sigma=3):
    filt = filter_not_bright_point(cat)

    μ_flux, med_flux, std_flux = sigma_clipped_stats(flux_param[filt], stdfunc="mad_std", sigma=3)

    return med_flux + sigma*std_flux


def filt_finite(col):
    filt = np.isfinite(col)
    if hasattr(col, "mask"):
        filt &= ~col.mask
        filt.mask &= ~col.mask

    return filt




def sample_point(poly):
    minx, miny, maxx, maxy = poly.bounds

    while True:
        x = random.uniform(minx, maxx)
        y = random.uniform(miny, maxy)
        p = Point(x, y)
        if poly.contains(p):
            return p

def to_polygon(flat_coords, coord0):
    ra = flat_coords[::2]
    dec = flat_coords[1::2]
    coords = SkyCoord(ra, dec, unit=u.degree)
    xi, eta = analysis.to_tangent(coords, coord0)

    points = [p for p in zip(xi*60, eta*60)]
    return Polygon(points)


def get_region_mask(obj):
    if obj == "yasone1":
        region_mask = [265.442,13.208, 265.485,13.206, 265.487,13.194, 265.518,13.195, 265.534,13.214, 265.524,13.228, 265.508,13.240, 265.491,13.228, 265.456,13.218, 265.442,13.213]

        region_poly = [265.446,13.118, 265.591,13.120, 265.592,13.250, 265.440,13.248]
    elif obj == "yasone2":
        region_poly = [262.276,6.370, 262.419,6.372, 262.420,6.502, 262.272,6.500]
        region_mask = None

    elif obj == "yasone3":
        region_poly = [292.8815,-26.5033, 293.0386,-26.5005, 293.0401,-26.3710, 292.8771,-26.3716]
        region_mask = None
    else:
        region_poly = None
        region_mask = None

    return region_mask, region_poly


def get_sample_region(radius, mask, poly, coord0):
    safe_region = to_polygon(poly, coord0)
    safe_region = safe_region.buffer(-radius)

    if mask is not None:
        mask_poly = to_polygon(mask, coord0)
        safe_region = safe_region.difference(mask_poly.buffer(radius))

    centre = Point(0,0).buffer(2*radius)
    safe_region = safe_region.difference(centre)
    return safe_region


def sample_from_poly(poly):
    minx, miny, maxx, maxy = poly.bounds

    while True:
        x = np.random.uniform(minx, maxx)
        y = np.random.uniform(miny, maxy)
        p = Point(x, y)
        if poly.contains(p):
            return p.x, p.y


def count_centre(cat, radius, xi0=0, eta0=0):
    return np.sum((cat["xi"] - xi0)**2 + (cat["eta"] - eta0)**2 < radius**2)


def count_not_centre(cat, radius, xymax=3, xi0=0, eta0=0):
    filt = (cat["xi"] - xi0)**2 + (cat["eta"] - eta0)**2 >= radius**2
    filt &= cat["xi"] < xymax
    filt &= cat["xi"] > -xymax
    filt &= cat["eta"] < xymax
    filt &= cat["eta"] > -xymax

    return np.sum(filt)



def counts_to_stats(Ntot, Ncen, area_scale):
    Ntot_normed = Ntot * area_scale
    Ncen_err = np.sqrt(Ncen)
    Ntot_normed_err = np.sqrt(Ntot) * area_scale

    num_err = np.sqrt(Ntot_normed_err**2 + Ncen_err**2)

    sigma = (Ncen - Ntot_normed) / Ntot_normed
    sigma_err = np.sqrt((Ncen_err / Ntot_normed)**2 + (Ncen/Ntot_normed**2 * Ntot_normed_err)**2)

    return sigma, sigma_err, Ntot, Ncen


def inner_overdensity(cat, radius, xymax=3):
    area_tot = (2*xymax)**2
    area_cen = np.pi * radius**2
    area_scale = area_cen / (area_tot - area_cen)

    Ntot = count_not_centre(cat, radius, xymax=xymax)
    Ncen = count_centre(cat, radius)

    return counts_to_stats(Ntot, Ncen, area_scale)





def make_polygon(iso, mag_err, dm, A_b, A_r, iso_width=0.05, b="SDSS_g", r="SDSS_r", err_scale=1):
    filt = iso["phase"] < 3

    x = iso[b][filt].data - iso[r][filt].data + A_b - A_r
    y = iso[b][filt].data + dm + A_b

    x_min = x - err_scale * mag_err(y) - iso_width
    x_max = x + err_scale * mag_err(y) + iso_width

    return np.hstack([x_min, x_max[::-1], x_min[0]]), np.hstack([y, y[::-1], y[0]])





def select_gr(cat, params):
    iso = isochrones[params.iso_fe_h][params.iso_log_age]
    gr_err = get_gr_err(cat, params)

    A_g, A_r, A_i = get_extinction(params)
    x_poly_gr, y_poly_gr = make_polygon(iso, gr_err, params.dm, A_b=A_g, A_r=A_r, r="SDSS_r", b="SDSS_g", iso_width=params.iso_width, err_scale=params.err_scale)

    cmd_filt_gr = analysis.is_in_poly(cat[params.gcol].data -
                                       cat[params.rcol].data, cat[params.gcol].data, x_poly_gr, y_poly_gr)

    return cmd_filt_gr



def select_ri(cat, params):
    iso = isochrones[params.iso_fe_h][params.iso_log_age]
    ri_err = get_ri_err(cat, params)

    A_g, A_r, A_i = get_extinction(params)
    x_poly_gr, y_poly_gr = make_polygon(iso, ri_err, params.dm, 
                                        A_b=A_r, A_r=A_i, r="SDSS_i", b="SDSS_r", iso_width=params.iso_width, err_scale=params.err_scale)

    cmd_filt = analysis.is_in_poly(cat[params.rcol].data - cat[params.icol].data, cat[params.rcol].data, x_poly_gr, y_poly_gr)

    return cmd_filt





def filter_nans(col):
    filt = np.isfinite(col)
    if hasattr(col, "mask"):
        filt &= ~col.mask

    return filt


def get_extinction(params):
    A_g, A_r, A_i = analysis.get_extinction(params.E_BV)
    return A_g, A_r, A_i + params.A_i_extra


def get_flags_filter_band(cat, params, band):
    filt = np.full(len(cat), True, )

    if params.flags_max is not None:
        filt &= cat[f"{band}_FLAGS"] <= params.flags_max
    if params.flags_weight_max is not None:
        filt &= cat[f"{band}_FLAGS_WEIGHT"] <= params.flags_weight_max
    if params.flags_iso_max is not None:
        filt &= cat[f"{band}_IMAFLAGS_ISO"] <= params.flags_iso_max

    if band == "G":
        filt &= filter_nans(cat[params.gcol])
    elif band == "R":
        filt &= filter_nans(cat[params.rcol])
    elif band == "I":
        filt &= filter_nans(cat[params.icol])

    return filt

def get_flags_filter(cat, params):
    filt = get_flags_filter_band(cat, params, "R")

    if params.mode == "gr":
        filt &= get_flags_filter_band(cat, params, "G")
    elif params.mode == "ri":
        filt &= get_flags_filter_band(cat, params, "I")
        filt &= get_flags_filter_band(cat, params, "G")
    elif params.mode == "gri":
        filt &= get_flags_filter_band(cat, params, "I")
    elif params.mode == "n":
        filt = np.full(len(cat), True, )
    elif params.mode == "r":
        pass
    else:
        raise Exception(f"mode not known {params.mode}")
        return None

    return filt


def get_star_galaxy_filter(cat, params):
    filt = np.full(len(cat), True, )

    if params.ellipticity_max is not None:
        filt &= cat["R_ELLIPTICITY"] <= params.ellipticity_max
    if params.fwhm_max is not None:
        filt &= cat["R_FWHM_IMAGE"] <= params.fwhm_max
    if params.r_half_max is not None:
        filt &= cat["R_FLUX_RADIUS"] <= params.r_half_max

    if params.r23_max_sigma is not None:
        thresh_23 = derive_threshhold(cat["MAG_23"], cat,
                                      sigma=params.r23_max_sigma)
        filt &= cat["MAG_23"] < thresh_23
    elif params.r23_max is not None:
        thresh_23 = params.r23_max
        filt &= cat["MAG_23"] < params.r23_max

    return filt

def filter_catalog(cat, params):    
    filt = get_flags_filter(cat, params)
    filt &= get_star_galaxy_filter(cat, params)

    if params.detection_sigma is not None:
        filt &= cat["R_SNR"] > params.detection_sigma

    if params.mode == "gr":
        filt &= select_gr(cat, params)
    elif params.mode == "ri":
        filt &= select_ri(cat, params)
    elif params.mode == "gri":
        filt &= select_gr(cat, params)
        filt &= select_ri(cat, params)

    return cat[filt]




def background_sample_counts(cat, radius, xymax=3, N=1000):
    points = cat[["xi", "eta"]].to_pandas().to_numpy()
    tree = cKDTree(points)

    # sampling radii from outside the central region but within xy max
    xi = np.random.uniform(-xymax+radius, -1.5*radius, N) * np.random.choice([1, -1], N)
    eta = np.random.uniform(-xymax+radius, -1.5*radius, N) * np.random.choice([1, -1], N)

    centres = np.column_stack([xi, eta])
    counts = tree.query_ball_point(centres, r=radius, return_length=True)
    return counts



def catalog_stats(cat, 
        params
    ):

    cat_filtered = filter_catalog(cat,
        params
    )
    sigma, sigma_err, Ntot, Ncen = inner_overdensity(cat_filtered, params.r_cen, xymax=params.xymax)



    N_bkgs = background_sample_counts(cat_filtered, params.r_cen, xymax=params.xymax)

    area_tot = (2*params.xymax)**2
    area_cen = np.pi * params.r_cen**2
    area_scale = area_cen / (area_tot - area_cen)
    Nbkg_med = np.median(N_bkgs)
    Nbkg_std = mad_std(N_bkgs)

    delta = (Ncen - Nbkg_med) / Nbkg_med
    delta_err = np.sqrt((np.sqrt(Ncen) / Nbkg_med)**2 + (Ncen/Nbkg_med**2 * Nbkg_std)**2)

    pval = np.mean(N_bkgs > Ncen)

    return pd.DataFrame(dict(
        pval = pval,
        sigma = sigma, 
        sigma_err = sigma_err,
        delta = delta,
        delta_err = delta_err,
        Ntot = Ntot,
        Ntot_scaled = Ntot * area_scale,
        area_scale = area_scale,
        Ncen = Ncen,
        Nbkg_med = Nbkg_med,
        Nbkg_std = Nbkg_std,
        **asdict(params)

    ), index=[0])




def rand_overdensity(cat, radius, safe_region=None):
    xi, eta = sample_from_poly(safe_region)    

    Ncen = count_centre(cat, radius, xi0=xi, eta0=eta)
    return Ncen


def plot_2d_counts(cat, radius, xymax=5, N=100, vmin=None, vmax=None, thresh=None):
    x = np.linspace(np.min(cat["xi"]), np.max(cat["xi"]), N)
    y = np.linspace(np.min(cat["eta"]), np.max(cat["eta"]), N)

    xgrid, ygrid = np.meshgrid(x, y)
    counts = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            counts[i, j] = count_centre(cat, radius, xi0=xgrid[i, j], eta0=ygrid[i, j])

    level_step = np.maximum(1, np.ceil(np.max(counts) / 15))

    p = plt.contourf(x, y, counts, levels=np.arange(np.min(counts), np.max(counts)+level_step, level_step), vmin=vmin, vmax=vmax)

    if thresh is not None:
        plt.contour(x, y, counts, levels=[thresh], color="k", lw=0.3)

    ax =  plt.gca()
    plt.colorbar(p, label="counts within circle", fraction=0.046, pad=0.04)
    plt.xlabel(r"$\xi$ / arcmin")
    plt.ylabel(r"$\eta$ / arcmin")

    circ = plt.Circle((0, 0), radius, color="k", lw=0.2, fill=False, alpha=1)
    ax.add_artist(circ)

    plt.scatter(0, 0, color="k", s=0.3, lw=0)
    ax.set_aspect(1)
    ax.invert_xaxis()
    return counts


def get_coord0(objname):
    props = analysis.get_obs_props(objname)
    return SkyCoord(props["ra"], props["dec"], unit=u.deg)


def plot_with_hist(cat, params, vmax=None, vmin=None, **kwargs):
    fig, axs = plt.subplots(2, 4, figsize=(8, 4))

    cat_filtered = filter_catalog(cat, params)
    results = catalog_stats(cat, params)
    

    coord0 = get_coord0(params.objname)
    region_mask, region_poly = get_region_mask(params.objname)
    safe_region = get_sample_region(params.r_cen, region_mask, region_poly,
                                    coord0)

    background_counts = [rand_overdensity(cat_filtered, params.r_cen, safe_region=safe_region) for i in range(10_000)]
    Ncen = count_centre(cat_filtered, params.r_cen)

    plot_detection(cat, params, results=results.iloc[0, :], fig=fig, axs=axs,
              **kwargs)


    plt.sca(axs[1][-1])
    grid_counts = plot_2d_counts(cat_filtered, params.r_cen, vmin=vmin, vmax=vmax, thresh=Ncen)



    plt.sca(axs[0][-1])
    plt.hist(background_counts, histtype="step", bins=-0.5 + np.arange(0, np.max(background_counts)))
    plt.hist(grid_counts.flatten(), histtype="step", bins=-0.5 + np.arange(0, np.max(grid_counts.flatten())))
    plt.gca().set_box_aspect(1)
    plt.xlabel("N within random circle")
    plt.ylabel("counts")
    p = np.mean(Ncen <= background_counts)
    p2 = np.mean(Ncen <= grid_counts)
    plt.text(0.9, 0.9, f"p = {p}", transform=plt.gca().transAxes, fontsize=8, ha="right", va="bottom", color="C0")
    plt.text(0.9, 0.85, f"p = {p2}", transform=plt.gca().transAxes, fontsize=8, ha="right", va="top", color="C1")

    plt.axvline(Ncen)


def get_gr_err(cat, params):
    if params.gr_err is not None:
        return params.gr_err
    x = cat[params.gcol]
    xerr = np.sqrt(cat[params.gcol + "_ERR"]**2 + cat[params.rcol + "_ERR"]**2)

    return analysis.fit_err(x, xerr)



def get_ri_err(cat, params):
    if params.ri_err is not None:
        return params.ri_err

    x = cat[params.rcol]
    xerr = np.sqrt(cat[params.rcol + "_ERR"]**2 + cat[params.icol + "_ERR"]**2)

    return analysis.fit_err(x, xerr)


def select_subsets(cat, params):
    obs_props = analysis.get_obs_props(params.objname)
    ra0 = obs_props["ra"]
    dec0 = obs_props["dec"]

    coord0 = SkyCoord(ra=ra0, dec=dec0, unit=u.deg)
    coord1 = SkyCoord(ra=ra0 + params.xi / 60/np.cos(np.deg2rad(dec0)), dec=dec0 + params.eta / 60, unit=u.deg)


    cat_filtered = filter_catalog(cat, params)
    cat_okay = cat[get_flags_filter(cat, params)]

    ri_err = get_ri_err(cat_okay, params)
    gr_err = get_gr_err(cat_okay, params)

    df_cen_all = analysis.select_centre(cat_okay, coord0, params.r_cen * u.arcmin)
    df_cen = analysis.select_centre(cat_filtered, coord0, params.r_cen * u.arcmin)

    df_bkg = analysis.select_centre(cat_filtered, coord1, params.r_cen * u.arcmin)
    df_bkg_all = analysis.select_centre(cat_okay, coord1, params.r_cen * u.arcmin)

    filt_no = ~np.isin(cat["INDEX"], cat_filtered["INDEX"])
    return {
        "ri_err": ri_err,
        "gr_err": gr_err,
        "centre_all": df_cen_all,
        "centre_selected": df_cen,
        "background_all": df_bkg_all,
        "background_selected": df_bkg,
        "selected": cat_filtered,
        "unselected": cat[filt_no],
        "not_flagged": cat_okay,
    }

def plot_circle(x, y, r, **kwargs):
    circ = plt.Circle((x, y), r, **kwargs)
    plt.gca().add_patch(circ)


def plot_unselected_points(subsets, params):
    cat_filtered = subsets["not_flagged"]
    plot_tangent(cat_filtered, **BKG_KWARGS)

    plot_circle(0, 0, params.r_cen, color="C0", lw=0.3, fill=False,
                      zorder=-1, alpha=1)
    plot_circle(params.xi, params.eta, params.r_cen, color="C1", lw=0.3,
                      fill=False, zorder=-1, alpha=1)

    plt.title("not flagged")



def plot_tangent(cat, **kwargs):
    tangent_axis()
    plt.scatter(cat["xi"], cat["eta"], **kwargs)


def plot_selected_points(subsets, params):

    cat_filtered = subsets["selected"]

    plot_tangent(cat_filtered, **SELECTED_KWARGS)

    plot_circle(0, 0, params.r_cen, color="C0", lw=0.3, fill=False,
                      zorder=-1, alpha=1)
    plot_circle(params.xi, params.eta, params.r_cen, color="C1", lw=0.3,
                      fill=False, zorder=-1, alpha=1)

    plt.title("selected")




def plot_gr_background(subsets, params):
    df_bkg_all = subsets["background_all"]
    df_bkg = subsets["background_selected"]
    gr_err = subsets["gr_err"]

    gr_axis()
    plot_iso_gr(params, gr_err)

    x = df_bkg_all[params.gcol] - df_bkg_all[params.rcol]
    y = df_bkg_all[params.gcol]
    plt.scatter(x, y, **BKG_KWARGS)

    x = df_bkg[params.gcol] - df_bkg[params.rcol]
    y = df_bkg[params.gcol]
    plt.scatter(x, y, **SELECTED_KWARGS)
    plt.title("background")


def plot_gr_centre(subsets, params):
    gr_axis()
    df_cen_all = subsets["centre_all"]
    df_cen = subsets["centre_selected"]
    gr_err = subsets["gr_err"]

    plot_iso_gr(params, gr_err)
    coord0 = get_coord0(params.objname)

    x = df_cen_all[params.gcol] - df_cen_all[params.rcol]
    y = df_cen_all[params.gcol]
    plt.scatter(x, y, **BKG_KWARGS)

    x = df_cen[params.gcol] - df_cen[params.rcol]
    y = df_cen[params.gcol]
    plt.scatter(x, y, **SELECTED_KWARGS)

    plt.title("centre")


def plot_ri_background(subsets, params):
    df_bkg_all = subsets["background_all"]
    df_bkg = subsets["background_selected"]
    ri_err = subsets["ri_err"]

    ri_axis()
    plot_iso_ri(params, ri_err)
    x = df_bkg_all[params.rcol] - df_bkg_all[params.icol]
    y = df_bkg_all[params.rcol]
    plt.scatter(x, y, **BKG_KWARGS)

    x = df_bkg[params.rcol] - df_bkg[params.icol]
    y = df_bkg[params.rcol]
    plt.scatter(x, y, **SELECTED_KWARGS)
    plt.title("background")



def plot_ri_centre(subsets, params):
    df_cen_all = subsets["centre_all"]
    df_cen = subsets["centre_selected"]
    ri_err = subsets["ri_err"]

    ri_axis()
    plot_iso_ri(params, ri_err)

    x = df_cen_all[params.rcol] - df_cen_all[params.icol]
    y = df_cen_all[params.rcol]
    plt.scatter(x, y, **BKG_KWARGS)

    x = df_cen[params.rcol] - df_cen[params.icol]
    y = df_cen[params.rcol]
    plt.scatter(x, y, **SELECTED_KWARGS)

    plt.title("centre")


def plot_detection(cat, params,  results=None, fig=None, axs=None):

    subsets = select_subsets(cat, params)

    if fig is None:
        fig, axs = plt.subplots(2, 3, figsize=(6, 4))

    plt.sca(axs[0][0])
    plot_unselected_points(subsets, params)

    plt.sca(axs[1][0])
    plot_selected_points(subsets, params)

    plt.sca(axs[1][1])
    plot_gr_background(subsets, params)
    if results is not None:
        plt.text(0.1, 0.9, f"${results.Nbkg_med:0.2f} \\pm {results.Nbkg_std:0.2f}$", transform=plt.gca().transAxes,
                 color="C1", fontsize=8)


    plt.sca(axs[0][1])
    plot_gr_centre(subsets, params)
    if results is not None:
        plt.text(0.1, 0.9, f"${results.Ncen}$", transform=plt.gca().transAxes,
                 color="C0", fontsize=8)

    plt.sca(axs[1][2])
    plot_ri_background(subsets, params)

    plt.sca(axs[0][2])
    plot_ri_centre(subsets, params)


    plt.tight_layout()





def tangent_axis(tangent_bounds=None):
    plt.ylabel(r"$\eta/\textrm{arcmin}$")
    plt.xlabel(r"$\xi/\textrm{arcmin}$")
    plt.xticks(np.arange(-5, 6))
    plt.yticks(np.arange(-5, 6))

    if tangent_bounds is not None:
        plt.xlim(tangent_bounds[0], tangent_bounds[2])
        plt.ylim(tangent_bounds[1], tangent_bounds[3])

    plt.gca().invert_xaxis()
    plt.gca().set_aspect(1)



def ri_axis():
    plt.xlabel(r"$r-i$ (mag)")
    plt.ylabel(r"$r$ (mag)")

    plt.xlim(-1.2, 1.8)
    plt.ylim(27, 17)
    plt.gca().set_box_aspect(1)


def gr_axis():
    plt.xlabel(r"$g-r$ (mag)")
    plt.ylabel(r"$g$ (mag)")
    plt.xlim(-1, 2)
    plt.ylim(27, 17)
    plt.gca().set_box_aspect(1)


def plot_iso(iso, dm, A_b, A_r, mag_err=None, phase_max=3, b="SDSS_g", r="SDSS_r", iso_width=0):
    filt = iso["phase"] < phase_max
    plt.plot(iso[b][filt] - iso[r][filt] + A_b-A_r, iso[b][filt] + dm+A_b, color="grey", zorder=-1)

    if mag_err is not None:
        x_poly, y_poly = make_polygon(iso, mag_err, dm=dm, A_b=A_b, A_r=A_r, b=b, r=r, iso_width=iso_width)
        plt.plot(x_poly, y_poly, color="grey", alpha=0.5, zorder=-1)


def plot_iso_gr(params, gr_err=None):
    A_g, A_r, A_i = get_extinction(params)
    iso = isochrones[params.iso_fe_h][params.iso_log_age]
    plot_iso(iso, dm=params.dm, A_b=A_g, A_r=A_r, mag_err=gr_err, iso_width=params.iso_width)


def plot_iso_ri(params, ri_err=None):
    A_g, A_r, A_i = get_extinction(params)
    iso = isochrones[params.iso_fe_h][params.iso_log_age]

    plot_iso(iso, dm=params.dm, A_b=A_r, A_r=A_i, mag_err=ri_err, b="SDSS_r", r="SDSS_i", iso_width=params.iso_width)


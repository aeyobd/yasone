import shapely
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

import tomllib
from scipy.optimize import curve_fit



def is_in_poly(x, y, x_poly, y_poly):
    N = len(x)
    poly = shapely.Polygon(np.vstack([x_poly, y_poly]).T)
    
    in_poly = np.full(N, False)
    for i in range(N):
        if np.isfinite(x[i]) and np.isfinite(y[i]):
            point = shapely.Point(x[i], y[i])
            in_poly[i] = point.within(poly)

    return in_poly


def make_polygon(iso, mag_err, dm, A_b, A_r, iso_width=0.05, b="SDSS_g", r="SDSS_r"):
    filt = iso["phase"] < 3

    x = iso[b][filt].data - iso[r][filt].data + A_b - A_r
    y = iso[b][filt].data + dm + A_b

    x_min = x - mag_err(y) - iso_width
    x_max = x + mag_err(y) + iso_width

    return np.hstack([x_min, x_max[::-1], x_min[0]]), np.hstack([y, y[::-1], y[0]])


def get_obs_props(objname):
    obs_props = {}

    with open("object_properties.toml", "rb") as f:
        obs_props = tomllib.load(f)[objname]

    return obs_props


def get_coord0(objname):
    obs_props = get_obs_props(objname)
    ra0 = obs_props["ra"]
    dec0 = obs_props["dec"]
    coord0 = SkyCoord(ra0, dec0, unit=u.degree)

    return coord0


def to_tangent(coord, coord0):
    alpha = coord.ra.to("rad").value
    delta = coord.dec.to("rad").value
    alpha0 = coord0.ra.to("rad").value
    delta0 = coord0.dec.to("rad").value


    denom = np.sin(delta) * np.sin(delta0) + np.cos(delta) * np.cos(delta0) * np.cos(alpha - alpha0)

    xi_num = np.cos(delta) * np.sin(alpha - alpha0)
    eta_num = np.sin(delta) * np.cos(delta0) - np.cos(delta) * np.sin(delta0) * np.cos(alpha-alpha0) 

    xi = np.rad2deg(xi_num/denom)
    eta = np.rad2deg(eta_num / denom)
    return xi, eta



def get_extinction(A_V):
    R_V = 3.1 # equal to A(V)/E(B-V)
    E_BV = A_V / R_V
    # values for R(V) = 3.1 from Schlafly & Finkbeiner (2011)
    R_g = 3.303 
    R_r = 2.285
    R_i = 1.698
    
    A_g = R_g * E_BV
    A_r = R_r * E_BV 
    A_i = R_i * E_BV 
    
    return A_g, A_r, A_i


def get_mag_shift(objname):
    shifts = {}
    with open("../photometry/" + objname + "/panstarrs_shift.toml", "rb") as f:
        shifts = tomllib.load(f)
    return shifts



def read_catalogue(objname, filter_bad=True):
    cat = Table.read(f"../photometry/{objname}/allcolours.cat")

    if filter_bad:
        cat = cat[(cat["RFLAGS"] < 4 ) & (cat["RFLAGS_WEIGHT"] == 0) ]

    add_columns(cat, objname)
    return cat



def add_columns(cat, objname):
    xi, eta = to_tangent(to_coords(cat), get_coord0(objname))
    cat["xi"] = xi * 60 * u.arcmin
    cat["eta"] = eta * 60 * u.arcmin

    shifts = get_mag_shift(objname)
    cat["G_MAG"] -= shifts["g"]
    cat["R_MAG"] -= shifts["r"]
    cat["I_MAG"] -= shifts["i"]

    cat["GR"] = cat["G_MAG"] - cat["R_MAG"]
    cat["RI"] = cat["R_MAG"] - cat["I_MAG"]
    
    cat["GR_ERR"] = np.sqrt(cat["G_MAG_ERR"]**2 + cat["R_MAG_ERR"]**2)
    cat["RI_ERR"] = np.sqrt(cat["R_MAG_ERR"]**2 + cat["I_MAG_ERR"]**2)

    # add pointsource depths
    cat["G_SNR"] = 1/cat["G_MAG_ERR"] * np.log(10) / 2.5
    cat["R_SNR"] = 1/cat["R_MAG_ERR"] * np.log(10) / 2.5
    cat["I_SNR"] = 1/cat["I_MAG_ERR"] * np.log(10) / 2.5




def fit_err(mag, magerr):
    filt = np.isfinite(magerr)
    if hasattr(magerr, "mask"):
        filt &= ~magerr.mask

    popt, covt = curve_fit(err_model, mag[filt], np.log10(magerr[filt]))

    def f(x):
        return 10**err_model(x, *popt)
        
    return f

def err_model(x, a, b, c):
    return np.log10(np.maximum(1e-12, 10**(x * a + b) + c))




def add_flux_param(cat):
    i = 0
    j = 1
    k = 3
    l = 1
    
    x = catalogue_matched[f"R_MAG_APER_{i}"] - catalogue_matched[f"R_MAG_APER_{j}"]
    y = -catalogue_matched[f"R_MAG_APER_{k}"] + catalogue_matched[f"R_MAG_APER_{l}"]

    flux_param = y - 0.5*x

    cat["FLUX_PARAM"] = flux_param


def select_centre(cat, coord0, radius):
    filt = coord0.separation(to_coords(cat)) < radius
    return cat[filt]


def get_upper_radius(radius, radius_bkg):
    return np.sqrt(radius**2 + radius_bkg**2)


def select_annuli(cat, coord0, radius, radius_bkg):
    radius_upper = get_upper_radius(radius, radius_bkg)
    # print(radius**2, radius_upper**2 - radius_bkg**2)
    seps = coord0.separation(to_coords(cat))

    filt = seps >= radius_bkg
    filt &= seps < radius_upper
    return cat[filt]


def to_coords(cat):
    return SkyCoord(cat["R_ALPHA_J2000"], cat["R_DELTA_J2000"], unit=u.degree)


import shapely
import numpy as np
import pathlib

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

import tomllib
from scipy.optimize import curve_fit
from dustmaps.sfd import SFDQuery

def correct_dust(cat):
    coords = to_coords(cat)
    # schlegel provides E(B-V) magnitudes
    E_BV = SFDQuery()(coords)
    A_g, A_r, A_i = get_extinction(E_BV)
    print("median extinction")
    print("A_g", np.median(A_g))
    print("A_r", np.median(A_r))
    print("A_i", np.median(A_i))

    cat["A_g"] = A_g
    cat["A_r"] = A_r
    cat["A_i"] = A_i
    cat["G_MAG"] -= A_g
    cat["R_MAG"] -= A_r
    cat["I_MAG"] -= A_i

    return cat


def get_extinction(E_BV):
    """
    Return the extinction A(g), A(r), and A(i) as a tuple
    given the E(B-V) extinctions.

    Uses the conversion from Schlafly & Finkbeiner (2011)
    for SDSS filters assuming R(V) = 3.1
    """
    R_g = 3.303 
    R_r = 2.285
    R_i = 1.698
    
    A_g = R_g * E_BV
    A_r = R_r * E_BV 
    A_i = R_i * E_BV 
    
    return A_g, A_r, A_i



def is_in_poly(x, y, x_poly, y_poly):
    """
    Given a set of coordinates (x, y) and a polyong (x_poly, y_poly)
    returns a bool vector stating whether each coordinate is fcontained within
    the polygon.
    """
    N = len(x)
    poly = shapely.Polygon(np.vstack([x_poly, y_poly]).T)
    
    in_poly = np.full(N, False)
    for i in range(N):
        if np.isfinite(x[i]) and np.isfinite(y[i]):
            point = shapely.Point(x[i], y[i])
            in_poly[i] = point.within(poly)

    return in_poly


def make_polygon(iso, mag_err, dm, A_b, A_r, iso_width=0.05, b="SDSS_g", r="SDSS_r"):
    """
    Given an isochrone, the magnitude colour error as a function of magnitude, 
    the distance modulus, the extinction in a "b" and "r" band, and optionally,
    a specified isochrone width and bands, returns a polygon selecting a regoin
    around the isochrone of the given width plus the magnitude colour error,
    shifted by the distance modulus and extinctions.
    """
    filt = iso["phase"] < 3

    x = iso[b][filt].data - iso[r][filt].data + A_b - A_r
    y = iso[b][filt].data + dm + A_b

    x_min = x - mag_err(y) - iso_width
    x_max = x + mag_err(y) + iso_width

    return np.hstack([x_min, x_max[::-1], x_min[0]]), np.hstack([y, y[::-1], y[0]])


def get_obs_props(objname):
    """
    Retrieve the information dictionary in the "object_properties.toml" file
    for the specified object (yasone1, yasone2, yasone3).
    """
    obs_props = {}

    pwd = pathlib.Path(__file__).parent.resolve()

    with open(pwd / "object_properties.toml", "rb") as f:
        obs_props = tomllib.load(f)[objname]

    return obs_props


def get_coord0(objname):
    """
    Retrieve the (arXiv) centre coordinate of the specified object as 
    a SkyCoord object.
    """
    obs_props = get_obs_props(objname)
    ra0 = obs_props["ra"]
    dec0 = obs_props["dec"]
    coord0 = SkyCoord(ra0, dec0, unit=u.degree)

    return coord0


def to_tangent(coord, coord0):
    """
    Converts the coordinate(s) to the tangent plane using the 
    specified reference coordinate. Expects both arguments to be SkyCoords.
    """
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



def get_mag_shift(objname, catname, shiftname):
    """
    Retrieve the magnitude shifts correcting the zeropoint of a catalogue as a
    dictionary,
    either using the file
    `../photometry/<objname>/<catname>_panstarrs_shift.toml`,
    the file `../photometry/<objname>/<shiftname>.toml` if shiftname is
    specified.

    """
    shifts = {}
    if shiftname is None:
        path = "../photometry/" + objname + f"/{catname}_panstarrs_shift.toml"
    else:
        path = "../photometry/" + objname + f"/{shiftname}.toml"

    with open(path, "rb") as f:
        shifts = tomllib.load(f)
    return shifts



def read_catalogue(objname, filter_bad=False, catname="", shiftname=None,
                   deredden=False):
    """Read in the catalogue with the specified name
    from the photometry/<objname> directory (suffix .cat automatically added). 

    Can optionally change the magnitude shift dictionary `shiftname` and filter
    out poor quality stars `filter_bad`, and deredden the magnitudes. 

    Adds columns for xi, eta, GR, RI, the colour errors, an index, possibly
    MAG_23, and SNR.
    """

    cat = Table.read(f"../photometry/{objname}/{catname}.cat")

    if filter_bad:
        cat = cat[(cat["R_FLAGS"] < 4 ) & (cat["R_FLAGS_WEIGHT"] == 0) ]

    add_columns(cat, objname, catname, shiftname=shiftname)
    if deredden:
        correct_dust(cat)
    return cat



def add_columns(cat, objname, catname, shiftname=None):
    xi, eta = to_tangent(to_coords(cat), get_coord0(objname))
    cat["xi"] = xi * 60 * u.arcmin
    cat["eta"] = eta * 60 * u.arcmin

    shifts = get_mag_shift(objname, catname, shiftname)
    print(shifts)
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

    cat["INDEX"] = np.arange(len(cat))
    if "R_MAG_APER_2" in cat.colnames:
        cat["MAG_23"] =  cat["R_MAG_APER_2"] - cat["R_MAG_APER_3"]



def flatten_arrays(cat):
    """
    For multidimensional array columns stored in the catalogue, flattens
    them into columns of floats by suffixing the column name with integers.
    """
    cat_out = Table()

    for col in cat.colnames:
        if len(cat[col].shape) == 2:
            for i in range(cat_r[col].shape[1]):
                cat_out[col] = cat[col][:, i]
            else:
                cat_out[col] = cat[col]
    return cat_out


def fit_err(mag, magerr):
    """
    Given a list of magnitudes and some uncertainties, fits a exponential +
    constant model to the errors, returning a function.
    """

    filt = np.isfinite(magerr)
    filt &= np.isfinite(mag)
    if hasattr(mag, "mask"):
        filt &= ~mag.mask
    if hasattr(magerr, "mask"):
        filt &= ~magerr.mask

    popt, covt = curve_fit(err_model, mag[filt], magerr[filt])

    def f(x):
        return err_model(x, *popt)
        
    return f

def err_model(x, a, b, c):
    return np.maximum(1e-12, 10**(x * a + b) + c)


def _add_flux_param(cat):
    """
    An old flux parameter to identify extended sources. 
    Deprecated in favour of MAG_23 (R_MAG_APER_2 - R_MAG_APER_3)
    """
    i = 0
    j = 1
    k = 3
    l = 1
    
    x = cat[f"R_MAG_APER_{i}"] - cat[f"R_MAG_APER_{j}"]
    y = -cat[f"R_MAG_APER_{k}"] + cat[f"R_MAG_APER_{l}"]

    flux_param = y - 0.5*x

    cat["FLUX_PARAM"] = flux_param


def select_centre(cat, coord0, radius):
    filt = coord0.separation(to_coords(cat)) < radius
    return cat[filt]


def get_upper_radius(radius, radius_bkg):
    """
    Calculates the outer radius forming an annulus with the given 
    inner radius (radius_bkg) with the area radius.
    """
    return np.sqrt(radius**2 + radius_bkg**2)


def select_annuli(cat, coord0, radius, radius_bkg):
    """
    Selects stars from `cat` within a background annulus with inner radius
    `radius_bkg` around `coord0` enclosing the same area of a circle of radius
    `radius`
    """
    radius_upper = get_upper_radius(radius, radius_bkg)

    assert np.isclose(radius**2, radius_upper**2 - radius_bkg**2)

    seps = coord0.separation(to_coords(cat))

    filt = seps >= radius_bkg
    filt &= seps < radius_upper
    return cat[filt]


def to_coords(cat, ra="ra", dec="dec"):
    """
    Returns SkyCoord of the positions of the cataluge, assuming
    the given ra and dec columns and units of degrees.
    """
    return SkyCoord(cat["ra"], cat["dec"], unit=u.degree)


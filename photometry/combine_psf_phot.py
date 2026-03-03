from astropy.stats import mad_std
from astropy.table import Table, join

import numpy as np


from pathlib import Path

import sys

from phot_utils import to_mag


def load_all_catalogues(objname, filtname, suffix=""):
    cats = []
    imgnames = []
    for path in Path(f"./{objname}/").glob(f"img_{filtname}_*"):
        catname = path / f"forced_psf_phot{suffix}.fits"
        if not catname.is_file():
            continue
        cat = Table.read(catname)
        cats.append(cat)
        imgnames.append(path.stem)

    return cats, imgnames


def load_ref_cat(obj, filt):
    cat = Table.read(f"./{obj}/coadd_median_{filt}/detection.cat", hdu=2)
    cat["index"] = np.arange(len(cat))

    return cat


def combine_all_catalogues(cat_ref, cats, catnames):
    cat_combined = cat_ref.copy()

    for i in range(0, len(cats)):
        cat = cats[i].copy()
        cat.rename_columns(cat.colnames, [c + "_" + catnames[i] for c in cat.colnames])

        cat_combined = join(cat_combined, cat, join_type="left", keys_left="index", keys_right="index_" + catnames[i])

    return cat_combined


def main(objname, filtname, suffix=""):

    cats, catnames = load_all_catalogues(objname, filtname, suffix)

    cat_ref = load_ref_cat(objname, filtname)

    cat_all = combine_all_catalogues(cat_ref, cats, catnames)

    colnames = get_safe_columns(cats[0])
    cat_reduced = reduce_cat(cat_ref, cat_all, catnames, colnames)
    add_mag_columns(cat_reduced)

    return cat_reduced


def add_mag_columns(cat_reduced):
    for name in ["flux_fit", "flux_fit_ana"]:
        med = cat_reduced[f"{name}_median"]
        err =  cat_reduced[f"{name}_err"]
        N =  cat_reduced[f"{name}_count"]

        mag, magerr = to_mag(med, err)

        if "ana" in name:
            cat_reduced[f"MED_MAG_PSF_ANA"] = mag
            cat_reduced[f"MED_MAG_PSF_ANA_ERR"] = magerr
        else:
            cat_reduced[f"MED_MAG_PSF"] = mag
            cat_reduced[f"MED_MAG_PSF_ERR"] = magerr



def get_safe_columns(cat):
    return [name for name in cat.colnames if np.issubdtype(cat[name].dtype, np.floating)]

def get_columns(cat_all, catnames, col):
    return cat_all[[col + "_" + name for name in catnames]]


def reduce_cat(cat_ref, cat_all, catnames, ap_colnames):

    cat_reduced = cat_ref.copy()

    for col in ap_colnames:
        data = get_columns(cat_all, catnames, col).to_pandas()


        cat_reduced[col + "_median"] = np.ma.median(data, axis=1)
        cat_reduced[col + "_err"] = mad_std(data, axis=1)
        cat_reduced[col + "_std"] = np.ma.std(data, axis=1)
        cat_reduced[col + "_count"] = np.ma.sum(np.isfinite(data), axis=1)
    return cat_reduced


if __name__ == "__main__":
    assert len(sys.argv) in [2, 3]
    
    if len(sys.argv) == 2:
        objname = sys.argv[1]
        suffix = ""
    else:
        objname, suffix = sys.argv[1:]

    for filtname in ["g", "r", "i"]:
        cat = main(objname, filtname, suffix)


        cat.write(f"./{objname}/forced_{filtname}/psf{suffix}.cat",
                  format="fits", overwrite=True)



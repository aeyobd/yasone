from astropy.stats import mad_std
from astropy.table import Table, join

import numpy as np


from pathlib import Path

import sys



def load_all_catalogues(objname, filtname, suffix=""):
    cats = []
    imgnames = []
    for path in Path(f"./{objname}/").glob(f"img_{filtname}_*"):
        cat = Table.read(path / f"forced_ap_phot{suffix}.fits")
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


def get_ap_colnames(cat):
    ap_colnames = [col for col in cat.colnames if col.startswith("aperture")]
    return ap_colnames


def main(objname, filtname, suffix=""):

    cats, catnames = load_all_catalogues(objname, filtname, suffix)

    cat_ref = load_ref_cat(objname, filtname)

    cat_all = combine_all_catalogues(cat_ref, cats, catnames)

    ap_colnames = get_ap_colnames(cats[0])
    cat_reduced = reduce_cat(cat_ref, cat_all, catnames, ap_colnames)
    add_mag_columns(cat_reduced, ap_colnames)

    return cat_reduced


def add_mag_columns(cat_reduced, ap_colnames):
    for name in ap_colnames:
        if not name.endswith("_err"):
            i = name[-1]
            med = cat_reduced[f"aperture_sum_{i}_median"]
            err =  cat_reduced[f"aperture_sum_{i}_err"]
            N =  cat_reduced[f"aperture_sum_{i}_count"]

            cat_reduced[f"MED_MAG_APER_{i}"] = -2.5*np.log10(med)

            cat_reduced[f"MED_MAG_APER_{i}_ERR"] = 2.5/np.log(10) * err / med * np.sqrt(np.pi/2) / np.sqrt(N-1)



def get_columns(cat_all, catnames, col):
    return cat_all[[col + "_" + name for name in catnames]]


def reduce_cat(cat_ref, cat_all, catnames, ap_colnames):

    cat_reduced = cat_ref.copy()

    for col in ap_colnames:
        data = get_columns(cat_all, catnames, col).to_pandas()


        cat_reduced[col + "_median"] = np.nanmedian(data, axis=1)
        cat_reduced[col + "_err"] = mad_std(data, axis=1)
        cat_reduced[col + "_std"] = np.nanstd(data, axis=1)
        cat_reduced[col + "_count"] = np.sum(np.isfinite(data), axis=1)
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


        cat.write(f"./{objname}/forced_{filtname}/aperture{suffix}.cat",
                  format="fits", overwrite=True)



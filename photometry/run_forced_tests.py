import sys
from forced_photometry import main, PhotometryParams

params_dict = {
    "": PhotometryParams(),
    "_small_bkg": PhotometryParams(back_size=64),
    "_large_bkg": PhotometryParams(back_size=256),
    "_small_bkg_filt": PhotometryParams(back_filter_size=3),
    "_large_bkg_filt": PhotometryParams(back_filter_size=11),
}


if __name__ == "__main__":
    objname = sys.argv[1]

    for filtname in ["g", "r", "i"]:
        for name, params in params_dict.items():
            main(objname, filtname, params, suffix=name)

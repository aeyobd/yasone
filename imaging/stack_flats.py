import os
from pathlib import Path
from astropy.nddata import CCDData

import numpy as np

import ccdproc
from ccdproc import ImageFileCollection
from convenience_functions import combine_images, mad_std

def get_flats(foldername, filt):
    imgfiles = ImageFileCollection(foldername, glob_include=f"flat_{filt}_*.fits")
    return imgfiles


def combine_images_err(imgfiles):
    err = mad_std(np.stack([x / np.median(x) for x in imgfiles.data()]), axis=0)
    return CCDData(err, unit="adu")

bad_flats = {
    "20230708": ["flat_i_05.fits"],
    "20230709": ["flat_g_06.fits", "flat_g_07.fits", "flat_i_01.fits",
    "flat_i_03.fits"],
    "20230723": ["flat_r_04.fits", "flat_i_04.fits"],
    "20230816": []
}



def flat_scale(A):
    return 1 / np.median(A)


def main():
    foldernames = ["20230708", "20230709", "20230723", "20230816"]

    # clean old files
    for folder in foldernames:
        for filt in ["g", "r", "i"]:
            for suffix in ["", "_all", "_filtered"]:
                outname = Path(folder, f"flat_{filt}_stacked{suffix}.fits")
                if outname.is_file():
                    os.remove(outname)

    for folder in foldernames:
        print(f"reading bias images in {folder}")

        for filt in ["g", "r", "i"]:
            imgfiles = get_flats(folder + "/unbiased/", filt)
            print(f"read in image files for band {filt} ", imgfiles.files)

            files = imgfiles.files_filtered(include_path=True)
            flat_stacked = combine_images(files, scale=flat_scale)
            flat_stacked.write(Path(folder, f"flat_{filt}_stacked.fits"))
            flat_stacked_err = combine_images_err(imgfiles)
            flat_stacked_err.write(Path(folder, f"flat_{filt}_stacked_err.fits"))


if __name__ == "__main__":
    main()

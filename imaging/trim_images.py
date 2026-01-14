# This script trims all of the raw images in `./raw/` and puts trimmed
# images (assuming ADU)_ into the `trimmed` directory in the calibration folders 
# 20230708, 20230709, 20230723, and 20230816
# and the object folders yasone1, yasone2, yasone3a and yasone3b

import os
import tomllib
from pathlib import Path
from glob import glob


import ccdproc
from ccdproc import ImageFileCollection
from astropy.nddata import CCDData

# taken from the CCD/fits header, also used by sausero.
IMG_SELECTION = "[28:2030,230:2026]" 

def read_img_keys():
    with open("image_keys.toml", "rb") as f:
        img_keys = tomllib.load(f)

    return img_keys


def clean_old_files(img_keys):
    for folder, file_keys in img_keys.items():
        newdir = Path(folder + "/trimmed")

        if not newdir.is_dir():
            newdir.mkdir()

        for f in glob(str(newdir / "*.fits")):
            os.remove(f)


"""Trim a CCDData image based on the IMG_SELECTION"""
def trim_image(img):
    return ccdproc.trim_image(img, fits_section=IMG_SELECTION)


"""
Trims images in the `raw` folder based on 
the keys in `file_keys` placing new images in `newdir`
"""
def process_folder(file_keys, newdir):
    for newname, oldname in file_keys.items():
        oldpath = Path("raw") / oldname
        newpath = newdir / newname
        print(f"trimming {oldpath} => {newpath}")

        img = CCDData.read(oldpath, unit="adu")
        img_trimmed = trim_image(img)
        img_trimmed.write(newpath)


def process_all_images(img_keys):
    for folder, file_keys in img_keys.items():
        newdir = Path(folder + "/trimmed")
        process_folder(file_keys, newdir)


def main():
    img_keys = read_img_keys()
    clean_old_files(img_keys)
    process_all_images(img_keys)


if __name__ == "__main__":
    main()

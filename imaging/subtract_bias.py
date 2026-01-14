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

from script_utils import calibration_folder, read_img_keys


def get_bias(foldername):
    return CCDData.read(calibration_folder(foldername) + "/bias_stacked.fits")

def clean(img_keys):
    for foldername, file_keys in img_keys.items():
        folder = Path(foldername) / "unbiased"
        if not folder.is_dir():
            folder.mkdir()

        for name in file_keys.keys():
            oldpath = folder / name
            if oldpath.is_file():
                os.remove(oldpath)

"""
Trims images in the `raw` folder based on 
the keys in `file_keys` placing new images in `newdir`
"""
def process_folder(file_keys, foldername):
    bias = get_bias(foldername)
    for name in file_keys.keys():
        oldpath = foldername + "/trimmed/" + name
        newpath = foldername + "/unbiased/" + name

        print("subtracting bias: ", oldpath, " => ", newpath)
        img = CCDData.read(oldpath)
        img_new = ccdproc.subtract_bias(img, bias)
        img_new.write(newpath)


def process_all(img_keys):
    for foldername, file_keys in img_keys.items():
        process_folder(file_keys, foldername)


def main():
    img_keys = read_img_keys()
    clean(img_keys)
    process_all(img_keys)


if __name__ == "__main__":
    main()

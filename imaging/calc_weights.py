# This script trims all of the raw images in `./raw/` and puts trimmed
# images (assuming ADU)_ into the `trimmed` directory in the calibration folders 
# 20230708, 20230709, 20230723, and 20230816
# and the object folders yasone1, yasone2, yasone3a and yasone3b

import os
import tomllib
from pathlib import Path
from glob import glob

import numpy as np

import ccdproc
from ccdproc import ImageFileCollection
from astropy.nddata import CCDData

from script_utils import calibration_folder, read_img_keys


def get_flat(foldername, filt):
    return CCDData.read(calibration_folder(foldername) + f"/flat_{filt}_stacked.fits")

def get_flat_err(foldername, filt):
    return CCDData.read(calibration_folder(foldername) +
                        f"/flat_{filt}_stacked_err.fits")

def get_bias_err(foldername):
    return CCDData.read(calibration_folder(foldername) + f"/bias_stacked_err.fits")



"""
Trims images in the `raw` folder based on 
the keys in `file_keys` placing new images in `newdir`
"""
def process_folder(file_keys, foldername):
    for name in file_keys.keys():
        if name.startswith("bias"):
            continue

        origpath = foldername + "/unbiased/" + name
        imgpath = foldername + "/flat_fielded/" + name
        newpath = foldername + "/flat_fielded/" + name.replace(".fits",
        "_err.fits")
        
        if Path(newpath).is_file():
            os.remove(newpath)

        img = CCDData.read(imgpath)
        filt = img.header["filter2"].split("Sloan_")[1]
        gain = img.header["gain"]

        counts = CCDData.read(origpath).data / gain

        flat = get_flat(foldername, filt).data
        flat_err = get_flat_err(foldername, filt).data
        bias_err = get_bias_err(foldername).data

        var_count = counts / flat**2
        var_bias = bias_err**2 / flat**2
        var_img = (flat_err / flat)**2 * (img.data)**2

        sigma2 = var_count + var_bias + var_img
        print("count: ", np.median(var_count))
        print("bias: ", np.median(var_bias))
        print("flat: ", np.median(var_img))
        CCDData(np.sqrt(sigma2), unit="adu").write(newpath)



def process_all(img_keys):
    for foldername, file_keys in img_keys.items():
        process_folder(file_keys, foldername)


def main():
    img_keys = read_img_keys()
    process_all(img_keys)


if __name__ == "__main__":
    main()

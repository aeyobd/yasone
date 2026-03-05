import tomllib
import os
from pathlib import Path

import numpy as np
from astropy.nddata import CCDData
from astropy.stats import mad_std
import ccdproc


def safe_rm(pathname):
    """Removes the pathname (string or pathlib.Path)
    only if the object exists as a file or symlink.
    """
    if isinstance(pathname, str):
        pathname = Path(pathname)
    if pathname.is_file() or pathname.is_symlink():
        os.remove(pathname)


def safe_mkdir(pathname):
    """Makes the specified directory 
    if not already present. 
    Takes a path as a string or pathlib.Path.
    """
    if isinstance(pathname, str):
        pathname = Path(pathname)

    if pathname.is_dir():
        return
    
    pathname.mkdir()

def read_img_keys():
    with open("image_keys.toml", "rb") as f:
        img_keys = tomllib.load(f)

    return img_keys



def calibration_folder(foldername):
    """retrieve the name of the folder 
    containing the appropriate calibration images
    for yasone 1, yasone2, yasone3a or b. 
    If the foldername is a calibration folder, returns the name
    """
    if foldername in ["20230708", "20230709", "20230723", "20230816"]:
        return foldername

    return {
        "yasone1": "20230708",
        "yasone2": "20230709",
        "yasone3a": "20230723",
        "yasone3b": "20230816",
    }[foldername]


def combine_images_err(imgfiles):
    err = mad_std(np.stack([x for x in imgfiles.data()]), axis=0)
    return CCDData(err, unit="adu")


"""
    combine_images(filenames, **kwargs)

Combines images using ccdproc. Expects CCDData and uses
median combination with median
"""
def combine_images(filenames, method="median", **kwargs):
    return ccdproc.combine(filenames, 
       method=method,
       sigma_clip=False, 
       mem_limit=1e9,
       **kwargs
        )


"""
    combine_images_average(filenames, **kwargs)

Combines images using ccdproc. Expects CCDData and uses
median combination with MAD_STD sigma clipping at 5 sigma.
"""
def combine_images_average(filenames, method="average", **kwargs):
    return ccdproc.combine(filenames, 
       method=method,
       sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, 
       sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, 
       mem_limit=1e9,
       **kwargs
        )

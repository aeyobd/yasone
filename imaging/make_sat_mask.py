from skimage import measure, morphology
from scipy import ndimage
from astropy.nddata import CCDData
import numpy as np

from script_utils import calibration_folder, read_img_keys
from pathlib import Path
import tomllib
import os

from trim_images import trim_image
FLAG_SAT = 4
FLAG_BRIGHTSAT = 8
FLAG_BADPIXEL = 16


def load_master_bpm():
    return trim_image(CCDData.read("master_calibration/BPM_OSIRIS_PLUS.fits",
                                   unit="adu"))

BPM = load_master_bpm()

def create_mask(img, sat_thresh=62_000):
    bpm_mask = BPM.data == 0.
    sat_mask = img.data > sat_thresh
    bright_sat_mask = create_bright_mask(sat_mask)

    flags = FLAG_SAT * sat_mask
    flags += FLAG_BRIGHTSAT * create_bright_mask(sat_mask) 
    flags += FLAG_BADPIXEL * (BPM.data == 0.)

    return flags


def create_bright_mask(sat_mask):
    # dialate to connect outline features
    mask = (ndimage.binary_dilation(sat_mask, iterations=4))
    
    mask = morphology.remove_small_holes(mask, max_size=20000)

    # dial back the dialation 
    mask = ndimage.binary_erosion(mask, iterations=3) 

    mask = mask | sat_mask # don't allow any previously masked pixels to remain

    mask = fill_small_convex_hulls(mask, 400)
    mask = remove_small_groups(mask, 16)

    return mask


def remove_small_groups(mask, min_size):
    labeled = measure.label(mask, connectivity=1)
    result = mask.copy()

    for region in measure.regionprops(labeled):
        if region.area < min_size:
            # Get convex hull of this region
            convex_hull = morphology.convex_hull_image(region.image)
            
            # Place convex hull back into result image
            minr, minc, maxr, maxc = region.bbox
            result[minr:maxr, minc:maxc][convex_hull] = False

    return result



def fill_small_convex_hulls(binary_image, max_size):
    """Fill convex hulls of True regions smaller than max_size"""
    # Label connected components
    labeled = measure.label(binary_image, connectivity=2)
    
    # Create output image
    result = binary_image.copy()
    
    # Process each region
    for region in measure.regionprops(labeled):
        if region.area < max_size:
            # Get convex hull of this region
            convex_hull = morphology.convex_hull_image(region.image)
            
            # Place convex hull back into result image
            minr, minc, maxr, maxc = region.bbox
            result[minr:maxr, minc:maxc][convex_hull] = True
    
    return result



"""
Trims images in the `raw` folder based on 
the keys in `file_keys` placing new images in `newdir`
"""
def process_folder(file_keys, foldername):
    for name in file_keys.keys():
        if ("std_" in name) or ("obj_" in name):
            oldpath = foldername + "/trimmed/" + name
            newpath = foldername + "/flat_fielded/" + name.replace(".fits",
            ".mask.fits")

            if Path(newpath).is_file():
                os.remove(newpath)

            print("creating mask ", oldpath, " => ", newpath)
            img = CCDData.read(oldpath)
            mask = np.asarray(create_mask(img), dtype=np.int16)
            img_new = CCDData(mask, unit="")
            img_new.write(newpath, overwrite=True, hdu_mask=None, hdu_uncertainty=None,)


def process_all(img_keys):
    for foldername, file_keys in img_keys.items():
        process_folder(file_keys, foldername)


def main():
    img_keys = read_img_keys()
    process_all(img_keys)


if __name__ == "__main__":
    main()

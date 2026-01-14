import os
from pathlib import Path

import ccdproc
from ccdproc import ImageFileCollection
from convenience_functions import combine_images

def get_biases(foldername):
    imgfiles = ImageFileCollection(foldername, glob_include="bias*.fits")
    return imgfiles


def main():
    foldernames = ["20230708", "20230709", "20230723", "20230816"]

    # clean old files
    for folder in foldernames:
        outname = Path(folder, "bias_stacked.fits")

        if outname.is_file():
            os.remove(outname)

    for folder in foldernames:
        print(f"reading bias images in {folder}")

        imgfiles = get_biases(folder + "/trimmed/")
        print("read in image files: ", imgfiles.files)

        files = imgfiles.files_filtered(include_path=True)
        bias_stacked = combine_images(files)

        bias_stacked.write(Path(folder, "bias_stacked.fits"))


if __name__ == "__main__":
    main()

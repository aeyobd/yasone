import sys
from pathlib import Path
from astropy.io import fits
import os

def main():
    if len(sys.argv) != 2:
        print("usage update_headers.py objectname")
        return

    objectname = sys.argv[1]

    for path in Path(objectname).glob(f"img_*/"):
        print("updating header in ", path)
        if (path / "detection.head").is_file():
            os.rename(path / "detection.head", path / "nobkg.head")

        header = fits.Header.fromtextfile(path / "nobkg.head")
        with fits.open(path / "nobkg.fits") as hdus:
            data = hdus[0].data
            hdus[0].header = header
            fits.writeto(path / "scamped.fits", data, header=header, overwrite=True)


if __name__ == "__main__":
    main()

from pathlib import Path
import os
import sys
import subprocess
from script_utils import safe_rm
import numpy as np
from astropy.nddata import CCDData 


def main():
    if len(sys.argv) != 3:
        print("usage run_swarp.py objectname filtername")
        return

    objectname, filtername = sys.argv[1:]

    indirs = [x for x in
              Path(objectname).glob(f"img_{filtername}_*/")]

    for (i, path) in enumerate(indirs):
        safe_rm(path / "nobkg.weight.fits")

    for (i, path) in enumerate(indirs):
        bkg_err = CCDData.read(path / "bkg_err.fits")
        weights = CCDData.read(path / "flat_fielded.weight.fits")
        img = CCDData.read(path / "nobkg.fits")
        
        sigma2 = bkg_err.data**2 + weights.data**2
        inv_var = 1 / (sigma2)
        sigma = np.sqrt(sigma2)
        tot_err = CCDData(inv_var, unit="")
        tot_err.write(path / "nobkg.weight.fits")

        f_bkg = np.median(bkg_err.data**2 / sigma2)
        print(path, f"frac bkg = {f_bkg:0.2f}")


if __name__ == "__main__":
    main()

from pathlib import Path
import shutil
import os
import sys
import subprocess
from script_utils import safe_rm


def main():
    if len(sys.argv) != 3:
        print("usage run_swarp.py objectname filtername")
        return

    objectname, filtername = sys.argv[1:]

    outdir = Path(objectname) / f"stacked_{filtername}"
    indirs = [x for x in
              Path(objectname).glob(f"img_{filtername}_*/flat_fielded-astrom.fits")]
    print("combining ", indirs, "into", outdir)

    for (i, path) in enumerate(indirs):
        safe_rm(path.parent / "with_bkg.head")
        safe_rm(path.parent / "with_bkg.fits")

        os.symlink("./nobkg.head", path.parent / "with_bkg.head")
        os.symlink("./flat_fielded-astrom.fits", path.parent / "with_bkg.fits")
        indirs[i] = "../.." / path.parent / "with_bkg.fits" # up two dirs since
        # we chdir
    if not outdir.is_dir():
        outdir.mkdir()
    else:
        [safe_rm(f) for f in outdir.glob("*")]

    os.chdir(outdir)

    subprocess.run(["swarp", "-c", "../../default.swarp"] + indirs)

if __name__ == "__main__":
    main()

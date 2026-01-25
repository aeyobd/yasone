from pathlib import Path
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
    indirs = ["../.." / x for x in
              Path(objectname).glob(f"img_{filtername}_*/flag.fits")]

    print("combining ", indirs, "into", outdir)

    os.chdir(outdir)

    for path in indirs:
        newpath = path.parent / "flag.head"
        if newpath.is_symlink():
            os.remove(newpath)
        os.symlink(path.parent / "nobkg.head", newpath)

    subprocess.run(["swarp", "-c", "../../flags.swarp"] + indirs)

if __name__ == "__main__":
    main()

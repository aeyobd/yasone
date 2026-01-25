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
              Path(objectname).glob(f"img_{filtername}_*/nobkg.fits")]
    print("combining ", indirs, "into", outdir)


    if not outdir.is_dir():
        outdir.mkdir()
    else:
        [os.remove(f) for f in outdir.glob("*")]

    os.chdir(outdir)

    #infiles = []
    #for (i, f) in enumerate(indirs):
    #    os.symlink(f, f"source_{i}.fits")
    #    os.symlink(f.parent / "nobkg.weight.fits", f"source_{i}.weight.fits")
    #    infiles.append(f"source_{i}.fits")


    subprocess.run(["swarp", "-c", "../../default.swarp"] + indirs)

if __name__ == "__main__":
    main()

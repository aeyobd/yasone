from pathlib import Path
import os
import sys
import subprocess


def main():
    if len(sys.argv) != 3:
        print("usage run_swarp.py objectname filtername")
        return

    objectname, filtername = sys.argv[1:]

    outdir = Path(objectname) / f"stacked_{filtername}"
    indirs = ["../.." / x for x in
              Path(objectname).glob(f"img_{filtername}_*/scamped.fits")]
    print("combining ", indirs, "into", outdir)

    if not outdir.is_dir():
        outdir.mkdir()
    else:
        [os.remove(f) for f in outdir.glob("*")]

    os.chdir(outdir)

    subprocess.run(["swarp", "-c", "../../default.swarp"] + indirs)

if __name__ == "__main__":
    main()

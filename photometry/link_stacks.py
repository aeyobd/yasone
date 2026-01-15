import os
from pathlib import Path


def get_newdir(file):
    program, prefix, system, filt, date, sky, suffix = file.stem.split("_")
    newname = f"julen_stack_{filt}_{prefix}"
    prefix = int(prefix)
    if prefix <= 3:
        destname = "yasone1"
    elif 3 < prefix <= 6:
        destname = "yasone2"
    elif 6 < prefix:
        destname = "yasone3"

    return Path(destname) / newname


def link_files():
    base_folder = Path(f"../imaging/julen_stacked")
    filenames = [f for f in (base_folder).glob(f"*.fits")]
    newdirs = [get_newdir(file) for file in filenames]

    for newdir in newdirs:
        newname = newdir / "nobkg.fits"

        if not newdir.is_dir():
            newdir.mkdir()

        if newname.exists():
            os.remove(newname)

    for file, newdir in zip(filenames, newdirs):
        newname = newdir / "stacked.fits"
        
        os.symlink("../.." / file, newname)
        print("linking", file, "=>", newname)

    return [file for file in filenames]



if __name__ == "__main__":
    link_files()

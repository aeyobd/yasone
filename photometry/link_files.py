import os
from pathlib import Path


def get_newdir(file, destname, increment):
    prefix, filt, num = file.stem.split("_")
    num = int(num) + increment
    newname = f"img_{filt}_{num:02d}"
    return Path(destname) / newname


def link_files(sourcename, destname=None, kind="obj", increment=0):
    if destname is None:
        destname = sourcename

    base_folder = Path(f"../imaging/{sourcename}")
    filenames = [f for f in (base_folder / "flat_fielded").glob(f"{kind}_*.fits")]
    newdirs = [get_newdir(file, destname, increment) for file in filenames]

    for newdir in newdirs:
        newname = newdir / "flat_fielded.fits"

        if not newdir.is_dir():
            newdir.mkdir()

        if newname.exists():
            os.remove(newname)

    for file, newdir in zip(filenames, newdirs):
        newname = newdir / "flat_fielded.fits"
        
        os.symlink("../.." / file, newname)
        print("linking", file, "=>", newname)

    return [file for file in filenames]


def main():
    link_files("yasone1")
    link_files("yasone2")
    link_files("yasone3a", "yasone3")
    link_files("yasone3b", "yasone3", increment=10)

    link_files("20230708", "std1", kind="std")
    link_files("20230709", "std1", kind="std", increment=10)
    link_files("20230723", "std2", kind="std")
    link_files("20230816", "std3", kind="std")



if __name__ == "__main__":
    main()

import tomllib
import os
from pathlib import Path


def safe_rm(pathname):
    if isinstance(pathname, str):
        pathname = Path(pathname)
    if pathname.is_file() or pathname.is_symlink():
        os.remove(pathname)




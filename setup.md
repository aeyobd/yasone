# Software
For the data processing, we require 

## Python
I am using conda. The current dependencies are:
``
conda env export --from-history

name: yasone
channels:
  - conda-forge
  - defaults
dependencies:
  - jupyter
  - scipy
  - numpy
  - matplotlib
  - astropy
  - pandas
  - seaborn
  - python=3.12
  - ccdproc
  - astroalign
  - fitsio
  - astroquery
  - tomli-w
prefix: /Users/daniel/miniconda3/envs/yasone
``
with a few pip dependencies
- arya-aeyobd==0.0.1
- astrometry-net-client==0.6.0
- astroquery==0.4.11
- sausero


## Astromatic.net
These software are tricker to install since the linear algebra package ATLAS may not be possible to compile on macos.

source extractor: via homebrew

SWARP configuration:
```
./configure --with-cfitsio-libdir=/opt/homebrew/Cellar/cfitsio/4.6.3/lib --with-cfitsio-incdir=/opt/homebrew/Cellar/cfitsio/4.6.3/include
```

SCAMP configuration:
```
./configure --enable-openblas --with-openblas-incdir=/opt/homebrew/Cellar/openblas/0.3.30/include --with-openblas-libdir=/opt/homebrew/Cellar/openblas/0.3.30/lib --with-curl-incdir=/opt/homebrew/opt/curl/include --with-curl-libdir=/opt/homebrew/opt/curl/lib
```

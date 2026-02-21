sex  coadd.fits -c ../../psf_calibration.sex -WEIGHT_IMAGE coadd.weight.fits -WEIGHT_TYPE MAP_WEIGHT
psfex psf_calibration.cat -c ../../default.psfex
rm psf_calibration.cat

sex coadd.fits -c ../../psf_phot.sex

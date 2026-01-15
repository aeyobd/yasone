#!/bin/bash
# default 128 with filtsize 6
params="flat_fielded-astrom.fits -CATALOG_NAME test.cat -CATALOG_TYPE NONE -CHECKIMAGE_TYPE BACKGROUND"

sex -c ../../default.sex $params -BACK_SIZE 64 -CHECKIMAGE_NAME bkg64.fits
sex -c ../../default.sex $params -BACK_SIZE 128 -CHECKIMAGE_NAME bkg128.fits
sex -c ../../default.sex $params -BACK_SIZE 256 -CHECKIMAGE_NAME bkg256.fits

sex -c ../../default.sex $params -BACK_FILTERSIZE 3 -CHECKIMAGE_NAME bkg_filt3.fits
sex -c ../../default.sex $params -BACK_FILTERSIZE 4 -CHECKIMAGE_NAME bkg_filt4.fits
sex -c ../../default.sex $params -BACK_FILTERSIZE 6 -BACK_SIZE 64 -CHECKIMAGE_NAME bkg64_filt6.fits
sex -c ../../default.sex $params -BACK_FILTERSIZE 10 -BACK_SIZE 64 -CHECKIMAGE_NAME bkg64_filt12.fits


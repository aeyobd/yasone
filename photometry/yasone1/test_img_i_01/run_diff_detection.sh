#!/bin/bash
params="flat_fielded-astrom.fits -CATALOG_NAME test.cat -CATALOG_TYPE NONE -CHECKIMAGE_TYPE SEGMENTATION,APERTURES"

run_se() {
  sex -c ../../default.sex $params -CHECKIMAGE_NAME seg_$1.fits,ap_$1.fits "${@:2}"
}


run_se min_area_1 -DETECT_MINAREA 1 -DETECT_THRESH 2
run_se min_area_6 -DETECT_MINAREA 6 -DETECT_THRESH 2
run_se min_area_12 -DETECT_MINAREA 12 -DETECT_THRESH 2

run_se 1sigma -DETECT_THRESH 1
run_se 2sigma -DETECT_THRESH 2
run_se 3sigma -DETECT_THRESH 3
run_se 5sigma -DETECT_THRESH 5

run_se filt_gaus2 -FILTER_NAME ../../filters/gaus_2.0_3x3.conv  -DETECT_THRESH 2
run_se filt_gaus3 -FILTER_NAME ../../filters/gaus_3.0_5x5.conv  -DETECT_THRESH 2
#run_se filt_tophat3 -FILTER_NAME ../../filters/tophat_3.0_3x3.conv  -DETECT_THRESH 2
run_se filt_none -FILTER N -DETECT_THRESH 2


run_se deblend_n32 -DEBLEND_NTHRESH 32 -DETECT_THRESH 2
run_se deblend_n128 -DEBLEND_NTHRESH 128 -DETECT_THRESH 2
run_se deblend_c0.0003 -DEBLEND_MINCONT 0.0003 -DETECT_THRESH 2
run_se deblend_c0.02 -DEBLEND_MINCONT 0.02 -DETECT_THRESH 2
run_se noclean -CLEAN N
run_se clean0.8 -CLEAN_PARAM 0.8
run_se clean1.2 -CLEAN_PARAM 1.2


#!/bin/bash
AstrometricSci="/home/bryn/UVIC/Research/Storage/Pointings/EXP_F140Core_drz_sci.fits" 
AstrometricWht="/home/bryn/UVIC/Research/Storage/Pointings/EXP_F140Core_drz_wht.fits"
Sci="/home/bryn/UVIC/Research/Storage/Pointings/EXP_F140Core_drz_sci.fits"
Wht="/home/bryn/UVIC/Research/Storage/Pointings/EXP_F140Core_drz_wht.fits"
Config="/home/bryn/UVIC/Research/Storage/Catalogues/CatConfig/Core-f140w.sex"
Param="/home/bryn/UVIC/Research/Storage/Catalogues/CatConfig/XLSSC-122.param"
Catname="/home/bryn/UVIC/Research/Storage/Catalogues/CatConfig/Core_DI_F140.cat"
Gain="216"
WT1="MAP_WEIGHT"
WT2="MAP_WEIGHT"

command="sex -c "$Config" "$AstrometricSci","$Sci" -WEIGHT_IMAGE "$AstrometricWht","$Wht" -CATALOG_NAME "$Catname" -GAIN $Gain -WEIGHT_TYPE $WT1,$WT2"

echo $command
eval "$command"


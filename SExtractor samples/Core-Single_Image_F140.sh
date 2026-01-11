
#!/bin/bash

Sci="/home/bryn/UVIC/Research/Storage/Pointings/EXP_F140Core_drz_sci.fits"
Wht="/home/bryn/UVIC/Research/Storage/Pointings/EXP_F140Core_drz_wht.fits"
Config="/home/bryn/UVIC/Research/Storage/Catalogues/CatConfig/Core-f140w.sex"
Param="/home/bryn/UVIC/Research/Storage/Catalogues/CatConfig/XLSSC-122.param"
Catname="/home/bryn/UVIC/Research/Storage/Catalogues/CatConfig/Core_SI_F140.cat"
Gain="216"
WT="MAP_WEIGHT"

command="sex -c "$Config" ""$Sci"" -WEIGHT_IMAGE "$Wht" -CATALOG_NAME "$Catname" -GAIN $Gain -WEIGHT_TYPE $WT"

echo $command
eval "$command"

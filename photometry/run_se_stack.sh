objectname=$1

cwd=$(pwd)
for folder in $objectname/coadd_median_*; do
  echo processing $folder
  cd $cwd
  cd $folder
  sex -c ../../stack.sex coadd.fits
done

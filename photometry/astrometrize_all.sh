
for f in $1/img_$2*/flat_fielded-astrom.fits; do
  rm $f
done

for f in $1/img_$2*/flat_fielded.fits; do
  echo \n\n
  echo astrometrizing $f
  python astrometrize.py $f
done



for f in $1/img_*/flat_fielded.fits; do
  python astrometrize.py $f
done

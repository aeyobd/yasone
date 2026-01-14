

for f in yasone2/img_*/flat_fielded.fits; do
  python astrometrize.py $f
done

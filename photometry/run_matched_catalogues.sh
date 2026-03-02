
suffix=""
#catname="psf_phot.cat"
#outname="allcolours_psf.cat"
#foldername="coadd_median"

#catname="aperture.fits"
#outname="allcolours_forced.cat"
#foldername="forced"

#catname="detection.cat"
#outname="allcolours.cat"
#foldername="coadd_median"

catname="photometry.cat"
outname="allcolours_julen_stack.cat"
foldername="julen_stack"


for obj in "yasone1" "yasone2" "yasone3"; do
  python match_catalogues.py $obj $catname $outname $foldername $suffix
done


# objname = "yasone3"
# catname = "photometry.cat"
# outname = "allcolours_julen.cat"
# foldername = "julen_stack"
# suffix = "_0007"

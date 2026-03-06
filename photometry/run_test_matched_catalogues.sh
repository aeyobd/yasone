


#foldername="forced"
#suffix=""
#catname="psf.cat"
#outname="allcolours_forced_psf"
#python match_catalogues.py yasone1 $catname $outname $foldername $suffix
#python match_catalogues.py yasone2 $catname $outname $foldername $suffix
#python match_catalogues.py yasone3 $catname $outname $foldername $suffix
#
#exit 
foldername="forced"
suffix=""
catname="aperture.cat"
outname="allcolours_forced"

python match_catalogues.py yasone1 $catname $outname $foldername $suffix

catname="aperture_large_bkg.cat"
outname="allcolours_forced_large_bkg"
python match_catalogues.py yasone1 $catname $outname $foldername $suffix

catname="aperture_large_bkg_filt.cat"
outname="allcolours_forced_large_bkg_filt"
python match_catalogues.py yasone1 $catname $outname $foldername $suffix

catname="aperture_small_bkg.cat"
outname="allcolours_forced_small_bkg"
python match_catalogues.py yasone1 $catname $outname $foldername $suffix

catname="aperture_small_bkg_filt.cat"
outname="allcolours_forced_small_bkg_filt"
python match_catalogues.py yasone1 $catname $outname $foldername $suffix

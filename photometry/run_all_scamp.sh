objectname=$1

cd $objectname;
rm -rf scamp
mkdir scamp
cd scamp
scamp -c ../../scamp.conf ../img_$2*/detection.cat

#for f in $objectname/img_$2*; do
#  bash run_scamp.sh $f/
#done


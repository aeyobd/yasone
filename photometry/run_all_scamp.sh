objectname=$1
header_dest=nobkg

cd $objectname;

rm -rf scamp
mkdir scamp

cd scamp
echo scamping $(echo ../img_$2*/detection.cat)
scamp -c ../../scamp.conf ../img_$2*/detection.cat
cd .. 
# we are in objectname directory
for folder in img_$2_*; do
  mv $folder/detection.head $folder/$header_dest.head
done


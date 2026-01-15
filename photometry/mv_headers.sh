
objectname=$1

for folder in $objectname/img_*; do
  mv $folder/detection.head $folder/nobkg.head
done

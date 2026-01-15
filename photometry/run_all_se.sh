objectname=$1

for folder in $objectname/img_*; do
  echo processing $folder
  bash run_se.sh $folder;
done

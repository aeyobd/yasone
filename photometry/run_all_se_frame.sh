objectname=$1

for folder in $objectname/img_*; do
  echo processing $folder
  bash run_se_frame.sh $folder;
done

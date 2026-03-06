if [ $# -ne 1 ]; then
  echo "usage $0 dirname"
  exit 1;
fi

echo "cd $1"
cwd=$(pwd)
cd $1

sex -c ../../frame.sex flat_fielded-astrom.fits

cd $cwd

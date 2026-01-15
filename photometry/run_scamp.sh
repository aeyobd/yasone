if [ $# -ne 1 ]; then
  echo "usage $0 dirname"
  exit 1;
fi

cwd=$(pwd)
cd $1

scamp -c ../../scamp.conf detection.cat
mv detection.head nobkg.head

cd $cwd

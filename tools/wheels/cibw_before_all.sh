set -xe

PROJECT_DIR="$1"

PLATFORM=$(PYTHONPATH=tools python3 -c "import get_plat; print(get_plat.get_plat())")

pip3 install scikit-build-core

mkdir /tmp/downloads
cd /tmp/downloads
curl https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.2.0.tar.gz --output cfitsio-4.2.0.tar.gz
tar -xvf cfitsio-4.2.0.tar.gz
cd cfitsio-4.2.0
mkdir build
cd build
cmake ../
make
make install

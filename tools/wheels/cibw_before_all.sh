set -xe

PROJECT_DIR="$1"

echo "Project directory: $PROJECT_DIR"

PLATFORM=$(PYTHONPATH=tools python -c "import get_plat; print(get_plat.get_plat())")

CI_DOWNLOAD_PATH=/tmp/downloads
CI_INSTALL_PREFIX=$CI_DOWNLOAD_PATH/local
LD_LIBRARY_PATH=$CI_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$CI_INSTALL_PREFIX/lib64:$LD_LIBRARY_PATH
DYLD_LIBRARY_PATH=$CI_INSTALL_PREFIX/lib:$DYLD_LIBRARY_PATH
DYLD_LIBRARY_PATH=$CI_INSTALL_PREFIX/lib64:$DYLD_LIBRARY_PATH

PHOTOSPLINE_COMMIT="1faf62b8ad7116fbfcdf1ac7ade763ab9e547402"
CFITSIO_VERSION="4.2.0"
ZLIB_VERSION="1.2.13"

mkdir -p $CI_INSTALL_PREFIX

if [[ $RUNNER_OS == "Linux" ]] ; then
    :
elif [[ $RUNNER_OS == "macOS" ]]; then
    :
elif [[ $RUNNER_OS == "Windows" ]]; then
    mkdir -p $CI_DOWNLOAD_PATH
    cd $CI_DOWNLOAD_PATH
    curl https://www.zlib.net/zlib-$ZLIB_VERSION.tar.gz --output zlib-$ZLIB_VERSION.tar.gz
    tar -xvf zlib-$ZLIB_VERSION.tar.gz
    mkdir -p zlib-$ZLIB_VERSION.build
    cd zlib-$ZLIB_VERSION.build
    cmake ../zlib-$ZLIB_VERSION -DCMAKE_INSTALL_PREFIX=$CI_INSTALL_PREFIX
    cmake --build . --config Release
    cmake --install .
fi

pip install scikit-build-core

if [[ $RUNNER_OS == "Linux" || $RUNNER_OS == "macOS" ]]; then
    mkdir -p $CI_DOWNLOAD_PATH
    cd $CI_DOWNLOAD_PATH
    curl https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-$CFITSIO_VERSION.tar.gz --output cfitsio-$CFITSIO_VERSION.tar.gz
    tar -xvf cfitsio-$CFITSIO_VERSION.tar.gz
    cd cfitsio-$CFITSIO_VERSION
    mkdir -p build
    cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=$CI_INSTALL_PREFIX
    make
    make install
elif [[ $RUNNER_OS == "Windows" ]]; then
    mkdir -p $CI_DOWNLOAD_PATH
    cd $CI_DOWNLOAD_PATH
    curl https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-$CFITSIO_VERSION.tar.gz --output cfitsio-$CFITSIO_VERSION.tar.gz
    tar -xvf cfitsio-$CFITSIO_VERSION.tar.gz
    cd cfitsio-$CFITSIO_VERSION
    mkdir -p cfitsio-$CFITSIO_VERSION.build
    cd cfitsio-$CFITSIO_VERSION.build
    cmake ../cfitsio-$CFITSIO_VERSION -DCMAKE_INSTALL_PREFIX=$CI_INSTALL_PREFIX
    cmake --build . --config Release
    cmake --install .
fi

#cd $CI_DOWNLOAD_PATH
#mkdir -p $CI_DOWNLOAD_PATH/photopline
#cd $CI_DOWNLOAD_PATH/photopline
#git init
#git remote add https://github.com/icecube/photospline.git
#git fetch origin $PHOTOSPLINE_COMMIT
#git reset --hard $PHOTOSPLINE_COMMIT
#mkdir -p build
#cd build
#cmake ../ -DCMAKE_INSTALL_PREFIX=$CI_INSTALL_PREFIX
#make
#make install

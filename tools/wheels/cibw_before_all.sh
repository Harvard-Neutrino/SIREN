set -xe

PROJECT_DIR="$1"

echo "Project directory: $PROJECT_DIR"

PLATFORM=$(PYTHONPATH=tools python -c "import get_plat; print(get_plat.get_plat())")

echo "Runner OS: $RUNNER_OS"

if [[ $RUNNER_OS == "Linux" ]] ; then
    CI_DOWNLOAD_PATH=/tmp/downloads
elif [[ $RUNNER_OS == "macOS" ]]; then
    CI_DOWNLOAD_PATH=/tmp/downloads
elif [[ $RUNNER_OS == "Windows" ]]; then
    CI_DOWNLOAD_PATH=C:/Users/AppData/Local/Temp/downloads
else
    echo "Unknown runner OS: $RUNNER_OS" 1>&2
    exit 1
fi

CI_INSTALL_PREFIX=$CI_DOWNLOAD_PATH/local
LD_LIBRARY_PATH=$CI_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=$CI_INSTALL_PREFIX/lib64:$LD_LIBRARY_PATH
DYLD_LIBRARY_PATH=$CI_INSTALL_PREFIX/lib:$DYLD_LIBRARY_PATH
DYLD_LIBRARY_PATH=$CI_INSTALL_PREFIX/lib64:$DYLD_LIBRARY_PATH

PHOTOSPLINE_COMMIT="1faf62b8ad7116fbfcdf1ac7ade763ab9e547402"
CFITSIO_VERSION="4.6.3"
ZLIB_VERSION="1.3.1"

mkdir -p $CI_INSTALL_PREFIX

if [[ $RUNNER_OS == "Linux" ]] ; then
    if command -v apk >/dev/null 2>&1; then
        apk add --no-cache gsl-dev || { echo "Failed to install GSL (apk)"; exit 1; }  # musllinux (Alpine)
    else
        (yum install -y gsl-devel || dnf install -y gsl-devel || microdnf install -y gsl-devel) || { echo "Failed to install GSL (yum/dnf/microdnf)"; exit 1; }  # manylinux
    fi
elif [[ $RUNNER_OS == "macOS" ]]; then
    export HOMEBREW_NO_AUTO_UPDATE=1
    brew list gsl >/dev/null 2>&1 || brew install gsl || { echo "Failed to install GSL (brew)"; exit 1; }
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

    pip install delvewheel
else
    echo "Unknown runner OS: $RUNNER_OS" 1>&2
    exit 1
fi

pip install scikit-build-core

pip install tomli-w

CFITSIO_TARBALL="cfitsio-$CFITSIO_VERSION.tar.gz"
# Cache dir lives under the project tree so the workflow's populate_cache
# job can populate it before this script runs and restore it for every
# matrix job. The fetch is delegated to fetch_cfitsio.sh (also called
# directly by the workflow) so the mirror list and retry logic live in
# one place.
CFITSIO_CACHE_DIR="$PROJECT_DIR/.cfitsio-cache"

mkdir -p "$CFITSIO_CACHE_DIR"
if [ -f "$CFITSIO_CACHE_DIR/$CFITSIO_TARBALL" ] && \
   tar -tzf "$CFITSIO_CACHE_DIR/$CFITSIO_TARBALL" >/dev/null 2>&1; then
    echo "Using cached $CFITSIO_CACHE_DIR/$CFITSIO_TARBALL"
else
    if [ -f "$CFITSIO_CACHE_DIR/$CFITSIO_TARBALL" ]; then
        echo "Cached tarball is corrupt; re-downloading"
        rm -f "$CFITSIO_CACHE_DIR/$CFITSIO_TARBALL"
    else
        echo "cfitsio tarball not in cache; downloading"
    fi
    bash "$PROJECT_DIR/tools/wheels/fetch_cfitsio.sh" \
        "$CFITSIO_VERSION" "$CFITSIO_CACHE_DIR/$CFITSIO_TARBALL"
fi

if [[ $RUNNER_OS == "Linux" || $RUNNER_OS == "macOS" ]]; then
    mkdir -p "$CI_DOWNLOAD_PATH/cfitsio-$CFITSIO_VERSION"
    tar -xzf "$CFITSIO_CACHE_DIR/$CFITSIO_TARBALL" \
        -C "$CI_DOWNLOAD_PATH/cfitsio-$CFITSIO_VERSION" --strip-components=1
    cd "$CI_DOWNLOAD_PATH/cfitsio-$CFITSIO_VERSION"
    mkdir -p build
    cd build
    cmake ../ -DCMAKE_INSTALL_PREFIX=$CI_INSTALL_PREFIX
    cmake --build . --config Release
    cmake --install .
elif [[ $RUNNER_OS == "Windows" ]]; then
    mkdir -p "$CI_DOWNLOAD_PATH/cfitsio-$CFITSIO_VERSION"
    tar -xzf "$CFITSIO_CACHE_DIR/$CFITSIO_TARBALL" \
        -C "$CI_DOWNLOAD_PATH/cfitsio-$CFITSIO_VERSION" --strip-components=1
    mkdir -p "$CI_DOWNLOAD_PATH/cfitsio-$CFITSIO_VERSION.build"
    cd "$CI_DOWNLOAD_PATH/cfitsio-$CFITSIO_VERSION.build"
    cmake "../cfitsio-$CFITSIO_VERSION" -DCMAKE_INSTALL_PREFIX=$CI_INSTALL_PREFIX
    cmake --build . --config Release
    cmake --install .
else
    echo "Unknown runner OS: $RUNNER_OS" 1>&2
    exit 1
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

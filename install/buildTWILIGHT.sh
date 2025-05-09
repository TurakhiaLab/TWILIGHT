
BUILD_TYPE=$1
# Default flags
CMAKE_OPTIONS=""

if [ "$BUILD_TYPE" == "cuda" ]; then
    # Build with CUDA
    CMAKE_OPTIONS="-DUSE_CUDA=ON -DUSE_HIP=OFF"
elif [ "$BUILD_TYPE" == "hip" ]; then
    # Build with HIP
    CMAKE_OPTIONS="-DUSE_CUDA=OFF -DUSE_HIP=ON"
else
    # Automatically detect and build
    CMAKE_OPTIONS="-DUSE_CUDA=OFF -DUSE_HIP=OFF"
fi

# Create build directory
startDir=$pwd
cd $(dirname "$0")
rm -rf ../build
mkdir -p ../build
cd ../build

# Build and install TBB 
wget https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.9.0.tar.gz
tar -xvzf v2021.9.0.tar.gz
cd oneTBB-2021.9.0
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${PWD}/../install
cmake --build . --config Release --parallel 
cmake --install .

# Build TWILIGHT
cd ../..
cmake -DTBB_DIR=${PWD}/oneTBB-2021.9.0/install/lib/cmake/TBB $CMAKE_OPTIONS ..
make -j

cd $startDir

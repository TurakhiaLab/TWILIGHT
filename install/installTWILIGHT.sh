
# Build
startDir=$pwd
cd $(dirname "$0")
mkdir -p ../build
cd ../build

wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake -DTBB_DIR=${PWD}/oneTBB-2019_U9 -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake ..
make -j

cd $startDir
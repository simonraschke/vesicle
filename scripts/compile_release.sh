#!/bin/bash

echo "rm -rf build/"
rm -rf build/
echo "mkdir build"
mkdir build
cd build
# echo "cmake --build ./build --target all --config Release -- -j"
# cmake --build ./build --target all --config Release -- -j
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
cd ..
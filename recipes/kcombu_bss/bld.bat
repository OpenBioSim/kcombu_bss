
# Build script for conda-forge package.

mkdir build
cd build
cmake ..\src
cmake --build . --config Release
cmake --install .

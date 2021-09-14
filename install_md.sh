#!/bin/bash

current=$(pwd)
echo $current

export INSTALL_DIR=${current}
export PATH=${INSTALL_DIR}/hdf5/bin:${PATH}
export LD_LIBRARY_PATH=${INSTALL_DIR}/hdf5/lib:${LD_LIBRARY_PATH}
export HDF5_ROOT_DIR=${INSTALL_DIR}/hdf5

echo "Installing MD program..."
cd src
make clean
make
cd ..

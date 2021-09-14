#!/bin/bash

current=$(pwd)
echo $current

echo "Building HDF5"
export INSTALL_DIR=${current}
export PATH=${INSTALL_DIR}/hdf5/bin:${PATH}
export LD_LIBRARY_PATH=${INSTALL_DIR}/hdf5/lib:${LD_LIBRARY_PATH}
export HDF5_ROOT_DIR=${INSTALL_DIR}/hdf5

cd hdf5
tar zxvf hdf5-*.*.*.tar.gz 
cd hdf5-*.*.* 

./configure --prefix=${INSTALL_DIR}/hdf5 --enable-parallel --enable-shared \
--enable-build-mode=production --disable-sharedlib-rpath --enable-static \
--enable-fortran CC=mpicc FC=mpif90

make -j 10

make install

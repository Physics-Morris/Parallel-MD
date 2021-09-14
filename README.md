# Parallel Molecular Dynamics

This is a parallel molecular dynamics code using MPI and OpenMP with Fortran 2003 standard. Currently under development.

## Install the dependency

Fist we need to build HDF5 library. Add following line in to your ~/.bashrc and modify your installation directory

```
export INSTALL_DIR=your_installation_directory
export PATH=${INSTALL_DIR}/hdf5/bin:${PATH}
export LD_LIBRARY_PATH=${INSTALL_DIR}/hdf5/lib:${LD_LIBRARY_PATH}
export HDF5_ROOT_DIR=${INSTALL_DIR}/hdf5
```

Then source ~/.bashrc file
```
source ~/.bashrc
```

Go into hdf5 directory and type
```
tar zxvf hdf5-*.*.*.tar.gz
```
```
cd hdf5-*.*.*
```
then use following configuration
```
./configure --prefix=${INSTALL_DIR}/hdf5 --enable-parallel --enable-shared --enable-build-mode=production --disable-sharedlib-rpath --enable-static --enable-fortran CC=mpicc FC=mpif90
```
finally preform make (you can use -j option for multicores)
```
make
```
```
make install
```

## Installation
Go to src/ directory and perform Makefile from there

```
cd src/
make
```

## Usage
To start the program, it needs a input file for simulation parameter.
Example input file is located in inp/ directory.

```
cd inp/
cp default new_inputfile
vi new_inputfile
```

Execute program located in bin/ directory, currently only support perfect cube number of 
processors. You must include -i option to specify your input file.

```
cd bin/
mpirun -np 8 ./MD -i put_your_inputfile_location_here
```

## Contributing
Pull requests are welcome.

## License
[MIT](https://choosealicense.com/licenses/mit/)

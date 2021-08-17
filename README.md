# Parallel Molecular Dynamics

This is a parallel molecular dynamics code using MPI and OpenMP with Fortran 2003 standard. Currently under development.

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

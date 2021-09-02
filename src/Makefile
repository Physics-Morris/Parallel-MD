FC=mpif90
FFLAGS=-O -fopenmp -std=f2003 -g -Wall -fbounds-check -fall-intrinsics -Wextra -pedantic-errors \
 	   -fcheck=all -fbacktrace -Wunderflow -ffpe-trap=zero,overflow,underflow

SRC = constants.f03 math.f03 shared_data.f03 error_handle.f03 particle.f03 helper.f03 \
      diagnostics.f03 mpi_routines.f03 md.f03

TARGET=MD
BINDIR=../bin
OBJ=${SRC:.f03=.o}
LIB=work/morris123/MD/hdf5/lib
INC=work/morris123/MD/hdf5/include

%.o: %.f03
	${FC} ${FFLAGS} -o $@ -I/$(INC) -c $<

${BINDIR}/${TARGET}: ${OBJ}
	${FC} ${FFLAGS} -L/$(LIB) -o $@ ${OBJ} -lhdf5 -lhdf5_fortran ; mv *.mod *.o ../obj

clean:
	@rm -f *.mod *.o ../obj/*.mod ../obj/*.o ../bin/${TARGET}


# dependency
constants.o: constants.f03
math.o: math.f03
shared_data.o: shared_data.f03
error_handle.o: error_handle.f03 constants.o shared_data.o
particle.o: particle.f03 constants.o shared_data.o math.o
helper.o: helper.f03 shared_data.o constants.o math.o
diagnostics.o: diagnostics.f03 shared_data.o constants.o math.o error_handle.o helper.o particle.o
mpi_routines.o: mpi_routines.f03 constants.o helper.o shared_data.o particle.o error_handle.o math.o
md.o: md.f03 constants.o math.o shared_data.o error_handle.o particle.o helper.o diagnostics.o \
	  mpi_routines.o diagnostics.o
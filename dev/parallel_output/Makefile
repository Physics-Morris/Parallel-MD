FC=mpif90
FFLAGS=-O3 -fopenmp -std=gnu -g -Wall -fbounds-check
SRC=hyperslab_by_col.f03
OBJ=${SRC:.f03=.o}
LIB=work/morris123/MD/hdf5/lib
INC=work/morris123/MD/hdf5/include

%.o: %.f03
	${FC} ${FFLAGS} -o $@ -I/$(INC) -c $<

main: ${OBJ}
	${FC} ${FFLAGS} -L/$(LIB) -o $@ ${OBJ} -lhdf5 -lhdf5_fortran

clean:
	@rm -f *.mod *.o main

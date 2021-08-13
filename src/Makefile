FC=mpif90
FFLAGS=-O3 -fopenmp -std=f2003 -g -Wall -fbounds-check -fall-intrinsics
SRC=constants.f03 shared_data.f03 particle.f03 math.f03 helper.f03 mpi_routines.f03 md.f03
TARGET=MD
BINDIR=../bin
OBJ=${SRC:.f03=.o}

%.o: %.f03
	${FC} ${FFLAGS} -o $@ -c $<

${BINDIR}/${TARGET}: ${OBJ}
	${FC} ${FFLAGS} -o $@ ${OBJ} ; mv *.mod *.o ../obj


clean:
	@rm -f *.mod *.o ../obj/*.mod ../obj/*.o ../bin/${TARGET}

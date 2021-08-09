FC=mpif90
FFLAGS=-O3 -fopenmp -std=f2003 -g -Wall -fbounds-check
SRC=helper.f90 variables.f90 md.f90
TARGET=MD
BINDIR=../bin
OBJ=${SRC:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<

${BINDIR}/${TARGET}: ${OBJ}
	${FC} ${FFLAGS} -o $@ ${OBJ} ; mv *.mod *.o ../obj


clean:
	@rm -f ../obj/*.mod ../obj/*.o ../bin/${TARGET}

# Compiler options
CC = mpic++
CFLAGS = -O3 

# Libraries
LIBS = -lboost_program_options -lblas

# Source files
SOURCES = swsolve.cpp swsolve_blas.cpp

# Object files
OBJS = $(SOURCES:.cpp=.o)

# Targets
all: swsolve swsolve_blas

swsolve: 
	$(CC) $(CFLAGS) -o swsolve cw.h swsolve.cpp cw_mpi.h $(LIBS)

swsolve_blas: 
	$(CC) $(CFLAGS) -o swsolve_blas cw_mpi.h swsolve_blas.cpp cw_blas.h $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

test1:
	mpirun -n 8 ./swsolve --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1

test2:
	mpirun -n 8 ./swsolve --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2

test3:
	mpirun -n 8 ./swsolve --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3

test4:
	mpirun -n 8 ./swsolve --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4

test1_blas:
	mpirun -n 8 ./swsolve_blas --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1

test2_blas:
	mpirun -n 8 ./swsolve_blas --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2

test3_blas:
	mpirun -n 8 ./swsolve_blas --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3

test4_blas:
	mpirun -n 8 ./swsolve_blas --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4
			
clean:
	rm -f swsolve swsolve_blas *.o

.PHONY: all clean test1


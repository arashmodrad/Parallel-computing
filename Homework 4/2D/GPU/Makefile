#apps
CUDA_APP=sw_2d
MPI_APP=sw_2d_mpi

#compilers
CUDA_COMPILE = nvcc
MPI_COMPILE = mpicc

#flags
CUDA_FLAGS = -lm
MPI_FLAGS = -lm -Wall

all: ${CUDA_APP} ${MPI_APP}

${CUDA_APP}: ${CUDA_APP}.cu
	${CUDA_COMPILE} ${CUDA_FLAGS} -o $@ $<
${MPI_APP}: ${MPI_APP}.c
	${MPI_COMPILE} ${MPI_FLAGS} -o $@ $<

clean:
	rm -f ${CUDA_APP} ${MPI_APP}

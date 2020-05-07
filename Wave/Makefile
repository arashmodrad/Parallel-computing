CC =gcc
CFLAGS= -Wall -lm
MPICC=mpicc

wave:
	$(MPICC) $(CFLAGS) -o wave wave.c

wave_mpi:
	$(MPICC) $(CFLAGS) -o wave_mpi wave_mpi.c

wave_debug:
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o wave wave.c

wave_mpi_debug:
	$(MPICC) $(CFLAGS) -DDEBUG=1 -o wave_mpi wave_mpi.c


clean:
	rm -f wave
	rm -f wave_mpi

CFLAGS= -Wall -lm
MPICC=mpicc

sw1d: sw_1d.c
	$(MPICC) $(CFLAGS) -o sw1d sw_1d.c

sw1d_mpi: sw_1d_mpi.c
	$(MPICC) $(CFLAGS) -o sw1d_mpi sw_1d_mpi.c

sw2d: sw_2d.c
	$(MPICC) $(CFLAGS) -o sw2d sw_2d.c

sw2d_mpi: sw_2d_mpi.c
	$(MPICC) $(CFLAGS) -o sw2d_mpi sw_2d_mpi.c

clean:
	rm -f sw1d
	rm -f sw1d_mpi
	rm -f sw2d 	rm -f sw2d_mpi

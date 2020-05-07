
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <functions.c>


// -----------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
int main(int argc, char* argv[]){
	//int matSz, subMatSz;
	int irank, nproc;

	MPI_Init(&argc, &argv); //initialize MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &irank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	// 1. Get the matrix size n from command line
	//MPI_Status status;
	int MatSize; // number of data in matrix A or B
	MatSize = getArgs_mpi(argc, argv, irank, MPI_COMM_WORLD); //get command-line argument
	int M = (int)sqrt((double)MatSize);



	int MatSize_Loc; // number of data in local matrixes
	int q = (int)sqrt((double)nproc);
	MatSize_Loc = (int)(MatSize / (double)nproc);
	int m = (int)sqrt((double)MatSize_Loc);
	printf("MatSize= %d MatSize_Loc= %d q= %d M= %d m= %d \n", MatSize, MatSize_Loc, q, M, m);
	// 2 & 3. Allocate matrices A, B, C of size n/q, where q = sqrt(p), and p is the number of processes, then populate them with numbers.
	
	double* A = get_data_debug((int)sqrt((double)MatSize));
	double Sum_A = mat_sum(A, M, M);
	double* B = get_data_debug((int)sqrt((double)MatSize));
	double Sum_B = mat_sum(B, M, M);
	double* C = allocate((int)sqrt((double)MatSize));
	//----------------------------------------------------
	// Compute serial time
	double* C_serial = allocate((int)sqrt((double)MatSize));
	struct timespec start, stop;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	C_serial = mat_mult(A, B, C_serial, M);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
	double serial_time = get_time(start, stop);

	//-----------------------------------------------------------
	//double Sum_C = mat_sum(C, M, M);
		// debugers
	double* C_Debug = allocate((int)sqrt((double)MatSize));
	C_Debug = mat_mult(A, B, C_Debug, M);
	double Sum_C_Debug = mat_sum(C_Debug, M, M);
	if (irank == 0)
	{
		printf("sum of C_Debug for rank %d is : %f \n", irank, Sum_C_Debug);
	}
	

	double* a = allocate((int)sqrt((double)MatSize_Loc));
	//double Sum_a = mat_sum(a, m, m);
	double* b = allocate((int)sqrt((double)MatSize_Loc));
	//double Sum_b = mat_sum(b, m, m);
	double* c = allocate((int)sqrt((double)MatSize_Loc));
	//double Sum_c = mat_sum(c, m, m);
	
	




#ifdef DEBUG
	if (irank == 0)
	{
		printf("sum of A for rank %d is : %f \n", irank, Sum_A);
		printf("sum of B for rank %d is : %f \n", irank, Sum_B);
		printf("sum of C for rank %d is : %f \n", irank, Sum_C);
		printf("sum of C_Debug for rank %d is : %f \n", irank, Sum_C_Debug);
		printf("sum of a for rank %d is : %f \n", irank, Sum_a);
		printf("sum of b for rank %d is : %f \n", irank, Sum_b);
		printf("sum of c for rank %d is : %f \n", irank, Sum_c);
	}
	else
	{
		printf("sum of a for rank %d is : %f \n", irank, Sum_a);
		printf("sum of b for rank %d is : %f \n", irank, Sum_b);
		printf("sum of c for rank %d is : %f \n", irank, Sum_c);
	}

#endif

	double* C_gather = allocate(10);
	double Sum_C = 0;
	double* time = allocate(10);
	for (int ii = 0; ii < 10; ii++){ //timing trick to avoid cache load overhead
		//Begin calculation
		//MPI_Barrier(MPI_COMM_WORLD);

		// 3a) Start the timer
		
		double starttime, stoptime;
		starttime = MPI_Wtime(); // start the race!

		// 4. Perform initial alignment

		int periods[] = { 1, 1 }; //both vertical and horizontal movement; 
		int dims[] = { q, q };
		int coords[2]; //2 Dimension topology so 2 coordinates 
		int right = 0, left = 0, down = 0, up = 0;    // neighbor ranks
		MPI_Comm Cart_Comm;
		MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &Cart_Comm);


		// Distribute A and B to a and b
		MPI_Scatter(A, MatSize_Loc, MPI_DOUBLE, a, MatSize_Loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter(B, MatSize_Loc, MPI_DOUBLE, b, MatSize_Loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//double Sum_a = mat_sum(a, m, m);
		//double Sum_b = mat_sum(b, m, m);
		//printf("sum of a for rank %d is : %f \n", irank, Sum_a);
		//printf("sum of b for rank %d is : %f \n", irank, Sum_b);

		// Determines process coords in cartesian topology given rank in group
		MPI_Cart_coords(Cart_Comm, irank, 2, coords);

		//    a) shift A left by the displacement equal to row number
		MPI_Cart_shift(Cart_Comm, 1, coords[0], &left, &right);
		MPI_Sendrecv_replace(a, MatSize_Loc, MPI_DOUBLE, left, 101, right, 101, Cart_Comm, MPI_STATUS_IGNORE);
		//double Sum_a = mat_sum(a, m, m);
		//printf("sum of a for rank %d is : %f \n", irank, Sum_a);

		//    b) shift B up by the displacement equal to column number
		MPI_Cart_shift(Cart_Comm, 0, coords[1], &up, &down);
		MPI_Sendrecv_replace(b, MatSize_Loc, MPI_DOUBLE, up, 101, down, 101, Cart_Comm, MPI_STATUS_IGNORE);
		//double Sum_b = mat_sum(b, m, m);
		//printf("sum of b for rank %d is : %f \n", irank, Sum_b);

		// 5. Compute C+=A*B
		//c = mult(c, a, b, m);
		c = mat_mult(a, b, c, m);
		//double Sum_c = mat_sum(c, m, m);
		//printf("Sum of c  for rank %d is = %.3f\n", irank, Sum_c);


		// 6. In a loop, repeat the following stages q-1 times:
		for (int i = 1; i < q; i++)
		{
			//	a) Shift A to the left by 1
			MPI_Cart_shift(Cart_Comm, 1, 1, &left, &right);
			MPI_Sendrecv_replace(a, MatSize_Loc, MPI_DOUBLE, left, 101, right, 101, Cart_Comm, MPI_STATUS_IGNORE);
			//    b) Shift B up by 1
			MPI_Cart_shift(Cart_Comm, 0, 1, &up, &down);
			MPI_Sendrecv_replace(b, MatSize_Loc, MPI_DOUBLE, up, 102, down, 102, Cart_Comm, MPI_STATUS_IGNORE);
			//    c) Update C+=A*B (Use our previous code)
			c = mat_mult(a, b, c, m);
			//c = mult(c, a, b, m);
		}
		//double Sum_c = mat_sum(c, m, m);
		//printf("Sum of c  for rank %d is = %.3f\n", irank, Sum_c);


		// 7. We are done
		// gather all c into C
		MPI_Gather(c, MatSize_Loc, MPI_DOUBLE, C, MatSize_Loc, MPI_DOUBLE, 0, Cart_Comm);
		Sum_C = mat_sum(C, M, M);
		C_gather[ii] = Sum_C;


		//    a) stop the timer - should do in a loop, and use fastest one
		stoptime = MPI_Wtime();
		time[ii] = stoptime - starttime;
		// reset 
		Sum_C = 0;
		a = zeros((int)sqrt((double)MatSize_Loc));
		b = zeros((int)sqrt((double)MatSize_Loc));
		c = zeros((int)sqrt((double)MatSize_Loc));
		C = zeros((int)sqrt((double)MatSize));
	}
	
	// calculate the minumum time
	//    b) check for correctness
	if (irank == 0)
	{
		double min = time[0];
		for (int i = 1; i < 10; i++)
		{
			if (time[i] < min)
			{
				min = time[i];
			}
		}
		printf("Dataset size: %d, number of processes = %d, parallel computation time = %e serial computation time = %e s\n", MatSize, nproc, min, serial_time);
		printf("Sum of C  for rank %d is = %.3f\n", irank, C_gather[1]);
	}


	//Deallocation

	free(A);
	free(B);
	free(C);
	free(C_Debug);

	free(a);
	free(b);
	free(c);

	free(time);

	MPI_Finalize();

	

	return 0;

}

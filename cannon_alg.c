
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

	// 2 & 3. Allocate matrices A, B, C of size n/q, where q = sqrt(p), and p is the number of processes, then populate them with numbers.
	if (irank == 0)
	{
		//create global
		double** A = get_data((int)sqrt((double)MatSize));
		double Sum_A = mat_sum(A, M, M);
		double** B = get_data((int)sqrt((double)MatSize));
		double Sum_B = mat_sum(B, M, M);
		double** C = zeros((int)sqrt((double)MatSize));
		// debugers
		double** C_Debug = zeros((int)sqrt((double)MatSize));
		//create local
		double** a = zeros((int)sqrt((double)MatSize_Loc));
		double Sum_a = mat_sum(a, m, m);
		double** b = zeros((int)sqrt((double)MatSize_Loc));
		double Sum_b = mat_sum(b, m, m);
		double** c = zeros((int)sqrt((double)MatSize_Loc));
	}
	else
	{
		//create local
		double** a = zeros((int)sqrt((double)MatSize_Loc));
		double Sum_a = mat_sum(a, m, m);
		double** b = zeros((int)sqrt((double)MatSize_Loc));
		double Sum_b = mat_sum(b, m, m);
		double** c = zeros((int)sqrt((double)MatSize_Loc));






#ifdef DEBUG
	if (irank == 0)
	{
		printf("sum of A for rank %f is : %f \n", Sum_A, irank);
		printf("sum of B for rank %f is : %f \n", Sum_B, irank);
		printf("sum of A for rank %f is : %f \n", Sum_a, irank);
		printf("sum of B for rank %f is : %f \n", Sum_b, irank);
	}
	else
	{
		printf("sum of A for rank %f is : %f \n", Sum_a, irank);
		printf("sum of B for rank %f is : %f \n", Sum_b, irank);
	}

#endif



  // 3a) Start the timer
	
	double starttime, stoptime;
	starttime = MPI_Wtime(); // start the race!

  // 4. Perform initial alignment

	int periods[] = { 1, 1 }; //both vertical and horizontal movement; 
	int dims[] = { m, m };
	int coords[2]; /* 2 Dimension topology so 2 coordinates */
	int right = 0, left = 0, down = 0, up = 0;    // neighbor ranks
	MPI_Comm Cart_Comm;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &Cart_Comm);

	
	// Distribute A and B to a and b
	MPI_Scatter(*A, MatSize_Loc, MPI_DOUBLE, &a[0], MatSize_Loc, MPI_DOUBLE, 0, Cart_Comm);
	MPI_Scatter(*B, MatSize_Loc, MPI_DOUBLE, &b[0], MatSize_Loc, MPI_DOUBLE, 0, Cart_Comm);


	// Determines process coords in cartesian topology given rank in group
	MPI_Cart_coords(Cart_Comm, irank, 2, coords);
 
	//    a) shift A left by the displacement equal to row number
	MPI_Cart_shift(Cart_Comm, 1, coords[0], &left, &right);
	MPI_Sendrecv_replace(&a, MatSize_Loc, MPI_DOUBLE, left, 101, right, 101, Cart_Comm, MPI_STATUS_IGNORE);



	//    b) shift B up by the displacement equal to column number
	MPI_Cart_shift(Cart_Comm, 0, coords[1], &up, &down);
	MPI_Sendrecv_replace(&b, MatSize_Loc, MPI_DOUBLE, up, 202, down, 202, Cart_Comm, MPI_STATUS_IGNORE);


	// 5. Compute C+=A*B
	c = mat_mult(a, b, c, MatSize_Loc);
	double Sum_c = mat_sum(c, m, m);
	printf("Sum of c  for rank %f is = %.3f\n", irank, Sum_c);


	// 6. In a loop, repeat the following stages q-1 times:
	for (i = 1; i<m; i++)
	{
		//	a) Shift A to the left by 1
		MPI_Cart_shift(Cart_Comm, 1, 1, &left, &right);
		MPI_Sendrecv_replace(&a, MatSize_Loc, MPI_DOUBLE, left, 101, right, 101, Cart_Comm, MPI_STATUS_IGNORE);
		//    b) Shift B up by 1
		MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
		MPI_Sendrecv_replace(&b, MatSize_Loc, MPI_DOUBLE, up, 202, down, 202, Cart_Comm, MPI_STATUS_IGNORE);
		//    c) Update C+=A*B (Use our previous code)
		c = mat_mult(a, b, c, MatSize_Loc);
	}


	// 7. We are done
	// gather all c into C
	MPI_Gather(&c, MatSize_Loc, MPI_DOUBLE, C, MatSize_Loc, MPI_DOUBLE, 0, Cart_Comm);
	//    a) stop the timer - should do in a loop, and use fastest one
	stoptime = MPI_Wtime();
	if (irank == 0)
	{
		double time[nproc];
	}
	MPI_Gather(&stoptime, 1, MPI_DOUBLE, time, 1, MPI_DOUBLE, 0, Cart_Comm);
	if (irank == 0)
	{
		double min = time[0];
		for (int i = 1; i < nproc; i++)
		{
			if (time[i] < min)
			{
				min = time[i];
			}
		}
		printf("Dataset size: %d, number of processes = %f, computation time = %e s\n", MatSize, nproc, min - starttime);
	}

	//    b) check for correctness
	if (irank == 0)
	{
		double Sum_C = mat_sum(C, M, M);
		printf("sum of C for rank %f is : %f \n", Sum_C, irank);
		C_Debug = mat_mult(A, B, C_Debug, MatSize);
		double Sum_C_Debug = mat_sum(C_Debug, M, M);
		printf("Sum of C: %d, Sum of C_Debug = %d s\n", Sum_C, Sum_C_Debug);
	}

	MPI_Finalize();

	return 0;

}
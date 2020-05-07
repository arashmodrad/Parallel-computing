#include <stdio.h>
#include <math.h>
// ----------------------------------------------------------------------------
// get_data creates an array of N numbers as output
double* get_data_debug(int N){
  
  
  double* x = (double*)malloc(N * N * sizeof(double));

  

  // Fill the matrix with data
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      *(x + i*N + j) = 1;
    }
  }

  return x;

}
// ----------------------------------------------------------------------------
// get_data creates an array of N numbers as output
double* allocate(int N){


	double* x = (double*)malloc(N * N * sizeof(double));



	return x;

}

// ----------------------------------------------------------------------------
// zeros(N) will produce a matrix of zeros
double* zeros(int N){

  double *x = (double*)malloc(N*N*sizeof(double));


  // Fill the matrix with zeros
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
		*(x + i*N + j) = 0;
    }
  }

  return x;

}

// ----------------------------------------------------------------------------
// get_data creates an array of N pseudo-random numbers as output
double* get_data(int N){

  double *x = (double*)malloc(N*N*sizeof(double));


  // Fill the matrix with random values
  srand(time(0)); // set random seed based on time
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      x[(i*N)+j] = (double)rand() / (double)RAND_MAX;
    }
  }

 return x;

}

// ----------------------------------------------------------------------------
// print_mat will print a matrix
void print_mat(double* mat, int m, int n){
  printf("\n");
  for (int row = 0; row < m; row++){
    printf("|");
    for (int col = 0; col < n; col++){
      printf("[%d][%d] = %f \n", row, col, (mat + row*n)[col]);
    }
    printf("|\n");
  }
  printf("\n");
        
  return;

}


// ----------------------------------------------------------------------------
// get_time will check wall clock time in seconds
double get_time(struct timespec start, struct timespec stop){
  return ((stop.tv_sec - start.tv_sec)*1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3)/1e6;
}
// ----------------------------------------------------------------------------
// get_gflops will return the number of Gigaflops computed
double get_gflops(double ops, struct timespec start, struct timespec stop){
  double time = get_time(start,stop);
  return ops/time/1e9;
}


// ----------------------------------------------------------------------------
int getArgs_mpi(int argc, char *argv[], int irank, MPI_Comm comm)
{
	
	int mtSz;

  if(irank==0){
    if ( argc != 2 ) /* argc should be 2 for correct execution */
      {
	//If not correct number of arguments, assume n=1000
	printf("Incorrect number of arguments. Usage: ./derivation N \n Assuming N=1024.\n");
	mtSz=1024;
      }
    else
      {
	//Use input argument
	mtSz = atoi(argv[1]);	
      }
  }
  MPI_Bcast(&mtSz,1,MPI_INT,0,comm);
  return mtSz;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
double* mat_mult(double* A, double* B, double* C, int matSz){
	double Abuff;
	int block, blockSz, i, j, k;
	if (matSz <= 256){
		blockSz = matSz;
	}
	else{
		blockSz = 256;
	}

	for (block = 0; block < matSz / blockSz; block++){
		for (i = 0; i < matSz; i++){
			for (k = 0; k < matSz; k++){
				Abuff = (A + i*matSz)[k];
				for (j = block*blockSz; j < (block + 1)*blockSz / 1; j++){
					//C[(i*matSz) + (1 * j)] += Abuff * B[(k*matSz) + (1 * j)];
					(C + i*matSz)[1 * j] += Abuff * (B + k*matSz)[1 * j];
					
				}
			}
		}
	}
	return C;
}

// -----------------------------------------------------------------------
double mat_sum(double* mat, int m, int n){
  double sum = 0.;
  for (int i = 0; i < m; i++){
    for (int j =0; j < n; j++){
      sum += (mat + i*n)[j];
    }
  }
  return sum;
}

// -----------------------------------------------------------------------
double* mult(double* c, double* a, double* b, int matSz){
	
	for (int i = 0; i < matSz; i++)
	{
		for (int j = 0; j < matSz; j++)
		{
			for (int k = 0; k < matSz; k++)
			{
				(c + i*matSz)[j] += (a + i*matSz)[k] * (b + k*matSz)[j];
				//c[(i*matSz) + j] += a[(i*matSz) + k] * b[(k*matSz) + j];
			}
		}

	}
	return c;
}
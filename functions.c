#include <stdio.h>
#include <math.h>
// ----------------------------------------------------------------------------
// get_data creates an array of N numbers as output
double** get_data_debug(int N){
  
  
  double **x = (double**)malloc(N*sizeof(double*));
  for (int i = 0; i <N; i++){
      x[i] = (double *)malloc(N*sizeof(double));
  }
  

  // Fill the matrix with data
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      x[i][j] = (double) (10*i + j);
    }
  }

  return x;

}

// ----------------------------------------------------------------------------
// get_data_rank creates an array of N numbers as output
double** get_data_rank(int N, int irank){

  // Preallocate the matrix
  double **x = (double**) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++){
    x[i] = (double *) malloc(N*sizeof(double));
  }

  // Fill the matrix with data
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      x[i][j] = ( (double) irank + (double) i/10 + (double) j/100);
    }
  }

  return x;

}

// ----------------------------------------------------------------------------
// zeros(N) will produce a matrix of zeros
double** zeros(int N){

  double **x = (double**)malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++){
    x[i] = (double *)malloc(N*sizeof(double));
  }

  // Fill the matrix with zeros
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      x[i][j] = 0;
    }
  }

  return x;

}

// ----------------------------------------------------------------------------
// get_data creates an array of N pseudo-random numbers as output
double** get_data(int N){

  double **x = (double**)malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++){
    x[i] = (double *)malloc(N*sizeof(double));
  }

  // Fill the matrix with random values
  srand(time(0)); // set random seed based on time
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      x[i][j] = (double)rand() / (double)RAND_MAX;
    }
  }

 return x;

}

// ----------------------------------------------------------------------------
// print_mat will print a matrix
void print_mat(double** mat, int m, int n){
  printf("\n");
  for (int row = 0; row < m; row++){
    double *x = *(mat + row);
    printf("|");
    for (int col = 0; col < n; col++){
      printf("%10.3f, ",*(x + col));
    }
    printf("|\n");
  }
  printf("\n");
        
  return;

}

// ----------------------------------------------------------------------------
// fprint_mat will print a matrix
void fprint_mat(FILE* fp, double** mat, int m, int n){
  fprintf(fp, "\n");
  for (int row = 0; row < m; row++){
    double *x = *(mat + row);
    fprintf(fp, "|");
    for (int col = 0; col < n; col++){
      fprintf(fp, "%10.3f, ",*(x + col));
    }
    fprintf(fp, "|\n");
  }
  fprintf(fp, "\n");
        
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
double** mat_mult(double** A, double** B, double** C, int matSz){
  double dmy;
  int block, blockSz, i, j, k;
  if (matSz <= 256){
    blockSz = matSz;
  }
  else{
    blockSz = 256;
  }

  for (block = 0; block < matSz/blockSz; block++){
    for (i = 0; i < matSz; i++){
      for (k = 0; k < matSz; k++){
        dmy = A[i][k];
        // for (j = block*blockSz; j < (block + 1)*blockSz/8; j++){
        //   C[ i ][ 8*j     ] += dmy * B[ k ][ 8*j     ];
        //   C[ i ][ 8*j + 1 ] += dmy * B[ k ][ 8*j + 1 ];
        //   C[ i ][ 8*j + 2 ] += dmy * B[ k ][ 8*j + 2 ];
        //   C[ i ][ 8*j + 3 ] += dmy * B[ k ][ 8*j + 3 ];
        //   C[ i ][ 8*j + 4 ] += dmy * B[ k ][ 8*j + 4 ];
        //   C[ i ][ 8*j + 5 ] += dmy * B[ k ][ 8*j + 5 ];
        //   C[ i ][ 8*j + 6 ] += dmy * B[ k ][ 8*j + 6 ];
        //   C[ i ][ 8*j + 7 ] += dmy * B[ k ][ 8*j + 7 ];
        for (j = block*blockSz; j < (block + 1)*blockSz/1; j++){
          C[ i ][ 1*j     ] += dmy * B[ k ][ 1*j     ];
          // C[ i ][ 1*j     ] += dmy * B[ k ][ 1*j     ];
        }
      }
    }
  }
  return C;
}

// -----------------------------------------------------------------------
double mat_sum(double** mat, int m, int n){
  double sum = 0.;
  for (int i = 0; i < m; i++){
    for (int j =0; j < n; j++){
      sum += mat[j][i];
    }
  }
  return sum;
}
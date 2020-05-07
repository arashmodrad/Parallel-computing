#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


double get_time(struct timespec start, struct timespec stop)
{
	return ((stop.tv_sec - start.tv_sec) * 1e6 + (stop.tv_nsec - start.tv_nsec) / 1e3) / 1e6;    // in seconds
}

double get_gflops(double ops, struct timespec start, struct timespec stop)
{
	double time = get_time(start, stop);
	return ops / time / 1e9;
}

/*int is_same(double *C, double *C_ref, int n)
{
	double sum = 0;
       	for (int i = 0; i < n; i++)
	{
	  sum += C[i] - C_ref[i];
	}
	if (sum / n>1e10) {
		printf(" Result does not match %d", sum / n);
		return 1;
	}
	else {
		return 0;
	}
	
}
*/
void reset(double **C, int n)
{
	//reset the result vector C
	for(int i = 0; i < n; i++)
	{
	  for(int j = 0; j < n; j++)
	    {		
	      C[i][j] = 0;
	    }
	}
     
}

void printC(double **C, int n, double alg)
{
  double sum;
	for(int i = 0; i < n; i++)
	{
	  for(int j = 0; j < n; j++)
	    {		
	     sum += C[i][j];
	    }
	}		
  printf("total sum of algorithem %f is: %f\n", alg, sum);
  sum = 0.0;        
}


// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
//main program
int main(int argc, char* argv[]){
    //driver

  int n, i, j, k;
  // double C[n][n], C1[n/4][n/4], C2[n/4][n/4], C3[n/4][n/4], C4[n/4][n/4];
  double gflops, NumOps, count, count2;
  struct timespec start, stop;
  // double C_ref[n][n];

  // Read command line arguments
  if ( argc != 2 ) /* argc should be 2 for correct execution */
    {
      printf("Incorrect number of arguments");
    }
  else
    {
      //Use input argument
      n = atoi(argv[1]);
    }


  NumOps = 2.0 * pow(n, 3);
  printf("ops %f\n",NumOps);
  

  double **A = (double**)malloc(n*sizeof(double*));
  double **A1 = (double**)malloc(n/4*sizeof(double*));
  double **A2 = (double**)malloc(n/4*sizeof(double*));
  double **A3 = (double**)malloc(n/4*sizeof(double*));
  double **A4 = (double**)malloc(n/4*sizeof(double*));
  double **B = (double**)malloc(n*sizeof(double*));
  double **B1 = (double**)malloc(n/4*sizeof(double*));
  double **B2 = (double**)malloc(n/4*sizeof(double*));
  double **B3 = (double**)malloc(n/4*sizeof(double*));
  double **B4 = (double**)malloc(n/4*sizeof(double*));
  double **C = (double**)malloc(n*sizeof(double*));
  double **A1B1 = (double**)malloc(n/4*sizeof(double*));
  double **A1B2 = (double**)malloc(n/4*sizeof(double*));
  double **A2B3 = (double**)malloc(n/4*sizeof(double*));
  double **A2B4 = (double**)malloc(n/4*sizeof(double*));
  double **A3B1 = (double**)malloc(n/4*sizeof(double*));
  double **A3B2 = (double**)malloc(n/4*sizeof(double*));
  double **A4B3 = (double**)malloc(n/4*sizeof(double*));
  double **A4B4 = (double**)malloc(n/4*sizeof(double*));
  double **C_ref = (double**)malloc(n*sizeof(double*));  
  for (i = 0; i < n; i++)
  {
	  A[i] = (double *)malloc(n * sizeof(double));
	  B[i] = (double *)malloc(n * sizeof(double));
	  C[i] = (double *)malloc(n * sizeof(double));
	  C_ref[i] = (double *)malloc(n * sizeof(double));
  }
  for (i = 0; i < n/4; i++)
  {
       	  A1[i] = (double *)malloc(n/4 * sizeof(double));
	  A2[i] = (double *)malloc(n/4 * sizeof(double));
	  A3[i] = (double *)malloc(n/4 * sizeof(double));
	  A4[i] = (double *)malloc(n/4 * sizeof(double));
	  B1[i] = (double *)malloc(n/4 * sizeof(double));
	  B2[i] = (double *)malloc(n/4 * sizeof(double));
	  B3[i] = (double *)malloc(n/4 * sizeof(double));
       	  A1B1[i] = (double *)malloc(n/4 * sizeof(double));
	  A1B2[i] = (double *)malloc(n/4 * sizeof(double));
	  A2B3[i] = (double *)malloc(n/4 * sizeof(double));
	  A2B4[i] = (double *)malloc(n/4 * sizeof(double));
       	  A3B1[i] = (double *)malloc(n/4 * sizeof(double));
	  A3B2[i] = (double *)malloc(n/4 * sizeof(double));
	  A4B3[i] = (double *)malloc(n/4 * sizeof(double));
	  A4B4[i] = (double *)malloc(n/4 * sizeof(double));
	  B4[i] = (double *)malloc(n/4 * sizeof(double));
  }
	  


  // create matrix --------------------------------------------------------------
  count = 0;
  count2 = 0;
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
	          A[i][j] = (double)1;
		  B[i][j] = (double)1; 
	  }

  }
  /*  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		  B[i][j] = (double)++count2; 
	  }

  }

  /*  
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		  printf("%d \n", &A[i][j]);
		  printf("B[%d][%d]= %d \n", i, j, B[i][j]); 
	  }

   } */
  double alg;
      //########################### Algorithem 1 #########################################
  alg = 1;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		  for (k = 0; k < n; k++)
		  {
			  C[i][j] += A[i][k] * B[k][j];
		  }
	  }
		  
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result1 = get_gflops(NumOps, start, stop);
  printf("%f \t", result1); 
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		  C_ref[i][j] = C[i][j];
	  }

  }
  
  /*  //print result
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		   printf("C[%f][%f] = %f\n", i, j, C[i][j]);
	  }

  }
  */
  printC(C, n, alg);
  reset(C, n); 
  	//########################### Algorithem 2 #########################################
  alg = 2;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (i = 0; i < n; i++)
  {
	  for (k = 0; k < n; k++)
	  {
		  for (j = 0; j < n; j++)
		  {
		    C[i][j] += A[i][k] * B[k][j];
		  }
	  }

  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result2 = get_gflops(NumOps, start, stop);
  printf("result2 = %f \t", result2);
  printC(C, n, alg);
  reset(C, n);
    //is_same(C, C_ref, n); //check whether the solution is close to the reference
	
  //########################### Algorithem 3 #########################################
  alg = 3;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (k = 0; k < n; k++)
  {
	  for (i = 0; i < n; i++)
	  {
		  for (j = 0; j < n; j++)
		  {
			  C[i][j] += A[i][k] * B[k][j];
		  }
	  }

  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result3 = get_gflops(NumOps, start, stop);
  printf("result3 = %f \t", result3);
  printC(C, n, alg);  
  reset(C, n);
  //is_same(C, C_ref, n); //check whether the solution is close to the reference
    	
  //########################### Algorithem 4 #########################################
  double Store;
  alg = 4;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		  Store = 0;
		  for (k = 0; k < n; k++)
		  {
			  Store += A[i][k] * B[k][j];
		  }
		  C[i][j] = Store;
	  }

  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result4 = get_gflops(NumOps, start, stop);
  printf("result4 = %f \t", result4);
  printC(C, n, alg);
  reset(C, n);
  //  is_same(C, C_ref, n); //check whether the solution is close to the reference
  		
 //########################### Algorithem 5A #########################################
  double S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16;
  alg = 5.1; 
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		  S1 = S2 = S3 = S4 = S5 = S6 = S7 = S8 = S9 = S10 = S11 = S12 = S13 = S14 = S15 = S16 = 0;
		  for (k = 0; k < n/16; k++)
		  {
			  S1 += A[i][16 * k] * B[16 * k][j];
			  S2 += A[i][16 * k + 1] * B[16 * k + 1][j];
			  S3 += A[i][16 * k + 2] * B[16 * k + 2][j];
			  S4 += A[i][16 * k + 3] * B[16 * k + 3][j];
			  S5 += A[i][16 * k + 4] * B[16 * k + 4][j];
			  S6 += A[i][16 * k + 5] * B[16 * k + 5][j];
			  S7 += A[i][16 * k + 6] * B[16 * k + 6][j];
			  S8 += A[i][16 * k + 7] * B[16 * k + 7][j];
			  S9 += A[i][16 * k + 8] * B[16 * k + 8][j];
			  S10 += A[i][16 * k + 9] * B[16 * k + 9][j];
			  S11 += A[i][16 * k + 10] * B[16 * k + 10][j];
			  S12 += A[i][16 * k + 11] * B[16 * k + 11][j];
			  S13 += A[i][16 * k + 12] * B[16 * k + 12][j];
			  S14 += A[i][16 * k + 13] * B[16 * k + 13][j];
			  S15 += A[i][16 * k + 14] * B[16 * k + 14][j];
			  S16 += A[i][16 * k + 15] * B[16 * k + 15][j];
		  }
		  C[i][j] = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8 + S9 + S10 + S11 + S12 + S13 + S14 + S15 + S16;
	  }

  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result5 = get_gflops(NumOps, start, stop);
  printf("result5A = %f \t", result5);
  printC(C, n, alg);
  reset(C, n);
  // is_same(C, C_ref, n); //check whether the solution is close to the reference
     	
//########################### Algorithem 5B #########################################
  alg = 5.2;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		  S1 = S2 = S3 = S4 = S5 = S6 = S7 = S8 = 0;
		  for (k = 0; k < n / 8; k++)
		  {
			  S1 += A[i][8 * k] * B[8 * k][j];
			  S2 += A[i][8 * k + 1] * B[8 * k + 1][j];
			  S3 += A[i][8 * k + 2] * B[8 * k + 2][j];
			  S4 += A[i][8 * k + 3] * B[8 * k + 3][j];
			  S5 += A[i][8 * k + 4] * B[8 * k + 4][j];
			  S6 += A[i][8 * k + 5] * B[8 * k + 5][j];
			  S7 += A[i][8 * k + 6] * B[8 * k + 6][j];
			  S8 += A[i][8 * k + 7] * B[8 * k + 7][j];
		  }
		  C[i][j] = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8;
	  }

  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result6 = get_gflops(NumOps, start, stop);
  printf("result 5B = %f \t", result6);
  printC(C, n, alg);
  reset(C, n);
  // is_same(C, C_ref, n); //check whether the solution is close to the reference
  	
//########################### Algorithem 5C #########################################
  alg = 5.3;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (i = 0; i < n; i++)
  {
	  for (j = 0; j < n; j++)
	  {
		  S1 = S2 = 0;
		  for (k = 0; k < n / 2; k++)
		  {
			  S1 += A[i][2 * k] * B[2 * k][j];
			  S2 += A[i][2 * k + 1] * B[2 * k + 1][j];
		  }
		  C[i][j] = S1 + S2;
	  }

  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result7 = get_gflops(NumOps, start, stop);
  printf("result 5C = %f \t", result7);
  printC(C, n, alg);
  reset(C, n);
  
  //  is_same(C, C_ref, n); //check whether the solution is close to the reference
/*	
  //########################### Algorithem 6 #########################################
  alg = 6.1;
  int cb = 16;
  int kk, jj;
  double sum1; // block size
  //int it = 0.0;
  S1 = S2 = S3 = S4 = S5 = S6 = S7 = S8 = 0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (kk = 0; kk < n / cb; kk++)
  {
    for (jj = 0; jj < n / cb; jj++)
      {
	  for (i = 0; i < n; i++)
	  {
	    for (j = jj * cb/8; j < (kk + 1)/8 *cb; j++)
		  {
			  S1 = S2 = S3 = S4 = S5 = S6 = S7 = S8 = 0;
		          sum1 = C[i][j];
			  for (k = kk * cb/8; k < (kk + 1)/8 * cb; k++)
			  {
				  S1 += A[i][8 * k] * B[8 * k][j];
				  S2 += A[i][(8 * k) + 1] * B[(8 * k) + 1][j];
				  S3 += A[i][(8 * k) + 2] * B[(8 * k) + 2][j];
				  S4 += A[i][(8 * k) + 3] * B[(8 * k) + 3][j];
				  S5 += A[i][(8 * k) + 4] * B[(8 * k) + 4][j];
				  S6 += A[i][(8 * k) + 5] * B[(8 * k) + 5][j];
				  S7 += A[i][(8 * k) + 6] * B[(8 * k) + 6][j];
				  S8 += A[i][(8 * k) + 7] * B[(8 * k) + 7][j];
				  
			  }
			  sum1 = S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8;
			  C[i][j] = sum1;
			  //printf("C[%d][%d] = %f \n", i, j, C[i][j]);
		  }

	  }
      }
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result10 = get_gflops(NumOps, start, stop);
  printf("result 6.1 = %f \t", result10);
  printC(C, n, alg);
  reset(C, n);
  
//  is_same(C, C_ref, n); //check whether the solution is close to the reference

*/	
  //########################### Algorithem 6 #########################################
  alg = 6.1;
  int cb = 4; // block size
  int it = 0.0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (int jj = 0; jj < n / cb; jj++)
  {
	  for (i = 0; i < n; i++)
	  {
		  for (j = 0; j < n; j++)
		  {
			  S1 = S2 = S3 = S4 = S5 = S6 = S7 = S8 = 0;
			  for (k = jj * cb/8; k < (jj + 1) * cb/8; k++)
			  {
				  S1 += A[i][8 * k] * B[8 * k][j];
				  S2 += A[i][(8 * k) + 1] * B[(8 * k) + 1][j];
				  S3 += A[i][(8 * k) + 2] * B[(8 * k) + 2][j];
				  S4 += A[i][(8 * k) + 3] * B[(8 * k) + 3][j];
				  S5 += A[i][(8 * k) + 4] * B[(8 * k) + 4][j];
				  S6 += A[i][(8 * k) + 5] * B[(8 * k) + 5][j];
				  S7 += A[i][(8 * k) + 6] * B[(8 * k) + 6][j];
				  S8 += A[i][(8 * k) + 7] * B[(8 * k) + 7][j];
				  //it++;
			  }
			  C[i][j] += S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8;
			  //printf("C[%d][%d] = %f \n", i, j, C[i][j]);
		  }

	  }
	  
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result81 = get_gflops(NumOps, start, stop);
  printf("result 6B = %f \t", result81);
  printC(C, n, alg);
  reset(C, n);
  
//  is_same(C, C_ref, n); //check whether the solution is close to the reference
  //########################### Algorithem 6 #########################################
  alg = 6.3;
  cb = 64; // block size
  //int it = 0.0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (int jj = 0; jj < n / cb; jj++)
  {
	  for (i = 0; i < n; i++)
	  {
		  for (j = 0; j < n; j++)
		  {
			  S1 = S2 = S3 = S4 = S5 = S6 = S7 = S8 = 0;
			  for (k = jj * cb/8; k < (jj + 1) * cb/8; k++)
			  {
				  S1 += A[i][8 * k] * B[8 * k][j];
				  S2 += A[i][(8 * k) + 1] * B[(8 * k) + 1][j];
				  S3 += A[i][(8 * k) + 2] * B[(8 * k) + 2][j];
				  S4 += A[i][(8 * k) + 3] * B[(8 * k) + 3][j];
				  S5 += A[i][(8 * k) + 4] * B[(8 * k) + 4][j];
				  S6 += A[i][(8 * k) + 5] * B[(8 * k) + 5][j];
				  S7 += A[i][(8 * k) + 6] * B[(8 * k) + 6][j];
				  S8 += A[i][(8 * k) + 7] * B[(8 * k) + 7][j];
				  //it++;
			  }
			  C[i][j] += S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8;
			  //printf("C[%d][%d] = %f \n", i, j, C[i][j]);
		  }

	  }
	  
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result83 = get_gflops(NumOps, start, stop);
  printf("result 6C = %f \t", result83);
  printC(C, n, alg);
  reset(C, n);
  
//  is_same(C, C_ref, n); //check whether the solution is close to the reference
  //########################### Algorithem 6 #########################################
  alg = 6;
  cb = 16; // block size
  // int it = 0.0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (int jj = 0; jj < n / cb; jj++)
  {
	  for (i = 0; i < n; i++)
	  {
		  for (j = 0; j < n; j++)
		  {
			  S1 = S2 = S3 = S4 = S5 = S6 = S7 = S8 = 0;
			  for (k = jj * cb/8; k < (jj + 1) * cb/8; k++)
			  {
				  S1 += A[i][8 * k] * B[8 * k][j];
				  S2 += A[i][(8 * k) + 1] * B[(8 * k) + 1][j];
				  S3 += A[i][(8 * k) + 2] * B[(8 * k) + 2][j];
				  S4 += A[i][(8 * k) + 3] * B[(8 * k) + 3][j];
				  S5 += A[i][(8 * k) + 4] * B[(8 * k) + 4][j];
				  S6 += A[i][(8 * k) + 5] * B[(8 * k) + 5][j];
				  S7 += A[i][(8 * k) + 6] * B[(8 * k) + 6][j];
				  S8 += A[i][(8 * k) + 7] * B[(8 * k) + 7][j];
				  //it++;
			  }
			  C[i][j] += S1 + S2 + S3 + S4 + S5 + S6 + S7 + S8;
			  //printf("C[%d][%d] = %f \n", i, j, C[i][j]);
		  }

	  }
	  
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result8 = get_gflops(NumOps, start, stop);
  printf("result 6B = %f \t", result8);
  printC(C, n, alg);
  reset(C, n);
  
//  is_same(C, C_ref, n); //check whether the solution is close to the reference

  //########################### Algorithem 6.2 #########################################
  alg = 6.2;
  // int cb = 16; // block size
  //int it = 0.0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  for (int jj = 0; jj < n / cb; jj++)
  {
	  for (i = 0; i < n; i++)
	  {
		  for (j = 0; j < n; j++)
		  {
			  S1 = 0;
			  for (k = jj * cb; k < (jj + 1) * cb; k++)
			  {
				  S1 += A[i][k] * B[k][j];
			  }
			  C[i][j] += S1;
			  //printf("C[%d][%d] = %f \n", i, j, C[i][j]);
		  }

	  }
	  
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result82 = get_gflops(NumOps, start, stop);
  printf("result 6D = %f \t", result82);
  printC(C, n, alg);
  reset(C, n);
  
//  is_same(C, C_ref, n); //check whether the solution is close to the reference

//########################### Extra 1 #########################################
  alg = 7;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
  //break A and B into four sub matrixes
  for (i = 0; i < n/4; i++)
  {
	  for (j = 0; j < n/4; j++)
	  {
		  A1[i][j] = A[i][j];
		  B1[i][j] = B[i][j];
	  }

  }
  for (i = 0; i < n/4; i++)
  {
	  for (j = 0; j < n/4; j++)
	  {
	    A2[i][j] = A[i+(n/4)][j+(n/4)];
	    B2[i][j] = B[i+(n/4)][j+(n/4)];
	  }

  }
  for (i = 0; i < n/4; i++)
  {
	  for (j = 0; j <  n/4; j++)
	  {
	    A3[i][j] = A[i+(n/2)][j+(n/2)];
	    B3[i][j] = B[i+(n/2)][j+(n/2)];
	  }

  }
  for (i = 0; i < n/4; i++)
  {
	  for (j =  0; j < n/4; j++)
	  {
	    A4[i][j] = A[i+(3*n/4)][j+(3*n/4)];
	    B4[i][j] = B[i+(3*n/4)][j+(3*n/4)];
	  }

  }
  // caclulate products
  for (i = 0; i < n/4; i++)
  {
	  for (j = 0; j < n/4; j++)
	  {
		  for (k = 0; k < n/4; k++)
		  {
			  A1B1[i][j] += A1[i][k] * B1[k][j];
			  A2B3[i][j] += A2[i][k] * B3[k][j];
			  A1B2[i][j] += A1[i][k] * B2[k][j];
			  A2B4[i][j] += A2[i][k] * B4[k][j];
			  A3B1[i][j] += A3[i][k] * B1[k][j];
			  A4B3[i][j] += A4[i][k] * B3[k][j];
			  A3B2[i][j] += A3[i][k] * B2[k][j];
			  A4B4[i][j] += A4[i][k] * B4[k][j];
		  }
	  }

  }
   // puting them back together
  for (i = 0; i < n/4; i++)
  {
	  for (j = 0; j < n/4; j++)
	  {
		  C[i][j] = A1B1[i][j] + A2B3[i][j];
	  }

  }
  for (i = 0; i < n/4; i++)
  {
	  for (j = 0; j < n / 4; j++)
	  {
	    C[i+(n/4)][j+(n/4)] = A1B2[i][j] + A2B4[i][j];
	  }

  }
  for (i = 0; i < n/4; i++)
  {
	  for (j = 0; j < n/4; j++)
	  {
	    C[i+(n/2)][j+(n/2)] = A3B1[i][j] + A4B3[i][j];
	  }

  }
  for (i = 0; i < n/4; i++)
  {
	  for (j = 0; j < n/4; j++)
	  {
	    C[i+(3*n/4)][j+(3*n/4)] = A3B2[i][j] + A4B4[i][j];
	  }

  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
  double result9 = get_gflops(NumOps, start, stop);
  printf("result extra1 = %f \t", result9);
  printC(C, n, alg);
  reset(C, n);
  //  is_same(C, C_ref, n); //check whether the solution is close to the reference
  
  //-------------------------------------------------------------------------------------
  
  //Deallocate memory
  for (i = 0; i < n; i++)
  {
	 free(A[i]);
         free(B[i]);
	 free(C[i]);
	 free(C_ref[i]);

  }
  free(A);
  free(B);
  free(C);
  free(C_ref);
  for (i = 0; i < n/4; i++)
  {
	 free(A1[i]);
         free(A2[i]);
	 free(A3[i]);
	 free(A4[i]);
	 free(B1[i]);
	 free(B2[i]);
	 free(B3[i]);
	 free(B4[i]);
	 free(A1B1[i]);
         free(A1B2[i]);
	 free(A2B3[i]);
	 free(A2B4[i]);
	 free(A3B1[i]);
	 free(A3B2[i]);
	 free(A4B3[i]);
	 free(A4B4[i]);

  }
  free(A1);
  free(A2);
  free(A3);
  free(A4);
  free(B1);
  free(B2);
  free(B3);
  free(B4);
  free(A1B1);
  free(A1B2);
  free(A2B3);
  free(A2B4);
  free(A3B1);
  free(A3B2);
  free(A4B3);
  free(A4B4);

  
  return 0;
   
}

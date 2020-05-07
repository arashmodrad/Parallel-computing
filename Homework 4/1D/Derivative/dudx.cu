#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#define THREADS_PER_BLOCK 15
#define N 15

void checkResult(double *hostRef, double dx, const int size)
{
    double epsilon = 1.0E-3;
    bool match = 1;
    double exact;
    double L2 = 0;
    for (int i = 1; i < size-1; i++) 
    {
      exact = -sin(i*dx);
      //      printf("%d exact %g solution %g, diff = %e\n",i,exact,hostRef[i],exact-hostRef[i]);

      L2+=(hostRef[i] - exact)*(hostRef[i] - exact);

      if (abs(hostRef[i] - exact) > epsilon)
        {
            match = 0;
            printf("Arrays do not match!\n");
            printf("host %10.8f exact %10.8f at location %d\n", hostRef[i],
                   exact, i);
            break;
	}
    }

    
    if (match) {
      printf("Arrays match.\n\n");
      printf("L2 error = %e\n",sqrt(L2)/size);
    }

    return;
}



void initialData(double *ip, double dx, int size)
{

    for (int i = 0; i < size; i++)
    {
      ip[i] = cos(i*dx);
    }

    return;
}

void printer (double *ip, int size)
{
	// generate the function values

	for (int i = 0; i < size; i++)
	{
		printf("matrix element %d is equall to %d\n", i, ip[i]);
	}

	return;
}

//COMPUTE_DUDX_GLOBAL computes derivative of function represented by u 
__global__ void compute_dudx_global(double *dudx, double *u, double dx)
{
	// CREATE GLOBAL INDEX
	int gidx = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (gidx<N + 1) printf("[block = %d, thread = %d] u[%d] = %g\n", blockIdx.x, threadIdx.x, gidx, u[gidx]);


	// COMPUTE DERIVATIVE
	double temp;
	if (gidx<N +1){
		double lam = 0.5 / dx; //compute coefficient
		temp = (u[gidx + 1] - u[gidx - 1])*lam;
	}

	//STORE RESULT BACK TO GLOBAL MEMORY
	dudx[gidx] = temp; //store result in a global memory

	__syncthreads();
	if (gidx<N + 1) printf("[block = %d, thread = %d] u[%d] = %10.8f\t dudx[%d] = %10.8f\n", blockIdx.x, threadIdx.x, gidx, u[gidx], gidx, dudx[gidx]);
}


// using finite difference method on points separated by dx
// kernel takes in the pointers to dudx, u, and a value of dx
// we are using double precision for accurate computations for large N
//__global__ void compute_dudx_global(double *dudx, double *u, double dx)
//{

  // K1) CREATE GLOBAL INDEX
  // use a combination of following variables to identify which memory location is a thread operating on
  // blockIdx.x - index of current block
  // blockDim.x - size of a block (number of threads per block)
  // threadIdx.x - index of a current thread


  // --- your code here


  // K1a) CHECK IF THE DATA ARRIVED TO THE KERNEL
  // if global index not greater than the size of data (N), print the value of the global index and  u to visually inspect it
  // you may also want to print the block and thread index 
  // remember to comment out this print when increasing the value of N


  // --- your code here


  // K2) COMPUTE DERIVATIVE
  // if global index is smaller than N-1 (we only compute derivative between 1 and N-2, and skip end points), 
  // compute the derivative using finite difference formula (u[i+1] - u[i-1])/(2*dx)
  // store the derivative in a temporary variable temp


  // --- your code here

    
  // K3) STORE RESULT BACK TO GLOBAL MEMORY
  // assign the temporary variable to an appropriate location in dudx stored in global memory


  // --- your code here


  // K3a) CHECK IF THE RESULT IS CORRECT
  // if global index less then N-1, print the value u, dudx, and the exact solution -sin(dx*i)
  // you may want to also print the thread index and block index
  // remember to comment out this print when increasing the value of N


  // --- your code here


//}


//this kernel uses shared memory to compute derivative of u
__global__ void compute_dudx_shared(double *dudx, double *u, double dx)
{

  // S1) declare shared array u_loc which has THREAD_PER_BLOCK+2 entries (to accommodate additional points we need for computation of derivatives using finite difference stencil)

  // S2) declare global index (the same as in the global kernel)

  // S3) declare local kernel which will point to a location in the shared memory
  // use threadIdx.x variable only, as the shared memory is accessible only within a block
  // and is shared only by threads within a block
  // remember that thread 0 has to work on memory location 1


  // S4) if global index is smaller than N, copu global memory entry u to local memory u_loc
  // remember to use appropriate indices

  // S4a) thread 0 has to copy a point to its left from global memory to a location in local memory which is also to its left
  
  // S4b) last thread in a block needs to copy a point to its right from global memory to a location in local memory which is also to its right
  // this should not happen if global index for that thread is greater than N, otherwise you are risking segfault by accessing invalid memory


  // S5) Synchronize threads to make sure the copy is complete

  // S5a) Check if the local memory holds the right values
  // thread 0 in each block prints the entire contents of u_loc
    /*
    if(threadIdx.x==0){
      for(int i=0; gidx+i<N && i<blockDim.x; i++){
	printf("[%d] i=%d, u = %g, u_loc = %g\n",blockIdx.x,i,u[gidx+i],u_loc[i]);
      }
      }*/


  // S6) Compute derivative using local memory 
  //almost identical like in global kernel, but index is nw local and variable is not u but u_loc
  // make sure the result is stored in temporary variable

  
  // S7) temporary variable is stored in dudx in the global memory 


}


int main(int argc, char **argv)
{
  printf("%s Starting...\n", argv[0]);
  
  // set up device
  int dev = 0;
  cudaDeviceProp deviceProp;
  CHECK(cudaGetDeviceProperties(&deviceProp, dev));
  printf("Using Device %d: %s\n", dev, deviceProp.name);
  CHECK(cudaSetDevice(dev));
  
  // set up data size of vectors
  int nElem = N + 2;
  double pi = acos(-1.0);
  double dx = 2*pi/(nElem-1);
  printf("Vector size %d dx= %d\n", nElem, dx);
  double iStart, iElaps;

  // 1) ALLOCATE HOST MEMORY
  // use variable nElem and function malloc to allocate host memory for u and dudx
  // it is not a bad idea to use h_u and h_dudx to indicate that those are host variables
  // to avoid writing sizeof() over and over, it is good to define the nBytes variable 
  size_t nBytes = nElem * sizeof(double);
  double *h_dudx, *h_u;
  double host_result;
  double gpu_result;
  h_u = (double*)malloc(nBytes);
  h_dudx = (double*)malloc(nBytes);
  int n = N;
  
  // --- your code here

  
  // 2) INITIALIZE HOST DATA
  // call function initialData (see definition above) to generate values for h_u (host u array)

  initialData(h_u, dx, n + 2);
 
  // printf("initialData Time elapsed %f sec\n", iElaps);
  host_result = 0;
  gpu_result = 0;
  

  // --- your code here


  // 2a) TEST WHETHER THE DATA IS CORRECTLY INITIALIZED
  // you can write a loop over all elements of h_u to inspect the values
  // selecting N to be an odd number will make sure you are hitting x=pi point where u=0
  printer(h_u, nElem);
  //for (int i=0; i<nElem; i++)
  //  printf("initial u[%d] = %g\n",i,h_u[i]); 


  // 3) ALLOCATE DEVICE MEMORY
  // use cudaMalloc to allocate device versions of u and dudx
  // again, the notation d_u and d_dudx is advised
  // use CHECK macro to make sure you are not getting any CUDA errors
  // malloc device global memory
  double *d_u, *d_dudx;
  CHECK(cudaMalloc((double **)&d_u, nBytes));
  CHECK(cudaMalloc((double **)&d_dudx, nBytes));

  
  // --- your code here

  
  // 4) COPY DATA TO THE DEVICE
  // use cudaMemcpy to transfer h_u to d_u
  // use cudaMemset to make sure d_dudx is set to 0
  // remember to use CHECK macro to ensure that errors are being reported
  CHECK(cudaMemcpy(d_u, h_u, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_dudx, h_dudx, nBytes, cudaMemcpyHostToDevice));


  // --- your code here


  // 5) INVOKE THE KERNEL 
  //set number of threads per block
  int nThreads = THREADS_PER_BLOCK;
  //make sure there is enough blocks to cover the entire u array
  int nBlocks =  ((nElem + nThreads - 1) / nThreads);

  //start timer
  iStart = seconds();

  //invoke the kernel
  compute_dudx_global<<<nBlocks, nThreads>>>(d_dudx, d_u, dx);
  // check kernel error
  CHECK(cudaGetLastError());

  //stop timer and report time
  CHECK(cudaDeviceSynchronize());
  iElaps = seconds() - iStart;
  printf("compute_dudx <<<  %d, %d  >>>  Time elapsed %f sec\n", nBlocks,
	 nThreads, iElaps);
  
  // check kernel error
  CHECK(cudaGetLastError());
  
  // 6) COPY THE RESULT BACK TO HOST
  // copy d_dudx to h_dudx using cudaMemcpy
  // use CHECK to make sure errors are reported
  CHECK(cudaMemcpy(h_dudx, d_dudx, nBytes, cudaMemcpyDeviceToHost));

  
  // --- your code here

  
  // 7) CHECK RESULTS
  checkResult(h_dudx,dx,N);
  

  // 8) FREE MEMORY
  // use cudaFree to free d_u and d_dudx
  // use CHECK to make sure there is no errors
  // use free function to free h_u and h_dudx

  // --- your code here
  CHECK(cudaFree(d_u));
  CHECK(cudaFree(d_dudx));
  
  // free host memory
  free(h_u);
  free(h_dudx);

  return(0);
}

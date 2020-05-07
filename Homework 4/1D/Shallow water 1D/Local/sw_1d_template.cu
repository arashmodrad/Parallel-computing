# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>
# include "common.h"

#define THREAD  1024


int main ( int argc, char *argv[] );
void initial_conditions ( int nx, double dx, double x_length, double x[], double h[], double uh[]);

//utilities
void getArgs(int *nx, double *dt, double *x_length, double *t_final, int *THREADS_PER_BLOCK, int argc, char *argv[]);
void write_results ( char *output_filename, int n, double x[], double h[], double uh[]);

__global__ void shallow_water_kernel(double *d_h, double *d_uh, double lambda, int nx, double g, int nBlocks)
{
  //declare global index
	__shared__ int temp_fh[THREAD + 2]; //define a shared array with entries for all threads in a block
	__shared__ int temp_fuh[THREAD + 2]; //define a shared array with entries for all threads in a block
  int gidx = blockIdx.x * blockDim.x + threadIdx.x + 1;
  int i, j;


  // K1) Compute interior fluxes - switch a loop for an if statement
  //for ( i = 1; i < nx+1; i++ )
  if(gidx<nx + 1) // want to do only interior nodes
    {
	  temp_fh[threadIdx.x + 1] = d_uh[gidx]; //flux for the height equation: u*h
	  temp_fuh[threadIdx.x + 1] = d_uh[gidx] * d_uh[gidx] / d_h[gidx] + 0.5*g*d_h[gidx] * d_h[gidx]; //flux for the momentum equation: u^2*h + 0.5*g*h^2
    }
  // compute local ghosts
  //most left
  if (threadIdx.x == 0 && blockIdx.x > 0)
  {
	  i = 0;
	  temp_fh[i] = d_uh[gidx - 1];
	  temp_fuh[i] = d_uh[gidx - 1] * d_uh[gidx - 1] / d_h[gidx - 1] + 0.5*g*d_h[gidx - 1] * d_h[gidx - 1];
  }
  //most right
  if (threadIdx.x == (THREAD - 1) && blockIdx.x < nBlocks)
  {
	  i = THREAD + 1;
      temp_fh[i] = d_uh[gidx + 1];
	  temp_fuh[i] = d_uh[gidx + 1] * d_uh[gidx + 1] / d_h[gidx + 1] + 0.5*g*d_h[gidx + 1] * d_h[gidx + 1];
  }
  // K2) Compute ghost fluxes (need ghost values) - switch loops for if statements
  //left ghost 
  if(gidx==1){
    i=0;
	temp_fh[i] = d_uh[i];
	temp_fuh[i] = d_uh[i] * d_uh[i] / d_h[i] + 0.5*g*d_h[i] * d_h[i];
  }
  if(gidx==nx){ 
    //right ghost
    j = nx+1;
	temp_fh[threadIdx.x + 1] = d_uh[j];
	temp_fuh[threadIdx.x + 1] = d_uh[j] * d_uh[j] / d_h[j] + 0.5*g*d_h[j] * d_h[j];
  }

  // K3) synchronize threads an print results
  __syncthreads();
  /*if (threadIdx.x == 0)
  {
	  for (int i = 0; gidx + i < nx && i < blockDim.x + 2; i++)
	  {
		  printf("I am alive");

		  //printf("[%d] i=%d, d_uh = %g, temp_fh = %g\n", blockIdx.x, i, d_uh[gidx + i], temp_fh[i]);
	  }
  }*/

  // K4) Compute updated variables - swtich a loop for an if statement
  // switch i index to global index, or store gidx as i if you want to save typing
  __shared__ int temp_hm[THREAD]; //define a shared array with entries for all threads in a block
  __shared__ int temp_uhm[THREAD]; //define a shared array with entries for all threads in a block
  int nloc = THREAD + 2;
  if (gidx > 0 && gidx  < nx + 1)
  {
	  i = gidx;
	  j = threadIdx.x + 1;
	  temp_hm[j - 1] = 0.5*(d_h[i + 1] + d_h[i - 1]) - lambda * (temp_fh[j + 1] - temp_fh[j - 1]);
	  temp_uhm[j - 1] = 0.5*(d_uh[i + 1] + d_uh[i - 1]) - lambda * (temp_fuh[j + 1] - temp_fuh[j - 1]);
  }

  /*for ( i = 1; i < nx; i++ )
    {
      hm[i] = 0.5*(d_h[i+1]+d_h[i-1]) - lambda * ( fh[i+1] - fh[i-1] );
      uhm[i] = 0.5*(d_uh[i+1]+d_uh[i-1]) - lambda * ( fuh[i+1] - fuh[i-1] );
    }
	*/
  // K5) Update the boundary conditions - only first and last thread (0 and nx)
  // will have to do this
  if (gidx == 1)
  {
	  d_h[0] = temp_hm[0];
	  d_uh[0] = -temp_uhm[0];
	  //temp_hm[0] = temp_hm[1];
	  //temp_uhm[0] = -temp_uhm[1];
  }
  if (gidx == nx + 1)
  {
	  d_h[nx + 1] = temp_hm[threadIdx.x];
	  d_uh[nx + 1] = -temp_uhm[threadIdx.x];
	  //temp_hm[nx + 1] = temp_hm[nx];
	  //temp_uhm[nx + 1] = -temp_uhm[nx];
  }
  
  
  // K6) Synchronize threads to make sure updates are computed everywhere
  __syncthreads();

  // K7) Update state variables - put the value from hm and uhm to solution arrays
  // replace loop with an if
  if (gidx > 0 && gidx < nx + 1)
  {
	  i = gidx;
	  j = threadIdx.x;
	  d_h[i] = temp_hm[j];
	  d_uh[i] = temp_uhm[j];
  }
  /*for (i = 1; i < nx+1; i++){
    d_h[i] = hm[i];
    d_uh[i] = uhm[i];
  }*/
  //store temporary to global

}

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:
    MAIN is the main program for SHALLOW_WATER_1D.

  Discussion:
    SHALLOW_WATER_1D approximates the 1D shallow water equations.
    The version of the shallow water equations being solved here is in
    conservative form, and omits the Coriolis force.  The state variables
    are H (the height) and UH (the mass velocity).

    The equations have the form
      dH/dt + d UH/dx = 0
      d UH/dt + d ( U^2 H + 1/2 g H^2 )/dx = 0

    Here U is the ordinary velocity, U = UH/H, and g is the gravitational
    acceleration.
    The initial conditions are used to specify ( H, UH ) at an equally
    spaced set of points, and then the Lax-Friedrichs method is used to advance
    the solution until a final time t_final, with
    boundary conditions supplying the first and last spatial values.
    Some input values will result in an unstable calculation that
    quickly blows up.  This is related to the Courant-Friedrichs-Levy
    condition, which requires that DT be small enough, relative to DX and
    the velocity, that information cannot cross an entire cell.

    A "reasonable" set of input quantities is
      shallow_water_1d 401 0.002 10.0 0.2

  Licensing:
    This code is distributed under the GNU LGPL license.

  Modified:
    26 March 2019 by Michal A. Kopera

  Parameters:
    Input, integer NX, the number of spatial nodes.
    Input, integer DT, the size of a time step.
    Input, real X_LENGTH, the length of the region.
    Input, real T_FINAL, the final time of simulation.

    Output, real X[NX], the X coordinates.
    Output, real H[NX], the height for all space points at time t_final.
    Output, real UH[NX], the mass velocity (discharge) for all space points at time t_final.
*/
{
printf("%s Starting...\n", argv[0]);

// set up device
int dev = 0;
cudaDeviceProp deviceProp;
  CHECK(cudaGetDeviceProperties(&deviceProp, dev));
  printf("Using Device %d: %s\n", dev, deviceProp.name);
  CHECK(cudaSetDevice(dev));
  int THREADS_PER_BLOCK;
  double dx;
  double dt;
  double g = 9.81; //[m^2/s] gravitational constant
  double *h;
  //double *fh;
  //double *hm;
  int nx;
  double t_final;
  double *uh;
 // double *fuh;
  //double *uhm;
  double *x;
  double x_length, time;


printf ( "\n" );
printf ( "SHALLOW_WATER_1D\n" );
printf ( "\n" );


  //get command line arguments
getArgs(&nx, &dt, &x_length, &t_final, &THREADS_PER_BLOCK, argc, argv);

    printf ( "  NX = %d\n", nx );
    printf ( "  DT = %g\n", dt );
    printf ( "  X_LENGTH = %g\n", x_length );
    printf ( "  T_FINAL = %g\n", t_final );
	printf("  THREADS_PER_BLOCK = %d\n", THREADS_PER_BLOCK);

  //M1) Allocate space (nx+2) long, to accound for ghosts
  //height array
    size_t nBytes = (nx+2)*sizeof(double);
  h = ( double * ) malloc ( nBytes );
  //discharge array
  uh = ( double * ) malloc ( nBytes);
  // location array
  x = ( double * ) malloc ( nx * sizeof ( double ) );

  //Define the locations of the nodes and time steps and the spacing.
  dx = x_length / ( double ) ( nx );

  // M2) Apply the initial conditions.
  initial_conditions ( nx, dx, x_length,  x, h, uh);

  // M3) Write initial condition to a file
  write_results((char *)"sw1d_cuda_init.dat",nx,x,h,uh);

  double lambda = 0.5*dt/dx;

  // M4) allocate device memory to hold h, uh, hm, uhm, fh, fuh
  // some memory will be moved between device and host (h, uh), so
  // I encourage you to use d_h and d_uh notation for those variables
  // variables which live only on the device do not need that
  // In principle, they can be allocated from the device, but since we
  // need to allocate them only once, I suggest to do it from the host
  // Use cudaMalloc and nBytes to allocate d_h, d_uh, hm, uhm, fh, fuh
  double *d_h, *d_uh;
  CHECK(cudaMalloc((double **)&d_h, nBytes));
  CHECK(cudaMalloc((double **)&d_uh, nBytes));
  //CHECK(cudaMalloc((double **)&hm, nBytes));
  //CHECK(cudaMalloc((double **)&uhm, nBytes));
  //CHECK(cudaMalloc((double **)&fh, nBytes));
  //CHECK(cudaMalloc((double **)&fuh, nBytes));

  // M5) transfer data from host to device
  // use cudaMemcpy to transfer h and uh to device
  CHECK(cudaMemcpy(d_h, h, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_uh, uh, nBytes, cudaMemcpyHostToDevice));

  // M6) initialize kernel data
  int nThreads = THREADS_PER_BLOCK;
  int nBlocks =  ((nx + nThreads - 1) / nThreads);

  //start timer
  double iStart, iElaps;
  iStart = seconds();

  time=0;
  while (time<t_final) //for ( it = 1; it <= nt; it++ )
    {
      //  Take a time step
      time=time+dt;
      //printf("time = %f\n",time);
	  //const int *THREAD = &THREADS_PER_BLOCK;
      // M7)  call computational kerne
	  shallow_water_kernel << <nBlocks, nThreads >> >(d_h, d_uh, lambda, nx, g, nBlocks);
      
      // check kernel error
      CHECK(cudaDeviceSynchronize());
      CHECK(cudaGetLastError());


      // PARALLELIZATION ENDS HERE
    }

  CHECK(cudaDeviceSynchronize());
  iElaps = seconds() - iStart;
  printf("shallow_water <<<  %d, %d  >>>  Time elapsed %f sec\n", nBlocks,
	 nThreads, iElaps);
    
  // M8) copy kernel result back to host side
  // use cudaMemcpy to get values from d_h and d_uh to the host
  CHECK(cudaMemcpy(h, d_h, nBytes, cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(uh, d_uh, nBytes, cudaMemcpyDeviceToHost));

  // M9) Write data to file
  write_results((char *)"sw1d_cuda_final.dat",nx,x,h,uh);
  
  // M10) Free host memory.
  free ( h );
  free ( uh );
  free ( x );

  // M11) Free device memory using cudaFree


 //Terminate.
  printf ( "\n" );
  printf ( "SHALLOW_WATER_1D:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );

  return 0;
}
/******************************************************************************/

void initial_conditions ( int nx, double dx, double x_length, double x[], double h[], 
			  double uh[])

/******************************************************************************/

{
  int i;
  
  for ( i = 1; i < nx+1; i++ )
    {
      x[i-1] = -x_length/2+dx/2+(i-1)*dx;
      double xx = x[i-1];
      h[i] = 1.0 + 0.4*exp ( -5 * ( xx*xx) );
    }
  for ( i = 1; i < nx+1; i++ )
    {
      uh[i] = 0.0;
    }
  h[0] = h[1];
  h[nx+1]=h[nx];
  uh[0] = 0.0;
  uh[nx+1] = 0.0;
  return;
}
/******************************************************************************/


void write_results ( char *output_filename, int n, double x[], double h[], double uh[])
/******************************************************************************/

{
  int j;
  FILE *output;
  
  //Open the file.
  output = fopen ( output_filename, "wt" );
    
  if ( !output ){
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "WRITE_RESULTS - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
    
  //Write the data.
  for ( j = 1; j < n+1; j++ )	{
    fprintf ( output, "  %24.16g\t %24.16g\t %24.16g\n", x[j-1], h[j], uh[j]);
  }
  
  //Close the file.
  fclose ( output );
  return;
}
/******************************************************************************/

void getArgs(int *nx, double *dt, double *x_length, double *t_final, int *THREADS_PER_BLOCK, int argc, char *argv[])
{

    /*
      Get the quadrature file root name:
    */
    if ( argc <= 1 ){
      *nx = 401;
    }else{
      *nx = atoi ( argv[1] );
    }
    
    if ( argc <= 2 ){
      *dt = 0.002;
    }else{
      *dt = atof ( argv[2] );
    }
    
    if ( argc <= 3 ){
      *x_length = 10.0;
    }else{
      *x_length = atof ( argv[3] );
    }
    
    if ( argc <= 4 ){
      *t_final = 0.5;
    }else{
      *t_final = atof ( argv[4] );
    }
	if (argc <= 5){
		*THREADS_PER_BLOCK = 32;
	}
	else{
		*THREADS_PER_BLOCK = atof(argv[5]);
	}
  
}

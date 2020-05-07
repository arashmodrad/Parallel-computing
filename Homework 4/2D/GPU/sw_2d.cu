# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>
# include <mpi.h>
# include "common.h"

#define ID_2D(i,j,nx) ((i)*(nx+2)+(j))
#define THREADS_PER_BLOCK_X 32
#define THREADS_PER_BLOCK_Y 32

int main ( int argc, char *argv[] );
void initial_conditions ( int nx, int ny, double dx, double dy,  double x_length, double y_length, double x[],double y[], double h[], double uh[] ,double vh[]);

//utilities
void getArgs(int *nx, double *dt, double *x_length, double *t_final, int argc, char *argv[]);
void write_results ( char *output_filename, int nx, int ny, double x[],double y[], double h[], double uh[], double vh[]);

__global__ void shallow_water_2d_kernel(double *d_h, double *d_uh, double *d_vh, 
					double *fh, double *fuh, double *fvh,
					double *gh, double *guh, double *gvh,
					double *hm, double *uhm, double *vhm,
					double lambda_x, double lambda_y, 
					int nx, int ny, double g)
{
  int gidx_x = blockIdx.x * blockDim.x + threadIdx.x + 1;//global index in x
  int gidx_y = blockIdx.y * blockDim.y + threadIdx.y + 1;//global index in y
  int gidx = ID_2D(gidx_y,gidx_x,nx); //global linear index
  int i,j, id;
  int id_left;
  int id_right;
  int id_bottom;
  int id_top;


  //if(gidx_x<nx+2 && gidx_y<ny+2) printf("block = (%d,%d), thread = (%d,%d), gidx_x = %d, gidx_y = %d, gidx = %d, x = %g, y = %d\n",blockIdx.x,blockIdx.y,threadIdx.x,threadIdx.y,gidx_x,gidx_y,gidx,d_h[gidx],d_uh[gidx]);

    //Compute interior fluxes 
    //for ( i = 1; i < nx_loc+1; i++ )
    if(gidx_x<nx+1 && gidx_y<ny+1) // want to do only interior nodes
	{
	  fh[gidx] = d_uh[gidx]; //flux for the height equation: u*h
	  fuh[gidx] = d_uh[gidx]*d_uh[gidx]/d_h[gidx] + 0.5*g*d_h[gidx]*d_h[gidx]; //flux for the momentum equation: u^2*h + 0.5*g*h^2
	  fvh[gidx] = d_uh[gidx]*d_vh[gidx]/d_h[gidx]; //flux for the momentum equation: u*v**h 
	  gh[gidx] = d_vh[gidx]; //flux for the height equation: v*h
	  guh[gidx] = d_uh[gidx]*d_vh[gidx]/d_h[gidx]; //flux for the momentum equation: u*v**h 
	  gvh[gidx] = d_vh[gidx]*d_vh[gidx]/d_h[gidx] + 0.5*g*d_h[gidx]*d_h[gidx]; //flux for the momentum equation: v^2*h + 
	}

      //Compute ghost fluxes (need ghost values)
      //left ghost
    if(gidx_x==1){
      j=0;
      id=ID_2D(gidx_y,j,nx);
      fh[id] = d_uh[id];
      fuh[id] = d_uh[id]*d_uh[id]/d_h[id] + 0.5*g*d_h[id]*d_h[id];
      fvh[id] = d_uh[id]*d_vh[id]/d_h[id]; //flux for the momentum equation: u*v**h 
      gh[id] = d_vh[id]; //flux for the height equation: v*h
      guh[id] = d_uh[id]*d_vh[id]/d_h[id]; //flux for the momentum equation: u*v**h 
      gvh[id] = d_vh[id]*d_vh[id]/d_h[id] + 0.5*g*d_h[id]*d_h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
      }

      //right ghost
    if(gidx_x==nx){
      j = nx+1;
      id=ID_2D(gidx_y,j,nx);
      fh[id] = d_uh[id];
      fuh[id] = d_uh[id]*d_uh[id]/d_h[id] + 0.5*g*d_h[id]*d_h[id];
      fvh[id] = d_uh[id]*d_vh[id]/d_h[id]; //flux for the momentum equation: u*v**h 
      gh[id] = d_vh[id]; //flux for the height equation: v*h
      guh[id] = d_uh[id]*d_vh[id]/d_h[id]; //flux for the momentum equation: u*v**h 
      gvh[id] = d_vh[id]*d_vh[id]/d_h[id] + 0.5*g*d_h[id]*d_h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
    }

    //bottom ghost
    if(gidx_y==1){
      i = 0;
      id=ID_2D(i,gidx_x,nx);
      fh[id] = d_uh[id];
      fuh[id] = d_uh[id]*d_uh[id]/d_h[id] + 0.5*g*d_h[id]*d_h[id];
      fvh[id] = d_uh[id]*d_vh[id]/d_h[id]; //flux for the momentum equation: u*v**h 
      gh[id] = d_vh[id]; //flux for the height equation: v*h
      guh[id] = d_uh[id]*d_vh[id]/d_h[id]; //flux for the momentum equation: u*v**h 
      gvh[id] = d_vh[id]*d_vh[id]/d_h[id] + 0.5*g*d_h[id]*d_h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
    }

      //top ghost
    if(gidx_y==ny){
      i = ny+1;
      id=ID_2D(i,gidx_x,nx);
	fh[id] = d_uh[id];
	fuh[id] = d_uh[id]*d_uh[id]/d_h[id] + 0.5*g*d_h[id]*d_h[id];
	fvh[id] = d_uh[id]*d_vh[id]/d_h[id]; //flux for the momentum equation: u*v**h 
	gh[id] = d_vh[id]; //flux for the height equation: v*h
	guh[id] = d_uh[id]*d_vh[id]/d_h[id]; //flux for the momentum equation: u*v**h 
	gvh[id] = d_vh[id]*d_vh[id]/d_h[id] + 0.5*g*d_h[id]*d_h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
      }

    __syncthreads();

    //Compute updated variables
    if(gidx_x<nx+1 && gidx_y<ny+1)
      {
	id_left=ID_2D(gidx_y,gidx_x-1,nx);
	id_right=ID_2D(gidx_y,gidx_x+1,nx);
	id_bottom=ID_2D(gidx_y-1,gidx_x,nx);
	id_top=ID_2D(gidx_y+1,gidx_x,nx);

	hm[gidx] = 0.25*(d_h[id_left]+d_h[id_right]+d_h[id_bottom]+d_h[id_top]) 
	  - lambda_x * ( fh[id_right] - fh[id_left] ) 
	  - lambda_y * ( gh[id_top] - gh[id_bottom] );
	uhm[gidx] = 0.25*(d_uh[id_left]+d_uh[id_right]+d_uh[id_bottom]+d_uh[id_top]) 
	  - lambda_x * ( fuh[id_right] - fuh[id_left] ) 
	  - lambda_y * ( guh[id_top] - guh[id_bottom] );
	vhm[gidx] = 0.25*(d_vh[id_left]+d_vh[id_right]+d_vh[id_bottom]+d_vh[id_top]) 
	  - lambda_x * ( fvh[id_right] - fvh[id_left] ) 
	  - lambda_y * ( gvh[id_top] - gvh[id_bottom] );
      }

    __syncthreads();

      //Update the boundary conditions
      //left
    if(gidx_x==1){
      j=1;
      i=gidx_y;
   
      id_left = ID_2D(i,j-1,nx);
      hm[id_left]  =   hm[gidx];
      uhm[id_left] = - uhm[gidx];
      vhm[id_left] =   vhm[gidx];
    }

      //right
    if(gidx_x==nx){
      j=nx;
      i=gidx_y;
      id_right = ID_2D(i,j+1,nx);
      hm[id_right]  =   hm[gidx];
      uhm[id_right] = - uhm[gidx];
      vhm[id_right] =   vhm[gidx];
    }

      //bottom
    if(gidx_y==1){
      i=1;
      j = gidx_x;
      id_bottom = ID_2D(i-1,j,nx);
      hm[id_bottom]  =   hm[gidx];
      uhm[id_bottom] =   uhm[gidx];
      vhm[id_bottom] = - vhm[gidx];
    }

      //top
    if(gidx_y==ny){
      i=ny;
      j = gidx_x;
      id_top = ID_2D(i+1,j,nx);
      hm[id_top]  =   hm[gidx];
      uhm[id_top] =   uhm[gidx];
      vhm[id_top] = - vhm[gidx];
      }


    __syncthreads();
      //update interior state variables
    if((gidx_x<nx+2) && (gidx_y<ny+2)){
	d_h[gidx] = hm[gidx];
	d_uh[gidx] = uhm[gidx];
	d_vh[gidx] = vhm[gidx];
      }

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

  double dx, dy;
  double dt;
  double g = 9.81; //[m^2/s] gravitational constant
  double *h;
  double *fh;
  double *gh;
  double *hm;
  int nx, ny;
  double t_final;
  double *uh;
  double *fuh;
  double *guh;
  double *uhm;
  double *vh;
  double *fvh;
  double *gvh;
  double *vhm;
  double *x;
  double *y;
  double x_length, y_length, time;


printf ( "\n" );
printf ( "SHALLOW_WATER_1D\n" );
printf ( "\n" );


  //get command line arguments
  getArgs(&nx, &dt, &x_length, &t_final, argc, argv);

    printf ( "  NX = %d\n", nx );
    printf ( "  DT = %g\n", dt );
    printf ( "  X_LENGTH = %g\n", x_length );
    printf ( "  T_FINAL = %g\n", t_final );

    ny = nx; //assuming perfect square 
    y_length = x_length; 

  //Allocate space (nx+2) long, to accound for ghosts
  //height array
    size_t nBytes = (nx+2)*(ny+2)*sizeof(double);
  h = ( double * ) malloc ( nBytes );
  //discharge array
  uh = ( double * ) malloc ( nBytes);
  //discharge array
  vh = ( double * ) malloc ( nBytes);
  // location array
  x = ( double * ) malloc ( nx * sizeof ( double ) );
  // location array
  y = ( double * ) malloc ( ny * sizeof ( double ) );

  //Define the locations of the nodes and time steps and the spacing.
  dx = x_length / ( double ) ( nx );
  dy = y_length / ( double ) ( ny );

  //Apply the initial conditions.
  initial_conditions ( nx, ny, dx, dy, x_length, y_length,  x, y, h, uh, vh);
  
  //Write initial condition to a file
  write_results((char *)"sw2d_cuda_init.dat",nx,ny,x,y,h,uh,vh);
  
  double lambda_x = 0.5*dt/dx; 
  double lambda_y = 0.5*dt/dy; 

  //allocate device memory
  double *d_h, *d_uh, *d_vh;
  CHECK(cudaMalloc((double **)&d_h, nBytes));
  CHECK(cudaMalloc((double **)&d_uh, nBytes));
  CHECK(cudaMalloc((double **)&d_vh, nBytes));
  CHECK(cudaMalloc((double **)&hm, nBytes));
  CHECK(cudaMalloc((double **)&uhm, nBytes));
  CHECK(cudaMalloc((double **)&vhm, nBytes));
  CHECK(cudaMalloc((double **)&fh, nBytes));
  CHECK(cudaMalloc((double **)&fuh, nBytes));
  CHECK(cudaMalloc((double **)&fvh, nBytes));
  CHECK(cudaMalloc((double **)&gh, nBytes));
  CHECK(cudaMalloc((double **)&guh, nBytes));
  CHECK(cudaMalloc((double **)&gvh, nBytes));

  //transfer data from host to device
  CHECK(cudaMemcpy(d_h, h, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_uh, uh, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_vh, vh, nBytes, cudaMemcpyHostToDevice));

  //initialize kernel data
  //  int nThreads = THREADS_PER_BLOCK;
  //int nBlocks =  ((nx + nThreads - 1) / nThreads);

  dim3 nThreads(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y);
  dim3 nBlocks((nx+nThreads.x-1)/nThreads.x, (ny+nThreads.y-1)/nThreads.y);

  //start timer
  double iStart, iElaps;
  iStart = seconds();

  time=0;
  while (time<t_final) //for ( it = 1; it <= nt; it++ )
    {
      //  Take a time step
      time=time+dt;
      //printf("time = %f\n",time);

      // call computational kernelsPARALLELIZATION STARTS HERE
      //      shallow_water_kernel<<<nBlocks,nThreads>>>(d_h, d_uh, fh, fuh, hm, uhm, lambda, nx, g);
      shallow_water_2d_kernel<<<nBlocks,nThreads>>>(d_h, d_uh, d_vh, 
						    fh, fuh, fvh, 
						    gh, guh, gvh,
						    hm, uhm, vhm, 
						    lambda_x, lambda_y,
						    nx, ny, g);

      // check kernel error
      CHECK(cudaDeviceSynchronize());
      CHECK(cudaGetLastError());


      // PARALLELIZATION ENDS HERE
    }

  CHECK(cudaDeviceSynchronize());
  iElaps = seconds() - iStart;
  printf("shallow_water <<< (%d, %d), (%d, %d)  >>>  Time elapsed %f sec\n", nBlocks.x, nBlocks.y,
	 nThreads.x, nThreads.y, iElaps);
    
  // copy kernel result back to host side
  CHECK(cudaMemcpy(h, d_h, nBytes, cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(uh, d_uh, nBytes, cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(vh, d_vh, nBytes, cudaMemcpyDeviceToHost));

  // Write data to file
  write_results((char *)"sw2d_cuda_final.dat",nx,ny,x,y,h,uh,vh);
  
  //Free memory.
  free ( h );
  free ( uh );
  free ( vh );
  free ( x );
  free ( y );

  CHECK(cudaFree(d_h));
  CHECK(cudaFree(d_uh));
  CHECK(cudaFree(d_vh));
  CHECK(cudaFree(hm));
  CHECK(cudaFree(uhm));
  CHECK(cudaFree(vhm));
  CHECK(cudaFree(fh));
  CHECK(cudaFree(fuh));
  CHECK(cudaFree(fvh));
  CHECK(cudaFree(gh));
  CHECK(cudaFree(guh));
  CHECK(cudaFree(gvh));


 //Terminate.
  printf ( "\n" );
  printf ( "SHALLOW_WATER_1D:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );

  return 0;
}
/******************************************************************************/

void initial_conditions ( int nx, int ny, double dx, double dy,  double x_length, double y_length, double x[],double y[], double h[], double uh[] ,double vh[])

/******************************************************************************/

{
  int i,j, id, id1;

  for ( i = 1; i < nx+1; i++ )
    {
      x[i-1] = -x_length/2+dx/2+(i-1)*dx;
      y[i-1] = -y_length/2+dy/2+(i-1)*dy;
    }

  for ( i = 1; i < nx+1; i++ )
    for( j = 1; j < ny+1; j++)
      {
	double xx = x[j-1];
	double yy = y[i-1];
	id=ID_2D(i,j,nx);
	h[id] = 1.0 + 0.4*exp ( -5 * ( xx*xx + yy*yy) );
      }
  
  for ( i = 1; i < nx+1; i++ )
    for( j = 1; j < ny+1; j++)
      {
	id=ID_2D(i,j,nx);
	uh[id] = 0.0;
	vh[id] = 0.0;
      }
  
  //set boundaries
  //bottom
  i=0;
  for( j = 1; j < nx+1; j++)
    {
      id=ID_2D(i,j,nx);
      id1=ID_2D(i+1,j,nx);
      
      h[id] = h[id1];
      uh[id] = 0.0;
      vh[id] = 0.0;
    }
  
  //top
  i=nx+1;
  for( j = 1; j < nx+1; j++)
    {
      id=ID_2D(i,j,nx);
      id1=ID_2D(i-1,j,nx);
      
      h[id] = h[id1];
      uh[id] = 0.0;
      vh[id] = 0.0;
    }
  
  //left
  j=0;
  for( i = 1; i < ny+1; i++)
    {
      id=ID_2D(i,j,nx);
      id1=ID_2D(i,j+1,nx);
      
      h[id] = h[id1];
      uh[id] = 0.0;
      vh[id] = 0.0;
    }
  
  //right
  j=nx+1;
  for( i = 1; i < ny+1; i++)
    {
      id=ID_2D(i,j,nx);
      id1=ID_2D(i,j-1,nx);
      
      h[id] = h[id1];
      uh[id] = 0.0;
      vh[id] = 0.0;
    }
  
  return;
}
/******************************************************************************/
 
 
 
 void write_results ( char *output_filename, int nx, int ny, double x[],double y[], double h[], double uh[], double vh[])
 /******************************************************************************/
   
 {
   int i, j, id;
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
   for ( i = 1; i < ny+1; i++ )	{
     for ( j = 1; j < nx+1; j++ )	{
       id=ID_2D(i,j,nx);
	fprintf ( output, "  %24.16g\t%24.16g\t%24.16g\t %24.16g\t %24.16g\n", x[j-1], y[i-1],h[id], uh[id], vh[id]);
     }
   }
   
  //Close the file.
  fclose ( output );
  return;
}
/******************************************************************************/

void getArgs(int *nx, double *dt, double *x_length, double *t_final, int argc, char *argv[])
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
  
}

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>
# include <mpi.h>

#define ID_2D(i,j,nx) ((i)*(nx+2)+(j))

int main ( int argc, char *argv[] );
void initial_conditions ( int nx, int ny, double dx, double dy,  double x_length, double x[],double y[], double h[], double uh[] ,double vh[], int irank, int q );
void compute_flux_2d_point(double fh, double fuh, double fvh, 
			   double gh, double guh, double gvh, 
			   double h, double uh, double vh, double g);

//utilities
void getArgs_mpi(int *nx, double *dt, double *x_length, double *t_final, int argc, char *argv[], int irank, MPI_Comm comm);

void write_results ( char *output_filename, int nx, int nx_loc, int ny, int ny_loc, double x[], double y[], double h[], double uh[], double vh[], int irank, int nproc );


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

  Author:
    John Burkardt

  Reference:
    Cleve Moler,
    "The Shallow Water Equations",
    Experiments with MATLAB.

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
  double dx;
  double dy;
  double dt;
  double g = 9.81; //[m^2/s] gravitational constant
  double *h;
  double *fh;
  double *gh;
  double *hm;
  int i,j, id, id_left, id_right, id_bottom, id_top;
  int nx, ny, nx_loc, ny_loc;
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
  double x_length, time;
  int irank, nproc;

  //initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  if(irank==0) {
    printf ( "\n" );
    printf ( "SHALLOW_WATER_2D\n" );
    printf ( "\n" );
  }

  //get command line arguments
  getArgs_mpi(&nx, &dt, &x_length, &t_final, argc, argv, irank, MPI_COMM_WORLD);
  if(irank==0){
    printf ( "  NX = %d\n", nx );
    printf ( "  DT = %g\n", dt );
    printf ( "  X_LENGTH = %g\n", x_length );
    printf ( "  T_FINAL = %g\n", t_final );
	printf("  nproc = %d\n", nproc);
  }

  //divide data among processors
  int q = (int) sqrt((double) nproc);
  ny=nx;
  nx_loc = nx/q;
  ny_loc = ny/q;

  //Allocate space (nx+2)((nx+2) long, to accound for ghosts
  //height array
  h = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  hm = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  fh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  gh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  //x momentum array
  uh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  uhm = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  fuh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  guh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  //y momentum array
  vh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  vhm = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  fvh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  gvh = ( double * ) malloc ( (nx_loc+2)*(nx_loc+2) * sizeof ( double ) );
  // location arrays
  x = ( double * ) malloc ( nx_loc * sizeof ( double ) );
  y = ( double * ) malloc ( nx_loc * sizeof ( double ) );

  //Define the locations of the nodes and time steps and the spacing.
  dx = x_length / ( double ) ( nx );
  dy = x_length / ( double ) ( nx );

  //Apply the initial conditions.
  initial_conditions ( nx_loc, ny_loc, dx, dy, x_length,  x, y, h, uh, vh, irank, q );

  //Write initial condition to a file
  //write_results("sw2d_init.dat",nx,nx_loc,nx,ny_loc,x,y,h,uh,vh, irank,nproc);
  
  double lambda_x = 0.5*dt/dx;
  double lambda_y = 0.5*dt/dy;

  MPI_Request request_recv_left;//, request_uh_recv_left;
  //  double h_left_recv, uh_left_recv;
  double buf_left_recv[3 * ny_loc], buf_right_recv[3 * ny_loc], buf_left_send[3 * ny_loc], buf_right_send[3 * ny_loc];
  double buf_down_recv[3 * nx_loc], buf_up_recv[3 * nx_loc], buf_down_send[3 * nx_loc], buf_up_send[3 * nx_loc];
  MPI_Request request_recv_right;//, request_uh_recv_right;
  //  double h_right_recv, uh_right_recv;
  MPI_Request request_recv_down;//, request_uh_recv_down;
  MPI_Request request_recv_up;//, request_uh_recv_up;
  MPI_Request request_send_left;//, request_uh_send_left;
  //  double h_left_send, uh_left_send;
  MPI_Request request_send_right;//, request_uh_send_right;
  //  double h_right_send, uh_right_send;
  MPI_Request request_send_down;//, request_uh_send_down;
  MPI_Request request_send_up;//, request_uh_send_up;
  MPI_Status status;

  //SET-UP CARTESIAN COMMUNICATOR
  int Snproc = (int)sqrt((double)nproc);
  MPI_Comm Comm_cart;
  int dims[] = { Snproc, Snproc };
  int rank_left, rank_right, rank_down, rank_up;
  int ndims = 2;
  int wrap_around[] = { 0, 0 };//no periodicity
  int coords[2]; //2 Dimension topology so 2 coordinates 

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, 0, &Comm_cart);


  //printf("previous rank: %d, left: %d, right: %d\n",irank,rank_left,rank_right);
  //printf("previous rank: %d, up: %d, down: %d\n", irank, rank_up, rank_down);
  
  // Determines process coords in cartesian topology given rank in group
  MPI_Cart_coords(Comm_cart, irank, 2, coords);

  //    a) shift A left by the displacement equal to row number
  MPI_Cart_shift(Comm_cart, 1, 1, &rank_left, &rank_right);
  MPI_Cart_shift(Comm_cart, 0, 1, &rank_down, &rank_up);

  //printf("rank: %d, left: %d, right: %d\n", irank, rank_left, rank_right);
  //printf("rank: %d, up: %d, down: %d\n", irank, rank_up, rank_down);
  //start timer
  double wtime;
  MPI_Barrier(MPI_COMM_WORLD);
  wtime = MPI_Wtime();


  time=0;
  while (time < t_final) //for ( it = 1; it <= nt; it++ )
  {
	  //  Take a time step
	  time = time + dt;
	  //printf("time = %f\n",time);
	  // PARALLELIZATION STARTS HERE
	  // 1)Initialize communication #################################################
	  if (rank_left >= 0){
		  //set up receiving from the left
		  MPI_Irecv(buf_left_recv, (3 * ny_loc), MPI_DOUBLE, irank - 1,
			  101, MPI_COMM_WORLD, &request_recv_left);


		  //set up sending to the left
		  for (int k = 1; k < (ny_loc + 1); k++)
		  {
			  buf_left_send[k - 1] = h[ID_2D(k, 1, ny_loc)];
			  buf_left_send[k - 1 + ny_loc] = uh[ID_2D(k, 1, ny_loc)];
			  buf_left_send[k - 1 + (2 * ny_loc)] = vh[ID_2D(k, 1, ny_loc)];
		  }
		  MPI_Isend(buf_left_send, (3 * ny_loc), MPI_DOUBLE, irank - 1,
			  101, MPI_COMM_WORLD, &request_send_left);

	  }

	  if (rank_right >= 0){
		  //set up receiving from the right
		  MPI_Irecv(buf_right_recv, (3 * ny_loc), MPI_DOUBLE, irank + 1,
			  101, MPI_COMM_WORLD, &request_recv_right);


		  //set up sending to the right
		  for (int k = 1; k < (ny_loc + 1); k++)
		  {
			  buf_left_send[k - 1] = h[ID_2D(k, nx_loc, ny_loc)];
			  buf_left_send[k - 1 + ny_loc] = uh[ID_2D(k, nx_loc, ny_loc)];
			  buf_left_send[k - 1 + (2 * ny_loc)] = vh[ID_2D(k, nx_loc, ny_loc)];
		  }
		  MPI_Isend(buf_right_send, (3 * ny_loc), MPI_DOUBLE, irank + 1,
			  101, MPI_COMM_WORLD, &request_send_right);

	  }

	  if (rank_down >= 0){
		  //set up receiving from the bottom
		  MPI_Irecv(buf_down_recv, (3 * nx_loc), MPI_DOUBLE, irank - q,
			  101, MPI_COMM_WORLD, &request_recv_down);

		  //set up sending to the bottom
		  for (int k = 1; k < (nx_loc + 1); k++)
		  {
			  buf_down_send[k - 1] = h[ID_2D(ny_loc, k, nx_loc)];
			  buf_down_send[k - 1 + nx_loc] = uh[ID_2D(ny_loc, k, nx_loc)];
			  buf_down_send[k - 1 + (2 * nx_loc)] = vh[ID_2D(ny_loc, k, nx_loc)];
		  }
		  MPI_Isend(buf_down_send, (3 * nx_loc), MPI_DOUBLE, irank - q,
			  101, MPI_COMM_WORLD, &request_send_down);

	  }

	  if (rank_up >= 0){
		  //set up receiving from the above
		  MPI_Irecv(buf_up_recv, (3 * nx_loc), MPI_DOUBLE, irank + q,
			  101, MPI_COMM_WORLD, &request_recv_up);


		  //set up sending to the above
		  for (int k = 1; k < (nx_loc + 1); k++)
		  {
			  buf_up_send[k - 1] = h[ID_2D(1, k, nx_loc)];
			  buf_up_send[k - 1 + nx_loc] = uh[ID_2D(1, k, nx_loc)];
			  buf_up_send[k - 1 + (2 * nx_loc)] = vh[ID_2D(1, k, nx_loc)];
		  }
		  MPI_Isend(buf_up_send, (3 * ny_loc), MPI_DOUBLE, irank + q,
			  101, MPI_COMM_WORLD, &request_send_up);

	  }
	  /*if (irank == 0)
	  {
	  printf("I am alive to this point 1 \n");
	  }*/
	  // 2)Compute what is possible without ghost elements ##################################
	  // a. Compute fluxes ##################################################################
	  // COMPUTATION GOES HERE
	  //Compute interior fluxes 
	  for (i = 1; i < ny_loc + 1; i++)
		  for (j = 1; j < nx_loc + 1; j++){
		  id = ID_2D(i, j, nx_loc);
		  //	  compute_flux_2d_point(fh[id], fuh[id], fvh[id], gh[id], guh[id], gvh[id], 
		  //			h[id], uh[id], vh[id], g);

		  fh[id] = uh[id]; //flux for the height equation: u*h
		  fuh[id] = uh[id] * uh[id] / h[id] + 0.5*g*h[id] * h[id]; //flux for the momentum equation: u^2*h + 0.5*g*h^2
		  fvh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gh[id] = vh[id]; //flux for the height equation: v*h
		  guh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gvh[id] = vh[id] * vh[id] / h[id] + 0.5*g*h[id] * h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
		  }

	  /*if (irank == 0)
	  {
	  printf("I am alive to this point 2.a \n");
	  }*/
	  // b.Compute updates ##################################################################
	  /*//Compute updated variables (interior)
	  for (i = 2; i < nx_loc; i++)
	  {
	  hm[i] = 0.5*(h[i + 1] + h[i - 1]) - lambda * (fh[i + 1] - fh[i - 1]);
	  uhm[i] = 0.5*(uh[i + 1] + uh[i - 1]) - lambda * (fuh[i + 1] - fuh[i - 1]);
	  vhm[i] = 0.5*(vh[i + 1] + vh[i - 1]) - lambda * (fvh[i + 1] - fvh[i - 1]);
	  }*/
	  //Compute updated variables
	  for (i = 2; i < ny_loc; i++)
		  for (j = 2; j < nx_loc; j++)
		  {
		  id = ID_2D(i, j, nx_loc);
		  id_left = ID_2D(i, j - 1, nx_loc);
		  id_right = ID_2D(i, j + 1, nx_loc);
		  id_bottom = ID_2D(i - 1, j, nx_loc);
		  id_top = ID_2D(i + 1, j, nx_loc);
		  hm[id] = 0.25*(h[id_left] + h[id_right] + h[id_bottom] + h[id_top])
			  - lambda_x * (fh[id_right] - fh[id_left])
			  - lambda_y * (gh[id_top] - gh[id_bottom]);
		  uhm[id] = 0.25*(uh[id_left] + uh[id_right] + uh[id_bottom] + uh[id_top])
			  - lambda_x * (fuh[id_right] - fuh[id_left])
			  - lambda_y * (guh[id_top] - guh[id_bottom]);
		  vhm[id] = 0.25*(vh[id_left] + vh[id_right] + vh[id_bottom] + vh[id_top])
			  - lambda_x * (fvh[id_right] - fvh[id_left])
			  - lambda_y * (gvh[id_top] - gvh[id_bottom]);
		  }
	  /*if (irank == 0)
	  {
	  printf("I am alive to this point 2.b \n");
	  }*/
	  // 3)Check whether communication has completed #######################################
	  //WAIT FOR COMMUNICATION TO COMPLETE
	  if (rank_left >= 0){
		  MPI_Wait(&request_recv_left, &status);
		  //MPI_Wait(&request_uh_recv_left,&status);
		  for (int k = 1; k < (ny_loc + 1); k++)
		  {
			  h[ID_2D(k, 0, nx_loc)] = buf_left_recv[k - 1];
			  uh[ID_2D(k, 0, nx_loc)] = buf_left_recv[k - 1 + ny_loc];
			  vh[ID_2D(k, 0, nx_loc)] = buf_left_recv[k - 1 + (2 * ny_loc)];
		  }
	  }
	  else{ //apply physical BC
		  for (int k = 1; k < (ny_loc + 1); k++)
		  {
			  h[ID_2D(k, 0, nx_loc)] = h[ID_2D(k, 1, nx_loc)];
			  uh[ID_2D(k, 0, nx_loc)] = -uh[ID_2D(k, 1, nx_loc)];
			  vh[ID_2D(k, 0, nx_loc)] = -vh[ID_2D(k, 1, nx_loc)];
		  }
	  }

	  if (rank_right >= 0){
		  MPI_Wait(&request_recv_right, &status);
		  //MPI_Wait(&request_uh_recv_right,&status);
		  for (int k = 1; k < (ny_loc + 1); k++)
		  {
			  h[ID_2D(k, nx_loc + 1, nx_loc)] = buf_right_recv[k - 1];
			  uh[ID_2D(k, nx_loc + 1, nx_loc)] = buf_right_recv[k - 1 + ny_loc];
			  vh[ID_2D(k, nx_loc + 1, nx_loc)] = buf_right_recv[k - 1 + (2 * ny_loc)];
		  }

	  }
	  else{ //apply physical BC
		  for (int k = 1; k < (ny_loc + 1); k++)
		  {
			  h[ID_2D(k, nx_loc + 1, nx_loc)] = h[ID_2D(k, nx_loc, nx_loc)];
			  uh[ID_2D(k, nx_loc + 1, nx_loc)] = -uh[ID_2D(k, nx_loc, nx_loc)];
			  vh[ID_2D(k, nx_loc + 1, nx_loc)] = -vh[ID_2D(k, nx_loc, nx_loc)];
		  }
	  }

	  if (rank_up >= 0){
		  MPI_Wait(&request_recv_up, &status);
		  //MPI_Wait(&request_uh_recv_left,&status);
		  for (int k = 1; k < (nx_loc + 1); k++)
		  {
			  h[ID_2D(ny_loc + 1, k, nx_loc)] = buf_up_recv[k - 1];
			  uh[ID_2D(ny_loc + 1, k, nx_loc)] = buf_up_recv[k - 1 + nx_loc];
			  vh[ID_2D(ny_loc + 1, k, nx_loc)] = buf_up_recv[k - 1 + (2 * nx_loc)];
		  }
	  }
	  else{ //apply physical BC
		  for (int k = 1; k < (nx_loc + 1); k++)
		  {
			  h[ID_2D(ny_loc + 1, k, nx_loc)] = h[ID_2D(ny_loc, k, nx_loc)];
			  uh[ID_2D(ny_loc + 1, k, nx_loc)] = -uh[ID_2D(ny_loc, k, nx_loc)];
			  vh[ID_2D(ny_loc + 1, k, nx_loc)] = -vh[ID_2D(ny_loc, k, nx_loc)];
		  }
	  }

	  if (rank_down >= 0){
		  MPI_Wait(&request_recv_down, &status);
		  //MPI_Wait(&request_uh_recv_right,&status);
		  for (int k = 1; k < (nx_loc + 1); k++)
		  {
			  h[ID_2D(0, k, nx_loc)] = buf_down_recv[k - 1];
			  uh[ID_2D(0, k, nx_loc)] = buf_down_recv[k - 1 + nx_loc];
			  vh[ID_2D(0, k, nx_loc)] = buf_down_recv[k - 1 + (2 * nx_loc)];
		  }
	  }
	  else{ //apply physical BC
		  for (int k = 1; k < (nx_loc + 1); k++)
		  {
			  h[ID_2D(0, k, nx_loc)] = h[ID_2D(1, k, nx_loc)];
			  uh[ID_2D(0, k, nx_loc)] = -uh[ID_2D(1, k, nx_loc)];
			  vh[ID_2D(0, k, nx_loc)] = -vh[ID_2D(1, k, nx_loc)];
		  }
	  }
	  /*if (irank == 0)
	  {
	  printf("I am alive to this point 3 \n");
	  }*/

	  // 4) Compute ghost fluxes (need ghost values) #################################################
	  //left ghost
	  j = 0;
	  for (i = 1; i < ny_loc + 1; i++){
		  id = ID_2D(i, j, nx_loc);
		  fh[id] = uh[id];
		  fuh[id] = uh[id] * uh[id] / h[id] + 0.5*g*h[id] * h[id];
		  fvh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gh[id] = vh[id]; //flux for the height equation: v*h
		  guh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gvh[id] = vh[id] * vh[id] / h[id] + 0.5*g*h[id] * h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
	  }

	  //right ghost
	  j = nx_loc + 1;
	  for (i = 1; i < ny_loc + 1; i++){
		  id = ID_2D(i, j, nx_loc);
		  fh[id] = uh[id];
		  fuh[id] = uh[id] * uh[id] / h[id] + 0.5*g*h[id] * h[id];
		  fvh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gh[id] = vh[id]; //flux for the height equation: v*h
		  guh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gvh[id] = vh[id] * vh[id] / h[id] + 0.5*g*h[id] * h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
	  }

	  //bottom ghost
	  i = 0;
	  for (j = 1; j < nx_loc + 1; j++){
		  id = ID_2D(i, j, nx_loc);
		  fh[id] = uh[id];
		  fuh[id] = uh[id] * uh[id] / h[id] + 0.5*g*h[id] * h[id];
		  fvh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gh[id] = vh[id]; //flux for the height equation: v*h
		  guh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gvh[id] = vh[id] * vh[id] / h[id] + 0.5*g*h[id] * h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
	  }

	  //top ghost
	  i = ny_loc + 1;
	  for (j = 1; j < nx_loc + 1; j++){
		  id = ID_2D(i, j, nx_loc);
		  fh[id] = uh[id];
		  fuh[id] = uh[id] * uh[id] / h[id] + 0.5*g*h[id] * h[id];
		  fvh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gh[id] = vh[id]; //flux for the height equation: v*h
		  guh[id] = uh[id] * vh[id] / h[id]; //flux for the momentum equation: u*v**h 
		  gvh[id] = vh[id] * vh[id] / h[id] + 0.5*g*h[id] * h[id]; //flux for the momentum equation: v^2*h + 0.5*g*h^2
	  }
	  /*if (irank == 0)
	  {
	  printf("I am alive to this point 4 \n");
	  }*/


	  // 5) Compute updated variables at the edges (need ghost values) ###################################################
	  //left boundary
	  j = 1;
	  for (i = 2; i < ny_loc; i++){
		  id = ID_2D(i, j, nx_loc);
		  id_left = ID_2D(i, j - 1, nx_loc);
		  id_right = ID_2D(i, j + 1, nx_loc);
		  id_bottom = ID_2D(i - 1, j, nx_loc);
		  id_top = ID_2D(i + 1, j, nx_loc);
		  hm[id] = 0.25*(h[id_left] + h[id_right] + h[id_bottom] + h[id_top])
			  - lambda_x * (fh[id_right] - fh[id_left])
			  - lambda_y * (gh[id_top] - gh[id_bottom]);
		  uhm[id] = 0.25*(uh[id_left] + uh[id_right] + uh[id_bottom] + uh[id_top])
			  - lambda_x * (fuh[id_right] - fuh[id_left])
			  - lambda_y * (guh[id_top] - guh[id_bottom]);
		  vhm[id] = 0.25*(vh[id_left] + vh[id_right] + vh[id_bottom] + vh[id_top])
			  - lambda_x * (fvh[id_right] - fvh[id_left])
			  - lambda_y * (gvh[id_top] - gvh[id_bottom]);
	  }
	  //right boundary
	  j = nx_loc;
	  for (i = 2; i < ny_loc; i++){
		  id = ID_2D(i, j, nx_loc);
		  id_left = ID_2D(i, j - 1, nx_loc);
		  id_right = ID_2D(i, j + 1, nx_loc);
		  id_bottom = ID_2D(i - 1, j, nx_loc);
		  id_top = ID_2D(i + 1, j, nx_loc);
		  hm[id] = 0.25*(h[id_left] + h[id_right] + h[id_bottom] + h[id_top])
			  - lambda_x * (fh[id_right] - fh[id_left])
			  - lambda_y * (gh[id_top] - gh[id_bottom]);
		  uhm[id] = 0.25*(uh[id_left] + uh[id_right] + uh[id_bottom] + uh[id_top])
			  - lambda_x * (fuh[id_right] - fuh[id_left])
			  - lambda_y * (guh[id_top] - guh[id_bottom]);
		  vhm[id] = 0.25*(vh[id_left] + vh[id_right] + vh[id_bottom] + vh[id_top])
			  - lambda_x * (fvh[id_right] - fvh[id_left])
			  - lambda_y * (gvh[id_top] - gvh[id_bottom]);
	  }

	  //bottom boundary
	  i = 1;
	  for (j = 1; j < nx_loc + 1; j++){
		  id = ID_2D(i, j, nx_loc);
		  id_left = ID_2D(i, j - 1, nx_loc);
		  id_right = ID_2D(i, j + 1, nx_loc);
		  id_bottom = ID_2D(i - 1, j, nx_loc);
		  id_top = ID_2D(i + 1, j, nx_loc);
		  hm[id] = 0.25*(h[id_left] + h[id_right] + h[id_bottom] + h[id_top])
			  - lambda_x * (fh[id_right] - fh[id_left])
			  - lambda_y * (gh[id_top] - gh[id_bottom]);
		  uhm[id] = 0.25*(uh[id_left] + uh[id_right] + uh[id_bottom] + uh[id_top])
			  - lambda_x * (fuh[id_right] - fuh[id_left])
			  - lambda_y * (guh[id_top] - guh[id_bottom]);
		  vhm[id] = 0.25*(vh[id_left] + vh[id_right] + vh[id_bottom] + vh[id_top])
			  - lambda_x * (fvh[id_right] - fvh[id_left])
			  - lambda_y * (gvh[id_top] - gvh[id_bottom]);
	  }

	  //top boundary
	  i = ny_loc;
	  for (j = 1; j < nx_loc + 1; j++){
		  id = ID_2D(i, j, nx_loc);
		  id_left = ID_2D(i, j - 1, nx_loc);
		  id_right = ID_2D(i, j + 1, nx_loc);
		  id_bottom = ID_2D(i - 1, j, nx_loc);
		  id_top = ID_2D(i + 1, j, nx_loc);
		  hm[id] = 0.25*(h[id_left] + h[id_right] + h[id_bottom] + h[id_top])
			  - lambda_x * (fh[id_right] - fh[id_left])
			  - lambda_y * (gh[id_top] - gh[id_bottom]);
		  uhm[id] = 0.25*(uh[id_left] + uh[id_right] + uh[id_bottom] + uh[id_top])
			  - lambda_x * (fuh[id_right] - fuh[id_left])
			  - lambda_y * (guh[id_top] - guh[id_bottom]);
		  vhm[id] = 0.25*(vh[id_left] + vh[id_right] + vh[id_bottom] + vh[id_top])
			  - lambda_x * (fvh[id_right] - fvh[id_left])
			  - lambda_y * (gvh[id_top] - gvh[id_bottom]);
	  }

	  /*if (irank == 0)
	  {
	  printf("I am alive to this point 5.a \n");
	  }
	  if (irank == 0)
	  {
	  printf("I am alive to this point 5.a1 \n");
	  }

	  }
	  if (irank == 0)
	  {
	  printf("I am alive to this point 5.b \n");
	  }
	  // Write data to file
	  write_results("sw2d_final.dat",nx,nx_loc,nx,ny_loc,x,y,h,uh,vh, irank,nproc);
	  */
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (irank == 0)
  {
	  wtime = MPI_Wtime() - wtime;
  }
  //Free memory.
  free ( h );
  free ( uh );
  free ( x );
  free(hm);
  free(fh);
  free(gh);
  free(uhm);
  free(fuh);
  free(guh);
  free(vh);
  free(vhm);
  free(fvh); 
  free(gvh);
  free(y);
  
 //Terminate.
  if(irank==0){
    printf ( "\n" );
    printf ( "SHALLOW_WATER_2D:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
	printf("time = %e\n", wtime);
  }
  if (irank == 0)
  {
	  printf("I am alive to this point 6 \n");
  }
  MPI_Finalize();
  return 0;
}
/******************************************************************************/

void compute_flux_2d_point(double fh, double fuh, double fvh, 
			   double gh, double guh, double gvh, 
			   double h, double uh, double vh, double g)
{
  fh = uh;
  fuh = uh*uh/h + 0.5*g*h*h;
  fvh = uh*vh/h; //flux for the momentum equation: u*v**h 
  gh = vh; //flux for the height equation: v*h
  guh = uh*vh/h; //flux for the momentum equation: u*v**h 
  gvh = vh*vh/h + 0.5*g*h*h; //flux for the momentum equation: v^2*h + 0.5*g*h^2
}

void initial_conditions ( int nx, int ny, double dx, double dy,  double x_length, double x[],double y[], double h[], double uh[] ,double vh[], int irank, int q )

/******************************************************************************/

{
  int i,j, id, id1;

  for ( i = 1; i < nx+1; i++ )
    {
      x[i-1] = -x_length/2+dx/2+(irank%q*nx+i-1)*dx;
      y[i-1] = -x_length/2+dy/2+(floor(irank/q)*ny+i-1)*dy;
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


void write_results ( char *output_filename, int nx, int nx_loc, int ny, int ny_loc, double x[], double y[], double h[], double uh[], double vh[], int irank, int nproc )
/******************************************************************************/

{
  int i,j, id;
  FILE *output;
  double *uh_print;
  double *vh_print;
  double *h_print;
  double *x_print;
  double *y_print;
  double *uh_send;
  double *vh_send;
  double *h_send;
  double *x_send;
  double *y_send;
  int recvcount[nproc];
  int displs[nproc];

  //allocate send arrays
    x_send = malloc(sizeof(double)*nx_loc*ny_loc);
    y_send = malloc(sizeof(double)*nx_loc*ny_loc);
    h_send = malloc(sizeof(double)*nx_loc*ny_loc);
    uh_send = malloc(sizeof(double)*nx_loc*ny_loc);
    vh_send = malloc(sizeof(double)*nx_loc*ny_loc);

  //allocate global array
  if(irank==0){ 
    x_print = malloc(sizeof(double)*nx*ny);
    y_print = malloc(sizeof(double)*nx*ny);
    h_print = malloc(sizeof(double)*nx*ny);
    uh_print = malloc(sizeof(double)*nx*ny);
    vh_print = malloc(sizeof(double)*nx*ny);
  }

  //prepare to send data
  for(i=0; i<ny_loc; i++)
    for(j=0; j<nx_loc; j++){
      id=ID_2D(i+1,j+1,nx_loc);
      x_send[i*nx_loc+j] = x[j];
      y_send[i*nx_loc+j] = y[i];
      h_send[i*nx_loc+j] = h[id];
      uh_send[i*nx_loc+j] = uh[id];
      vh_send[i*nx_loc+j] = vh[id];
  }
  
  //gather data to one process 
  // mind that I am not accounting for extra points this time
  // so could have used MPI_Gather instead of MPI_Gatherv
  for(i=0; i<nproc; i++)
    recvcount[i] = nx*ny/nproc;
  
  displs[0] = 0;
  for(i=1; i<nproc; i++)
    displs[i]=displs[i-1]+recvcount[i];
  
  MPI_Gatherv(uh_send,nx_loc*ny_loc,MPI_DOUBLE,uh_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(vh_send,nx_loc*ny_loc,MPI_DOUBLE,vh_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(h_send,nx_loc*ny_loc,MPI_DOUBLE,h_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(x_send,nx_loc*ny_loc,MPI_DOUBLE,x_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Gatherv(y_send,nx_loc*ny_loc,MPI_DOUBLE,y_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  if(irank==0){  
    //Open the file.
    output = fopen ( output_filename, "wt" );
    
    if ( !output ){
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "WRITE_RESULTS - Fatal error!\n" );
      fprintf ( stderr, "  Could not open the output file.\n" );
      exit ( 1 );
    }
    
    //Write the data.
    for ( i = 0; i < ny; i++ ) 
      for ( j = 0; j < nx; j++ ){
	fprintf ( output, "  %24.16g\t%24.16g\t%24.16g\t %24.16g\t %24.16g\n", x_print[i*nx+j], y_print[i*nx+j],h_print[i*nx+j], uh_print[i*nx+j], vh_print[i*nx+j]);
      }
    
    //Close the file.
    fclose ( output );
    
    //Clean-up
    free(vh_print);
    free(uh_print);
    free(h_print);
    free(y_print);
    free(x_print);
    free(vh_send);
    free(uh_send);
    free(h_send);
    free(x_send);
    free(y_send);
  }
  return;
}
/******************************************************************************/

void getArgs_mpi(int *nx, double *dt, double *x_length, double *t_final, int argc, char *argv[], int irank, MPI_Comm comm)
{

  if(irank==0){  
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
  
  MPI_Bcast(nx,1,MPI_INT,0,comm);
  MPI_Bcast(dt,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(x_length,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(t_final,1,MPI_DOUBLE,0,comm);
}

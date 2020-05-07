#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

// mpirun -np... wave npts dt t-point (400 0.1 5) 
//function headers (definitions at the end of the file)

// getArgs_mpi reads input parameters npts, dt, time_final 
// from command line
void getArgs_mpi(int *npts, double *dt, double *time_final, int argc, char *argv[], int irank, MPI_Comm comm);

// computes the parameter alpha and checks whether user provided stable dt
double get_alpha( double c, double dt, double dx);

//reference solutions
double exact_solution(double x, double t, double c);
double exact_solution_derivative(double x, double t, double c);

//utility functions:
double maxValue(double myArray[], int size); 
double minValue(double myArray[], int size); 
double L2_error_MPI(int npts_loc, double dx, double* data, double* reference);
void print_results(double *u, double *x, int npts, int npts_loc, int irank, int nproc, double time, double c);


/********************************************************************************
 WAVE solves a wave equation u_tt = c^2 u_xx in 1D 
 using 2nd order finite difference formula for both space and time derivative
 *******************************************************************************/

int main(int argc, char * argv[]){

  double *u0; // u^{n-1} solution at previous time-step
  double *u1; // u^n     solution at current time-step
  double *u2; // u^{n+1} solution at next time-step
  double *x;  //         spatial locations of points

  double *u_init; //     initial solution
  double *u_t_init; //   derivative of initial solution

  double dt;  //         time-step 
  double time_final; //  final time for the simulation
  double time; //        current time of the simulation
  double dx;  //         spatial resolution between points

  double c = 1.0; //     wave speed (set to 1.0 by default)
  double alpha; //       CFL non-dimensional parameter parameter alpha = c*dt/dx 

  int npts, npts_loc; // global and local number of points

  int i, irank, nproc; 
  double wtime, wtime_global;


  MPI_Init(&argc,&argv); //initialize MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  //******
  // INITIALIZE SIMULATION
  //******

  //get command-line arguments
  getArgs_mpi(&npts, &dt, &time_final, argc,argv,irank,MPI_COMM_WORLD); 

  //set dx based on the number of points
  dx = 1.0/(npts-1); 

  // compute parameter alpha = c*dt/dx
  alpha = get_alpha(c,dt,dx);
  
  //decide how many points is owned by local rank
  npts_loc = npts/nproc;
  //deal with npts non-divisible by nproc:
  int extra_pts = npts%nproc;

  //add one extra point to first extra_pts processes
  if(irank<extra_pts) npts_loc++ ;

#ifdef DEBUG
  printf("rank=%d, npts =%d, npts_loc=%d, alpha = %f\n",irank,npts,npts_loc,alpha);
#endif

  //set-up coordinates and reference solutions
  x = malloc(sizeof(double)*npts_loc);
  u0 = malloc(sizeof(double)*(npts_loc+2));
  u1 = malloc(sizeof(double)*(npts_loc + 2));
  u2 = malloc(sizeof(double)*(npts_loc + 2));

  u_init = malloc(sizeof(double)*(npts_loc + 2));
  u_t_init = malloc(sizeof(double)*(npts_loc + 2));

  //******  
  // INITIAL CONDITIONS
  //******

  //offset used to compute physical location when dealing with extra_pts
  int offset = 0; 
  if(irank>=extra_pts) offset = extra_pts;

  for(i=0;i<npts_loc;i++){
    // create x coordinates 
    x[i] = (irank*npts_loc+i + offset)*dx;
    //use exact solution as initial condition for u 
    u_init[i] = exact_solution(x[i],0,c); 
    //store initial time deriative of initial condition
    u_t_init[i] = exact_solution_derivative(x[i],0,c); 
    //set initial condition
    u1[i] = u_init[i]; 
    //set initial condition for u0 based on initial dudt
    u0[i] = u1[i] - dt*u_t_init[i]; 
  }  

#ifdef DEBUG
  printf("rank=%d, xmin = %f, xmax = %f\n",irank,minValue(x,npts_loc),maxValue(x,npts_loc));
#endif  


  for(int test=0; test<2; test++){ //timing trick

    //start timer
    MPI_Barrier(MPI_COMM_WORLD);
    wtime = MPI_Wtime();

    
    /**********
     BEGIN SIMULATION - parallelize code below this line
    **********/

    time = 0;
	double ghost;
	MPI_Status status;
    //compute alpha^2
    double alpha2 = alpha*alpha;

    //begin time loop
    while(time<time_final){
      
      //compute interior points
      for(i=1;i<npts_loc-1;i++)
		  u2[i] = alpha2*(u1[i-1] - 2*u1[i] + u1[i+1]) + 2*u1[i] - u0[i];


	  if (irank == 0)
		  //compute right end points
		  //speficy boundary condition for domain left end
		  u2[0] = exact_solution(x[0], time + dt, c);

	  if (irank == (nproc - 1))
		  //compute right end points
		  //specify boundary condition for domain right end
		  u2[npts_loc - 1] = exact_solution(x[npts_loc - 1], time + dt, c);

      //send right boundary, receive ghost
	  //Inter-process boundary
	  //left
	  if (irank == 0)
	  {
		  ghost = u1[npts_loc - 1]; //right ghost point

		  //send ghost to right neighbor and receive ghost from it
		  MPI_Sendrecv_replace(&ghost, 1, MPI_DOUBLE, 1, 101, 1, 102, MPI_COMM_WORLD, &status);

		  //compute u2 at right rank boundary
		  u2[npts_loc - 1] = alpha2*(u1[npts_loc - 2] - 2 * u1[npts_loc - 1] + ghost) + 2 * u1[npts_loc - 1] - u0[npts_loc - 1];
	  }
	  else if (irank == (nproc - 1)){ //Last rank communicates only to the left

		  ghost = u1[0]; //left ghost point
		  //send ghost to left and receive ghost from it
		  MPI_Sendrecv_replace(&ghost, 1, MPI_DOUBLE, irank - 1, 102, irank - 1, 101, MPI_COMM_WORLD, &status);

		  //compute u2 at left rank-boundary
		  u2[0] = alpha2*(ghost - 2 * u1[0] + u1[1]) + 2 * u1[0] - u0[0];
	  } else{

		//left ghost
		ghost = u1[0]; 
		//send ghost to left rank receive ghost from left rank
		MPI_Sendrecv_replace(&ghost,1,MPI_DOUBLE,irank-1,102,irank-1,101,MPI_COMM_WORLD, &status);
		//compute u2 at left rank-boundary
		u2[0] = alpha2*(ghost - 2 * u1[0] + u1[1]) + 2 * u1[0] - u0[0];


		//right ghost 
		ghost = u1[npts_loc - 1];
		//send ghost to right and revceive ghost from right
		MPI_Sendrecv_replace(&ghost,1,MPI_DOUBLE,irank+1,101,irank+1,102,MPI_COMM_WORLD,&status);
		//compute right rank-boundary
		u2[npts_loc - 1] = alpha2*(u1[npts_loc - 2] - 2 * u1[npts_loc - 1] + ghost) + 2 * u1[npts_loc - 1] - u0[npts_loc - 1];
	  }

      
      
      //send left boundary, receive ghost
      
      //increment time
      time = time+dt;
      
      //update u0,u1
      for(i=0;i<npts_loc;i++){
	u0[i] = u1[i];
	u1[i] = u2[i];
      }
      
      
    }//end time loop
    wtime = MPI_Wtime()-wtime;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&wtime, &wtime_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  }

  /**************
    END SIMULATION - do not modify below this line, unless you have to
  **************/

  // Computation of the L2 error, assuming we have made a full revolution 
  
  double u_error = L2_error_MPI(npts_loc, dx, u1, u_init);
  if(irank==0)
    printf("%d\t%d\t%e\t%e\n",nproc,npts,u_error,wtime_global);

  //write output data to file
  print_results(u1,x,npts,npts_loc,irank,nproc,time,c);

  
  //Deallocation
  free(u0);
  free(u1);
  free(u2);
  free(x);


  MPI_Finalize();

  return 0;
} 

void getArgs_mpi(int *npts, double *dt, double *time_final, int argc, char *argv[], int irank, MPI_Comm comm)
{

  if(irank==0){
    if ( argc != 4 ) /* argc should be 2 for correct execution */
      {
	//If not correct number of arguments, assume n=1000
	printf("Incorrect number of arguments. Usage: ./wave npts dt time_final \n ABORTING!\n");
	MPI_Finalize ( );
	exit ( 1 ); // exit with error code 1
      }
    else
      {
	//Use input argument
	*npts = atoi(argv[1]);	
	*dt = atof(argv[2]);
	*time_final = atof(argv[3]);
      }
  }

  MPI_Bcast(npts,1,MPI_INT,0,comm);
  MPI_Bcast(dt,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(time_final,1,MPI_DOUBLE,0,comm);
}

double get_alpha( double c, double dt, double dx)
{
  double alpha = c*dt/dx;

    if ( 1.0 <= fabs ( alpha ) )
    {
      int irank;
      MPI_Comm_rank(MPI_COMM_WORLD,&irank);
      if ( irank == 0 )
	{
	  printf ("\n" );
	  printf ("  CFL condition not met!");
	  printf ("  c = %g\n", c );
	  printf ("  dt = %g\n", dt );
	  printf ("  dx = %g\n", dx );
	  printf ("  alpha = d*dt/dx %g - abs value exceeds 1\n", alpha );
	  printf ("  Computation will not be stable!\n" );
	}
      MPI_Finalize ( );
      exit ( 1 );
    }

    return alpha;

}

double exact_solution(double x, double t, double c)
{
  double pi = acos(-1.0);
  
  double u = sin(2.0 * pi*(x - c*t));

  return u;
}

double exact_solution_derivative(double x, double t, double c)
{
  double pi = acos(-1.0);

  double dudt = -2.0*pi*c*cos (2.0*pi*(x-c*t));

  return dudt;

}

double maxValue(double myArray[], int size) {

  double maxV = myArray[0];

  for (int i = 1; i < size; ++i) {
    if ( myArray[i] > maxV ) {
      maxV = myArray[i];
    }
  }
  return maxV;
}

double minValue(double myArray[], int size) {

  double minV = myArray[0];

  for (int i = 1; i < size; ++i) {
    if ( myArray[i] < minV ) {
      minV = myArray[i];
    }
  }
  return minV;
}

double L2_error_MPI(int npts_loc, double dx, double* data, double* reference){
    double error=0.0; 
    double error_total = 0.0;

    for(int j=0;j<npts_loc;j++){
      error += dx*pow((data[j]-reference[j]),2);
    }
    MPI_Allreduce(&error,&error_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return sqrt(error_total);
}

void print_results(double *u, double *x, int npts, int npts_loc, int irank, int nproc, double time, double c)
  {
    double *u_print;
    double *x_print;

    int recvcount[nproc];
    int displs[nproc];
    int i;

    int extra_pts = npts%nproc;

    for(i=0; i<nproc; i++)
      {
	recvcount[i] = npts/nproc;
	if(i<extra_pts) recvcount[i]++;	
      }

    displs[0] = 0;
    for(i=1; i<nproc; i++)
      displs[i]=displs[i-1]+recvcount[i];

    if(irank==0){ 
      x_print = malloc(sizeof(double)*npts);
      u_print = malloc(sizeof(double)*npts);
    }

    //gather data to one process
    MPI_Gatherv(u,npts_loc,MPI_DOUBLE,u_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gatherv(x,npts_loc,MPI_DOUBLE,x_print,recvcount,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);


    if(irank==0){
      FILE *f;
      f = fopen("results_wave.dat","a"); //open file
      fprintf(f,"x\t u\t u_exact\n");
      for(i=0; i<npts; i++){
	fprintf(f,"%f\t%f\t%f\n",x_print[i],u_print[i],exact_solution(x_print[i],time,c));
      }
      fclose(f);
      free(x_print);
      free(u_print);
    }

  }

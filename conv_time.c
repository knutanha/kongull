#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>

// Parallel toolboxes
#include <mpi.h>
#include <omp.h>

typedef double Real;
// GLOBAL VARS
#define pi 3.141592653589793238462643383279502884197169399375105820974944592307816406286


// ALLOCATION FUNCTIONS
Real *createRealArray (int n);
//Real **createReal2DArray (int m, int n);
Real **createReal2DArray (int n1, int n2);

// TRANSPOSER
void transpose_paral(Real** bt, Real** b, int* sdispl,int* scount, int m, int dist, int* distribution, int size,  MPI_Comm communicator);

// FORTRAN FUNCTIONS
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

// I/O FUNCTIONS (Only to be used on very small systems)
void printMatrix(int* distribution, int m, Real** b, int dist, int clusters, int globalRank, MPI_Comm mpicom, MPI_Status *status);
void printRealArray(Real* input, int inputlength, int rank);
void print2DArray(Real** input, int rows, int cols, int rank);
void free2DArray(Real **input);

Real f(int j, int i, Real h)
{
    Real x,y;
    x = h + h*(double)i;
    y = h + h*(double)j;
    //printf("f(%i,%i) %e\n",j,i,5.0*pi*pi * sin(pi*x) * sin(2.0*pi*y));
    return 5.0*pi*pi * sin(pi*x) * sin(2.0*pi*y);
}

Real evaluate(int j, int i, Real u_ji, Real h)
{
    Real x,y,err;
    x = h + h*(double)i;
    y = h + h*(double)j;
    err = fabs(sin(pi*x) * sin(2.0*pi*y) - u_ji);
    //printf("f(%i,%i) u_ji = %e, trueVal = %e, err = %e\n",j,i,u_ji,sin(pi*x) * sin(2.0*pi*y), err);
    return err;
    
}


int main(int argc, char **argv)
{
    //printf("\ntest\n");
    Real *diag, **b, **bt, **z;
    Real h, umax, globalMax, globalMin;
    int i, j, k, l, m, n, nn, omp_threads;
    double wtime, time;
    
    // Parallel variables
    int size, rank, dist, leftovers, globalRank, distDispl;
    int* distribution,* scount, *sdispl;
    
    if( argc < 2 ) {
        printf("need a problem size\n");
        return 0;
    }
	
    // Problem size
    n = atoi( argv[argc - 1] );
    m  = n-1;
    nn = 4*n;
    h    = 1./(Real)n;
    
    // PI
    //pi   = 4.*atan(1.);
    
    // MPI INIT:
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // INIT TIME
    wtime = MPI_Wtime(); 
    time = MPI_Wtime();
    
    
    
    // openMP init:
    #pragma omp parallel
    {
        omp_threads = omp_get_num_threads();
    }
    
    // -----------------------------------
    
    distribution = (int*) malloc (size * sizeof(int));
    scount 		 = (int*) malloc (size * sizeof(int));
    sdispl 	     = (int*) malloc (size * sizeof(int));
	
    // Distribution of data
    dist =  m / size;
    leftovers = m - dist * size;
    for (i = 0; i < size; i++)
    {
	    distribution[i] = (dist + (i < leftovers ? 1:0) );
    }
    
    for (i = 0, distDispl = 0; i < rank; i++)
    {
        distDispl += distribution[i];
    }
    dist = distribution[rank];
	
	
    diag = createRealArray(m);
    b = createReal2DArray(dist,m);
    bt = createReal2DArray(dist,m);
    z = createReal2DArray(omp_threads,nn);
    
    
    // FILL MATRIX 
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < dist; j++)
        {
            b[j][i] = h*h * f(j+distDispl,i,h);
        }
    }
    
    
	// ---- STATUS PRINT ----- //
    if(rank == 0)
    {
        printf("\n------ RANK 0 PRINTING ------\nSize = %i, omp_threads = %i, n = %i\nDISTRIBUTION = {",size, omp_threads, n);
        for (i = 0; i < size; i++)
        {
            printf("%i",distribution[i]);
            if ( i < size-1) printf(",");
        }
        printf("}\n");
    }    
    
    // ------- POISSON SOLVER ------- //
    
    // PRINTTIME
    if( rank==0 ) printf("Init_time \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
    
    
    // Set eigenvalues
    for (i=0; i < m; i++) {
        diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    }
    
    // PRINTTIME
    if( rank==0 ) printf("Eigen_time \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
    
    #pragma omp parallel for schedule(dynamic)
    for (j=0; j < dist; j++) 
    {
        fst_(b[j], &n, z[omp_get_thread_num()], &nn);
    }
    
    // PRINTTIME
    if( rank==0 ) printf("fst_time_1 \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
    
    transpose_paral(bt,b,sdispl, scount, m, dist, distribution, size, MPI_COMM_WORLD);
    
    // PRINTTIME
    if( rank==0 )  printf("trans_time_1 \t= \t%1.8e\n", MPI_Wtime() - wtime);
    time = MPI_Wtime();
    
    #pragma omp parallel for schedule(dynamic)
    for (j=0; j < dist; j++) 
    {
        fstinv_(bt[j], &n, z[omp_get_thread_num()], &nn);
    }
    
    // PRINTTIME
    if( rank==0 ) printf("fstinv_time_1 \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
    
    #pragma omp parallel for schedule(dynamic) private(i)
    for (j=0; j < dist; j++) 
    {
        for (i=0; i < m; i++) 
        {
            bt[j][i] = bt[j][i]/(diag[j+distDispl]+diag[i]);
        }
    }
    
    // PRINTTIME
    if( rank==0 ) printf("eigdiv_time_1 \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
    
    
    #pragma omp parallel for schedule(dynamic)
    for (j=0; j < dist; j++) 
    {
        fst_(bt[j], &n, z[omp_get_thread_num()], &nn);
    }
    
    // PRINTTIME
    if( rank==0 ) printf("fst_time_2 \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
    
    transpose_paral(b,bt,sdispl, scount, m, dist, distribution, size,  MPI_COMM_WORLD);
    
    // PRINTTIME
    if( rank==0 ) printf("trans_time_2 \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
    
    #pragma omp parallel for schedule(dynamic)
    for (j=0; j < dist; j++) 
    { 
        fstinv_(b[j], &n, z[omp_get_thread_num()], &nn);
    }
    
    // PRINTTIME
    if( rank==0 ) printf("fstinv_time_2 \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
    
    umax = 0.0;
    #pragma omp parallel for schedule(dynamic) private(i) shared(umax)
    for (j=0; j < dist; j++) 
    {
        Real temp;
        for (i=0; i < m; i++) 
        {
            temp = evaluate(j+distDispl,i,b[j][i],h);
            if (temp > umax) umax = temp;
        }
    }
    
    // PRINTTIME
    if( rank==0 ) printf("umax_time_local \t= \t%1.8e\n", MPI_Wtime() - time);
    time = MPI_Wtime();
   
    
    MPI_Reduce(&umax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //PRINT TIME AND MAX GLOBAL ERROR
    if (rank == 0) printf ("globalMaxErr = %e, uerr(0) = %e \n umax_time_total \t= \t%1.8e\n",globalMax,umax, MPI_Wtime() - time);

    
    umax = 0.0;
    #pragma omp parallel for schedule(dynamic) private(i) shared(umax)
    for (j=0; j < dist; j++) 
    {
        for (i=0; i < m; i++) 
        {
            if (b[j][i] > umax) umax = b[j][i];
        }
    }
    MPI_Reduce(&umax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("globalMax = %e, umax(0) = %e \n",globalMax,umax);
   
    
    free2DArray(b);
    free2DArray(bt);
    free(z);
    free(diag);
    
    
    wtime = MPI_Wtime() - wtime; 
    MPI_Reduce(&wtime, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&wtime, &globalMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if( rank==0 ) printf("wtime(Rank 0) = %e\nmax(wtime) is %E\nmin(wtime) = %E\n\n", wtime, globalMax, globalMin);
    if (rank!=0) printf("wtime(Rank %i) = %e\n", rank, wtime);
    //printf("Finalizing Rank %i!\n", rank);
    MPI_Finalize();

}



// ------ PARALLELL TRANSPOSE FUNCTION ------ //
void transpose_paral(Real** bt, Real** b, int* sdispl,int* scount, int m, int dist, int* distribution, int size,  MPI_Comm communicator)
{
    int i,j;
	sdispl[0] = 0;
    scount[0] = distribution[0]*dist;
    for (i=1; i < size; i++)
    {
        scount[i] = distribution[i] * dist;
	    sdispl[i] = scount[i-1] + sdispl[i-1];
    }
    
    #pragma omp parallel for schedule(dynamic) private(i)
	for (j = 0; j < size; j++)
	{
		for (i = 0; i < dist; i++)
		{
			memcpy(bt[0] + sdispl[j] + distribution[j]*i, b[i] + sdispl[j]/dist, distribution[j]*sizeof(double));
		}
	}
	
	MPI_Alltoallv(bt[0], scount, sdispl, MPI_DOUBLE, b[0], scount, sdispl, MPI_DOUBLE, communicator);
	
	#pragma omp parallell for schedule(dynamic) private(i)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < dist; j++)
		{
			bt[j][i] = *(b[0] + i*dist + j);
		}
	}
}

//-------------- UTILITY -------------------//
        
Real *createRealArray (int n){
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
    int i, n;
    Real **a;
    a    = (Real **)malloc(n1   *sizeof(Real *));
    a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
    for (i=1; i < n1; i++) 
    {
        a[i] = a[i-1] + n2;
    }
    n = n1*n2;
    memset(a[0],0 ,n*sizeof(Real));
    return (a);
}

void free2DArray(Real **input) {
    free(input[0]);
    free(input);
}



void printRealArray(Real* input, int inputlength, int rank)
{
    int i;
	printf("Rank %i = {", rank);
	for (i=0; i<inputlength; i++)
	{
		printf("%e",input[i]);
		if (i < inputlength-1)
			printf(", ");
		else
			printf("}\n");
	}
}      


void print2DArray(Real** input, int rows, int cols, int rank)
{
    int i,j;
	printf("Rank %i = [\n", rank);
	for (i = 0; i < rows; i++)
	{
	    for (j = 0; j < cols; j++)
	    {
            printf( "(r%i,%i,%i)%e",rank,j,i,input[j][i] );
            if(j < cols-1) printf(" ");
        }
        printf("]\n");
    }
    printf("]\n\n");
}
        
void printMatrix(int* distribution, int m, Real** b, int dist, int size, int rank, MPI_Comm mpicom, MPI_Status *status)
{
    int i,j,k,l,count;

    Real * buffer = createRealArray(m*dist);
    if (rank == 0)
    {
        Real ** A = createReal2DArray (m,m);
        for (i = 0; i < m; i++)
        {
            for (j = 0; j<dist; j++)
            {
                A[j][i] = b[j][i];
            }
        }
        
        l = dist;
        for (k = 1; k < size; k++)
        {
            MPI_Recv(buffer, count, MPI_DOUBLE, k, 100, mpicom, status);
            for (i = 0; i < m; i++)
            {
                for ( j = 0; j < distribution[k]; j++)
                { 
                    A[l + j][i] = buffer[i + j*m];
                }
            }
            l = l+distribution[k];
        }
        
        
        // PRINT //
        printf("[\n");
        for (i = 0; i < m; i++)
        {   
            printf("[");
            for( j = 0; j < m; j++)
            {
                printf( "%e",A[j][i] );
                if(j < m-1) printf(" ");
            }
            printf("]\n");
        }
        printf("]\n\n");
    }
    else
    {
        count = dist*m;
        MPI_Send(b[0], count, MPI_DOUBLE, 0,100, mpicom);
    }
}    



/* PRINTFUNKSJONER

MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i<size; i++)
    {
        if(rank == i) print2DArray(bt,m,dist,rank);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0)printf("\n ---------------------------------------------- \n\n");
    MPI_Barrier(MPI_COMM_WORLD);
    */



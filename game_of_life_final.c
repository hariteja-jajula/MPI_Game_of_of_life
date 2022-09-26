#include "mpi.h"
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
//Declare Global Variables
int DIES=0;
int ALIVE=1;

//compute Function
void compute(int *life, int *temp, int M, int N) {
  int i, j, value;
  for (i = 1; i < M+1; i++) {
    for (j = 1; j < N+1; j++) {
      /* find out the value of the current cell */
      value = life[(i-1)*(N+2) + (j-1)] + life[(i-1)*(N+2) + j] + 
              life[(i-1)*(N+2) + (j+1)] + life[i*(N+2) + (j-1)] + 
              life[i*(N+2) + (j+1)] + life[(i+1)*(N+2) + (j-1)] + 
              life[(i+1)*(N+2) + j] + life[(i+1)*(N+2) + (j+1)] ;
      
      /* check if the cell dies or life is born */
      if (life[i*(N+2) + j]) { // cell was alive in the earlier iteration
	if (value < 2 || value > 3) {
	  temp[i*(N+2) + j] = DIES ;
	}
	else // value must be 2 or 3, so no need to check explicitly
	  temp[i*(N+2) + j] = ALIVE ; // no change
      } 
      else { // cell was dead in the earlier iteration
	if (value == 3) {
	  temp[i*(N+2) + j] = ALIVE;
	}
	else
	  temp[i*(N+2) + j] = DIES; // no change
      }
    }
  }

}

int main(int argc, char *argv[])  
{
    int numtasks, rank, sendcount, recvcount, source;
    int N,NTIMES, *life=NULL, *TEMP1=NULL,*GATHER=NULL,*TEMP2=NULL;
    int i, j, k;
    double t1, t2;
    N = atoi(argv[1]);
    NTIMES = atoi(argv[2]);
   int n = (N+2)*(N+2);
    recvcount=n/(numtasks);
    /**** START MPI PROCESS ****/
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    life = malloc(n*sizeof(int));
    GATHER = malloc(n*sizeof(int));
 
    int block_size=sqrt(recvcount);
    MPI_Request reqs[N];   // required variable for non-blocking calls
    MPI_Status stats[N];   // required variable for Waitall routine

    
    /** INITITALIZE THE ARRAY **/
  if(rank==0){
    for (i = 0; i < N+2; i++) {
    life[i*(N+2)] = life[i*(N+2) + (N+1)] = DIES ;
    GATHER[i*(N+2)] = GATHER[i*(N+2) + (N+1)] = DIES ;
    }
    for (j = 0; j < N+2; j++) {
    life[j] = life[(N+1)*(N+2) + j] = DIES ;
    GATHER[j] = GATHER[(N+1)*(N+2) + j] = DIES ;
    }
    
    for (i = 1; i < N+1; i++) 
    {
    for (j = 1; j< N+1; j++) {
        if (drand48() < 0.5) 
	        life[i*(N+2) + j] = ALIVE ;
        else
	        life[i*(N+2) + j] = DIES ;
    }
    }
  }
    /****************************/


    /****************************/
    TEMP1 = malloc(recvcount*sizeof(int));
    TEMP2 = malloc(recvcount*sizeof(int));
    source=0;
    MPI_Scatter(life,n,MPI_INT,TEMP1,recvcount,MPI_INT,source,MPI_COMM_WORLD);

    /*PROCESSING WITH RECVD BLOCKS*/
    int top_block =1;
    int bottom_block=2;
    int prev = rank-1;
    int next = rank+1;
    if (rank==0)
    {
      prev=MPI_PROC_NULL;
    }
    if (rank==numtasks-1)
    {
      next=MPI_PROC_NULL;
    }
    MPI_Irecv(&TEMP1, N+2, MPI_INT, prev, top_block, MPI_COMM_WORLD, &reqs[0]);
    MPI_Irecv(&TEMP1, N+2, MPI_INT, next, bottom_block, MPI_COMM_WORLD, &reqs[1]);
    MPI_Isend(&TEMP1, N+2, MPI_INT, prev, top_block, MPI_COMM_WORLD, &reqs[2]);
    MPI_Isend(&TEMP1, N+2, MPI_INT, next, bottom_block, MPI_COMM_WORLD, &reqs[3]);
    MPI_Waitall(4, reqs, stats); 
    compute(TEMP1,TEMP2,(N+2)/numtasks ,N+2);
    compute(TEMP2,TEMP1,(N+2)/numtasks ,N+2);
       
MPI_Finalize();

}

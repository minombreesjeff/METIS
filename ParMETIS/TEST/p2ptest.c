/*
 * main.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define cleartimer(tmr) (tmr = 0.0)
#define starttimer(tmr) (tmr -= MPI_Wtime())
#define stoptimer(tmr) (tmr += MPI_Wtime())
#define gettimer(tmr) (tmr)

#define MAXLEN  10000

/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int npes, mype;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);


  Point2PointTest(mype, npes);

  MPI_Finalize();
}



/*************************************************************************
* Let the game begin
**************************************************************************/
Point2PointTest(int mype, int npes)
{
  int i, j, k, min, max, len;
  int *send, *recv;
  double tmr;

  send = (int *)malloc(sizeof(int)*MAXLEN*npes);
  recv = (int *)malloc(sizeof(int)*MAXLEN*npes);

  for (len = 0; len < MAXLEN; len += 1000) {
    tmr = 0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    starttimer(tmr);
    for (i=0; i<10; i++)
      TotalExchange(mype, npes, len, send, recv);
    MPI_Barrier(MPI_COMM_WORLD);
    stoptimer(tmr);
    if (mype == 0)
      printf("%5d, %7.5lf\n", len, tmr/10.0);
  }

}


/*************************************************************************
* This function performs the gather/scatter for the boundary vertices
**************************************************************************/
TotalExchange(int mype, int npes, int len, int *send, int *recv)
{
  int i, j, k; 
  MPI_Request sreq[100], rreq[100];
  MPI_Status statuses[100];

  /* Issue the receives first */
  for (i=0; i<npes; i++) 
    MPI_Irecv((void *)(recv+i*len), len, MPI_INT, i, 1, MPI_COMM_WORLD, rreq+i);

  MPI_Barrier(MPI_COMM_WORLD);

  /* Issue the sends next */
  for (i=0; i<npes; i++) 
    MPI_Isend((void *)(send+i*len), len, MPI_INT, i, 1, MPI_COMM_WORLD, sreq+i);

  /* OK, now get into the loop waiting for the operations to finish */
  MPI_Waitall(npes, rreq, statuses); 
  MPI_Waitall(npes, sreq, statuses); 

}


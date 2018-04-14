/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * io.c
 *
 * This file contains routines related to I/O
 *
 * Started 10/19/94
 * George
 *
 * $Id: io.c,v 1.2 1997/07/18 00:32:07 karypis Exp $
 *
 */

#include <par_kmetis.h>


/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
void ReadGraph(GraphType *graph, char *filename, MPI_Comm comm)
{
  int i, j, k, l, npes, mype;
  int nvtxs, nedges, penum, snedges, snvtxs;
  FILE *fpin;
  idxtype *vtxdist, *sxadj, *sadjncy, *ssize, scratch;
  MPI_Status status;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist = idxsmalloc(npes+1, 0, "ReadGraph: vtxdist");

  if (mype == 0) {
    ssize = idxsmalloc(npes, 0, "ReadGraph: ssize");

    if ((fpin = fopen(filename, "rb")) == NULL) 
      errexit("Failed to open file %s", filename);

    /* Read the nvtxs and nedges */
    fread(&scratch, sizeof(idxtype), 1, fpin);
    nvtxs = scratch;
    fread(&scratch, sizeof(idxtype), 1, fpin);
    nedges = scratch;

    printf("Nvtxs: %d, Nedges: %d\n", nvtxs, nedges);

    /* Construct vtxdist and send it to all the processors */
    vtxdist[0] = 0;
    for (i=0,k=nvtxs; i<npes; i++) {
      l = k/(npes-i);
      vtxdist[i+1] = vtxdist[i]+l;
      k -= l;
    }
  }

  MPI_Bcast((void *)vtxdist, npes+1, IDX_DATATYPE, 0, comm);

  graph->gnvtxs = vtxdist[npes];
  graph->nvtxs = vtxdist[mype+1]-vtxdist[mype];
  graph->xadj = idxmalloc(graph->nvtxs+1, "ReadGraph: xadj");

  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
      snvtxs = vtxdist[penum+1]-vtxdist[penum];
      sxadj = idxmalloc(snvtxs+1, "ReadGraph: sxadj");
      if (penum == 0) {
        fread(sxadj, sizeof(idxtype), snvtxs+1, fpin);
        scratch = sxadj[snvtxs];
      }
      else {
        fread(sxadj+1, sizeof(idxtype), snvtxs, fpin);
        sxadj[0] = scratch;
        scratch = sxadj[snvtxs];
      }

      for (i=snvtxs; i>=0; i--)
        sxadj[i] -= sxadj[0];

      ssize[penum] = sxadj[snvtxs];

      if (penum == mype) 
        idxcopy(snvtxs+1, sxadj, graph->xadj);
      else
        MPI_Send((void *)sxadj, snvtxs+1, IDX_DATATYPE, penum, 1, comm); 

      free(sxadj);
    }
  }
  else 
    MPI_Recv((void *)graph->xadj, graph->nvtxs+1, IDX_DATATYPE, 0, 1, comm, &status);


  graph->nedges = graph->xadj[graph->nvtxs];
  graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: graph->adjncy");


  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
      sadjncy = idxmalloc(ssize[penum], "ReadGraph: sadjncy");
      fread(sadjncy, sizeof(idxtype), ssize[penum], fpin);

      if (penum == mype) 
        idxcopy(ssize[penum], sadjncy, graph->adjncy);
      else
        MPI_Send((void *)sadjncy, ssize[penum], IDX_DATATYPE, penum, 1, comm); 

      free(sadjncy);
    }

    free(ssize);
    fclose(fpin);
  }
  else 
    MPI_Recv((void *)graph->adjncy, graph->nedges, IDX_DATATYPE, 0, 1, comm, &status);
}


/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
void ReadPartitionedGraph(GraphType *graph, char *filename, MPI_Comm comm)
{
  int i, j, k, l, npes, mype;
  int nvtxs, nedges, penum, snedges, snvtxs;
  FILE *fpin;
  idxtype *vtxdist, *sxadj, *sadjncy, *ssize, scratch;
  MPI_Status status;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist = idxsmalloc(npes+1, 0, "ReadGraph: vtxdist");

  if (mype == 0) {
    ssize = idxsmalloc(npes, 0, "ReadGraph: ssize");

    if ((fpin = fopen(filename, "rb")) == NULL) 
      errexit("Failed to open file %s", filename);

    /* Read the nvtxs and nedges */
    fread(&scratch, sizeof(idxtype), 1, fpin);
    nvtxs = scratch;
    fread(&scratch, sizeof(idxtype), 1, fpin);
    nedges = scratch;
    fread(&scratch, sizeof(idxtype), 1, fpin);
    if (scratch != npes)
      printf("The graph does not match with the number of processors!\n");

    printf("Nvtxs: %d, Nedges: %d\n", nvtxs, nedges);

    /* Read vtxdist from the file as well */
    fread(vtxdist, sizeof(idxtype), npes+1, fpin);
  }

  MPI_Bcast((void *)vtxdist, npes+1, IDX_DATATYPE, 0, comm);

  graph->gnvtxs = vtxdist[npes];
  graph->nvtxs = vtxdist[mype+1]-vtxdist[mype];
  graph->xadj = idxmalloc(graph->nvtxs+1, "ReadGraph: xadj");

  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
      snvtxs = vtxdist[penum+1]-vtxdist[penum];
      sxadj = idxmalloc(snvtxs+1, "ReadGraph: sxadj");
      if (penum == 0) {
        fread(sxadj, sizeof(idxtype), snvtxs+1, fpin);
        scratch = sxadj[snvtxs];
      }
      else {
        fread(sxadj+1, sizeof(idxtype), snvtxs, fpin);
        sxadj[0] = scratch;
        scratch = sxadj[snvtxs];
      }

      for (i=snvtxs; i>=0; i--)
        sxadj[i] -= sxadj[0];

      ssize[penum] = sxadj[snvtxs];

      if (penum == mype) 
        idxcopy(snvtxs+1, sxadj, graph->xadj);
      else
        MPI_Send((void *)sxadj, snvtxs+1, IDX_DATATYPE, penum, 1, comm); 

      free(sxadj);
    }
  }
  else 
    MPI_Recv((void *)graph->xadj, graph->nvtxs+1, IDX_DATATYPE, 0, 1, comm, &status);


  graph->nedges = graph->xadj[graph->nvtxs];
  graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: graph->adjncy");


  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
      sadjncy = idxmalloc(ssize[penum], "ReadGraph: sadjncy");
      fread(sadjncy, sizeof(idxtype), ssize[penum], fpin);

      if (penum == mype) 
        idxcopy(ssize[penum], sadjncy, graph->adjncy);
      else
        MPI_Send((void *)sadjncy, ssize[penum], IDX_DATATYPE, penum, 1, comm); 

      free(sadjncy);
    }

    free(ssize);
    fclose(fpin);
  }
  else 
    MPI_Recv((void *)graph->adjncy, graph->nedges, IDX_DATATYPE, 0, 1, comm, &status);
}


/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
float *ReadCoordinates(GraphType *graph, char *filename, MPI_Comm comm)
{
  int i, j, k, l, npes, mype;
  int nvtxs, nedges, penum, snedges, snvtxs;
  float *xyz, *txyz;
  FILE *fpin;
  idxtype *vtxdist, *sxadj, *sadjncy, *ssize, scratch;
  MPI_Status status;
  char xyzfile[256];

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist;

  xyz = fmalloc(graph->nvtxs*3, "io");

  if (mype == 0) {
    sprintf(xyzfile, "%s.xyz", filename);
    if ((fpin = fopen(xyzfile, "rb")) == NULL) 
      errexit("Failed to open file %s", filename);
  }

  if (mype == 0) {
    txyz = fmalloc(2*graph->nvtxs*3, "io");

    for (penum=0; penum<npes; penum++) {
      if (penum == mype) 
        fread(xyz, sizeof(float), 3*(vtxdist[penum+1]-vtxdist[penum]), fpin);
      else {
        fread(txyz, sizeof(float), 3*(vtxdist[penum+1]-vtxdist[penum]), fpin);
        MPI_Send((void *)txyz, 3*(vtxdist[penum+1]-vtxdist[penum]), MPI_FLOAT, penum, 1, comm); 
      }
    }
    free(txyz);
    fclose(fpin);
  }
  else 
    MPI_Recv((void *)xyz, 3*graph->nvtxs, MPI_FLOAT, 0, 1, comm, &status);

  return xyz;
}



/*************************************************************************
* This function writes out a partition vector
**************************************************************************/
void WritePVector(char *gname, idxtype *vtxdist, idxtype *part, MPI_Comm comm)
{
  int i, j, k, l, rnvtxs, npes, mype, penum;
  FILE *fpin;
  idxtype *rpart;
  char partfile[256];
  MPI_Status status;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (mype == 0) {
    sprintf(partfile, "%s.part", gname);
    if ((fpin = fopen(partfile, "wb")) == NULL) 
      errexit("Failed to open file %s", partfile);

    fwrite(part, sizeof(idxtype), vtxdist[1], fpin);

    for (penum=1; penum<npes; penum++) {
      rnvtxs = vtxdist[penum+1]-vtxdist[penum];
      rpart = idxmalloc(rnvtxs, "rpart");
      MPI_Recv((void *)rpart, rnvtxs, IDX_DATATYPE, penum, 1, comm, &status);
      fwrite(rpart, sizeof(idxtype), rnvtxs, fpin);
      for (i=0; i<rnvtxs; i++)
        if (rpart[i] >= npes)
          printf("part[%d] = %d\n", vtxdist[penum]+i, rpart[i]);
      free(rpart);
    }
    fclose(fpin);
  }
  else
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_DATATYPE, 0, 1, comm); 

}


/*************************************************************************
* This function writes out a partition vector
**************************************************************************/
void WriteOVector(char *gname, idxtype *vtxdist, idxtype *part, MPI_Comm comm)
{
  int i, j, k, l, rnvtxs, npes, mype, penum;
  FILE *fpin;
  idxtype *rpart;
  char partfile[256];
  MPI_Status status;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (mype == 0) {
    sprintf(partfile, "%s.order.%d", gname, npes);
    if ((fpin = fopen(partfile, "wb")) == NULL) 
      errexit("Failed to open file %s", partfile);

    fwrite(part, sizeof(idxtype), vtxdist[1], fpin);

    for (penum=1; penum<npes; penum++) {
      rnvtxs = vtxdist[penum+1]-vtxdist[penum];
      rpart = idxmalloc(rnvtxs, "rpart");
      MPI_Recv((void *)rpart, rnvtxs, IDX_DATATYPE, penum, 1, comm, &status);
      fwrite(rpart, sizeof(idxtype), rnvtxs, fpin);
      free(rpart);
    }
    fclose(fpin);
  }
  else
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_DATATYPE, 0, 1, comm); 

}

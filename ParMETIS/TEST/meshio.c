/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * meshio.c
 *
 * This file contains routines related to I/O
 *
 * Started 10/19/94
 * George
 *
 * $Id: meshio.c,v 1.2 1998/09/22 23:17:19 karypis Exp $
 *
 */

#include <parmetis.h>
#define	MAXLINE	8192

/*************************************************************************
* This function reads a mesh from a file
**************************************************************************/
void ParallelReadMesh(MeshType *mesh, char *filename, MPI_Comm comm)
{
  int i, j, k, pe;
  int npes, mype, ier;
  int gnelms, nelms, your_nelms, etype, maxnelms;
  int maxnode, gmaxnode, minnode, gminnode;
  idxtype *elmdist, *elements;
  idxtype *your_elements;
  MPI_Status stat;
  char *line = NULL, *oldstr, *newstr;
  FILE *fpin = NULL;
  int esize, esizes[5] = {-1, 3, 4, 8, 4};
  int mgcnum, mgcnums[5] = {-1, 2, 3, 4, 2};

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  elmdist = mesh->elmdist = idxsmalloc(npes+1, 0, "ReadGraph: elmdist");

  if (mype == npes-1) {
    ier = 0;
    fpin = fopen(filename, "r");

    if (fpin == NULL){
      printf("COULD NOT OPEN FILE '%s' FOR SOME REASON!\n", filename);
      ier++;
    }

    MPI_Bcast(&ier, 1, MPI_INT, npes-1, comm);
    if (ier > 0){
      fclose(fpin);
      MPI_Finalize();
      exit(0);
    }

    line = (char *)GKmalloc(sizeof(char)*(MAXLINE+1), "line");

    fgets(line, MAXLINE, fpin);
    sscanf(line, "%d %d", &gnelms, &etype);

    /* Construct elmdist and send it to all the processors */
    elmdist[0] = 0;
    for (i=0,j=gnelms; i<npes; i++) {
      k = j/(npes-i);
      elmdist[i+1] = elmdist[i]+k;
      j -= k;
    }

    MPI_Bcast((void *)elmdist, npes+1, IDX_DATATYPE, npes-1, comm);
  }
  else {
    MPI_Bcast(&ier, 1, MPI_INT, npes-1, comm);
    if (ier > 0){
      MPI_Finalize();
      exit(0);
    }

    MPI_Bcast((void *)elmdist, npes+1, IDX_DATATYPE, npes-1, comm);
  }

  MPI_Bcast((void *)(&etype), 1, MPI_INT, npes-1, comm);

  gnelms = mesh->gnelms = elmdist[npes];
  nelms = mesh->nelms = elmdist[mype+1]-elmdist[mype];
  mesh->etype = etype;
  esize = esizes[etype];
  mgcnum = mgcnums[etype];

  elements = mesh->elements = idxmalloc(nelms*esize, "ParallelReadMesh: elements");

  if (mype == npes-1) {
    maxnelms = 0;
    for (i=0; i<npes; i++) {
      maxnelms = (maxnelms > elmdist[i+1]-elmdist[i]) ?
      maxnelms : elmdist[i+1]-elmdist[i];
    }

    your_elements = idxmalloc(maxnelms*esize, "your_elements");

    for (pe=0; pe<npes; pe++) {
      your_nelms = elmdist[pe+1]-elmdist[pe];
      for (i=0; i<your_nelms; i++) {

        fgets(line, MAXLINE, fpin);
        oldstr = line;
        newstr = NULL;

        /*************************************/
        /* could get element weigts here too */
        /*************************************/

        for (j=0; j<esize; j++) {
          your_elements[i*esize+j] = (int)strtol(oldstr, &newstr, 10);
          oldstr = newstr;
        }
      }

      if (pe < npes-1) {
        MPI_Send((void *)your_elements, your_nelms*esize, IDX_DATATYPE, pe, 0, comm);
      }
      else {
        for (i=0; i<your_nelms*esize; i++)
          elements[i] = your_elements[i];
      }
    }
    fclose(fpin);
    free(your_elements);
  }
  else {
    MPI_Recv((void *)elements, nelms*esize, IDX_DATATYPE, npes-1, 0, comm, &stat);
  }

  /*********************************/
  /* now check for number of nodes */
  /*********************************/
  minnode = elements[idxamin(nelms*esize, elements)];
  MPI_Allreduce((void *)&minnode, (void *)&gminnode, 1, MPI_INT, MPI_MIN, comm);
  for (i=0; i<nelms*esize; i++)
    elements[i] -= gminnode;

  maxnode = elements[idxamax(nelms*esize, elements)];
  MPI_Allreduce((void *)&maxnode, (void *)&gmaxnode, 1, MPI_INT, MPI_MAX, comm);
  mesh->gnns = gmaxnode+1;

  if (mype==0) printf("Nelements: %d, Nnodes: %d, EType: %d\n", gnelms, mesh->gnns, etype);
}



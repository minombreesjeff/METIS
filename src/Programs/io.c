/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * io.c
 *
 * This file contains routines related to I/O
 *
 * Started 8/28/94
 * George
 *
 * $Id: io.c,v 1.2 1998/01/28 17:23:11 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void ReadGraph(GraphType *graph, char *filename)
{
  int i, j, k, fmt, readew, readvw, edge, ewgt;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  char line[MAXLINE+1], *oldstr, *newstr;
  FILE *fpin;

  InitGraph(graph);

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  if (feof(fpin)) {
    graph->nvtxs = 0;
    return;
  }

  if (sscanf(line,"%d %d %d",&(graph->nvtxs), &(graph->nedges), &fmt) == 2)
    fmt = 0;

  readew = (fmt%10 > 0);
  readvw = ((fmt/10)%10 > 0);
  if (fmt >= 100)
    errexit("Cannot read this type of file format!");

  graph->nedges *=2;

  if (graph->nvtxs > MAXIDX) 
    errexit("\nThe matrix is too big: %d", graph->nvtxs);

  xadj = graph->xadj = idxsmalloc(graph->nvtxs+1, 0, "ReadGraph: xadj");
  adjncy = graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: adjncy");

  vwgt = graph->vwgt = (readvw ? idxmalloc(graph->nvtxs, "ReadGraph: vwgt") : idxsmalloc(graph->nvtxs, 1, "ReadGraph: vwgt"));
  adjwgt = graph->adjwgt = (readew ? idxmalloc(graph->nedges, "ReadGraph: adjwgt") : idxsmalloc(graph->nedges, 1, "ReadGraph: adjwgt"));

  /* Start reading the graph file */
  for (xadj[0]=0, k=0, i=0; i<graph->nvtxs; i++) {
    do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%' && !feof(fpin));
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

    if (readvw) {
      vwgt[i] = (int)strtol(oldstr, &newstr, 10);
      oldstr = newstr;
    }

    for (;;) {
      edge = (int)strtol(oldstr, &newstr, 10) -1;
      oldstr = newstr;

      if (readew) {
        ewgt = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }

      if (edge < 0)
        break;

      adjncy[k] = edge;
      if (readew) 
        adjwgt[k] = ewgt;
      k++;
    } 
    xadj[i+1] = k;
  }

  fclose(fpin);

  if (k != graph->nedges)
    errexit("ReadGraph: Something wrong with the edges from input file %d %d",graph->nedges, k);

}



/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePartition(char *fname, idxtype *part, int n, int nparts)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.part.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n",part[i]);

  fclose(fpout);

}


/*************************************************************************
* This function writes out the partition vectors for a mesh
**************************************************************************/
void WriteMeshPartition(char *fname, int nparts, int ne, idxtype *epart, int nn, idxtype *npart)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.epart.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<ne; i++)
    fprintf(fpout,"%d\n", epart[i]);

  fclose(fpout);

  sprintf(filename,"%s.npart.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<nn; i++)
    fprintf(fpout,"%d\n", npart[i]);

  fclose(fpout);


}



/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePermutation(char *fname, idxtype *iperm, int n)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.iperm",fname);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the permutation file: %s", filename);

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n", iperm[i]);

  fclose(fpout);

}



/*************************************************************************
* This function checks if a graph is valid
**************************************************************************/
int CheckGraph(GraphType *graph)
{
  int i, j, k, l, nvtxs, err=0;
  idxtype *xadj, *adjncy, *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;


  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];

      if (i == k) {
        printf("Vertex %d contains a self-loop (i.e., diagonal entry in the matrix)!\n", i);
        err++;
      }
      else {
        for (l=xadj[k]; l<xadj[k+1]; l++) {
          if (adjncy[l] == i) {
            if (adjwgt[l] != adjwgt[j]) {
              printf("Edges (%d %d) and (%d %d) do not have the same weight! %d %d\n", i,k,k,i, adjwgt[l], adjwgt[adjncy[j]]);
              err++;
            }
            break;
          }
        }
        if (l == xadj[k+1]) {
          printf("Missing edge: (%d %d)!\n", k, i);
          err++;
        }
      }
    }
  }

  if (err > 0) 
    printf("A total of %d errors exist in the input file. Correct them, and run again!\n", err);

  return (err == 0 ? 1 : 0);
}


/*************************************************************************
* This function reads the element node array of a mesh
**************************************************************************/
idxtype *ReadMesh(char *filename, int *ne, int *nn, int *etype)
{
  int i, j, k, esize;
  idxtype *elmnts;
  FILE *fpin;

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fscanf(fpin, "%d %d", ne, etype);

  switch (*etype) {
    case 1:
      esize = 3;
      break;
    case 2:
      esize = 4;
      break;
    case 3:
      esize = 8;
      break;
    case 4:
      esize = 4;
      break;
    default:
      errexit("Unknown mesh-element type: %d\n", *etype);
  }

  elmnts = idxmalloc(esize*(*ne), "ReadMesh: elmnts");

  for (j=esize*(*ne), i=0; i<j; i++) {
    fscanf(fpin, "%d", elmnts+i);
    elmnts[i]--;
  }

  fclose(fpin);

  *nn = elmnts[idxamax(j, elmnts)]+1;

  return elmnts;
}


/*************************************************************************
* This function reads the element node array of a mesh
**************************************************************************/
void WriteGraph(char *filename, int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, j;
  FILE *fpout;

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d", nvtxs, xadj[nvtxs]/2);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n");
    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fpout, " %d", adjncy[j]+1);
  }

  fclose(fpout);
}

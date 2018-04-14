/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * main.c
 *
 * This file contains the driving routine for multilevel method
 *
 * Started 8/28/94
 * George
 *
 * $Id: main.c,v 1.1 1995/08/18 14:02:07 karypis Exp $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>


#define MAXLINE	16384

/*************************************************************************
* The following data structure is a node in a doubly-link node of graphs
* used to represent the coarsion sequence.
* Each node in this list will contain the current graph, and mapping
* information to obtain the coarser graph.
**************************************************************************/
struct graphdef {
  int nvtxs;			/* The number of vertices */
  int nedges;			/* The total number of edges */
  int weightflag;
  int *xadj;			/* CRS storage scheme */
  int *adjncy;
  int *vwgts;
  int *ewgts;
};

typedef struct graphdef GraphType;


int *imalloc(int, char *);



/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int i;
  GraphType graph;
  char filename[256];
  int *partition, options[5], edgecut, nparts, numbering;

  if (argc != 3) {
    printf("Usage: %s <GraphFile> <nparts>\n",argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);
  nparts = atoi(argv[2]);

  ReadGraph(&graph, filename);

  printf("PMETIS: Nvtxs: %6d, Nedges: %8d\n", graph.nvtxs, graph.nedges);

  partition = imalloc(graph.nvtxs, "main: partition");
  options[0] = 0;
  numbering = 1;

  PMETIS(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgts, graph.ewgts,
                  &graph.weightflag, &nparts, options, &numbering, &edgecut, partition);

  printf("EdgeCut: %d\n", edgecut);

}



/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
ReadGraph(GraphType *graph, char *filename)
{
  char line[MAXLINE+1], *oldstr, *newstr;
  int i, j, k, fmt, edge;
  int readew, readvw;
  FILE *fpin;
  int *xadj, *adjncy, *vwgts, *ewgts;


  if ((fpin = fopen(filename, "r")) == NULL)  {
    printf("Failed to open file %s", filename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%');

  if (sscanf(line,"%d %d %d",&(graph->nvtxs), &(graph->nedges), &fmt) == 2)
    fmt = 0;

  readew = (fmt%10 > 0);
  readvw = ((fmt/10)%10 > 0);
  if (fmt >= 100) {
    printf("Cannot read this type of file format!");
    exit(0);
  }

  graph->weightflag = 0;
  if (readvw)
    graph->weightflag += 2;
  if (readew)
    graph->weightflag += 1;

  graph->nedges *= 2;

  xadj = graph->xadj = imalloc(graph->nvtxs+1, "readgraph: xadj");
  adjncy = graph->adjncy = imalloc(graph->nedges, "readgraph: adjncy");
  if (readvw)
    vwgts = graph->vwgts = imalloc(graph->nvtxs, "readgraph: vwgts");
  else
    vwgts = NULL;
  if (readew)
    ewgts = graph->ewgts = imalloc(graph->nedges, "readgraph: ewgts");
  else
    ewgts = NULL;


  /* Start reading the graph file */
  xadj[0] = 1;
  k = 0;
  for (i=1; i<=graph->nvtxs; i++) {
    do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%');
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE)  {
      printf("\nBuffer for fgets not big enough!\n");
      exit(0);
    }

    if (readvw) {
      vwgts[i-1] = (int)strtol(oldstr, &newstr, 10);
      oldstr = newstr;
    }

    for (;;) {
      edge = (int)strtol(oldstr, &newstr, 10);
      oldstr = newstr;

      if (edge <= 0)
        break;

      adjncy[k] = edge;

      if (readew) {
        ewgts[k] = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }

      k++;

    } 
    xadj[i] = k+1;
  }

  if (k != graph->nedges) {
    printf("readgraph: Something wrong with the edges from input file %d %d",graph->nedges, k);
    exit(0);
  }

  fclose(fpin);

}








/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
int *imalloc(int n, char *msg)
{
  int *ptr;

  if (n == 0)
    return NULL;

  ptr = (int *)malloc(sizeof(int)*n);

  if (ptr == NULL) {
    printf("Memory allocation failed: %s\n", msg);
    exit(0);
  }

  return ptr;

}

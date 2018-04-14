/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * pwtest.c
 *
 * This file contains code for testing METIS_WPartGraphRecursive 
 *
 * Started 9/15/97
 * George
 *
 * $Id: pwtest.c,v 1.1 1997/11/04 23:19:54 karypis Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int i, nparts, options[10], pwgts[256];
  idxtype *part;
  GraphType graph;
  char filename[256];
  float sum, tpwgts[256];
  int edgecut, numflag=0, wgtflag=3;


  if (argc != 3) {
    printf("Usage: %s <GraphFile> <nparts>\n",argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);
  nparts = atoi(argv[2]);

  ReadGraph(&graph, filename);

  printf("Enter the tpwgts fractions for %d partitions: ", nparts);
  for (sum=0.0, i=0; i<nparts-1; i++) {
    scanf("%f", tpwgts+i);
    sum += tpwgts[i];
  }
  tpwgts[nparts-1] = 1.0 - sum;

  printf("You entered the following weights: ");
  for (i=0; i<nparts; i++)
    printf("%5.4f ", tpwgts[i]);
  printf("\n");

  part = idxmalloc(graph.nvtxs, "main: part");

  options[0] = 0;
  METIS_WPartGraphRecursive(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, &wgtflag, &numflag, &nparts, tpwgts, options, &edgecut, part);

  printf("%d-way cut: %d\n", nparts, edgecut);
  for (i=0; i<nparts; i++)
    pwgts[i] = 0;
  for (i=0; i<graph.nvtxs; i++)
    pwgts[part[i]]++;
  for (i=0; i<nparts; i++)
    printf("%5d ", pwgts[i]);
  printf("\n");


  GKfree(&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, &part, -1);
}  



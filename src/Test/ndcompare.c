/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ndcompare.c
 *
 * This file contains for comparing the various orderings routines.
 *
 * Started 8/12/97
 * George
 *
 * $Id: ndcompare.c,v 1.1 1997/11/04 23:19:53 karypis Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int i, ii, options[10];
  idxtype *perm, *iperm, fillin[10];
  GraphType graph;
  char filename[256];
  int numflag=0;


  if (argc != 2) {
    printf("Usage: %s <GraphFile>\n",argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);

  ReadGraph(&graph, filename);

  perm = idxmalloc(graph.nvtxs, "main: order");
  iperm = idxmalloc(graph.nvtxs, "main: order");

  options[0] = 1;
  options[OPTION_RTYPE] = 2;
  options[OPTION_ITYPE] = 1;
  options[OPTION_DBGLVL] = 0;
  options[OPTION_OFLAGS] = 1;
  options[OPTION_PFACTOR] = 0;
  options[OPTION_NSEPS] = 1;

  options[OPTION_CTYPE] = 1;
  METIS_NodeND(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, options, perm, iperm);
  fillin[0] = ComputeFillIn2(&graph, iperm);

  options[OPTION_CTYPE] = 2;
  METIS_NodeND(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, options, perm, iperm);
  fillin[1] = ComputeFillIn2(&graph, iperm);

  options[OPTION_CTYPE] = 3;
  METIS_NodeND(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, options, perm, iperm);
  fillin[2] = ComputeFillIn2(&graph, iperm);

  options[0] = 0;
  METIS_EdgeND(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, options, perm, iperm);
  fillin[3] = ComputeFillIn2(&graph, iperm);


  printf("%15s: ", argv[1]);
  for (i=0; i<4; i++)
    printf("%10d ", fillin[i]);

  printf("\tMin: %d\n", idxamin(4, fillin));
  fflush(stdout);

  GKfree(&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, &perm, &iperm, -1);
}  



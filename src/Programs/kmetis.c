/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This file contains the driving routine for kmetis
 *
 * Started 8/28/94
 * George
 *
 * $Id: kmetis.c,v 1.6 2003/07/31 16:15:50 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(idxtype argc, char *argv[])
{
  idxtype i, nparts, options[10];
  idxtype *part;
  float rubvec[MAXNCON], lbvec[MAXNCON];
  GraphType graph;
  char filename[256];
  idxtype numflag = 0, wgtflag = 0, edgecut;
  timer TOTALTmr, METISTmr, IOTmr;

  if (argc != 3) {
    printf("Usage: %s <GraphFile> <Nparts>\n",argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);
  nparts = atoi(argv[2]);

  if (nparts < 2) {
    printf("The number of partitions should be greater than 1!\n");
    exit(0);
  }

  cleartimer(TOTALTmr);
  cleartimer(METISTmr);
  cleartimer(IOTmr);

  starttimer(TOTALTmr);
  starttimer(IOTmr);
  ReadGraph(&graph, filename, &wgtflag);

  /* The following is for debuging empty graphs... 
  graph.nedges = 0;
  idxset(graph.nvtxs+1, 0, graph.xadj);
  */

  if (graph.nvtxs <= 0) {
    printf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  stoptimer(IOTmr);

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Graph Information ---------------------------------------------------\n");
  printf("  Name: %s, #Vertices: %d, #Edges: %d, #Parts: %d\n", filename, graph.nvtxs, graph.nedges/2, nparts);
  if (graph.ncon > 1)
    printf("  Balancing Constraints: %d\n", graph.ncon);
  printf("\nK-way Partitioning... -----------------------------------------------\n");

  part = idxmalloc(graph.nvtxs, "main: part");
  options[0] = 0;

  starttimer(METISTmr);
  if (graph.ncon == 1) {
    METIS_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
          &wgtflag, &numflag, &nparts, options, &edgecut, part);
  }
  else {
    for (i=0; i<graph.ncon; i++)
      rubvec[i] = HORIZONTAL_IMBALANCE;

    METIS_mCPartGraphKway(&graph.nvtxs, &graph.ncon, graph.xadj, graph.adjncy, graph.vwgt, 
          graph.adjwgt, &wgtflag, &numflag, &nparts, rubvec, options, &edgecut, part);
  }
  stoptimer(METISTmr);

  ComputePartitionBalance(&graph, nparts, part, lbvec);

  printf("  %d-way Edge-Cut: %7d, Balance: ", nparts, edgecut);
  for (i=0; i<graph.ncon; i++)
    printf("%5.2f ", lbvec[i]);
  printf("\n");

  starttimer(IOTmr);
  WritePartition(filename, part, graph.nvtxs, nparts); 
  stoptimer(IOTmr);
  stoptimer(TOTALTmr);

  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Partitioning: \t\t %7.3f   (KMETIS time)\n", gettimer(METISTmr));
  printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
  printf("**********************************************************************\n");


  GKfree((void *)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, &part, LTERM);
}  



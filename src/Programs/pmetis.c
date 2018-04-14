/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * pmetis.c
 *
 * This file contains the driving routine for multilevel method
 *
 * Started 8/28/94
 * George
 *
 * $Id: pmetis.c,v 1.3 2002/08/10 06:57:51 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(idxtype argc, char *argv[])
{
  idxtype i, options[10];
  idxtype *part;
  float lbvec[MAXNCON];
  GraphType graph;
  idxtype numflag = 0, wgtflag = 0, edgecut;
  ParamType params;
  timer TOTALTmr, METISTmr, IOTmr;


  parse_cmdline(&params, argc, argv);


  if (params.nparts < 2) {
    printf("The number of partitions should be greater than 1!\n");
    exit(0);
  }

  cleartimer(TOTALTmr);
  cleartimer(METISTmr);
  cleartimer(IOTmr);

  starttimer(TOTALTmr);
  starttimer(IOTmr);
  ReadGraph(&graph, params.filename, &wgtflag);
  if (graph.nvtxs <= 0) {
    printf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  stoptimer(IOTmr);

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Graph Information ---------------------------------------------------\n");
  printf("  Name: %s, #Vertices: %d, #Edges: %d, #Parts: %d\n", params.filename, graph.nvtxs, graph.nedges/2, params.nparts);
  if (graph.ncon > 1)
    printf("  Balancing Constraints: %d\n", graph.ncon);
  printf("\nRecursive Partitioning... -------------------------------------------\n");

  part = idxmalloc(graph.nvtxs, "main: part");
  options[0] = 0;

  starttimer(METISTmr);
  if (graph.ncon == 1) {
    METIS_PartGraphRecursive(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
          &wgtflag, &numflag, &(params.nparts), options, &edgecut, part);
  }
  else {
    METIS_mCPartGraphRecursive(&graph.nvtxs, &graph.ncon, graph.xadj, graph.adjncy, graph.vwgt, 
          graph.adjwgt, &wgtflag, &numflag, &(params.nparts), options, &edgecut, part);
  }

  stoptimer(METISTmr);

  ComputePartitionBalance(&graph, params.nparts, part, lbvec);

  printf("  %d-way Edge-Cut: %7d, Balance: ", params.nparts, edgecut);
  for (i=0; i<graph.ncon; i++)
    printf("%5.2f ", lbvec[i]);
  printf("\n");


  starttimer(IOTmr);
  WritePartition(params.filename, part, graph.nvtxs, params.nparts); 
  stoptimer(IOTmr);
  stoptimer(TOTALTmr);

  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Partitioning: \t\t %7.3f   (PMETIS time)\n", gettimer(METISTmr));
  printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
  printf("**********************************************************************\n");


  GKfree((void *)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, &part, LTERM);
}  



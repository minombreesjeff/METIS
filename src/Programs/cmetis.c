/*
 * Copyright 2003, Regents of the University of Minnesota
 *
 * cmetis.c
 *
 * This file contains the driving routine for partitioning for 
 * sub-domain direct factorization.
 *
 * Started 3/25/03
 * George
 *
 * $Id: cmetis.c,v 1.5 2003/04/08 12:52:59 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(idxtype argc, char *argv[])
{
  idxtype i, options[10], nnodes, nclean, naclean, ndirty, maxdepth;
  idxtype *part, *sflag;
  float lbvec[MAXNCON];
  GraphType graph;
  idxtype numflag = 0, wgtflag = 0, edgecut;
  ParamType params;
  void *cinfo;
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
  ReadCoordinates(&graph, params.xyzfilename);
  stoptimer(IOTmr);

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Graph Information ---------------------------------------------------\n");
  printf("  Name: %s, #Vertices: %d, #Edges: %d, #Parts: %d\n", params.filename, graph.nvtxs, graph.nedges/2, params.nparts);
  if (graph.ncon > 1)
    printf("  Balancing Constraints: %d\n", graph.ncon);
  printf("\nRecursive Partitioning... -------------------------------------------\n");

  part  = idxmalloc(graph.nvtxs, "main: part");
  sflag = idxsmalloc(graph.nvtxs, 1, "main: sflag");
  options[0] = 0;

  starttimer(METISTmr);
  cinfo = METIS_PartGraphForContact(&graph.nvtxs, graph.xadj, graph.adjncy, graph.coords, sflag, 
                    &numflag, &(params.nparts), options, &edgecut, part);

  METIS_UpdateContactInfo(cinfo, &graph.nvtxs, graph.coords, sflag);
  stoptimer(METISTmr);

  ComputePartitionBalance(&graph, params.nparts, part, lbvec);

  graph.vwgt = ismalloc(graph.nvtxs, 1, "main: graph->vwgt");

  printf("  %d-way Edge-Cut: %7d, Volume: %7d, Balance: ", params.nparts, edgecut, 
                                                           ComputeVolume(&graph, part));
  for (i=0; i<graph.ncon; i++)
    printf("%5.2f ", lbvec[i]);
  printf("\n");

  starttimer(IOTmr);
  WritePartition(params.filename, part, graph.nvtxs, params.nparts); 
  stoptimer(IOTmr);
  stoptimer(TOTALTmr);

  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Partitioning: \t\t %7.3f   (CMETIS time)\n", gettimer(METISTmr));
  printf("  Total:        \t\t %7.3f\n", gettimer(TOTALTmr));
  printf("**********************************************************************\n");

  GKfree((void *)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, &graph.coords, &part, &sflag, LTERM);
}  


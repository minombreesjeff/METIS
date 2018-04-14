/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * partnmesh.c
 *
 * This file reads in the element node connectivity array of a mesh and 
 * partitions both the elements and the nodes using KMETIS on the dual graph.
 *
 * Started 9/29/97
 * George
 *
 * $Id: partnmesh.c,v 1.2 2002/08/10 06:02:54 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(idxtype argc, char *argv[])
{
  idxtype i, j, ne, nn, etype, mtype, numflag=0, nparts, edgecut;
  idxtype *elmnts, *epart, *npart, *metype, *weights;
  timer IOTmr, DUALTmr;
  char etypestr[5][5] = {"TRI", "TET", "HEX", "QUAD", "LINE"};
  GraphType graph;

  if (argc != 3) {
    printf("Usage: %s <meshfile> <nparts>\n",argv[0]);
    exit(0);
  }



  nparts = atoi(argv[2]);
  if (nparts < 2) {
    printf("nparts must be greater than one.\n");
    exit(0);
  }
   
  cleartimer(IOTmr);
  cleartimer(DUALTmr);
 
  mtype=MeshType(argv[1]);
  ne=MixedElements(argv[1]);
  metype = idxmalloc(ne, "main: metype");
  weights = idxmalloc(ne, "main: weights");

  starttimer(IOTmr);
 
  if(mtype==1)
       elmnts = ReadMesh(argv[1], &ne, &nn, &etype);
  else if(mtype==3)
       elmnts = ReadMeshWgt(argv[1], &ne, &nn, &etype, weights);
  else if(mtype==0)
       elmnts = ReadMixedMesh(argv[1], &ne, &nn, metype);
  else
       elmnts = ReadMixedMeshWgt(argv[1], &ne, &nn, metype, weights);

  stoptimer(IOTmr);

  epart = idxmalloc(ne, "main: epart");
  npart = idxmalloc(nn, "main: npart");

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Mesh Information ----------------------------------------------------\n");
  if (mtype==1)
  printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n", argv[1], ne, nn, etypestr[etype-1]);
  else
  printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n", argv[1], ne, nn, "Mixed");
  printf("Partitioning Nodal Graph... -----------------------------------------\n");


  starttimer(DUALTmr);
  
  if (mtype==1 || mtype==3)
  METIS_PartMeshNodal(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
  else 
  METIS_PartMixedMeshNodal(&ne, &nn, elmnts, metype, &numflag, &nparts, &edgecut, epart, npart);
  
  stoptimer(DUALTmr);

  printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, edgecut, ComputeElementBalance(ne, nparts, epart));

  starttimer(IOTmr);
  WriteMeshPartition(argv[1], nparts, ne, epart, nn, npart);
  stoptimer(IOTmr);


  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Partitioning: \t\t %7.3f\n", gettimer(DUALTmr));
  printf("**********************************************************************\n");

/*
  graph.nvtxs = ne;
  graph.xadj = idxmalloc(ne+1, "xadj");
  graph.vwgt = idxsmalloc(ne, 1, "vwgt");
  graph.adjncy = idxmalloc(10*ne, "adjncy");
  graph.adjwgt = idxsmalloc(10*ne, 1, "adjncy");

  METIS_MeshToDual(&ne, &nn, elmnts, &etype, &numflag, graph.xadj, graph.adjncy);

  ComputePartitionInfo(&graph, nparts, epart);

  GKfree((void *)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
*/

  GKfree((void *)&elmnts, &epart, &npart, &metype, &weights,  LTERM);

}



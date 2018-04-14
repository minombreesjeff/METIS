/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh2dual.c
 *
 * This file reads in the element node connectivity array of a mesh and writes
 * out its dual in the format suitable for Metis.
 *
 * Started 9/29/97
 * George
 *
 * $Id: mesh2dual.c,v 1.2 2002/08/10 06:02:53 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(idxtype argc, char *argv[])
{
  idxtype i, j, ne, nn, etype, mtype, cnt, numflag=0;
  idxtype *elmnts, *xadj, *adjncy, *metype;
  idxtype *conmat, *elms, *weights; 
  timer IOTmr, DUALTmr;
  char fileout[256], etypestr[5][5] = {"TRI", "TET", "HEX", "QUAD", "LINE"};

  if (argc <2) {
    printf("Usage: %s <meshfile> [confile]\n",argv[0]);
    exit(0);
  }


  mtype=MeshType(argv[1]);
  ne=MixedElements(argv[1]);
  metype = idxmalloc(ne, "main: metype");
  weights = idxmalloc(ne, "main: weights");
  
  if (mtype==1 || mtype==3){

  cleartimer(IOTmr);
  cleartimer(DUALTmr);

  starttimer(IOTmr);
  if (mtype==1)
     elmnts = ReadMesh(argv[1], &ne, &nn, &etype);
  else
     elmnts = ReadMeshWgt(argv[1], &ne, &nn, &etype, weights);
  stoptimer(IOTmr);

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Mesh Information ----------------------------------------------------\n");
  printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n", argv[1], ne, nn, etypestr[etype-1]);
  printf("Forming Dual Graph... -----------------------------------------------\n");

  xadj = idxmalloc(ne+1, "main: xadj");
  elms = idxsmalloc(ne+1, 0, "main: elms");


  starttimer(DUALTmr);
  cnt=METIS_MeshToDualCount(&ne, &nn, elmnts, elms, &etype, &numflag);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MeshToDual(&ne, &nn, elmnts, elms, &etype, &numflag, xadj, adjncy);
  stoptimer(DUALTmr);

  printf("  Dual Information: #Vertices: %d, #Edges: %d\n", ne, xadj[ne]/2);

  sprintf(fileout, "%s.dgraph", argv[1]);
  starttimer(IOTmr);
  if (mtype==1)
     WriteGraph(fileout, ne, xadj, adjncy);
  else
     WriteWgtGraph(fileout, ne, xadj, adjncy, weights);
     
  stoptimer(IOTmr);


  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Dual Creation:\t\t %7.3f\n", gettimer(DUALTmr));
  printf("**********************************************************************\n");

  }

  else {

  
  cleartimer(IOTmr);
  cleartimer(DUALTmr);

  starttimer(IOTmr);

  
  if(mtype==0)
     elmnts = ReadMixedMesh(argv[1], &ne, &nn, metype);
  else
     elmnts = ReadMixedMeshWgt(argv[1], &ne, &nn, metype, weights);

  if (argc==3)  
  conmat = ReadMgcnums(argv[2]);
  stoptimer(IOTmr);

   

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Mesh Information ----------------------------------------------------\n");
  printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n", argv[1], ne, nn, "Mixed");
  printf("Forming Dual Graph... ----------------------------------------------\n");

  xadj = idxmalloc(ne+1, "main: xadj");
  elms = idxsmalloc(ne+1, 0, "main: elms");

  
  starttimer(DUALTmr);

  if (argc==3){
  cnt=METIS_MixedMeshToDualCount(&ne, &nn, elmnts, elms, metype, &numflag, 
conmat, 1);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MixedMeshToDual(&ne, &nn, elmnts, elms, metype, &numflag, xadj, adjncy, 
 conmat, 1);
  } 
  else{
  cnt=METIS_MixedMeshToDualCount(&ne, &nn, elmnts, elms, metype, &numflag, 
conmat, 0);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MixedMeshToDual(&ne, &nn, elmnts, elms, metype, &numflag, xadj, adjncy, conmat, 0);
  } 
  stoptimer(DUALTmr);

  printf("  Dual Information: #Vertices: %d, #Edges: %d\n", ne, xadj[ne]/2);

  sprintf(fileout, "%s.dgraph", argv[1]);
  starttimer(IOTmr);

  if (mtype==0)
     WriteGraph(fileout, ne, xadj, adjncy);
  else
     WriteWgtGraph(fileout, ne, xadj, adjncy, weights);
  stoptimer(IOTmr);


  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Dual Creation:\t\t %7.3f\n", gettimer(DUALTmr));
  printf("**********************************************************************\n");

  }

  GKfree((void *)&elmnts, &xadj, &adjncy, &metype, &weights, &elms,  LTERM);
}



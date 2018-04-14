/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh2nodal.c
 *
 * This file reads in the element node connectivity array of a mesh and writes
 * out its dual in the format suitable for Metis.
 *
 * Started 9/29/97
 * George
 *
 * $Id: mesh2nodal.c,v 1.2 2002/08/10 06:02:53 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(idxtype argc, char *argv[])
{
  idxtype i, j, ne, nn, etype, mtype, numflag=0;
  idxtype *elmnts, *xadj, *adjncy, *metype, *weights;
  timer IOTmr, DUALTmr;
  char fileout[256], etypestr[5][5] = {"TRI", "TET", "HEX", "QUAD", "LINE"};

  if (argc != 2) {
    printf("Usage: %s <meshfile>\n",argv[0]);
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
  printf("Forming Nodal Graph... ----------------------------------------------\n");

  xadj = idxmalloc(nn+1, "main: xadj");
  adjncy = idxmalloc(20*nn, "main: adjncy");

  starttimer(DUALTmr);
  METIS_MeshToNodal(&ne, &nn, elmnts, &etype, &numflag, xadj, adjncy);
  stoptimer(DUALTmr);

  printf("  Nodal Information: #Vertices: %d, #Edges: %d\n", nn, xadj[nn]/2);

  sprintf(fileout, "%s.ngraph", argv[1]);
  starttimer(IOTmr);
  WriteGraph(fileout, nn, xadj, adjncy);
  stoptimer(IOTmr);


  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Nodal Creation:\t\t %7.3f\n", gettimer(DUALTmr));
  printf("**********************************************************************\n");


  }

  else{

  cleartimer(IOTmr);
  cleartimer(DUALTmr);

  starttimer(IOTmr);
  if(mtype==0)
     elmnts = ReadMixedMesh(argv[1], &ne, &nn, metype);
  else
     elmnts = ReadMixedMeshWgt(argv[1], &ne, &nn, metype, weights);
  stoptimer(IOTmr);


  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Mesh Information ----------------------------------------------------\n");
  printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n", argv[1], ne, nn, "Mixed");
  printf("Forming Nodal Graph... ----------------------------------------------\n");

  xadj = idxmalloc(nn+1, "main: xadj");
  adjncy = idxmalloc(20*nn, "main: adjncy");

  starttimer(DUALTmr);
  METIS_MixedMeshToNodal(&ne, &nn, elmnts, metype, &numflag, xadj, adjncy);
  stoptimer(DUALTmr);

  printf("  Nodal Information: #Vertices: %d, #Edges: %d\n", nn, xadj[nn]/2);

  sprintf(fileout, "%s.ngraph", argv[1]);
  starttimer(IOTmr);
  WriteGraph(fileout, nn, xadj, adjncy);
  stoptimer(IOTmr);


  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Nodal Creation:\t\t %7.3f\n", gettimer(DUALTmr));
  printf("**********************************************************************\n");


  }
  
  GKfree((void *)&elmnts, &xadj, &adjncy, &metype, &weights, LTERM);

}



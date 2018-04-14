/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * meshtest.c
 *
 * This file contains code for testing the mesh conversion routines
 *
 * Started 8/18/97
 * George
 *
 * $Id: meshtest.c,v 1.1 1997/11/04 23:19:52 karypis Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int i, j, nelmnts, nvtxs, ndims;
  idxtype *elmnts, *xadj, *adjncy, *part;
  FILE *fp;
  timer tmr1, tmr2;
  int nparts=32, edgecut, etype, numflag=1, wgtflag=0, options[10];

  if (argc != 5) {
    printf("Usage: %s <meshfile> <nn> <ne> <ndims>\n",argv[0]);
    exit(0);
  }

  nvtxs = atoi(argv[2]);
  nelmnts = atoi(argv[3]);
  ndims = atoi(argv[4]);

  if (ndims == 2)
    etype = 1;
  else
    etype = 2;
    
  /* Work on the dual graph */   
  elmnts = idxmalloc((ndims+1)*nelmnts, "main: elmnts");
  fp = fopen(argv[1], "rb");
  fread(elmnts, sizeof(int), (ndims+1)*nelmnts, fp);
  fclose(fp);

  printf("Read a mesh with %d elements\n", nelmnts);

  xadj = idxmalloc(nelmnts+1, "main: xadj");
  adjncy = idxmalloc((ndims+1)*nelmnts, "main: adjncy");

  cleartimer(tmr1); starttimer(tmr1);
  METIS_MeshToDual(&nelmnts, &nvtxs, elmnts, &etype, &numflag, xadj, adjncy);
  stoptimer(tmr1);
  printf("Tmr1: %7.3f\n", gettimer(tmr1));

  printf("Dual contains %d edges!\n", xadj[nelmnts]-1);
  free(elmnts);

  part = idxmalloc(nelmnts, "main: part");

  /* Call KMETIS */
  options[0] = 0;

  cleartimer(tmr1); starttimer(tmr1);
  METIS_PartGraphKway(&nelmnts, xadj, adjncy, NULL, NULL, &wgtflag, &numflag, &nparts, options, &edgecut, part);
  stoptimer(tmr1);
  printf("Partitioning time: %7.3f\n", gettimer(tmr1));
  printf("Cut: %d\n", edgecut);

  GKfree(&xadj, &adjncy, &part, -1);


  /* Work on the nodal graph */   
  elmnts = idxmalloc((ndims+1)*nelmnts, "main: elmnts");
  fp = fopen(argv[1], "rb");
  fread(elmnts, sizeof(int), (ndims+1)*nelmnts, fp);
  fclose(fp);

  printf("Read a mesh with %d elements\n", nelmnts);

  xadj = idxmalloc(nvtxs+1, "main: xadj");
  adjncy = idxmalloc(16*nvtxs, "main: adjncy");

  cleartimer(tmr1);
  starttimer(tmr1);
  METIS_MeshToNodal(&nelmnts, &nvtxs, elmnts, &etype, &numflag, xadj, adjncy);
  stoptimer(tmr1);
  printf("Tmr1: %7.3f\n", gettimer(tmr1));

  printf("Nodal contains %d edges!\n", xadj[nvtxs]-1);
  free(elmnts);

  part = idxmalloc(nvtxs, "main: part");

  /* Call KMETIS */
  options[0] = 0;

  cleartimer(tmr1); starttimer(tmr1);
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, NULL, NULL, &wgtflag, &numflag, &nparts, options, &edgecut, part);
  stoptimer(tmr1);
  printf("Partitioning time: %7.3f\n", gettimer(tmr1));
  printf("Cut: %d\n", edgecut);

  GKfree(&xadj, &adjncy, &part, -1);


}



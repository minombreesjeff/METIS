/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * main.c
 *
 * This file contains code for testing ParMetis
 *
 * Started 10/19/95
 * George
 *
 * $Id: main.c,v 1.8 1997/07/18 00:32:09 karypis Exp $
 *
 */

#include <par_kmetis.h>

/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int npes, mype;
  MPI_Comm comm;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  if (argc != 2) {
    if (mype == 0)
      printf("Usage: %s <graph-file>\n", argv[0]);
    exit(0);
  }

  TestParMetis(argv[1], MPI_COMM_WORLD);

  MPI_Finalize();
}



/***********************************************************************************
* This function is the testing routine for the adaptive multilevel partitioning code.
* It computes a partition from scratch, it then moves the graph and changes some
* of the vertex weights and then call the adaptive code.
************************************************************************************/
void TestParMetis(char *filename, MPI_Comm comm)
{
  int i, j, k, l, npes, mype, realcut;
  int opt1, opt2;
  GraphType graph, mgraph;
  idxtype *part, *mpart, *order, *sizes;
  int options[5];
  float *xyz;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);


  ReadGraph(&graph, filename, comm);
  xyz = ReadCoordinates(&graph, filename, comm);

  part = idxmalloc(graph.nvtxs, "TestParMetis: part");


  /*======================================================================
  / PARKMETIS 
  /=======================================================================*/
  options[3] = 0;
  options[4] = 3;
  for (opt1=0; opt1<=300; opt1+=100) {
    for (opt2=1; opt2<=2; opt2++) {
      options[1] = opt1;
      options[2] = opt2;

      if (mype == 0)
        printf("\nTesting PARKMETIS with options[1-4] = {%d %d %d %d}\n", options[1], options[2], options[3], options[4]);

      PARKMETIS(graph.vtxdist, graph.xadj, NULL, graph.adjncy, NULL, part, options, comm);

      realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
      if (mype == 0) {
        if (realcut == options[0])
          printf("PARKMETIS reported a cut of %d which is correct!\n", options[0]);
        else
          printf("PARKMETIS reported a cut of %d which is incorrect (realcut = %d)!\n", options[0], realcut);
      }
    }
  }


  /*======================================================================
  / PARGKMETIS 
  /=======================================================================*/
  options[3] = 0;
  options[4] = 3;
  for (opt1=0; opt1<=300; opt1+=100) {
    for (opt2=1; opt2<=2; opt2++) {
      options[1] = opt1;
      options[2] = opt2;

      if (mype == 0)
        printf("\nTesting PARGKMETIS with options[1-4] = {%d %d %d %d}\n", options[1], options[2], options[3], options[4]);

      PARGKMETIS(graph.vtxdist, graph.xadj, NULL, graph.adjncy, NULL, 3, xyz, part, options, comm);

      realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
      if (mype == 0) {
        if (realcut == options[0])
          printf("PARGKMETIS reported a cut of %d which is correct!\n", options[0]);
        else
          printf("PARGKMETIS reported a cut of %d which is incorrect (realcut = %d)!\n", options[0], realcut);
      }
    }
  }


  /*======================================================================
  / PARGRMETIS 
  /=======================================================================*/
  options[3] = 0;
  options[4] = 3;

  if (mype == 0)
    printf("\nTesting PARGRMETIS with options[3-4] = {%d %d}\n", options[3], options[4]);

  PARGRMETIS(graph.vtxdist, graph.xadj, NULL, graph.adjncy, NULL, 3, xyz, part, options, comm);

  realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
  if (mype == 0) {
    if (realcut == options[0])
      printf("PARGRMETIS reported a cut of %d which is correct!\n", options[0]);
    else
      printf("PARGRMETIS reported a cut of %d which is incorrect (realcut = %d)!\n", options[0], realcut);
  }


  /*======================================================================
  / PARGMETIS 
  /=======================================================================*/
  options[3] = 0;
  options[4] = 3;

  if (mype == 0)
    printf("\nTesting PARGMETIS with options[3-4] = {%d %d}\n", options[3], options[4]);

  PARGMETIS(graph.vtxdist, graph.xadj, graph.adjncy, 3, xyz, part, options, comm);

  realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
  if (mype == 0) {
    if (realcut == options[0])
      printf("PARGMETIS reported a cut of %d which is correct!\n", options[0]);
    else
      printf("PARGMETIS reported a cut of %d which is incorrect (realcut = %d)!\n", options[0], realcut);
  }


  /*======================================================================
  / PARRMETIS 
  /=======================================================================*/
  options[3] = 0;
  options[4] = 3;

  if (mype == 0)
    printf("\nTesting PARRMETIS with options[3-4] = {%d %d}\n", options[3], options[4]);

  PARRMETIS(graph.vtxdist, graph.xadj, NULL, graph.adjncy, NULL, part, options, comm);

  realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
  if (mype == 0) {
    if (realcut == options[0])
      printf("PARRMETIS reported a cut of %d which is correct!\n", options[0]);
    else
      printf("PARRMETIS reported a cut of %d which is incorrect (realcut = %d)!\n", options[0], realcut);
  }


  /* Compute a good partition and move the graph. Do so quietly! */
  options[1] = 0;
  options[2] = 1;
  options[3] = 0;
  options[4] = 0;
  PARKMETIS(graph.vtxdist, graph.xadj, NULL, graph.adjncy, NULL, part, options, comm);
  TestMoveGraph(&graph, &mgraph, part, comm);
  mpart = idxmalloc(mgraph.nvtxs, "TestParMetis: mpart");


  /* Test #2 */
  options[3] = 0;
  options[4] = 3;

  if (mype == 0)
    printf("\nTesting PARRMETIS with options[3-4] = {%d %d} (after move)\n", options[3], options[4]);

  PARRMETIS(mgraph.vtxdist, mgraph.xadj, NULL, mgraph.adjncy, NULL, mpart, options, comm);

  realcut = ComputeRealCut2(graph.vtxdist, mgraph.vtxdist, part, mpart, filename, comm);
  if (mype == 0) {
    if (realcut == options[0])
      printf("PARRMETIS reported a cut of %d which is correct!\n", options[0]);
    else
      printf("PARRMETIS reported a cut of %d which is incorrect (realcut = %d)!\n", options[0], realcut);
  }


  /*======================================================================
  / PARUAMETIS, PARDAMETIS 
  /=======================================================================*/
  mgraph.vwgt = idxsmalloc(mgraph.nvtxs, 1, "TestParKMetis: mgraph.vwgt");
  AdaptGraph(&mgraph, 2, comm); 

  options[3] = 0;
  options[4] = 3;

  if (mype == 0)
    printf("\nTesting PARUAMETIS with options[3-4] = {%d %d}\n", options[3], options[4]);

  PARUAMETIS(mgraph.vtxdist, mgraph.xadj, mgraph.vwgt, mgraph.adjncy, NULL, mpart, options, comm);

  realcut = ComputeRealCut2(graph.vtxdist, mgraph.vtxdist, part, mpart, filename, comm);
  if (mype == 0) {
    if (realcut == options[0])
      printf("PARUAMETIS reported a cut of %d which is correct!\n", options[0]);
    else
      printf("PARUAMETIS reported a cut of %d which is incorrect (realcut = %d)!\n", options[0], realcut);
  }

  if (mype == 0)
    printf("\nTesting PARDAMETIS with options[3-4] = {%d %d}\n", options[3], options[4]);

  PARDAMETIS(mgraph.vtxdist, mgraph.xadj, mgraph.vwgt, mgraph.adjncy, NULL, mpart, options, comm);

  realcut = ComputeRealCut2(graph.vtxdist, mgraph.vtxdist, part, mpart, filename, comm);
  if (mype == 0) {
    if (realcut == options[0])
      printf("PARUAMETIS reported a cut of %d which is correct!\n", options[0]);
    else
      printf("PARUAMETIS reported a cut of %d which is incorrect (realcut = %d)!\n", options[0], realcut);
  }


  /*======================================================================
  / PAOMETIS 
  /=======================================================================*/
  sizes = idxmalloc(2*npes, "TestParMetis: sizes");
  order = idxmalloc(graph.nvtxs, "TestParMetis: sizes");

  options[1] = 0;
  options[3] = 0;
  options[4] = 3;

  for (opt2=1; opt2<=2; opt2++) {
    options[2] = opt2;

    if (mype == 0)
      printf("\nTesting PAROMETIS with options[1-4] = {%d %d %d %d}\n", options[1], options[2], options[3], options[4]);

    PAROMETIS(graph.vtxdist, graph.xadj, NULL, graph.adjncy, NULL, order, sizes, options, comm);
  }


  GKfree(&part, &mpart, &order, &sizes, -1);

}





/*************************************************************************
* This function implements a simple graph adaption strategy.
**************************************************************************/
void AdaptGraph(GraphType *graph, int afactor, MPI_Comm comm)
{
  int i, j, k, nvtxs, nadapt;
  int npes, mype, mypwgt, max, min, sum;
  idxtype *vwgt, *perm;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  srand(mype);
  srand48(mype);

  nvtxs = graph->nvtxs;
  vwgt = graph->vwgt;

  perm = idxmalloc(nvtxs, "AdaptGraph: perm");
  FastRandomPermute(nvtxs, perm, 1);

  nadapt = RandomInRange(nvtxs);

  for (i=0; i<nadapt; i++)
    vwgt[perm[i]] = afactor*vwgt[perm[i]];

  mypwgt = idxsum(nvtxs, vwgt);

  MPI_Allreduce((void *)&mypwgt, (void *)&max, 1, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce((void *)&mypwgt, (void *)&min, 1, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce((void *)&mypwgt, (void *)&sum, 1, MPI_INT, MPI_SUM, comm);

  if (mype == 0)
    printf("Initial Load Imbalance: %5.4f, [%5d %5d %5d]\n", (1.0*max*npes)/(1.0*sum), min, max, sum);

  free(perm);
}


/******************************************************************************
* This function takes a partition vector that is distributed and reads in
* the original graph and computes the edgecut
*******************************************************************************/
int ComputeRealCut(idxtype *vtxdist, idxtype *part, char *filename, MPI_Comm comm)
{
  int i, j, k, nvtxs, nedges, mype, npes, cut;
  idxtype *xadj, *adjncy, *gpart, scratch;
  MPI_Status status;
  FILE *fpin;


  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (mype != 0) {
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_DATATYPE, 0, 1, comm);
  }
  else {  /* Processor 0 does all the rest */
    gpart = idxmalloc(vtxdist[npes], "ComputeRealCut: gpart");
    idxcopy(vtxdist[1], part, gpart);

    for (i=1; i<npes; i++) 
      MPI_Recv((void *)(gpart+vtxdist[i]), vtxdist[i+1]-vtxdist[i], IDX_DATATYPE, i, 1, comm, &status);

    /* Ok, now read the graph from the file */
    if ((fpin = fopen(filename, "rb")) == NULL) 
      errexit("Failed to open file %s", filename);

    fread(&scratch, sizeof(idxtype), 1, fpin);
    nvtxs = scratch;
    fread(&scratch, sizeof(idxtype), 1, fpin);
    nedges = scratch;

    xadj = idxmalloc(nvtxs+1, "ComputeRealCut: xadj");
    adjncy = idxmalloc(nedges, "ComputeRealCut: adjncy");
    fread(xadj, sizeof(idxtype), nvtxs+1, fpin);
    fread(adjncy, sizeof(idxtype), nedges, fpin);
    fclose(fpin);

    /* OK, now compute the cut */
    for (cut=0, i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (gpart[i] != gpart[adjncy[j]])
          cut++;
      }
    }
    cut = cut/2;

    GKfree(&gpart, &xadj, &adjncy, -1);

    return cut;
  }
}


/******************************************************************************
* This function takes a partition vector that is distributed and reads in
* the original graph and computes the edgecut
*******************************************************************************/
int ComputeRealCut2(idxtype *vtxdist, idxtype *mvtxdist, idxtype *part, idxtype *mpart, char *filename, MPI_Comm comm)
{
  int i, j, k, nvtxs, nedges, mype, npes, cut;
  idxtype *xadj, *adjncy, *gpart, *gmpart, *perm, *sizes, scratch;
  MPI_Status status;
  FILE *fpin;


  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (mype != 0) {
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_DATATYPE, 0, 1, comm);
    MPI_Send((void *)mpart, mvtxdist[mype+1]-mvtxdist[mype], IDX_DATATYPE, 0, 1, comm);
  }
  else {  /* Processor 0 does all the rest */
    gpart = idxmalloc(vtxdist[npes], "ComputeRealCut: gpart");
    idxcopy(vtxdist[1], part, gpart);
    gmpart = idxmalloc(mvtxdist[npes], "ComputeRealCut: gmpart");
    idxcopy(mvtxdist[1], mpart, gmpart);

    for (i=1; i<npes; i++) {
      MPI_Recv((void *)(gpart+vtxdist[i]), vtxdist[i+1]-vtxdist[i], IDX_DATATYPE, i, 1, comm, &status);
      MPI_Recv((void *)(gmpart+mvtxdist[i]), mvtxdist[i+1]-mvtxdist[i], IDX_DATATYPE, i, 1, comm, &status);
    }

    /* OK, now go and reconstruct the permutation to go from the graph to mgraph */
    perm = idxmalloc(vtxdist[npes], "ComputeRealCut: perm");
    sizes = idxsmalloc(npes+1, 0, "ComputeRealCut: sizes");

    for (i=0; i<vtxdist[npes]; i++)
      sizes[gpart[i]]++;
    MAKECSR(i, npes, sizes);
    for (i=0; i<vtxdist[npes]; i++)
      perm[i] = sizes[gpart[i]]++;


    /* Ok, now read the graph from the file */
    if ((fpin = fopen(filename, "rb")) == NULL) 
      errexit("Failed to open file %s", filename);

    fread(&scratch, sizeof(idxtype), 1, fpin);
    nvtxs = scratch;
    fread(&scratch, sizeof(idxtype), 1, fpin);
    nedges = scratch;

    xadj = idxmalloc(nvtxs+1, "ComputeRealCut: xadj");
    adjncy = idxmalloc(nedges, "ComputeRealCut: adjncy");
    fread(xadj, sizeof(idxtype), nvtxs+1, fpin);
    fread(adjncy, sizeof(idxtype), nedges, fpin);
    fclose(fpin);

    /* OK, now compute the cut */
    for (cut=0, i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (gmpart[perm[i]] != gmpart[perm[adjncy[j]]])
          cut++;
      }
    }
    cut = cut/2;

    GKfree(&gpart, &gmpart, &perm, &sizes, &xadj, &adjncy, -1);

    return cut;
  }
}



/******************************************************************************
* This function takes a graph and its partition vector and creates a new
* graph corresponding to the one after the movement
*******************************************************************************/
void TestMoveGraph(GraphType *ograph, GraphType *omgraph, idxtype *part, MPI_Comm comm)
{
  int i, j, k;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;
  int options[5] = {0, 0, 1, 0, 0};

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  SetUpCtrl(&ctrl, options, comm); 
  graph = SetUpGraph(&ctrl, ograph->vtxdist, ograph->xadj, NULL, ograph->adjncy, NULL);
  PreAllocateMemory(&ctrl, graph, &wspace);

  SetUp(&ctrl, graph, &wspace);
  graph->where = part;
  mgraph = MoveGraph(&ctrl, graph, &wspace);

  omgraph->nvtxs = mgraph->nvtxs;
  omgraph->vtxdist = mgraph->vtxdist;
  omgraph->xadj = mgraph->xadj;
  omgraph->adjncy = mgraph->adjncy;
  mgraph->vtxdist = NULL;
  mgraph->xadj = NULL;
  mgraph->adjncy = NULL;
  FreeGraph(mgraph);

  graph->where = NULL;
  FreeInitialGraph(graph, 1, 1);
  FreeWSpace(&wspace);
}  


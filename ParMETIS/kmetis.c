/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This is the entry point of parallel kmetis
 *
 * Started 10/19/96
 * George
 *
 * $Id: kmetis.c,v 1.10 1997/07/18 00:32:08 karypis Exp $
 *
 */

#include <par_kmetis.h>


/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partitionioner. 
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
void PARKMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
               idxtype *part, int *options, MPI_Comm comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    options[OPTION_CUT] = 0;
    return;
  }

  SetUpCtrl(&ctrl, options, comm);

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 25*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt);

  tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl.CoarsenTo;

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  if (npes%2 == 0 && npes > 2)
    FldGlobal_Partition(&ctrl, graph, &wspace, 0);
  else
    Global_Partition(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  options[OPTION_CUT] = graph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f [%d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt), 
          graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs));

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  FreeInitialGraph(graph, vwgt == NULL, adjwgt == NULL);
  FreeWSpace(&wspace);

}



/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partition 
* refinement algorithm.  It does not compute an initial partition but refines
* an existing partition using multilevel $k$-way refinement. 
* This function utilizes local coarsening.
************************************************************************************/
void PARRMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
               idxtype *part, int *options, MPI_Comm comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    options[OPTION_CUT] = 0;
    return;
  }

  SetUpCtrl(&ctrl, options, comm);

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 50*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt);

  tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl.CoarsenTo;

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Refine_Partition(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  options[OPTION_CUT] = graph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f [%d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt), 
          graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs));

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  FreeInitialGraph(graph, vwgt == NULL, adjwgt == NULL);
  FreeWSpace(&wspace);
}


/***********************************************************************************
* This function is the entry point of the parallel multilevel undirected diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
void PARUAMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                idxtype *part, int *options, MPI_Comm comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;
  int nmoved, maxin, maxout;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    options[OPTION_CUT] = 0;
    return;
  }

  SetUpCtrl(&ctrl, options, comm);

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 70*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt);
  graph->vsize = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vsize");

  tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl.CoarsenTo;

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  AdaptiveUndirected_Partition(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  options[OPTION_CUT] = graph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_INFO, ComputeMoveStatistics(&ctrl, graph, &nmoved, &maxin, &maxout));
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f \nNMoved: %d %d %d %d [%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt),
          nmoved, maxin, maxout, maxin+maxout, npes, graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs));

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  GKfree(&graph->vsize, -1);
  FreeInitialGraph(graph, vwgt == NULL, adjwgt == NULL);
  FreeWSpace(&wspace);
}


/***********************************************************************************
* This function is the entry point of the parallel multilevel directed diffusion
* algorithm. It uses directed diffusion at the coarsest graph followed by adaptive 
* k-way refinement. This function utilizes local coarsening.
************************************************************************************/
void PARDAMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                idxtype *part, int *options, MPI_Comm comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;
  int nmoved, maxin, maxout;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    options[OPTION_CUT] = 0;
    return;
  }

  SetUpCtrl(&ctrl, options, comm);

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 50*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt);
  graph->vsize = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vsize");

  tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl.CoarsenTo;

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  AdaptiveDirected_Partition(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  options[OPTION_CUT] = graph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_INFO, ComputeMoveStatistics(&ctrl, graph, &nmoved, &maxin, &maxout));
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f \nNMoved: %d %d %d %d [%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt),
          nmoved, maxin, maxout, maxin+maxout, npes, graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs));

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  GKfree(&graph->vsize, -1);
  FreeInitialGraph(graph, vwgt == NULL, adjwgt == NULL);
  FreeWSpace(&wspace);
}




/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the 
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void PAROMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                idxtype *order, idxtype *sizes, int *options, MPI_Comm comm)
{
  int i, j, k, min, max, tvwgt, nparts;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;
  idxtype *morder;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (!ispow2(npes)) {
    if (mype == 0)
      printf("Error: The number of processors must be a power of 2!\n");
    return;
  }

  SetUpCtrl(&ctrl, options, comm);

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 1);

  ctrl.ipart = IPART_RB;
  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 50*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt);

  tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl.CoarsenTo;

  PreAllocateMemory(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial k-way partitioning 
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  if (npes > 2)
    FldGlobal_Partition(&ctrl, graph, &wspace, 0);
  else
    Global_Partition(&ctrl, graph, &wspace);


  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  /*=======================================================
   * Move the graph according to the partitioning
   =======================================================*/

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.MoveTmr));

  mgraph = MoveGraph(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.MoveTmr));

  /*=======================================================
   * Now compute an ordering of the moved graph
   =======================================================*/

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  FreeWSpace(&wspace);
  PreAllocateMemory(&ctrl, mgraph, &wspace);

  ctrl.ipart = options[OPTION_IPART];
  ctrl.CoarsenTo = amin(vtxdist[npes]-1, amax(20*npes, 1000));
  mgraph->maxvwgt = graph->maxvwgt;

  morder = idxmalloc(mgraph->nvtxs, "PAROMETIS: morder");
  MultilevelOrder(&ctrl, mgraph, morder, sizes, &wspace);

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(&ctrl, graph, order, morder, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 0);

  free(morder);
  FreeGraph(mgraph);
  FreeInitialGraph(graph, vwgt == NULL, adjwgt == NULL);
  FreeWSpace(&wspace);

}


/***********************************************************************************
* This function is the entry point of the parallel kmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
void PARGKMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm)
{
  int i, j, k, min, max, tvwgt, nparts;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    options[OPTION_CUT] = 0;
    return;
  }

  SetUpCtrl(&ctrl, options, comm);

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 25*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt);

  tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl.CoarsenTo;

  PreAllocateMemory(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial geometric partitioning 
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, ndims, xyz, 1, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  /*=======================================================
   * Move the graph according to the partitioning
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.MoveTmr));

  mgraph = MoveGraph(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.MoveTmr));

  if (ctrl.dbglvl&DBG_INFO) {
    ComputePartitionParams(&ctrl, graph, &wspace);
    rprintf(&ctrl, "XYZ Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt), 
          graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs,
          GlobalSEMax(&ctrl, graph->nrecv), GlobalSESum(&ctrl, graph->nrecv), 
          GlobalSEMax(&ctrl, graph->nsend), GlobalSESum(&ctrl, graph->nsend));
  }

  /*=======================================================
   * Now compute the partition of the moved graph
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  FreeWSpace(&wspace);
  PreAllocateMemory(&ctrl, mgraph, &wspace);

  mgraph->maxvwgt = graph->maxvwgt;

  if (npes%2 == 0 && npes > 2)
    FldGlobal_Partition(&ctrl, mgraph, &wspace, 0);
  else
    Global_Partition(&ctrl, mgraph, &wspace);

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(&ctrl, graph, part, mgraph->where, &wspace);

  options[OPTION_CUT] = mgraph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));

  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          mgraph->mincut, 1.0*npes*mgraph->gpwgts[idxamax(npes, mgraph->gpwgts)]/(1.0*tvwgt), 
          mgraph->gpwgts[idxamax(npes, mgraph->gpwgts)], tvwgt, mgraph->gnvtxs, 
          GlobalSEMax(&ctrl, mgraph->nrecv), GlobalSESum(&ctrl, mgraph->nrecv), 
          GlobalSEMax(&ctrl, mgraph->nsend), GlobalSESum(&ctrl, mgraph->nsend)));

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  FreeGraph(mgraph);
  FreeInitialGraph(graph, vwgt == NULL, adjwgt == NULL);
  FreeWSpace(&wspace);

}




/***********************************************************************************
* This function is the entry point of the parallel rmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
void PARGRMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm)
{
  int i, j, k, min, max, tvwgt, nparts;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    options[OPTION_CUT] = 0;
    return;
  }

  SetUpCtrl(&ctrl, options, comm);

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 25*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt);

  tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl.CoarsenTo;

  PreAllocateMemory(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial geometric partitioning 
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, ndims, xyz, 1, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  /*=======================================================
   * Move the graph according to the partitioning
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.MoveTmr));

  mgraph = MoveGraph(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.MoveTmr));

  if (ctrl.dbglvl&DBG_INFO) {
    ComputePartitionParams(&ctrl, graph, &wspace);
    rprintf(&ctrl, "XYZ Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt), 
          graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs,
          GlobalSEMax(&ctrl, graph->nrecv), GlobalSESum(&ctrl, graph->nrecv), 
          GlobalSEMax(&ctrl, graph->nsend), GlobalSESum(&ctrl, graph->nsend));
  }

  /*=======================================================
   * Now compute the partition of the moved graph
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  FreeWSpace(&wspace);
  PreAllocateMemory(&ctrl, mgraph, &wspace);

  mgraph->maxvwgt = graph->maxvwgt;

  Refine_Partition(&ctrl, mgraph, &wspace);

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(&ctrl, graph, part, mgraph->where, &wspace);

  options[OPTION_CUT] = mgraph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));

  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          mgraph->mincut, 1.0*npes*mgraph->gpwgts[idxamax(npes, mgraph->gpwgts)]/(1.0*tvwgt), 
          mgraph->gpwgts[idxamax(npes, mgraph->gpwgts)], tvwgt, mgraph->gnvtxs, 
          GlobalSEMax(&ctrl, mgraph->nrecv), GlobalSESum(&ctrl, mgraph->nrecv), 
          GlobalSEMax(&ctrl, mgraph->nsend), GlobalSESum(&ctrl, mgraph->nsend)));

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  FreeGraph(mgraph);
  FreeInitialGraph(graph, vwgt == NULL, adjwgt == NULL);
  FreeWSpace(&wspace);

}




/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the 
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void PARGMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int ndims, float *xyz, 
               idxtype *part, int *options, MPI_Comm comm)
{
  int i, j, k, min, max, tvwgt, nparts;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    options[OPTION_CUT] = 0;
    return;
  }

  SetUpCtrl(&ctrl, options, comm);

  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  ctrl.CoarsenTo = amin(vtxdist[npes]-1, 25*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, NULL, adjncy, NULL);

  tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl.CoarsenTo;

  PreAllocateMemory(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial geometric partitioning 
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, ndims, xyz, 0, &wspace);

  options[OPTION_CUT] = -1;
  idxcopy(graph->nvtxs, graph->where, part);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(comm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  if (ctrl.dbglvl&DBG_INFO) {
    SetUp(&ctrl, graph, &wspace);
    graph->where = (idxtype *)realloc(graph->where, sizeof(idxtype)*(graph->nvtxs+graph->nrecv));
    ComputePartitionParams(&ctrl, graph, &wspace);
    options[0] = graph->mincut;

    rprintf(&ctrl, "XYZ Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt), 
          graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs,
          GlobalSEMax(&ctrl, graph->nrecv), GlobalSESum(&ctrl, graph->nrecv), 
          GlobalSEMax(&ctrl, graph->nsend), GlobalSESum(&ctrl, graph->nsend));
  }


  if (options[OPTION_NUMBERING] == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  FreeInitialGraph(graph, 1, 1);
  FreeWSpace(&wspace);
}


/*************************************************************************
* This function setsup the CtrlType structure
**************************************************************************/
GraphType *SetUpGraph(CtrlType *ctrl, idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt)
{
  GraphType *graph;

  graph = CreateGraph();
  graph->level = 0;
  graph->gnvtxs = vtxdist[ctrl->npes];
  graph->nvtxs = vtxdist[ctrl->mype+1]-vtxdist[ctrl->mype];
  graph->nedges = xadj[graph->nvtxs];
  graph->xadj = xadj;
  graph->vwgt = vwgt;
  graph->adjncy = adjncy;
  graph->adjwgt = adjwgt;
  graph->vtxdist = vtxdist;

  if (graph->vwgt == NULL)
    graph->vwgt = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vwgt");
  if (graph->adjwgt == NULL)
    graph->adjwgt = idxsmalloc(graph->nedges, 1, "Par_KMetis: vwgt");

  return graph;
}


/*************************************************************************
* This function setsup the CtrlType structure
**************************************************************************/
void SetUpCtrl(CtrlType *ctrl, int *options, MPI_Comm comm)
{
  int j;

  MPI_Comm_rank(comm, &ctrl->mype);
  MPI_Comm_size(comm, &ctrl->npes);

  ctrl->nparts = ctrl->npes;  /* Set the # of partitions equal to the # of PEs */
  ctrl->gcomm = comm; 
  ctrl->comm = comm;
  ctrl->dbglvl = options[OPTION_DBGLVL];
  ctrl->foldf = options[OPTION_FOLDF];
  ctrl->ipart = options[OPTION_IPART];
  ctrl->xyztype = XYZ_SPFILL;

  /* Take care the ipart type for non-power of two */
  if (!ispow2(ctrl->npes))
    ctrl->ipart = IPART_SER;

  srand(ctrl->mype);
  srand48(ctrl->mype);
}


/*************************************************************************
* This function changes the numbering from 1 to 0 or 0 to 1
**************************************************************************/
void ChangeNumbering(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *part, int npes, int mype, int from)
{
  int i, nvtxs, nedges;

  if (from == 1) {  /* Change it from 1 to 0 */
    for (i=0; i<npes+1; i++)
      vtxdist[i]--;

    nvtxs = vtxdist[mype+1]-vtxdist[mype];
    for (i=0; i<nvtxs+1; i++)
      xadj[i]--;

    nedges = xadj[nvtxs];
    for (i=0; i<nedges; i++)
      adjncy[i]--;
  }
  else {  /* Change it from 0 to 1 */
    for (i=0; i<npes+1; i++)
      vtxdist[i]++;

    nvtxs = vtxdist[mype+1]-vtxdist[mype];
    for (i=0; i<nvtxs+1; i++)
      xadj[i]++;

    nedges = xadj[nvtxs];
    for (i=0; i<nedges; i++)
      adjncy[i]++;

    for (i=0; i<nvtxs; i++)
      part[i]++;
  }
}


/*************************************************************************
* This function initializes the various timers
**************************************************************************/
void InitTimers(CtrlType *ctrl)
{
  cleartimer(ctrl->TotalTmr); 
  cleartimer(ctrl->InitPartTmr);
  cleartimer(ctrl->MatchTmr); 
  cleartimer(ctrl->ContractTmr); 
  cleartimer(ctrl->CoarsenTmr); 
  cleartimer(ctrl->RefTmr);
  cleartimer(ctrl->SetupTmr); 
  cleartimer(ctrl->ProjectTmr); 
  cleartimer(ctrl->KWayInitTmr); 
  cleartimer(ctrl->KWayTmr);
  cleartimer(ctrl->MoveTmr);

  cleartimer(ctrl->AuxTmr1); 
  cleartimer(ctrl->AuxTmr2); 
  cleartimer(ctrl->AuxTmr3);
  cleartimer(ctrl->AuxTmr4); 
  cleartimer(ctrl->AuxTmr5); 
  cleartimer(ctrl->AuxTmr6);
}


/*************************************************************************
* This function prints timing information about KMETIS
**************************************************************************/
void PrintTimingInfo(CtrlType *ctrl)
{
/*  PrintTimer(ctrl, ctrl->CoarsenTmr,  " Coarsening"); */
  PrintTimer(ctrl, ctrl->SetupTmr,    "      Setup");
  PrintTimer(ctrl, ctrl->MatchTmr,    "   Matching");
  PrintTimer(ctrl, ctrl->ContractTmr, "Contraction");
  PrintTimer(ctrl, ctrl->InitPartTmr, "   InitPart");
/*  PrintTimer(ctrl, ctrl->RefTmr,      " Refinement"); */
  PrintTimer(ctrl, ctrl->ProjectTmr,  "    Project");
  PrintTimer(ctrl, ctrl->KWayInitTmr, " Initialize");
  PrintTimer(ctrl, ctrl->KWayTmr,     "      K-way");
  PrintTimer(ctrl, ctrl->MoveTmr,     "       Move");
  PrintTimer(ctrl, ctrl->TotalTmr,    "      Total");
/*
  PrintTimer(ctrl, ctrl->AuxTmr1,     "       Aux1");
  PrintTimer(ctrl, ctrl->AuxTmr2,     "       Aux2");
  PrintTimer(ctrl, ctrl->AuxTmr3,     "       Aux3");
  PrintTimer(ctrl, ctrl->AuxTmr4,     "       Aux4");
  PrintTimer(ctrl, ctrl->AuxTmr5,     "       Aux5");
*/
  PrintTimer(ctrl, ctrl->AuxTmr6,     "       Aux6");
}


/*************************************************************************
* This function prints timer stat
**************************************************************************/
void PrintTimer(CtrlType *ctrl, timer tmr, char *msg)
{
  double sum, max, tsec;

  tsec = gettimer(tmr);
  MPI_Reduce((void *)&tsec, (void *)&sum, 1, MPI_DOUBLE, MPI_SUM, 0, ctrl->comm);

  tsec = gettimer(tmr);
  MPI_Reduce((void *)&tsec, (void *)&max, 1, MPI_DOUBLE, MPI_MAX, 0, ctrl->comm);

  if (ctrl->mype == 0 && sum != 0.0)
    printf("%s: Max: %7.3lf, Sum: %7.3lf, Balance: %7.3lf\n", 
            msg, max, sum, max*ctrl->npes/sum);
}


/*************************************************************************
* This function randomly permutes the locally stored adjacency lists
**************************************************************************/
void GraphRandomPermute(GraphType *graph) 
{
  int i, j, k, tmp;

  for (i=0; i<graph->nvtxs; i++) {
    for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++) {
      k = graph->xadj[i] + RandomInRangeFast(graph->xadj[i+1]-graph->xadj[i]);
      SWAP(graph->adjncy[j], graph->adjncy[k], tmp);
      SWAP(graph->adjwgt[j], graph->adjwgt[k], tmp);
    }
  }
}

/*************************************************************************
* This function computes movement statistics for adaptive refinement
* schemes
**************************************************************************/
void ComputeMoveStatistics(CtrlType *ctrl, GraphType *graph, int *nmoved, int *maxin, int *maxout)
{
  int i, j, k, nvtxs;
  idxtype *vwgt, *where;
  idxtype *lpvtxs, *gpvtxs;

  nvtxs = graph->nvtxs;
  vwgt = graph->vwgt;
  where = graph->where;

  lpvtxs = idxsmalloc(ctrl->nparts, 0, "ComputeMoveStatistics: lpvtxs");
  gpvtxs = idxsmalloc(ctrl->nparts, 0, "ComputeMoveStatistics: gpvtxs");

  for (j=i=0; i<nvtxs; i++) {
    lpvtxs[where[i]]++;
    if (where[i] != ctrl->mype)
      j++;
  }
  MPI_Allreduce((void *)lpvtxs, (void *)gpvtxs, ctrl->nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);

  *nmoved = GlobalSESum(ctrl, j);
  *maxout = GlobalSEMax(ctrl, j);
  *maxin = GlobalSEMax(ctrl, gpvtxs[ctrl->mype]-(nvtxs-j));

  GKfree(&lpvtxs, &gpvtxs, -1);
}

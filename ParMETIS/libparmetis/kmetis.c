/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This is the entry point of ParMETIS_PartKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: kmetis.c 10391 2011-06-23 19:00:08Z karypis $
 *
 */

#include <parmetislib.h>

/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partitionioner. 
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
void ParMETIS_V3_PartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
              idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	      real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	      MPI_Comm *comm)
{
  idx_t h, i;
  idx_t nvtxs = -1, npes, mype;
  ctrl_t ctrl;
  graph_t *graph;
  real_t avg, maximb;
  idx_t seed, dbglvl = 0;
  idx_t iwgtflag, inumflag, incon, inparts, ioptions[10], moptions[METIS_NOPTIONS];
  real_t *itpwgts, iubvec[MAXNCON];

  gkMPI_Comm_size(*comm, &npes);
  gkMPI_Comm_rank(*comm, &mype);

  /* Deal with poor vertex distributions */
  ctrl.comm = *comm;
  if (GlobalSEMin(&ctrl, vtxdist[mype+1]-vtxdist[mype]) < 1) {
    if (mype == 0)
      printf("Error: Poor vertex distribution (processor with no vertices).\n");
    return;
  }


  /********************************/
  /* Try and take care bad inputs */
  /********************************/
  if (options != NULL && options[0] == 1)
    dbglvl = options[PMV3_OPTION_DBGLVL];

  CheckInputs(STATIC_PARTITION, npes, dbglvl, wgtflag, &iwgtflag, numflag, &inumflag, 
              ncon, &incon, nparts, &inparts, tpwgts, &itpwgts, ubvec, iubvec, NULL, 
              NULL, options, ioptions, part, comm);


  /**********************************/
  /* Take care the nparts == 1 case */
  /**********************************/
  if (inparts <= 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], 0, part); 
    *edgecut = 0;
    return;
  }

  /*******************************/
  /* Take care of npes == 1 case */
  /*******************************/
  if (npes == 1) {
    nvtxs = vtxdist[1] - vtxdist[0];
    METIS_SetDefaultOptions(moptions);
    moptions[METIS_OPTION_NUMBERING] = inumflag;

    METIS_PartGraphKway(&nvtxs, &incon, xadj, adjncy, vwgt, NULL, adjwgt, 
          &inparts, itpwgts, iubvec, moptions, edgecut, part);
 
    return;
  }


  if (inumflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  /*****************************/
  /* Set up control structures */
  /*****************************/
  if (ioptions[0] == 1) {
    dbglvl = ioptions[PMV3_OPTION_DBGLVL];
    seed   = ioptions[PMV3_OPTION_SEED];
  }
  else {
    dbglvl = GLOBAL_DBGLVL;
    seed   = GLOBAL_SEED;
  }
  SetUpCtrl(&ctrl, inparts, dbglvl, *comm);
  ctrl.CoarsenTo   = gk_min(vtxdist[npes]+1, 25*incon*gk_max(npes, inparts));
  ctrl.seed        = (seed == 0) ? mype : seed*mype;
  ctrl.sync        = GlobalSEMax(&ctrl, seed);
  ctrl.partType    = STATIC_PARTITION;
  ctrl.ps_relation = -1;
  ctrl.tpwgts      = itpwgts;
  rcopy(incon, iubvec, ctrl.ubvec);

  graph = SetUpGraph(&ctrl, incon, vtxdist, xadj, vwgt, adjncy, adjwgt, &iwgtflag);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  AllocateWSpace(&ctrl, graph);

  /*******************************************/
  /* Check for funny cases                   */
  /* 	-graph with no edges                 */
  /* 	-graph with self edges               */
  /* 	-graph with less than 20*npe nodes   */
  /*******************************************/
  if (vtxdist[npes] < SMALLGRAPH || 
      vtxdist[npes] < npes*20 || 
      GlobalSESum(&ctrl, graph->nedges) == 0) {
    IFSET(ctrl.dbglvl, DBG_INFO, 
        rprintf(&ctrl, "Partitioning a graph of size %"PRIDX" serially\n", vtxdist[npes]));
    PartitionSmallGraph(&ctrl, graph);
    rprintf(&ctrl, "Small graph\n");
  }
  else {
    /***********************/
    /* Partition the graph */
    /***********************/
    Global_Partition(&ctrl, graph);
  }
  ParallelReMapGraph(&ctrl, graph);

  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  icopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  /*******************/
  /* Print out stats */
  /*******************/
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));

  if (ctrl.dbglvl&DBG_INFO) {
    rprintf(&ctrl, "Final %"PRIDX"-way CUT: %6"PRIDX" \tBalance: ", inparts, graph->mincut);
    avg = 0.0;
    for (h=0; h<incon; h++) {
      maximb = 0.0;
      for (i=0; i<inparts; i++)
        maximb =gk_max(maximb, graph->gnpwgts[i*incon+h]/itpwgts[i*incon+h]);
      avg += maximb;
      rprintf(&ctrl, "%.3"PRREAL" ", maximb);
    }
    rprintf(&ctrl, "  avg: %.3"PRREAL"\n", avg/(real_t)incon);
  }

  gk_free((void **)&itpwgts, (void **)&graph->lnpwgts, (void **)&graph->gnpwgts, 
         (void **)&graph->nvwgt, LTERM);
  FreeInitialGraphAndRemap(graph, iwgtflag, 1);
  FreeCtrl(&ctrl);

  if (inumflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}



/*************************************************************************
* This function is the driver to the multi-constraint partitioning algorithm.
**************************************************************************/
void Global_Partition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ncon, nparts;
  real_t ftmp, ubavg, lbavg, lbvec[MAXNCON];
 
  ncon   = graph->ncon;
  nparts = ctrl->nparts;
  ubavg  = ravg(graph->ncon, ctrl->ubvec);

  CommSetup(ctrl, graph);

  if (ctrl->dbglvl&DBG_PROGRESS) {
    rprintf(ctrl, "[%6"PRIDX" %8"PRIDX" %5"PRIDX" %5"PRIDX"] [%"PRIDX"] [", graph->gnvtxs, GlobalSESum(ctrl, graph->nedges),
	    GlobalSEMin(ctrl, graph->nvtxs), GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo);
    for (i=0; i<ncon; i++)
      rprintf(ctrl, " %.3"PRREAL"", GlobalSEMinFloat(ctrl,graph->nvwgt[rargmin_strd(graph->nvtxs, graph->nvwgt+i, ncon)*ncon+i]));  
    rprintf(ctrl, "] [");
    for (i=0; i<ncon; i++)
      rprintf(ctrl, " %.3"PRREAL"", GlobalSEMaxFloat(ctrl, graph->nvwgt[rargmax_strd(graph->nvtxs, graph->nvwgt+i, ncon)*ncon+i]));  
    rprintf(ctrl, "]\n");
  }

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo ||
	(graph->finer != NULL &&
	graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {

    /* Done with coarsening. Find a partition */
    AllocateRefinementWorkSpace(ctrl, 2*graph->nedges);
    graph->where = imalloc(graph->nvtxs+graph->nrecv, "graph->where");

    InitPartition(ctrl, graph);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputePartitionParams(ctrl, graph);
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
          graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");
    }

    /* In case no coarsening took place */
    if (graph->finer == NULL) {
      ComputePartitionParams(ctrl, graph);
      KWayFM(ctrl, graph, NGR_PASSES);
    }
  }
  else {
    Match_Global(ctrl, graph);

    Global_Partition(ctrl, graph->coarser);

    ProjectPartition(ctrl, graph);

    ComputePartitionParams(ctrl, graph);

    if (graph->ncon > 1 && graph->level < 3) {
      for (i=0; i<ncon; i++) {
        ftmp = rsum(nparts, graph->gnpwgts+i, ncon);
        if (ftmp != 0.0)
          lbvec[i] = (real_t)(nparts) *
          graph->gnpwgts[rargmax_strd(nparts, graph->gnpwgts+i, ncon)*ncon+i]/ftmp;
        else
          lbvec[i] = 1.0;
      }
      lbavg = ravg(graph->ncon, lbvec);

      if (lbavg > ubavg + 0.035) {
        if (ctrl->dbglvl&DBG_PROGRESS) {
          ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
          rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
              graph->gnvtxs, graph->mincut);
          for (i=0; i<graph->ncon; i++) 
            rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
          rprintf(ctrl, " [b]\n");
	}

        KWayBalance(ctrl, graph, graph->ncon);
      }
    }

    KWayFM(ctrl, graph, NGR_PASSES);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
          graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");
    }

    if (graph->level != 0)
      gk_free((void **)&graph->lnpwgts, (void **)&graph->gnpwgts, LTERM);
  }

  return;
}



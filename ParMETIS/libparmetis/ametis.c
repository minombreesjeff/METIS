/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ametis.c
 *
 * This is the entry point of parallel difussive repartitioning routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: ametis.c 10391 2011-06-23 19:00:08Z karypis $
 *
 */

#include <parmetislib.h>


/***********************************************************************************
* This function is the entry point of the parallel multilevel local diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
void ParMETIS_V3_AdaptiveRepart(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy,
  idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag,
  idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, real_t *ipc2redist,
  idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm)
{
  idx_t h, i, npes, mype;
  ctrl_t ctrl;
  graph_t *graph;
  idx_t tewgt, tvsize, nmoved, maxin, maxout, vtx_factor;
  real_t gtewgt, gtvsize, avg, maximb;
  idx_t ps_relation, seed, dbglvl = 0;
  idx_t iwgtflag, inumflag, incon, inparts, ioptions[10];
  real_t iipc2redist, *itpwgts, iubvec[MAXNCON];

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
  CheckInputs(ADAPTIVE_PARTITION, npes, dbglvl, wgtflag, &iwgtflag, numflag, &inumflag,
              ncon, &incon, nparts, &inparts, tpwgts, &itpwgts, ubvec, iubvec, 
	      ipc2redist, &iipc2redist, options, ioptions, part, comm);

  /* ADD: take care of disconnected graph */

  /*********************************/
  /* Take care the nparts = 1 case */
  /*********************************/
  if (inparts == 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], 0, part); 
    *edgecut = 0;
    return;
  }

  /**************************/
  /* Set up data structures */
  /**************************/
  if (inumflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  /*****************************/
  /* Set up control structures */
  /*****************************/
  if (ioptions[0] == 1) {
    dbglvl      = ioptions[PMV3_OPTION_DBGLVL];
    seed        = ioptions[PMV3_OPTION_SEED];
    ps_relation = (npes == inparts ? ioptions[PMV3_OPTION_PSR] : PARMETIS_PSR_UNCOUPLED);
  }
  else {
    dbglvl      = GLOBAL_DBGLVL;
    seed        = GLOBAL_SEED;
    ps_relation = (npes == inparts ? PARMETIS_PSR_COUPLED : PARMETIS_PSR_UNCOUPLED);
  }

  SetUpCtrl(&ctrl, inparts, dbglvl, *comm);
  vtx_factor         = (gk_max(npes, inparts) > 256) ? 20 : 50;
  ctrl.CoarsenTo     = gk_min(vtxdist[npes]+1, vtx_factor*incon*gk_max(npes, inparts));
  ctrl.ipc_factor    = iipc2redist;
  ctrl.redist_factor = 1.0;
  ctrl.redist_base   = 1.0;
  ctrl.seed          = (seed == 0 ? mype : seed*mype);
  ctrl.sync          = GlobalSEMax(&ctrl, seed);
  ctrl.partType      = ADAPTIVE_PARTITION;
  ctrl.ps_relation   = ps_relation;
  ctrl.tpwgts        = itpwgts;

  graph = SetUpGraph(&ctrl, incon, vtxdist, xadj, vwgt, adjncy, adjwgt, &iwgtflag);
  graph->vsize = (vsize == NULL ? ismalloc(graph->nvtxs, 1, "vsize") : vsize);

  graph->home = imalloc(graph->nvtxs, "home");
  if (ctrl.ps_relation == PARMETIS_PSR_COUPLED)
    iset(graph->nvtxs, mype, graph->home);
  else {
    /* Downgrade the partition numbers if part[] has more partitions that nparts */
    for (i=0; i<graph->nvtxs; i++)
      part[i] = (part[i] >= ctrl.nparts ? 0 : part[i]);

    icopy(graph->nvtxs, part, graph->home);
  }

  tewgt   = isum(graph->nedges, graph->adjwgt, 1);
  tvsize  = isum(graph->nvtxs, graph->vsize, 1);
  gtewgt  = (real_t) GlobalSESum(&ctrl, tewgt) + 1.0/graph->gnvtxs;  /* The +1/graph->gnvtxs were added to remove any FPE */
  gtvsize = (real_t) GlobalSESum(&ctrl, tvsize) + 1.0/graph->gnvtxs;
  ctrl.edge_size_ratio = gtewgt/gtvsize;
  rcopy(incon, iubvec, ctrl.ubvec);

  AllocateWSpace(&ctrl, graph);

  /***********************/
  /* Partition and Remap */
  /***********************/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Adaptive_Partition(&ctrl, graph);
  ParallelReMapGraph(&ctrl, graph);

  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  icopy(graph->nvtxs, graph->where, part);
  if (edgecut != NULL)
    *edgecut = graph->mincut;

  /***********************/
  /* Take care of output */
  /***********************/
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));

  if (ctrl.dbglvl&DBG_INFO) {
    Mc_ComputeMoveStatistics(&ctrl, graph, &nmoved, &maxin, &maxout);
    rprintf(&ctrl, "Final %3"PRIDX"-way Cut: %6"PRIDX" \tBalance: ", inparts, graph->mincut);
    avg = 0.0;
    for (h=0; h<incon; h++) {
      maximb = 0.0;
      for (i=0; i<inparts; i++)
        maximb =gk_max(maximb, graph->gnpwgts[i*incon+h]/itpwgts[i*incon+h]);
      avg += maximb;
      rprintf(&ctrl, "%.3"PRREAL" ", maximb);
    }
    rprintf(&ctrl, "\nNMoved: %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", nmoved, maxin, maxout, maxin+maxout);
  }

  /*************************************/
  /* Free memory, renumber, and return */
  /*************************************/
  gk_free((void **)&graph->lnpwgts, &graph->gnpwgts, &graph->nvwgt, &graph->home, 
      &itpwgts, LTERM);
      
  FreeInitialGraphAndRemap(graph, iwgtflag, vsize == NULL);
  FreeCtrl(&ctrl);

  if (inumflag == 1)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  return;
}


/*************************************************************************
* This function is the driver for the adaptive refinement mode of ParMETIS
**************************************************************************/
void Adaptive_Partition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i;
  idx_t tewgt, tvsize;
  real_t gtewgt, gtvsize;
  real_t ubavg, lbavg, lbvec[MAXNCON];

  /************************************/
  /* Set up important data structures */
  /************************************/
  CommSetup(ctrl, graph);

  ubavg   = ravg(graph->ncon, ctrl->ubvec);
  tewgt   = isum(graph->nedges, graph->adjwgt, 1);
  tvsize  = isum(graph->nvtxs, graph->vsize, 1);
  gtewgt  = (real_t) GlobalSESum(ctrl, tewgt) + 1.0/graph->gnvtxs;  /* The +1/graph->gnvtxs were added to remove any FPE */
  gtvsize = (real_t) GlobalSESum(ctrl, tvsize) + 1.0/graph->gnvtxs;
  ctrl->redist_factor = ctrl->redist_base * ((gtewgt/gtvsize)/ ctrl->edge_size_ratio);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6"PRIDX" %8"PRIDX" %5"PRIDX" %5"PRIDX"][%"PRIDX"]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs), GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo));

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo ||
     (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {

    AllocateRefinementWorkSpace(ctrl, 2*graph->nedges);

    /***********************************************/
    /* Balance the partition on the coarsest graph */
    /***********************************************/
    graph->where = ismalloc(graph->nvtxs+graph->nrecv, -1, "graph->where");
    icopy(graph->nvtxs, graph->home, graph->where);

    ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
    lbavg = ravg(graph->ncon, lbvec);

    if (lbavg > ubavg + 0.035 && ctrl->partType != REFINE_PARTITION)
      Balance_Partition(ctrl, graph);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputePartitionParams(ctrl, graph);
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");
    }

    /* check if no coarsening took place */
    if (graph->finer == NULL) {
      ComputePartitionParams(ctrl, graph);
      KWayBalance(ctrl, graph, graph->ncon);
      KWayAdaptiveRefine(ctrl, graph, NGR_PASSES);
    }
  }
  else {
    /*******************************/
    /* Coarsen it and partition it */
    /*******************************/
    switch (ctrl->ps_relation) {
      case PARMETIS_PSR_COUPLED:
        Match_Local(ctrl, graph);
        break;
      case PARMETIS_PSR_UNCOUPLED:
      default:
        Match_Global(ctrl, graph);
        break;
    }

    Adaptive_Partition(ctrl, graph->coarser);

    /********************************/
    /* project partition and refine */
    /********************************/
    ProjectPartition(ctrl, graph);
    ComputePartitionParams(ctrl, graph);

    if (graph->ncon > 1 && graph->level < 4) {
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      lbavg = ravg(graph->ncon, lbvec);

      if (lbavg > ubavg + 0.025) {
        KWayBalance(ctrl, graph, graph->ncon);
      }
    }

    KWayAdaptiveRefine(ctrl, graph, NGR_PASSES);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");
    }
  }
}


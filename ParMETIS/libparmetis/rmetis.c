/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rmetis.c
 *
 * This is the entry point of the partitioning refinement routine
 *
 * Started 10/19/96
 * George
 *
 * $Id: rmetis.c 10391 2011-06-23 19:00:08Z karypis $
 *
 */

#include <parmetislib.h>



/***********************************************************************************
* This function is the entry point of the parallel multilevel local diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
void ParMETIS_V3_RefineKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy,
              idx_t *vwgt, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, 
	      idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, 
	      idx_t *part, MPI_Comm *comm)
{
  idx_t h, i;
  idx_t npes, mype;
  ctrl_t ctrl;
  graph_t *graph;
  idx_t tewgt, tvsize, nmoved, maxin, maxout;
  real_t gtewgt, gtvsize, avg, maximb;
  idx_t ps_relation, seed, dbglvl = 0;
  idx_t iwgtflag, inumflag, incon, inparts, ioptions[10];
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

  CheckInputs(REFINE_PARTITION, npes, dbglvl, wgtflag, &iwgtflag, numflag, &inumflag,
              ncon, &incon, nparts, &inparts, tpwgts, &itpwgts, ubvec, iubvec, 
              NULL, NULL, options, ioptions, part, comm);

  /* ADD: take care of disconnected graph */
  /* ADD: take care of highly unbalanced vtxdist */
  /*********************************/
  /* Take care the nparts = 1 case */
  /*********************************/
  if (inparts <= 1) {
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
    ps_relation = (npes == inparts) ? ioptions[PMV3_OPTION_PSR] : PARMETIS_PSR_UNCOUPLED;
  }
  else {
    dbglvl      = GLOBAL_DBGLVL;
    seed        = GLOBAL_SEED;
    ps_relation = (npes == inparts) ? PARMETIS_PSR_COUPLED : PARMETIS_PSR_UNCOUPLED;
  }

  SetUpCtrl(&ctrl, inparts, dbglvl, *comm);
  ctrl.CoarsenTo     = gk_min(vtxdist[npes]+1, 50*incon*gk_max(npes, inparts));
  ctrl.ipc_factor    = 1000.0;
  ctrl.redist_factor = 1.0;
  ctrl.redist_base   = 1.0;
  ctrl.seed          = (seed == 0) ? mype : seed*mype;
  ctrl.sync          = GlobalSEMax(&ctrl, seed);
  ctrl.partType      = REFINE_PARTITION;
  ctrl.ps_relation   = ps_relation;
  ctrl.tpwgts        = itpwgts;

  graph = SetUpGraph(&ctrl, incon, vtxdist, xadj, vwgt, adjncy, adjwgt, &iwgtflag);
  graph->vsize = ismalloc(graph->nvtxs, 1, "vsize");

  graph->home = imalloc(graph->nvtxs, "home");
  if (ctrl.ps_relation == PARMETIS_PSR_COUPLED)
    iset(graph->nvtxs, mype, graph->home);
  else
    icopy(graph->nvtxs, part, graph->home);

  tewgt   = isum(graph->nedges, graph->adjwgt, 1);
  tvsize  = isum(graph->nvtxs, graph->vsize, 1);
  gtewgt  = (real_t) GlobalSESum(&ctrl, tewgt) + 1.0/graph->gnvtxs;
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
      &graph->vsize, &itpwgts, LTERM);

  FreeInitialGraphAndRemap(graph, iwgtflag, 1);
  FreeCtrl(&ctrl);

  if (inumflag == 1)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

  return;
}



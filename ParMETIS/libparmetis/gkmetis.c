/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * gkmetis.c
 *
 * This is the entry point of parallel geometry based partitioning
 * routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: gkmetis.c 10385 2011-06-22 23:07:13Z karypis $
 *
 */

#include <parmetislib.h>




/***********************************************************************************
* This function is the entry point of the parallel kmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
void ParMETIS_V3_PartGeomKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy,
              idx_t *vwgt, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, 
	      real_t *xyz, idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, 
	      idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm)
{
  idx_t h, i, j, npes, mype;
  idx_t nvtxs = -1;
  idx_t uwgtflag, cut, gcut, maxnvtxs;
  idx_t ltvwgts[MAXNCON];
  idx_t moptions[METIS_NOPTIONS];
  ctrl_t ctrl;
  idx_t *uvwgt;
  graph_t *graph, *mgraph;
  real_t avg, maximb, balance;
  idx_t seed, dbglvl = 0;
  idx_t iwgtflag, inumflag, incon, inparts, ioptions[10];
  real_t *itpwgts, iubvec[MAXNCON];

  gkMPI_Comm_size(*comm, &npes);
  gkMPI_Comm_rank(*comm, &mype);

  /* If too many processors switch to non-geometric partitioning.
     This is to take care the current poor implementation of sorting
     that has an npes*npes memory complexity. The following fix assumes
     that the machine can allocate about 128MB of memory per node just 
     for sorting alone. 
     Also, if each processor does not have npes vertices, switch to the 
     non-geometric version of the code.
  */
  ctrl.comm = *comm;
  if (npes > 4096 || GlobalSEMin(&ctrl, vtxdist[mype+1]-vtxdist[mype]) < npes) {
    ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, 
         numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, part, comm);
    return;
  }


  /********************************/
  /* Try and take care bad inputs */
  /********************************/
  if (options != NULL && options[0] == 1)
    dbglvl = options[PMV3_OPTION_DBGLVL];

  CheckInputs(STATIC_PARTITION, npes, dbglvl, wgtflag, &iwgtflag, numflag, &inumflag,
              ncon, &incon, nparts, &inparts, tpwgts, &itpwgts, ubvec, iubvec, 
	      NULL, NULL, options, ioptions, part, comm);


  /*********************************/
  /* Take care the nparts = 1 case */
  /*********************************/
  if (inparts <= 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], 0, part);
    *edgecut = 0;
    return;
  }

  /******************************/
  /* Take care of npes = 1 case */
  /******************************/
  if (npes == 1 && inparts > 1) {
    nvtxs = vtxdist[1];

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
  SetUpCtrl(&ctrl, npes, dbglvl, *comm);
  ctrl.CoarsenTo   = gk_min(vtxdist[npes]+1, 25*incon*gk_max(npes, inparts));
  ctrl.seed        = (seed == 0) ? mype : seed*mype;
  ctrl.sync        = GlobalSEMax(&ctrl, seed);
  ctrl.partType    = STATIC_PARTITION;
  ctrl.ps_relation = -1;
  ctrl.tpwgts      = itpwgts;
  rcopy(incon, iubvec, ctrl.ubvec);

  uwgtflag = iwgtflag|2;
  uvwgt = ismalloc(vtxdist[mype+1]-vtxdist[mype], 1, "uvwgt");
  graph = SetUpGraph(&ctrl, 1, vtxdist, xadj, uvwgt, adjncy, adjwgt, &uwgtflag);
  gk_free((void **)&graph->nvwgt, &uvwgt, LTERM); 

  AllocateWSpace(&ctrl, graph);

  /*=================================================================
   * Compute the initial npes-way partitioning geometric partitioning
   =================================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, *ndims, xyz, 1);

  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  /*=================================================================
   * Move the graph according to the partitioning
   =================================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.MoveTmr));

  graph->vwgt = ((iwgtflag&2) != 0) ? vwgt : ismalloc(graph->nvtxs*incon, 1, "vwgt");
  graph->ncon = incon;
  j = ctrl.nparts;
  ctrl.nparts = ctrl.npes;
  mgraph = MoveGraph(&ctrl, graph);
  ctrl.nparts = j;

  /**********************************************************/
  /* Do the same functionality as SetUpGraph for mgraph */
  /**********************************************************/
  /* compute tvwgts */
  for (j=0; j<incon; j++)
    ltvwgts[j] = 0;

  for (i=0; i<graph->nvtxs; i++)
    for (j=0; j<incon; j++)
      ltvwgts[j] += mgraph->vwgt[i*incon+j];

  for (j=0; j<incon; j++) {
    ctrl.tvwgts[j] = GlobalSESum(&ctrl, ltvwgts[j]);
    ctrl.invtvwgts[j] = 1.0/ctrl.tvwgts[j];
  }

  /* check for zero wgt constraints */
  for (i=0; i<incon; i++) {
    /* ADD: take care of the case in which tvwgts is zero */
    if (ctrl.tvwgts[i] == 0) {
      if (ctrl.mype == 0) printf("ERROR: sum weight for constraint %"PRIDX" is zero\n", i);
      gkMPI_Finalize();
      exit(-1);
    }
  }

  /* compute nvwgt */
  mgraph->nvwgt = rmalloc(mgraph->nvtxs*incon, "mgraph->nvwgt");
  for (i=0; i<mgraph->nvtxs; i++) {
    for (j=0; j<incon; j++)
      mgraph->nvwgt[i*incon+j] = ctrl.invtvwgts[j]*mgraph->vwgt[i*incon+j];
  }


  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.MoveTmr));

  if (ctrl.dbglvl&DBG_INFO) {
    cut = 0;
    for (i=0; i<graph->nvtxs; i++)
      for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
        if (graph->where[i] != graph->where[graph->adjncy[j]])
          cut += graph->adjwgt[j];
    gcut = GlobalSESum(&ctrl, cut)/2;
    maxnvtxs = GlobalSEMax(&ctrl, mgraph->nvtxs);
    balance = (real_t)(maxnvtxs)/((real_t)(graph->gnvtxs)/(real_t)(npes));
    rprintf(&ctrl, "XYZ Cut: %6"PRIDX" \tBalance: %6.3"PRREAL" [%"PRIDX" %"PRIDX" %"PRIDX"]\n",
      gcut, balance, maxnvtxs, graph->gnvtxs, npes);
  }

  /*=================================================================
   * Set up the newly moved graph
   =================================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  ctrl.nparts = inparts;

  /*=======================================================
   * Now compute the partition of the moved graph
   =======================================================*/
  if (vtxdist[npes] < SMALLGRAPH || 
      vtxdist[npes] < npes*20 || 
      GlobalSESum(&ctrl, mgraph->nedges) == 0) {
    IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Partitioning a graph of size %"PRIDX" serially\n", vtxdist[npes]));
    PartitionSmallGraph(&ctrl, mgraph);
  }
  else 
    Global_Partition(&ctrl, mgraph);

  ParallelReMapGraph(&ctrl, mgraph);

  /* Invert the ordering back to the original graph */
  ctrl.nparts = npes;
  ProjectInfoBack(&ctrl, graph, part, mgraph->where);

  *edgecut = mgraph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  /*******************/
  /* Print out stats */
  /*******************/
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));

  if (ctrl.dbglvl&DBG_INFO) {
    rprintf(&ctrl, "Final %"PRIDX"-way CUT: %6"PRIDX" \tBalance: ", inparts, mgraph->mincut);
    avg = 0.0;
    for (h=0; h<incon; h++) {
      maximb = 0.0;
      for (i=0; i<inparts; i++)
        maximb =gk_max(maximb, mgraph->gnpwgts[i*incon+h]/itpwgts[i*incon+h]);
      avg += maximb;
      rprintf(&ctrl, "%.3"PRREAL" ", maximb);
    }
    rprintf(&ctrl, "  avg: %.3"PRREAL"\n", avg/(real_t)incon);
  }

  gk_free((void **)&itpwgts, LTERM);
  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph, iwgtflag, 1);
  FreeCtrl(&ctrl);

  if (inumflag == 1)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}



/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void ParMETIS_V3_PartGeom(idx_t *vtxdist, idx_t *ndims, real_t *xyz, idx_t *part, MPI_Comm *comm)
{
  idx_t i, nvtxs, firstvtx, dbglvl, npes, mype;
  idx_t *xadj, *adjncy;
  ctrl_t ctrl;
  graph_t *graph;
  idx_t zeroflg = 0;

  gkMPI_Comm_size(*comm, &npes);
  gkMPI_Comm_rank(*comm, &mype);

  if (npes == 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], 0, part);
    return;
  }

  /* Return without computing a partitioning under the following cases: 
     - The number of processors is greater than 4096 (due to npes^2
       memory complexity of the sorting algorithm implemented).
     - When each processor does not have at least npes elements.
     These retrictions will be fixed in 4.0.
  */
  ctrl.comm = *comm;
  if (npes > 4096 || GlobalSEMin(&ctrl, vtxdist[mype+1]-vtxdist[mype]) < npes) {
    if (mype == 1)
      printf("ParMETIS_V3_PartGeom can only be used for less than 4096 processors "
             "and when each processor has at least npes elements.\n");
    return;
  }


  /* Setup a fake graph to allow the rest of the code to work unchanged */
  dbglvl = 0;

  nvtxs = vtxdist[mype+1]-vtxdist[mype];
  firstvtx = vtxdist[mype];
  xadj = imalloc(nvtxs+1, "ParMETIS_PartGeom: xadj");
  adjncy = imalloc(nvtxs, "ParMETIS_PartGeom: adjncy");
  for (i=0; i<nvtxs; i++) {
    xadj[i] = i;
    adjncy[i] = firstvtx + (i+1)%nvtxs;
  }
  xadj[nvtxs] = nvtxs;

  /* Proceed with the rest of the code */
  SetUpCtrl(&ctrl, npes, dbglvl, *comm);
  ctrl.seed      = mype;
  ctrl.CoarsenTo = gk_min(vtxdist[npes]+1, 25*npes);

  graph = SetUpGraph(&ctrl, 1, vtxdist, xadj, NULL, adjncy, NULL, &zeroflg);

  AllocateWSpace(&ctrl, graph);

  /*=======================================================
   * Compute the initial geometric partitioning
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, *ndims, xyz, 0);

  icopy(graph->nvtxs, graph->where, part);

  IFSET(ctrl.dbglvl, DBG_TIME, gkMPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  FreeInitialGraphAndRemap(graph, 0, 1);
  FreeCtrl(&ctrl);

  gk_free((void **)&xadj, (void **)&adjncy, LTERM);
}





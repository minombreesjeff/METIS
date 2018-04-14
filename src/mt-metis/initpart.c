/**
 * @file initpart.c
 * @brief Parallel initial partitioning routines (see pmetis and kmetis)
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2012, Regents of the University of Minnesota
 * @version 1
 * @date 2012-07-06
 */


#include "includes.h"

idx_t ParInitKWayPartitioning(dctrl_t * const dctrl, dgraph_t * const graph)
{
  INIT_PARALLEL();

  idx_t i,j,idx, cut;
  idx_t *cuts, *cutidx;

  const idx_t nvtxs = graph->nvtxs; 
  const idx_t nedges = graph->nedges;
  const idx_t nparts = dctrl->ctrl->nparts;
  const idx_t ndist = graph->ndist;

  ParAllocateKWayPartitionMemory(dctrl->ctrl, graph);

  const idx_t tcuts = gk_max(NSOLUTIONS,nthreads/TPPRATIO);

  idx_t myncuts = (tcuts/ nthreads);
  idx_t a = tcuts % nthreads, b = nthreads, c = a;

  /* factor */
  while (c > 0) {
    if (a % c == 0 && b % c == 0) {
      a /= c;
      b /= c;
    }
    --c;
  }
  myncuts += (myid % b < a) ? 1 : 0;

  idx_t *where;
  startwctimer(dctrl->convertTmr);
  graph_t * mgraph = 
      ParConvertDGraph(dctrl,graph,(graph_t**)&dctrl->tptr_void[0]);

  ASSERT(mgraph->nvtxs == nvtxs);
  ASSERT(mgraph->nedges == nedges);

  stopwctimer(dctrl->convertTmr);

  if (myncuts > 0) {
    startwctimer(dctrl->ipTmr);

    idx_t options[METIS_NOPTIONS];

    where = dctrl->tptr_idx[myid] = imalloc(nvtxs,"X");

    METIS_SetDefaultOptions(options);

    options[METIS_OPTION_NITER] = 10;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_SEED] = dctrl->ctrl->seed + myid;
    options[METIS_OPTION_NCUTS] = myncuts;

    ASSERT(options[METIS_OPTION_NCUTS] > 0);

    ctrl_t * ctrl = SetupCtrl(METIS_OP_PMETIS, options, 1, nparts, 
        dctrl->ctrl->tpwgts, NULL);

    idx_t * myxadj, *myadjncy,*myvwgt,*myadjwgt;

    #ifdef LGCOPY
    myxadj = imalloc(nvtxs+1,"X");
    myvwgt = imalloc(nvtxs,"X");
    myadjncy = imalloc(mgraph->nedges,"X");
    myadjwgt = imalloc(mgraph->nedges,"X");

    startwctimer(dctrl->lCopyTmr);
    icopy(nvtxs+1,mgraph->xadj,myxadj);
    icopy(nvtxs,mgraph->vwgt,myvwgt);
    icopy(mgraph->nedges,mgraph->adjncy,myadjncy);
    icopy(mgraph->nedges,mgraph->adjwgt,myadjwgt);
    stopwctimer(dctrl->lCopyTmr);
    #else
    myxadj = mgraph->xadj;
    myvwgt = mgraph->vwgt;
    myadjncy = mgraph->adjncy;
    myadjwgt = mgraph->adjwgt;
    #endif

    graph_t * mygraph = SetupGraph(ctrl,nvtxs,1,myxadj,
        myadjncy,myvwgt,NULL,myadjwgt);

    ASSERT(mgraph->nvtxs == mygraph->nvtxs);
    ASSERT(mgraph->nedges == mygraph->nedges);

    ctrl->ubfactors[0] = 
      (real_t)pow(dctrl->ctrl->ubfactors[0], 1.0/log(ctrl->nparts));

    AllocateWorkSpace(ctrl,mygraph);


    startwctimer(dctrl->parBisTmr);
    cut = MlevelRecursiveBisection(ctrl,mygraph,ctrl->nparts,where,
        ctrl->tpwgts,0);
    stopwctimer(dctrl->parBisTmr);
    
    #ifdef LGCOPY
    gk_free((void**)&myxadj,&myvwgt,&myadjncy,&myadjwgt,LTERM);
    #endif

    FreeWorkSpace(ctrl);
    FreeCtrl(&ctrl);
  } else {
    cut = 1<<((sizeof(idx_t)*8)-2);
    ASSERT(cut > 0);
    where = dctrl->tptr_idx[myid] = NULL;
  }

  dl_omp_minreduce_index(myid,idx,cut,dctrl->buffer1,nthreads);

  if (myid == idx) {
    ASSERT(where != NULL);
    graph->mincut = cut;
  }

  startwctimer(dctrl->renParTmr);

  idx_t mynvtxs;
  mynvtxs = graph->mynvtxs[myid];
  for (i=0;i<mynvtxs;++i) {
    graph->where[myid][i] = dctrl->tptr_idx[idx][graph->rename[myid][i]];
  }
  #pragma omp barrier

  stopwctimer(dctrl->renParTmr);

  gk_free((void**)&graph->rename[myid],&where,LTERM);

  stopwctimer(dctrl->ipTmr);

  #pragma omp master
  {
    FreeGraph(&mgraph);
  }

  ASSERT(ParComputeCut(dctrl,graph,(const idx_t **)graph->where) 
      == graph->mincut);

  return graph->mincut;
}


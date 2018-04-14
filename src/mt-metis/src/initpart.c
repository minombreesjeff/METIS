/**
 * @file initpart.c
 * @brief Parallel initial partitioning routines (see pmetis and kmetis)
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2012-07-06
 */




#ifndef MTMETIS_INITPART_C
#define MTMETIS_INITPART_C




#include "initpart.h"

#undef real_t
#include <metis.h>
#define real_t mtmetis_real_t




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static wgt_t __initpart_metis_kway(
    ctrl_t * const ctrl,
    size_t const ncuts,
    vtx_t nvtxs,
    adj_t * const xadj,
    vtx_t * const adjncy,
    wgt_t * const vwgt,
    wgt_t * const adjwgt,
    pid_t * const where)
{
  idx_t nparts;
  idx_t cut, ncon;
  idx_t options[METIS_NOPTIONS];
  real_t ubf;

  tid_t const myid = omp_get_thread_num();

  METIS_SetDefaultOptions(options);

  ncon = 1;

  options[METIS_OPTION_NITER] = 10;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_SEED] = ctrl->seed + myid;
  options[METIS_OPTION_NCUTS] = ncuts;
  options[METIS_OPTION_DBGLVL] = 0;
  nparts = ctrl->nparts;
  ubf = ctrl->ubfactor;

  #ifdef KWAY_INIT
  METIS_PartGraphKway((idx_t*)&nvtxs,&ncon,(idx_t*)xadj,(idx_t*)adjncy, \
      (idx_t*)vwgt,NULL,(idx_t*)adjwgt,&nparts,NULL,&ubf,options,&cut, \
      (idx_t*)where);
  #else
  METIS_PartGraphRecursive((idx_t*)&nvtxs,&ncon,(idx_t*)xadj,(idx_t*)adjncy, \
      (idx_t*)vwgt,NULL,(idx_t*)adjwgt,&nparts,NULL,&ubf,options,&cut, \
      (idx_t*)where);
  #endif

  return (wgt_t)cut;
}



/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/

static pid_t * __ik_where;
vtx_t initpart_kway(
    ctrl_t * const ctrl,
    graph_t * const graph)
{
  vtx_t voff, idx;
  wgt_t cut;
  adj_t * xadj;
  vtx_t * adjncy;
  wgt_t * adjwgt, * vwgt;
  pid_t * where = NULL;

  tid_t const nthreads = omp_get_num_threads();
  tid_t const myid = omp_get_thread_num();

  vtx_t const nvtxs = graph->nvtxs; 

  size_t const tcuts = dl_max(NSOLUTIONS,nthreads/TPPRATIO);

  size_t myncuts = (tcuts / nthreads);

  size_t a = tcuts % nthreads, b = nthreads, c = a;

  #pragma omp master
  {
    dl_start_timer(&ctrl->timers.initpart);
  }

  graph_gather(graph,&xadj,&adjncy,&vwgt,&adjwgt,&voff);

  /* factor */
  while (c > 0) {
    if (a % c == 0 && b % c == 0) {
      a /= c;
      b /= c;
    }
    --c;
  }
  myncuts += (myid % b < a) ? 1 : 0;

  if (myncuts > 0) {
    where = pid_alloc(nvtxs);

    cut = __initpart_metis_kway(ctrl,myncuts,nvtxs,xadj,adjncy,vwgt,adjwgt, \
        where);
  } else {
    cut = graph->tadjwgt+1;
  }

  idx = wgt_omp_minreduce_index(cut);

  if (myid == idx) {
    DL_ASSERT(where != NULL,"Non-cutting thread chosen");
    graph->mincut = cut;
    __ik_where = where;
  }
  #pragma omp barrier

  graph_alloc_partmemory(ctrl,graph);

  par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_MEDIUM,"Selected initial " \
      "partition with cut of %"PF_WGT_T"\n",graph->mincut);

  /* save the best where */
  pid_copy(graph->where[myid],__ik_where+voff,graph->mynvtxs[myid]);

  #pragma omp barrier
  #pragma omp master
  {
    /* free the gathered graph */
    dl_free(xadj);
    dl_free(adjncy);
    dl_free(vwgt);
    dl_free(adjwgt);
  }

  if (where) {
    dl_free(where);
  }

  #pragma omp master
  {
    dl_stop_timer(&ctrl->timers.initpart);
  }

  return graph->mincut;
}




#endif

/**
 * @file mtmetis.c
 * @brief Library entry points
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-01
 */




#ifndef MTMETIS_C
#define MTMETIS_C




#include "base.h"
#include "ctrl.h"
#include "graph.h"
#include "partition.h"




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __partition_unify_where(
    pid_t * const * const where, 
    graph_t const * const graph, 
    pid_t * const uwhere)
{
  vtx_t i;

  tid_t const myid = omp_get_thread_num();

  for (i=0;i<graph->mynvtxs[myid];++i) {
    uwhere[graph->label[myid][i]] = where[myid][i];
  }
}



/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


/**
 * @brief Innermost function for spawning parallelism.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * @param where The output partition IDs for each vertex.
 */
void mtmetis_partition_internal(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * where);




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void mtmetis_partition_internal(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const where)
{
  pid_t ** dwhere;
  tid_t const nthreads = ctrl->nthreads;

  dwhere = r_vtx_alloc(ctrl->nthreads);

  #pragma omp parallel shared(dwhere) num_threads(nthreads) \
    default(none)
  {
    tid_t const myid = omp_get_thread_num();

    vtx_omp_reduction_init(nthreads);
    wgt_omp_reduction_init(nthreads);
    adj_omp_reduction_init(nthreads);
    double_omp_reduction_init(nthreads);

    dwhere[myid] = pid_alloc(graph->mynvtxs[myid]);

    #pragma omp master
    {
      dl_start_timer(&ctrl->timers.partitioning);
    }

    partition_kway(ctrl,graph,dwhere);

    #pragma omp master
    {
      dl_stop_timer(&ctrl->timers.partitioning);
    }

    if (where) {
      __partition_unify_where(dwhere,graph,where);
    }

    if (ctrl->verbosity >= MTMETIS_VERBOSITY_LOW) {
      #pragma omp barrier
      #pragma omp master
      {
        partition_print_info(ctrl,graph,(pid_t const **)dwhere);
      }
      #pragma omp barrier
    }

    dl_free(dwhere[myid]);

    vtx_omp_reduction_clear();
    adj_omp_reduction_clear();
    wgt_omp_reduction_clear();
    double_omp_reduction_clear();
  }

  dl_free(dwhere);
}


double * mtmetis_init_options(void)
{
  double * options = double_init_alloc(MTMETIS_VAL_OFF,MTMETIS_NOPTIONS);
  return options; 
}


int mtmetis_partition_explicit(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const vwgt,
    wgt_t const * const adjwgt,
    double const * const options,
    pid_t * const where,
    wgt_t * const r_edgecut)
{
  int rv;
  ctrl_t * ctrl = NULL;
  graph_t * graph = NULL;
  vtx_t ** dwhere = NULL;
  
  if ((rv = ctrl_parse(options,&ctrl)) != MTMETIS_SUCCESS) {
    goto CLEANUP;
  }
  
  ctrl_setup(ctrl,nvtxs);

  graph = graph_distribute(MTMETIS_DISTRIBUTION_BLOCKCYCLIC,nvtxs,xadj, \
      adjncy,vwgt,adjwgt,ctrl->nthreads);

  mtmetis_partition_internal(ctrl,graph,where);

  if (r_edgecut) {
    *r_edgecut = graph->mincut;
  }

  CLEANUP:

  if (graph) {
    graph_free(graph);
  }
  if (ctrl) {
    ctrl_free(ctrl);
  }
  if (dwhere) {
    r_pid_free(dwhere,ctrl->nthreads);
  }

  return rv;
}


int mtmetis_partkway(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const vwgt,
    wgt_t const * const adjwgt,
    pid_t const nparts,
    pid_t * const where,
    wgt_t * const r_edgecut)
{
  int rv;
  tid_t nthreads;
  double options[MTMETIS_NOPTIONS];
  
  if (omp_in_parallel()) {
    nthreads = 1;
  } else {
    nthreads = omp_get_max_threads();
  }

  options[MTMETIS_OPTION_NTHREADS] = (double)nthreads;
  options[MTMETIS_OPTION_NPARTS] = (double)nparts;

  rv = mtmetis_partition_explicit(nvtxs,xadj,adjncy,vwgt,adjwgt,options, \
      where,r_edgecut);

  return rv;
}




#endif

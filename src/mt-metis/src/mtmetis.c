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
#include "order.h"
#include "internal.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


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
  pid_t ** dwhere = NULL;
  
  if ((rv = ctrl_parse(options,&ctrl)) != MTMETIS_SUCCESS) {
    goto CLEANUP;
  }
  
  ctrl_setup(ctrl,NULL,nvtxs);

  graph = graph_distribute(ctrl->dist,nvtxs,xadj, \
      adjncy,vwgt,adjwgt,ctrl->nthreads);

  switch (ctrl->ptype) {
    case MTMETIS_PTYPE_RB:
      mtmetis_partition_rb_int(ctrl,graph,where);
      break;
    case MTMETIS_PTYPE_KWAY:
      mtmetis_partition_kway_int(ctrl,graph,where);
      break;
    case MTMETIS_PTYPE_ESEP:
      mtmetis_partition_esep_int(ctrl,graph,where);
      break;
    case MTMETIS_PTYPE_VSEP:
      mtmetis_partition_vsep_int(ctrl,graph,where);
      break;
    case MTMETIS_PTYPE_ND:
      order_nd(ctrl,graph,where);
      break;
    default:
      dl_error("Unknown ptype '%d'",ctrl->ptype);
  }

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

  double_set(options,MTMETIS_VAL_OFF,MTMETIS_NOPTIONS);
  options[MTMETIS_OPTION_NTHREADS] = (double)nthreads;
  options[MTMETIS_OPTION_NPARTS] = (double)nparts;

  rv = mtmetis_partition_explicit(nvtxs,xadj,adjncy,vwgt,adjwgt,options, \
      where,r_edgecut);

  return rv;
}



int mtmetis_nd(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const vwgt,
    wgt_t const * const adjwgt,
    pid_t * const perm)
{
  int rv;
  tid_t nthreads;
  ctrl_t * ctrl = NULL;
  graph_t * graph = NULL;
  pid_t ** dperm = NULL;

  double options[MTMETIS_NOPTIONS];
  
  if (omp_in_parallel()) {
    nthreads = 1;
  } else {
    nthreads = omp_get_max_threads();
  }

  double_set(options,MTMETIS_VAL_OFF,MTMETIS_NOPTIONS);
  options[MTMETIS_OPTION_NTHREADS] = (double)nthreads;
  options[MTMETIS_OPTION_PTYPE] = MTMETIS_PTYPE_ND;
  options[MTMETIS_OPTION_NPARTS] = 3;

  if ((rv = ctrl_parse(options,&ctrl)) != MTMETIS_SUCCESS) {
    goto CLEANUP;
  }
  
  ctrl_setup(ctrl,NULL,nvtxs);

  graph = graph_distribute(ctrl->dist,nvtxs,xadj,adjncy,vwgt,adjwgt,nthreads);

  order_nd(ctrl,graph,perm);

  CLEANUP:

  if (dperm) {
    r_pid_free(dperm,nthreads);
  }
  if (graph) {
    graph_free(graph);
  }
  if (ctrl) {
    ctrl_free(ctrl);
  }

  return MTMETIS_SUCCESS;
}




#endif

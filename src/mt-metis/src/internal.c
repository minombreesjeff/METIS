/**
 * @file internal.c
 * @brief Internal mtmetis toplevel functions. 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 * @date 2015-01-17
 */




#ifndef MTMETIS_INTERNAL_C
#define MTMETIS_INTERNAL_C




#include "internal.h"
#include "partition.h"
#include "order.h"
#include "check.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct arg_t {
  ctrl_t * ctrl;
  graph_t * graph;
  pid_t * where;
  pid_t ** dwhere;
} arg_t;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __unify_where(
    pid_t * const * const where, 
    graph_t const * const graph, 
    pid_t * const uwhere)
{
  vtx_t i;

  tid_t const myid = dlthread_get_id(graph->comm);

  for (i=0;i<graph->mynvtxs[myid];++i) {
    uwhere[graph->label[myid][i]] = where[myid][i];
  }
}


static void __kway_func(
    void * const ptr)
{
  tid_t myid;
  arg_t * arg;
  ctrl_t * ctrl;
  graph_t * graph;
  pid_t * where;
  pid_t ** dwhere;

  arg = ptr;

  ctrl = arg->ctrl;
  graph = arg->graph;
  where = arg->where;
  dwhere = arg->dwhere;

  myid = dlthread_get_id(ctrl->comm);

  dlthread_barrier(ctrl->comm);

  dwhere[myid] = pid_alloc(graph->mynvtxs[myid]);

  par_partition_kway(ctrl,graph,dwhere);

  if (where) {
    __unify_where(dwhere,graph,where);
  }

  if (ctrl->verbosity >= MTMETIS_VERBOSITY_LOW) {
    dlthread_barrier(ctrl->comm);
    if (myid == 0) {
      partition_print_info(ctrl,graph,(pid_t const **)dwhere);
    }
    dlthread_barrier(ctrl->comm);
  }

  dl_free(dwhere[myid]);
}


static void __rb_func(
    void * const ptr)
{
  tid_t myid;
  arg_t * arg;
  ctrl_t * ctrl;
  graph_t * graph;
  pid_t * where;
  pid_t ** dwhere;

  arg = ptr;

  ctrl = arg->ctrl;
  graph = arg->graph;
  where = arg->where;
  dwhere = arg->dwhere;

  myid = dlthread_get_id(ctrl->comm);

  dlthread_barrier(ctrl->comm);

  dwhere[myid] = pid_alloc(graph->mynvtxs[myid]);

  par_partition_rb(ctrl,graph,dwhere);

  if (where) {
    __unify_where(dwhere,graph,where);
  }

  if (ctrl->verbosity >= MTMETIS_VERBOSITY_LOW) {
    dlthread_barrier(ctrl->comm);
    if (myid == 0) {
      partition_print_info(ctrl,graph,(pid_t const **)dwhere);
    }
    dlthread_barrier(ctrl->comm);
  }

  dl_free(dwhere[myid]);
}


static void __esep_func(
    void * const ptr)
{
  tid_t myid;
  arg_t * arg;
  ctrl_t * ctrl;
  graph_t * graph;
  pid_t * where;
  pid_t ** dwhere;

  arg = ptr;

  ctrl = arg->ctrl;
  graph = arg->graph;
  where = arg->where;
  dwhere = arg->dwhere;

  myid = dlthread_get_id(ctrl->comm);

  dlthread_barrier(ctrl->comm);

  dwhere[myid] = pid_alloc(graph->mynvtxs[myid]);

  par_partition_edgeseparator(ctrl,graph,dwhere);

  if (where) {
    __unify_where(dwhere,graph,where);
  }

  if (ctrl->verbosity >= MTMETIS_VERBOSITY_LOW) {
    dlthread_barrier(ctrl->comm);
    if (myid == 0) {
      partition_print_info(ctrl,graph,(pid_t const **)dwhere);
    }
    dlthread_barrier(ctrl->comm);
  }

  dl_free(dwhere[myid]);
}


static void __vsep_func(
    void * const ptr)
{
  tid_t myid, nthreads;
  arg_t * arg;
  ctrl_t * ctrl;
  graph_t * graph;
  pid_t * where;
  pid_t ** dwhere;

  arg = ptr;

  ctrl = arg->ctrl;
  graph = arg->graph;
  where = arg->where;
  dwhere = arg->dwhere;

  myid = dlthread_get_id(ctrl->comm);
  nthreads = dlthread_get_nthreads(ctrl->comm);

  dlthread_barrier(ctrl->comm);

  if (nthreads > 1) {
    par_partition_pre(ctrl,graph);
  }

  dwhere[myid] = pid_alloc(graph->mynvtxs[myid]);

  par_partition_vertexseparator(ctrl,graph,dwhere);

  if (where) {
    __unify_where(dwhere,graph,where);
  }

  if (ctrl->verbosity >= MTMETIS_VERBOSITY_LOW) {
    dlthread_barrier(ctrl->comm);
    if (myid == 0) {
      partition_print_info(ctrl,graph,(pid_t const **)dwhere);
    }
    dlthread_barrier(ctrl->comm);
  }

  dl_free(dwhere[myid]);
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void mtmetis_partition_kway_int(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const where)
{
  arg_t arg;
  pid_t ** dwhere;
  tid_t const nthreads = ctrl->nthreads;

  dwhere = r_pid_alloc(ctrl->nthreads);

  ctrl->comm = DLTHREAD_COMM_ROOT;
  graph->comm = DLTHREAD_COMM_ROOT;

  arg.ctrl = ctrl;
  arg.graph = graph;
  arg.where = where;
  arg.dwhere = dwhere;

  dlthread_launch(nthreads,&__kway_func,&arg);

  dl_free(dwhere);
}


void mtmetis_partition_rb_int(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const where)
{
  arg_t arg;
  pid_t ** dwhere;
  tid_t const nthreads = ctrl->nthreads;

  dwhere = r_pid_alloc(ctrl->nthreads);

  ctrl->comm = DLTHREAD_COMM_ROOT;
  graph->comm = DLTHREAD_COMM_ROOT;

  arg.ctrl = ctrl;
  arg.graph = graph;
  arg.where = where;
  arg.dwhere = dwhere;

  dlthread_launch(nthreads,&__rb_func,&arg);

  dl_free(dwhere);
}


void mtmetis_partition_esep_int(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const where)
{
  arg_t arg;
  pid_t ** dwhere;

  tid_t const nthreads = ctrl->nthreads;

  dwhere = r_pid_alloc(ctrl->nthreads);

  graph->comm = DLTHREAD_COMM_ROOT;
  ctrl->comm = DLTHREAD_COMM_ROOT;

  arg.ctrl = ctrl;
  arg.graph = graph;
  arg.where = where;
  arg.dwhere = dwhere;

  dlthread_launch(nthreads,&__esep_func,&arg);

  dl_free(dwhere);
}


void mtmetis_partition_vsep_int(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const where)
{
  arg_t arg;
  pid_t ** dwhere;

  tid_t const nthreads = ctrl->nthreads;

  dwhere = r_pid_alloc(ctrl->nthreads);

  graph->comm = DLTHREAD_COMM_ROOT;
  ctrl->comm = DLTHREAD_COMM_ROOT;

  arg.ctrl = ctrl;
  arg.graph = graph;
  arg.where = where;
  arg.dwhere = dwhere;

  dlthread_launch(nthreads,&__vsep_func,&arg);

  dl_free(dwhere);
}




#endif

/**
 * @file order.c
 * @brief Functions for creating a partition induced ordering (such as nested
 * dissection).
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-10-13
 */




#ifndef MTMETIS_ORDER_C
#define MTMETIS_ORDER_C




#include "order.h"
#include "partition.h"

#undef real_t
#include <metis.h>
#define real_t mtmetis_real_t




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct arg_t {
  ctrl_t * ctrl;
  graph_t * graph;
  pid_t * perm;
  pid_t ** dperm;
} arg_t;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __unify_perm(
    pid_t const * const * const dperm,
    graph_t * const graph,
    pid_t * const perm)
{
  vtx_t i;

  tid_t const myid = dlthread_get_id(graph->comm);

  for (i=0;i<graph->mynvtxs[myid];++i) {
    perm[graph->label[myid][i]] = dperm[myid][i];
  }
}


static void __order_nd(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t const offset,
    pid_t ** const perm) 
{
  vtx_t i;
  pid_t * fperm;
  idx_t options[METIS_NOPTIONS];

  tid_t const myid = dlthread_get_id(ctrl->comm);

  DL_ASSERT_EQUALS(graph->dist.nthreads,1,"%"PF_TID_T);

  fperm = pid_alloc(graph->mynvtxs[0]);

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.metis));
  }

  dlthread_pool_init(omp_get_num_threads());

  __METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NITER] = ctrl->nrefpass;
  options[METIS_OPTION_DBGLVL] = 0;
  options[METIS_OPTION_SEED] = ctrl->seed + myid;
  options[METIS_OPTION_NSEPS] = ctrl->ncuts;
  options[METIS_OPTION_UFACTOR] = 1000*(ctrl->ubfactor - 1.0);

  __METIS_NodeND((idx_t*)&graph->mynvtxs[0],(idx_t*)graph->xadj[0], \
      (idx_t*)graph->adjncy[0],(idx_t*)graph->vwgt[0],options,(idx_t*)fperm, \
      (idx_t*)perm[0]);
  
  dlthread_pool_finalize();

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.metis));
  }

  dl_free(fperm);

  /* add my offset */
  for (i=0;i<graph->mynvtxs[0];++i) {
    perm[0][i] += offset;
  }
}


static void __par_order_nd(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t const offset,
    pid_t ** const perm) 
{
  vtx_t v, g, mynvtxs, lvtx, i, nsep, lastvtx;
  tid_t hmyid, mygroup, lid;
  dlthread_comm_t lcomm;
  pid_t hoff;
  vtx_t * prefix, * sep;
  ctrl_t * myctrl;

  graph_t ** hgraphs;
  pid_t *** hperm;

  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);
  tid_t const myid = dlthread_get_id(ctrl->comm);

  /* handle the serial case */
  if (nthreads == 1) {
    __order_nd(ctrl,graph,offset,perm);
    return;
  }

  /* initial bisection */
  par_partition_vertexseparator(ctrl,graph,perm);

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.recursion));
  }

  hgraphs = dlthread_get_shmem((sizeof(graph_t*)*2) + (sizeof(pid_t**)*2) + \
      (sizeof(vtx_t)*(nthreads+1)),ctrl->comm);

  hperm = (pid_t***)(hgraphs+2);
  prefix = (vtx_t*)(hperm+2);

  sep = vtx_alloc(graph->mynvtxs[myid]);
  nsep = 0;
  for (i=0;i<graph->mynvtxs[myid];++i) {
    if (perm[myid][i] == MTMETIS_VSEP_SEP) {
      sep[nsep++] = i;
    }
  }
  prefix[myid] = nsep;
  dlthread_barrier(graph->comm);
  
  if (myid == 0) {
    vtx_prefixsum_exc(prefix,nthreads);
  }

  /* extract subgraphs and structure based on number of calling threads */
  mygroup = par_graph_extract_halves(graph,(pid_t const **)perm,hgraphs);

  if (mygroup == 0) {
    hoff = offset;
  } else {
    hoff = hgraphs[0]->nvtxs + offset;
  }
  dlthread_barrier(ctrl->comm);

  /* order my portion of the separator */
  lastvtx = graph->nvtxs + offset - prefix[myid];
  for (i=0;i<nsep;++i) {
    perm[myid][sep[i]] = --lastvtx;
  }
  dl_free(sep);

  lcomm = hgraphs[mygroup]->comm;

  hmyid = dlthread_get_id(lcomm);

  myctrl = par_ctrl_split(ctrl,hgraphs[mygroup]->nvtxs,MTMETIS_VSEP_NPARTS, \
      lcomm);

  DL_ASSERT_EQUALS((size_t)myctrl->comm,(size_t)hgraphs[mygroup]->comm,"%zu");

  hperm[mygroup] = dlthread_get_shmem(sizeof(pid_t*)*nthreads,lcomm);
  hperm[mygroup][hmyid] = pid_alloc(hgraphs[mygroup]->mynvtxs[hmyid]);

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.recursion));
  }

  __par_order_nd(myctrl,hgraphs[mygroup],hoff,hperm[mygroup]);

  if (myid == 0) {
    ctrl_combine_timers(ctrl,myctrl);
  }

  par_ctrl_free(myctrl);

  /* project my newly ordered vertices */
  mynvtxs = hgraphs[mygroup]->mynvtxs[hmyid];
  for (v=0;v<mynvtxs;++v) {
    g = hgraphs[mygroup]->label[hmyid][v];
    lid = gvtx_to_tid(g,graph->dist);
    lvtx = gvtx_to_lvtx(g,graph->dist);
    perm[lid][lvtx] = hperm[mygroup][hmyid][v];
  }

  dl_free(hperm[mygroup][hmyid]);
  dlthread_free_shmem(hperm[mygroup],lcomm);

  par_graph_free(hgraphs[mygroup]);

  dlthread_comm_finalize(lcomm);
  dlthread_free_shmem(hgraphs,ctrl->comm);
}


static void __nd_func(
    void * const ptr)
{
  tid_t myid, nthreads;
  arg_t * arg;
  ctrl_t * ctrl;
  graph_t * graph;
  pid_t * perm;
  pid_t ** dperm;

  arg = ptr;

  ctrl = arg->ctrl;
  graph = arg->graph;
  perm = arg->perm;
  dperm = arg->dperm;

  myid = dlthread_get_id(ctrl->comm);
  nthreads = dlthread_get_nthreads(ctrl->comm);

  dlthread_barrier(ctrl->comm);

  if (nthreads > 1) {
    par_partition_pre(ctrl,graph);
  }

  dperm[myid] = pid_alloc(graph->mynvtxs[myid]);

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.ordering));
  }

  __par_order_nd(ctrl,graph,0,dperm);

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.ordering));
  }

  if (perm) {
    __unify_perm((pid_t const **)dperm,graph,perm);
  }

  dl_free(dperm[myid]);
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void order_nd(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const perm) 
{
  pid_t ** dperm;
  arg_t arg;

  tid_t const nthreads = graph->dist.nthreads;

  dperm = r_pid_alloc(nthreads);

  graph->comm = DLTHREAD_COMM_ROOT;
  ctrl->comm = DLTHREAD_COMM_ROOT;

  arg.ctrl = ctrl;
  arg.graph = graph;
  arg.perm = perm;
  arg.dperm = dperm;

  dlthread_launch(nthreads,&__nd_func,&arg);

  dl_free(dperm);
}




#endif

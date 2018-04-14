/**
 * @file imetis.c
 * @brief Metis wrappers 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 * @date 2015-06-08
 */




#ifndef MTMETIS_IMETIS_C
#define MTMETIS_IMETIS_C




#include "imetis.h"




#undef real_t
#include "../metis/include/metis.h"
#define real_t mtmetis_real_t





/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __create_arrays(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const vwgt,
    wgt_t const * const adjwgt,
    pid_t * const where,
    idx_t ** r_xadj,
    idx_t ** r_adjncy,
    idx_t ** r_vwgt,
    idx_t ** r_adjwgt,
    idx_t ** r_where)
{
  vtx_t i;
  adj_t j;

  /* copy graph if need-be */
  if (sizeof(idx_t) != sizeof(vtx_t)) {
    *r_adjncy = malloc(sizeof(idx_t)*xadj[nvtxs]);
    for (j=0;j<xadj[nvtxs];++j) {
      (*r_adjncy)[j] = (idx_t)adjncy[j];
    }
  } else {
    *r_adjncy = (idx_t*)adjncy;
  }
  if (sizeof(idx_t) != sizeof(adj_t)) {
    *r_xadj = malloc(sizeof(idx_t)*(nvtxs+1));
    for (i=0;i<nvtxs+1;++i) {
      (*r_xadj)[i] = (idx_t)xadj[i];
    }
  } else {
    *r_xadj = (idx_t*)xadj; 
  }
  if (sizeof(idx_t) != sizeof(wgt_t)) {
    *r_vwgt = malloc(sizeof(idx_t)*nvtxs);
    for (i=0;i<nvtxs;++i) {
      (*r_vwgt)[i] = (idx_t)vwgt[i];
    }
    if (adjwgt && r_adjwgt) {
      *r_adjwgt = malloc(sizeof(idx_t)*xadj[nvtxs]);
      for (j=0;j<xadj[nvtxs];++j) {
        (*r_adjwgt)[j] = (idx_t)adjwgt[j];
      }
    }
  } else {
    *r_vwgt = (idx_t*)vwgt;
    if (r_adjwgt) {
      *r_adjwgt = (idx_t*)adjwgt;
    }
  }
  if (sizeof(pid_t) != sizeof(idx_t)) {
    *r_where = malloc(sizeof(idx_t)*nvtxs);
  } else {
    *r_where = (idx_t*)where;
  }
}


static void __destroy_arrays(
    vtx_t const nvtxs,
    pid_t * const where,
    idx_t * m_xadj,
    idx_t * m_adjncy,
    idx_t * m_vwgt,
    idx_t * m_adjwgt,
    idx_t * m_where)
{
  vtx_t i;

  /* free and copy idx_t junk */
  if (sizeof(idx_t) != sizeof(vtx_t)) {
    dl_free(m_adjncy);
  }
  if (sizeof(idx_t) != sizeof(adj_t)) {
    dl_free(m_xadj);
  }
  if (sizeof(idx_t) != sizeof(wgt_t)) {
    dl_free(m_vwgt);
    if (m_adjwgt) {
      dl_free(m_adjwgt);
    }
  }
  if (sizeof(idx_t) != sizeof(pid_t)) {
    for (i=0;i<nvtxs;++i) {
      where[i] = (pid_t)m_where[i];
    }
    dl_free(m_where);
  }
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


wgt_t metis_initcut(
    ctrl_t * const ctrl,
    pid_t const nparts,
    real_t * tpwgts,
    size_t const ncuts,
    int const rb,
    vtx_t const nvtxs,
    adj_t * const xadj,
    vtx_t * const adjncy,
    wgt_t * const vwgt,
    wgt_t * const adjwgt,
    pid_t * const where)
{
  idx_t m_nvtxs, m_nparts, cut, m_ncon, status;
  idx_t options[METIS_NOPTIONS];
  real_t ubf;
  idx_t * m_xadj, * m_adjncy, * m_vwgt, * m_adjwgt, * m_where;

  tid_t const myid = dlthread_get_id(ctrl->comm);

  __METIS_SetDefaultOptions(options);

  m_ncon = 1;

  options[METIS_OPTION_NITER] = 10;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_SEED] = ctrl->seed + myid;
  options[METIS_OPTION_NCUTS] = ncuts;
  options[METIS_OPTION_DBGLVL] = 0;
  options[METIS_OPTION_NO2HOP] = !ctrl->leafmatch;

  m_nparts = (idx_t)nparts;
  m_nvtxs = (idx_t)nvtxs;
  ubf = pow(ctrl->ubfactor,1.0/log(nparts));

  __create_arrays(nvtxs,xadj,adjncy,vwgt,adjwgt,where,&m_xadj,&m_adjncy, \
      &m_vwgt,&m_adjwgt,&m_where);

  status = METIS_OK;
  if (rb || nparts == 2) {
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
    status = __METIS_PartGraphRecursive(&m_nvtxs,&m_ncon,m_xadj, \
        m_adjncy,m_vwgt,NULL,m_adjwgt,&m_nparts,tpwgts, \
        &ubf,options,&cut,m_where);
  } else {
    status = __METIS_PartGraphKway(&m_nvtxs,&m_ncon,m_xadj, \
        m_adjncy,m_vwgt,NULL,m_adjwgt,&m_nparts,tpwgts, \
        &ubf,options,&cut,m_where);
  }

  __destroy_arrays(nvtxs,where,m_xadj,m_adjncy,m_vwgt,m_adjwgt,m_where);

  if (status != METIS_OK) {
    dl_error("Metis returned '%"PRIDX"' during initial partitioning\n",status);
  }

  return (wgt_t)cut;
}


wgt_t metis_initsep(
    ctrl_t * const ctrl,
    size_t const nseps,
    vtx_t const nvtxs,
    adj_t * const xadj,
    vtx_t * const adjncy,
    wgt_t * const vwgt,
    wgt_t * const adjwgt,
    pid_t * const where)
{
  idx_t sep, m_nvtxs;
  idx_t options[METIS_NOPTIONS];
  idx_t * m_xadj, * m_adjncy, * m_vwgt, * m_where;

  tid_t const myid = dlthread_get_id(ctrl->comm);

  __METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NITER] = 10;
  options[METIS_OPTION_DBGLVL] = 0;
  options[METIS_OPTION_SEED] = ctrl->seed + myid;
  options[METIS_OPTION_NSEPS] = nseps;
  options[METIS_OPTION_NCUTS] = nseps;
  options[METIS_OPTION_UFACTOR] = 1000*(ctrl->ubfactor - 1.0);

  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_NODE;

  m_nvtxs = nvtxs;

  __create_arrays(nvtxs,xadj,adjncy,vwgt,NULL,where,&m_xadj,&m_adjncy, \
      &m_vwgt,NULL,&m_where);

  __METIS_ComputeVertexSeparator(&m_nvtxs,m_xadj,m_adjncy, \
      m_vwgt,options,&sep,m_where);

  __destroy_arrays(nvtxs,where,m_xadj,m_adjncy,m_vwgt,NULL,m_where);

  return (wgt_t)sep;
}


wgt_t metis_kway(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const * const where,
    int const rb)
{
  idx_t m_nparts, m_nvtxs, cut, m_ncon;
  idx_t options[METIS_NOPTIONS];
  real_t ubf;
  idx_t * m_xadj, * m_adjncy, * m_vwgt, * m_adjwgt, * m_where;

  tid_t const myid = dlthread_get_id(ctrl->comm);

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.metis));
  }

  __METIS_SetDefaultOptions(options);

  m_ncon = 1;

  options[METIS_OPTION_NITER] = ctrl->nrefpass;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_SEED] = ctrl->seed + myid;
  options[METIS_OPTION_NCUTS] = ctrl->ncuts;
  options[METIS_OPTION_NO2HOP] = !ctrl->leafmatch;
  options[METIS_OPTION_DBGLVL] = 0;
  m_nparts = ctrl->nparts;
  ubf = ctrl->ubfactor;

  __create_arrays(graph->mynvtxs[0],graph->xadj[0],graph->adjncy[0], \
      graph->vwgt[0],graph->adjwgt[0],where[0],&m_xadj,&m_adjncy, \
      &m_vwgt,&m_adjwgt,&m_where);

  m_nvtxs = graph->mynvtxs[0];

  if (rb) {
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
    __METIS_PartGraphRecursive(&m_nvtxs,&m_ncon,m_xadj,m_adjncy,m_vwgt,NULL, \
        m_adjwgt,&m_nparts,NULL,&ubf,options,&cut,m_where);
  } else {
    __METIS_PartGraphKway(&m_nvtxs,&m_ncon,m_xadj,m_adjncy,m_vwgt,NULL, \
        m_adjwgt,&m_nparts,NULL,&ubf,options,&cut,m_where);
  }

  __destroy_arrays(graph->mynvtxs[0],where[0],m_xadj,m_adjncy,m_vwgt, \
      m_adjwgt,m_where);

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.metis));
  }

  graph->mincut = cut;

  return cut;
}


wgt_t metis_esep(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const * const where)
{
  vtx_t i;
  pid_t me;
  idx_t curobj;
  idx_t options[METIS_NOPTIONS];
  idx_t ncon = 1, nparts = 2;
  real_t ubf = 1.0;
  idx_t * m_xadj, * m_adjncy, * m_vwgt, * m_adjwgt, * m_where;

  tid_t const myid = dlthread_get_id(ctrl->comm);

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.metis));
  }

  __METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NITER] = 10;
  if (ctrl->verbosity == MTMETIS_VERBOSITY_MAXIMUM) {
    options[METIS_OPTION_DBGLVL] = 15;
  } else {
    options[METIS_OPTION_DBGLVL] = 0;
  }
  options[METIS_OPTION_NSEPS] = 1;
  options[METIS_OPTION_NCUTS] = 1;
  options[METIS_OPTION_NITER] = ctrl->nrefpass;
  options[METIS_OPTION_UFACTOR] = 1000*(ctrl->ubfactor - 1.0);
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

  ubf = ctrl->ubfactor;

  options[METIS_OPTION_SEED] = ctrl->seed;

  graph->pwgts = wgt_init_alloc(0,2);
  graph->where = r_pid_alloc(1);
  graph->where[0] = pid_alloc(graph->mynvtxs[0]);

  __create_arrays(graph->mynvtxs[0],graph->xadj[0],graph->adjncy[0], \
      graph->vwgt[0],graph->adjwgt[0],where[0],&m_xadj,&m_adjncy, \
      &m_vwgt,&m_adjwgt,&m_where);

  __METIS_PartGraphRecursive((idx_t*)graph->mynvtxs,&ncon, \
      (idx_t*)graph->xadj[0],(idx_t*)graph->adjncy[0], \
      (idx_t*)graph->vwgt[0],NULL,(idx_t*)graph->adjwgt[0],&nparts,NULL, \
      &ubf,options,&curobj,(idx_t*)graph->where[0]);

  __destroy_arrays(graph->mynvtxs[0],where[0],m_xadj,m_adjncy,m_vwgt, \
      m_adjwgt,m_where);

  for (i=0;i<graph->mynvtxs[0];++i) {
    me = graph->where[0][i];
    graph->pwgts[me] += graph->vwgt[0][i];
  }

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.metis));
  }

  return (wgt_t)curobj;
}


wgt_t metis_vsep(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const * const where)
{
  vtx_t i;
  pid_t me;
  idx_t curobj, m_nvtxs;
  idx_t options[METIS_NOPTIONS];
  idx_t * m_xadj, * m_adjncy, * m_vwgt, * m_where;

  tid_t const myid = dlthread_get_id(ctrl->comm);

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.metis));
  }

  __METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NITER] = 10;
  if (ctrl->verbosity == MTMETIS_VERBOSITY_MAXIMUM) {
    options[METIS_OPTION_DBGLVL] = 15;
  } else {
    options[METIS_OPTION_DBGLVL] = 0;
  }
  options[METIS_OPTION_NSEPS] = 1;
  options[METIS_OPTION_NCUTS] = 1;
  options[METIS_OPTION_NITER] = ctrl->nrefpass;
  options[METIS_OPTION_UFACTOR] = 1000*(ctrl->ubfactor - 1.0);
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_NODE;

  options[METIS_OPTION_SEED] = ctrl->seed;

  graph->pwgts = wgt_init_alloc(0,3);
  graph->where = r_pid_alloc(1);
  graph->where[0] = pid_alloc(graph->mynvtxs[0]);

  m_nvtxs = graph->mynvtxs[0];

  __create_arrays(graph->mynvtxs[0],graph->xadj[0],graph->adjncy[0], \
      graph->vwgt[0],NULL,where[0],&m_xadj,&m_adjncy, \
      &m_vwgt,NULL,&m_where);

  __METIS_ComputeVertexSeparator(&m_nvtxs,m_xadj,m_adjncy, \
      m_vwgt,options,&curobj,m_where);

  __destroy_arrays(graph->mynvtxs[0],where[0],m_xadj,m_adjncy,m_vwgt, \
      NULL,m_where);

  for (i=0;i<graph->mynvtxs[0];++i) {
    me = graph->where[0][i];
    graph->pwgts[me] += graph->vwgt[0][i];
  }

  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.metis));
  }

  return (wgt_t)curobj;
}





#endif







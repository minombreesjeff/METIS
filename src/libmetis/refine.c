/*
\file
\brief This file contains the driving routines for multilevel refinement

\date   Started 7/24/1997
\author George  
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version\verbatim $Id: refine.c 9987 2011-05-26 16:35:25Z karypis $ \endverbatim
*/

#include "metislib.h"


/*************************************************************************/
/*! This function is the entry point of refinement */
/*************************************************************************/
void Refine2Way(ctrl_t *ctrl, graph_t *orggraph, graph_t *graph, real_t *tpwgts)
{

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->UncoarsenTmr));

  /* Compute the parameters of the coarsest graph */
  Compute2WayPartitionParams(ctrl, graph);

  for (;;) {
    ASSERT(CheckBnd(graph));

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->RefTmr));

    Balance2Way(ctrl, graph, tpwgts);

    FM_2WayRefine(ctrl, graph, tpwgts, ctrl->niter); 

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->RefTmr));

    if (graph == orggraph)
      break;

    graph = graph->finer;
    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->ProjectTmr));
    Project2WayPartition(ctrl, graph);
    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->ProjectTmr));
  }

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->UncoarsenTmr));
}


/*************************************************************************/
/*! This function allocates memory for 2-way edge refinement */
/*************************************************************************/
void Allocate2WayPartitionMemory(ctrl_t *ctrl, graph_t *graph)
{
  idx_t nvtxs, ncon;

  nvtxs = graph->nvtxs;
  ncon  = graph->ncon;

  graph->pwgts  = imalloc(2*ncon, "Allocate2WayPartitionMemory: pwgts");
  graph->where  = imalloc(nvtxs, "Allocate2WayPartitionMemory: where");
  graph->bndptr = imalloc(nvtxs, "Allocate2WayPartitionMemory: bndptr");
  graph->bndind = imalloc(nvtxs, "Allocate2WayPartitionMemory: bndind");
  graph->id     = imalloc(nvtxs, "Allocate2WayPartitionMemory: id");
  graph->ed     = imalloc(nvtxs, "Allocate2WayPartitionMemory: ed");
}


/*************************************************************************/
/*! This function computes the initial id/ed */
/*************************************************************************/
void Compute2WayPartitionParams(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, l, nvtxs, ncon, nbnd, mincut;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *pwgts;
  idx_t *where, *bndptr, *bndind, *id, *ed;
  idx_t me, other;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where  = graph->where;
  pwgts  = iset(2*ncon, 0, graph->pwgts);
  bndptr = iset(nvtxs, -1, graph->bndptr);
  bndind = graph->bndind;

  nbnd = 0;


  /* Compute pwgts */
  if (ncon == 1) {
    for (i=0; i<nvtxs; i++) {
      ASSERT(where[i] >= 0 && where[i] <= 1);
      pwgts[where[i]] += vwgt[i];
    }
    ASSERT(pwgts[0]+pwgts[1] == graph->tvwgt[0]);
  }
  else {
    for (i=0; i<nvtxs; i++) {
      me = where[i];
      for (j=0; j<ncon; j++)
        pwgts[me*ncon+j] += vwgt[i*ncon+j];
    }
  }


  /* Compute the required info for refinement  */
  id = iset(nvtxs, 0, graph->id);
  ed = iset(nvtxs, 0, graph->ed);
  
  for (mincut=0, i=0; i<nvtxs; i++) {
    me = where[i];
  
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]])
        id[i] += adjwgt[j];
      else
        ed[i] += adjwgt[j];
    }
  
    if (ed[i] > 0 || xadj[i] == xadj[i+1]) {
      BNDInsert(nbnd, bndind, bndptr, i);
      mincut += ed[i];
    }
  }

  graph->mincut = mincut/2;
  graph->nbnd   = nbnd;

}


/*************************************************************************/
/*! Projects a partition and computes the refinement params. */
/*************************************************************************/
void Project2WayPartition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, nvtxs, nbnd, me;
  idx_t *xadj, *adjncy, *adjwgt, *adjrsum;
  idx_t *cmap, *where, *bndptr, *bndind;
  idx_t *cwhere, *cbndptr;
  idx_t *id, *ed, *cid, *ced;
  graph_t *cgraph;

  Allocate2WayPartitionMemory(ctrl, graph);

  cgraph  = graph->coarser;
  cwhere  = cgraph->where;
  cbndptr = cgraph->bndptr;

  nvtxs   = graph->nvtxs;
  cmap    = graph->cmap;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  adjrsum = graph->adjrsum;

  where  = graph->where;
  bndptr = iset(nvtxs, -1, graph->bndptr);
  bndind = graph->bndind;

  nbnd = 0;

  /* Project the partition and record which of these nodes came from the
     coarser boundary */
  for (i=0; i<nvtxs; i++) {
    k = cmap[i];
    where[i] = cwhere[k];
    cmap[i]  = cbndptr[k];
  }


  /* Compute the refinement information of the nodes */
  id  = iset(nvtxs, 0, graph->id);
  ed  = iset(nvtxs, 0, graph->ed);
  cid = cgraph->id;
  ced = cgraph->ed;
  
  for (i=0; i<nvtxs; i++) {
    me = where[i];
  
    id[i] = adjrsum[i];
  
    if (xadj[i] == xadj[i+1]) {
      BNDInsert(nbnd, bndind, bndptr, i);
    }
    else {
      if (cmap[i] != -1) { /* If it is an interface node. Note that cmap[i] = cbndptr[cmap[i]] */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          if (me != where[adjncy[j]])
            ed[i] += adjwgt[j];
        }
        id[i] -= ed[i];
  
        if (ed[i] > 0 || xadj[i] == xadj[i+1]) {
          BNDInsert(nbnd, bndind, bndptr, i);
        }
      }
    }
  }
  graph->mincut = cgraph->mincut;
  graph->nbnd   = nbnd;

  /* copy pwgts */
  icopy(2*graph->ncon, cgraph->pwgts, graph->pwgts);

  FreeGraph(&graph->coarser);
  graph->coarser = NULL;

}


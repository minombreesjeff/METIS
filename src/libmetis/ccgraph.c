/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ccgraph.c
 *
 * This file contains the functions that create the coarse graph
 *
 * Started 8/11/97
 * George
 *
 */

#include "metislib.h"



/*************************************************************************/
/*! This function creates the coarser graph */
/*************************************************************************/
void CreateCoarseGraph(ctrl_t *ctrl, graph_t *graph, idx_t cnvtxs, 
         idx_t *match, idx_t *perm)
{
  idx_t i, j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon, cnedges, 
        v, u, mask, dovsize;
  idx_t *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjrsum, *auxadj;
  idx_t *cmap, *htable;
  idx_t *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjrsum;
  graph_t *cgraph;

  dovsize = (ctrl->objtype == METIS_OBJTYPE_VOL ? 1 : 0);

  mask = HTLENGTH;
  if (cnvtxs < 8*mask || graph->nedges/graph->nvtxs > 15) { 
    CreateCoarseGraphNoMask(ctrl, graph, cnvtxs, match, perm);
    return;
  }

  WCOREPUSH;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->ContractTmr));

  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  xadj    = graph->xadj;
  vwgt    = graph->vwgt;
  vsize   = graph->vsize;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  adjrsum = graph->adjrsum;
  cmap    = graph->cmap;

  /* Initialize the coarser graph */
  cgraph   = SetupCoarseGraph(graph, cnvtxs, dovsize);
  cxadj    = cgraph->xadj;
  cvwgt    = cgraph->vwgt;
  cvsize   = cgraph->vsize;
  cadjrsum = cgraph->adjrsum;
  cadjncy  = cgraph->adjncy;
  cadjwgt  = cgraph->adjwgt;


  iend = xadj[nvtxs];
  auxadj = icopy(iend, adjncy, iwspacemalloc(ctrl, iend));
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  htable = iset(mask+1, -1, iwspacemalloc(ctrl, mask+1)); 

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      icopy(ncon, vwgt+v*ncon, cvwgt+cnvtxs*ncon);

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    cadjrsum[cnvtxs] = adjrsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m] += adjwgt[j];
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj] += adjwgt[j];
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges] = k;
          cadjwgt[nedges++] = adjwgt[j];
        }
      }
    }

    if (v != u) { 
      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        iaxpy(ncon, 1, vwgt+u*ncon, 1, cvwgt+cnvtxs*ncon, 1);

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      cadjrsum[cnvtxs] += adjrsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[kk] = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m] += adjwgt[j];
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj] += adjwgt[j];
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges] = k;
            cadjwgt[nedges++] = adjwgt[j];
          }
        }
      }

      /* Remove the contracted adjacency weight */
      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs) 
            break;
        }
      }
      if (jj >= 0 && cadjncy[jj] == cnvtxs) { /* This 2nd check is needed for non-adjacent matchings */
        cadjrsum[cnvtxs] -= cadjwgt[jj];
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      }
    }

    ASSERTP(cadjrsum[cnvtxs] == isum(nedges, cadjwgt, 1), 
        ("%"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", 
         cnvtxs, cadjrsum[cnvtxs], isum(nedges, cadjwgt, 1), adjrsum[u], adjrsum[v]));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* Zero out the htable */
    htable[cnvtxs&mask] = -1;

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  for (i=0; i<ncon; i++) {
    cgraph->tvwgt[i]    = isum(cgraph->nvtxs, cgraph->vwgt+i, ncon);
    cgraph->invtvwgt[i] = 1.0/(cgraph->tvwgt[i] > 0 ? cgraph->tvwgt[i] : 1);
  }


  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->ContractTmr));

  WCOREPOP;
}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraphNoMask(ctrl_t *ctrl, graph_t *graph, idx_t cnvtxs, 
         idx_t *match, idx_t *perm)
{
  idx_t i, j, k, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, dovsize;
  idx_t *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjrsum, *auxadj;
  idx_t *cmap, *htable;
  idx_t *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjrsum;
  graph_t *cgraph;

  WCOREPUSH;

  dovsize = (ctrl->objtype == METIS_OBJTYPE_VOL ? 1 : 0);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->ContractTmr));

  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  xadj    = graph->xadj;
  vwgt    = graph->vwgt;
  vsize   = graph->vsize;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  adjrsum = graph->adjrsum;
  cmap    = graph->cmap;


  /* Initialize the coarser graph */
  cgraph = SetupCoarseGraph(graph, cnvtxs, dovsize);
  cxadj    = cgraph->xadj;
  cvwgt    = cgraph->vwgt;
  cvsize   = cgraph->vsize;
  cadjrsum = cgraph->adjrsum;
  cadjncy  = cgraph->adjncy;
  cadjwgt  = cgraph->adjwgt;

  htable = iset(cnvtxs, -1, iwspacemalloc(ctrl, cnvtxs));

  iend = xadj[nvtxs];
  auxadj = icopy(iend, adjncy, iwspacemalloc(ctrl, iend));
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      icopy(ncon, vwgt+v*ncon, cvwgt+cnvtxs*ncon);

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    cadjrsum[cnvtxs] = adjrsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      if ((m = htable[k]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[k] = nedges++;
      }
      else {
        cadjwgt[m] += adjwgt[j];
      }
    }

    if (v != u) { 
      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        iaxpy(ncon, 1, vwgt+u*ncon, 1, cvwgt+cnvtxs*ncon, 1);

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      cadjrsum[cnvtxs] += adjrsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        if ((m = htable[k]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[k] = nedges++;
        }
        else {
          cadjwgt[m] += adjwgt[j];
        }
      }

      /* Remove the contracted adjacency weight */
      if ((j = htable[cnvtxs]) != -1) {
        ASSERT(cadjncy[j] == cnvtxs);
        cadjrsum[cnvtxs] -= cadjwgt[j];
        cadjncy[j] = cadjncy[--nedges];
        cadjwgt[j] = cadjwgt[nedges];
        htable[cnvtxs] = -1;
      }
    }

    ASSERTP(cadjrsum[cnvtxs] == isum(nedges, cadjwgt, 1), 
        ("%"PRIDX" %"PRIDX"\n", cadjrsum[cnvtxs], isum(nedges, cadjwgt, 1)));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]] = -1;  /* Zero out the htable */

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  for (i=0; i<ncon; i++) {
    cgraph->tvwgt[i]    = isum(cgraph->nvtxs, cgraph->vwgt+i, ncon);
    cgraph->invtvwgt[i] = 1.0/(cgraph->tvwgt[i] > 0 ? cgraph->tvwgt[i] : 1);
  }

  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->ContractTmr));

  WCOREPOP;
}


/*************************************************************************
* Setup the various arrays for the coarse graph
**************************************************************************/
graph_t *SetupCoarseGraph(graph_t *graph, idx_t cnvtxs, idx_t dovsize)
{
  graph_t *cgraph;

  cgraph = CreateGraph();

  cgraph->nvtxs = cnvtxs;
  cgraph->ncon  = graph->ncon;

  cgraph->finer  = graph;
  graph->coarser = cgraph;


  /* Allocate memory for the coarser graph */
  cgraph->xadj     = imalloc(cnvtxs+1, "SetupCoarseGraph: xadj");
  cgraph->adjrsum  = imalloc(cnvtxs,   "SetupCoarseGraph: adjrsum");
  cgraph->cmap     = imalloc(cnvtxs,   "SetupCoarseGraph: cmap");
  cgraph->adjncy   = imalloc(graph->nedges,   "SetupCoarseGraph: adjncy");
  cgraph->adjwgt   = imalloc(graph->nedges,   "SetupCoarseGraph: adjwgt");
  cgraph->vwgt     = imalloc(cgraph->ncon*cnvtxs, "SetupCoarseGraph: vwgt");
  cgraph->tvwgt    = imalloc(cgraph->ncon, "SetupCoarseGraph: tvwgt");
  cgraph->invtvwgt = rmalloc(cgraph->ncon, "SetupCoarseGraph: invtvwgt");

  if (dovsize)
    cgraph->vsize = imalloc(cnvtxs,   "SetupCoarseGraph: vsize");

  return cgraph;
}


/*************************************************************************
* This function re-adjusts the amount of memory that was allocated if
* it will lead to significant savings
**************************************************************************/
void ReAdjustMemory(graph_t *graph, graph_t *cgraph, idx_t dovsize) 
{
  if (cgraph->nedges > 10000 && cgraph->nedges < 0.8*graph->nedges) {
    cgraph->adjncy = irealloc(cgraph->adjncy, cgraph->nedges, "ReAdjustMemory: adjncy");
    cgraph->adjwgt = irealloc(cgraph->adjwgt, cgraph->nedges, "ReAdjustMemory: adjwgt");
  }
}

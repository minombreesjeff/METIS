/**
 * @file flow.c
 * @brief Functions for computing maximum-flows / minimum-cuts
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Dominique LaSalle
 * @version 1
 * @date 2014-12-11
 */




#ifndef BOWSTRING_FLOW_C
#define BOWSTRING_FLOW_C




#include "flow.h"
#include "graph.h"




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static inline vtx_t __inc_cycle(
    vtx_t const i,
    vtx_t const n)
{
  if (i == n-1) {
    return 0;
  } else {
    return i+1;
  }
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


wgt_t maxflow_edge(
    vtx_t src,
    vtx_t dst,
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt,
    wgt_t * flow)
{
  vtx_t v, j, k, h, sq, nq, minh;
  wgt_t cap, maxflow;
  vtx_t * q, * height;
  wgt_t * excess;
  int * active;
  adj_t * radj;

  adj_t const nedges = xadj[nvtxs];
  wgt_t const tadjwgt = wgt_sum(adjwgt,nedges);

  /* flow structures */
  q = vtx_alloc(nvtxs);
  active = int_init_alloc(0,nvtxs);
  excess = wgt_init_alloc(0,nvtxs);
  height = vtx_init_alloc(0,nvtxs);

  height[src] = nvtxs;
  excess[src] = tadjwgt; 
  active[src] = 1;

  radj = adj_alloc(nedges);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  wgt_set(flow,0,nedges);

  sq = 0;
  nq = 0;

  /* saturate source edges */
  for (j=xadj[src];j<xadj[src+1];++j) {
    k = adjncy[j];
    /* perform a push */
    flow[j] = adjwgt[j];
    flow[radj[j]] = -adjwgt[j];
    excess[k] += adjwgt[j];
    active[k] = 1;
    q[nq] = k;
    nq = __inc_cycle(nq,nvtxs);
  }

  while (sq != nq) {
    /* this will loop until a steady state is reached -- vertices are only
     * added to the queue as their flows change */
    v = q[sq];
    sq = __inc_cycle(sq,nvtxs); 
    active[v] = 0;

    DL_ASSERT(excess[v] > 0,"Inspecting vertex without excess\n");

    if (v == dst) {
      /* skip destination vertex */
      continue;
    }

    h = height[v];

    minh = nvtxs;
    for (j=xadj[v];j<xadj[v+1];++j) {
      k = adjncy[j];
      if (height[k] < h) {
        /* perform a push */
        cap = dl_min(adjwgt[j]-flow[j],excess[v]);
        if (cap > 0) {
          /* check to see which way this edge is flowing */
          excess[v] -= cap;
          excess[k] += cap;
          flow[j] += cap;
          flow[radj[j]] -= cap;
          /* mark this vertex as active */
          if (!active[k]) {
            active[k] = 1;
            q[nq] = k;
            nq = __inc_cycle(nq,nvtxs);
          }
        }
      } else if (height[k] < minh) {
        minh = height[k];
      }
      if (excess[v] <= 0) {
        /* if we have no more exces, exit */
        break;
      }
    }
    if (minh < nvtxs && excess[v] > 0) {
      /* perform a relabel operation */
      height[v] = minh+1;
      /* put this vertex back in teh queue */
      active[v] = 1;
      q[nq] = v;
      nq = __inc_cycle(nq,nvtxs);
    }
  }

  maxflow = excess[dst];

  dl_free(height);
  dl_free(excess);
  dl_free(active);
  dl_free(q);
  dl_free(radj);

  return maxflow;
}


wgt_t maxflow_vertex(
    vtx_t const src,
    vtx_t const dst,
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const vwgt,
    wgt_t * const flow)
{
  vtx_t v, j, k, h, sq, nq, minh;
  wgt_t cap, maxflow;
  vtx_t * q, * height;
  wgt_t * excess, * in;
  int * active;
  adj_t * radj;

  adj_t const nedges = xadj[nvtxs];
  wgt_t const tvwgt = wgt_sum(vwgt,nvtxs);

  /* flow structures */
  q = vtx_alloc(nvtxs);
  active = int_init_alloc(0,nvtxs);
  excess = wgt_init_alloc(0,nvtxs);
  in = wgt_init_alloc(0,nvtxs);
  height = vtx_init_alloc(0,nvtxs);

  height[src] = nvtxs;
  excess[src] = tvwgt; 
  active[src] = 1;
  in[dst] = -tvwgt;

  radj = adj_alloc(nedges);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  wgt_set(flow,0,nedges);

  sq = 0;
  nq = 0;

  /* saturate source edges */
  for (j=xadj[src];j<xadj[src+1];++j) {
    k = adjncy[j];
    /* perform a push */
    flow[j] = vwgt[j];
    flow[radj[j]] = -vwgt[k];
    excess[k] += vwgt[k];
    in[k] += vwgt[k];
    active[k] = 1;
    q[nq] = k;
    nq = __inc_cycle(nq,nvtxs);
  }

  while (sq != nq) {
    /* this will loop until a steady state is reached -- vertices are only
     * added to the queue as their flows change */
    v = q[sq];
    sq = __inc_cycle(sq,nvtxs); 
    active[v] = 0;

    DL_ASSERT(excess[v] > 0,"Inspecting vertex without excess\n");

    if (v == dst) {
      /* skip destination vertex */
      continue;
    }

    h = height[v];

    minh = nvtxs;
    for (j=xadj[v];j<xadj[v+1];++j) {
      k = adjncy[j];
      if (height[k] < h) {
        /* if there is no longer a valid flow on this edge, cancel it */
        if (flow[j] < 0) {
          in[v] += flow[j];
          excess[v] += flow[j];
          excess[k] -= flow[j];
          flow[j] = flow[radj[j]] = 0;
          if (!active[k]) {
            /* activate the vertex */
            active[k] = 1;
            q[nq] = k;
            nq = __inc_cycle(nq,nvtxs);
          }
        } 
        /* perform a push */
        cap = dl_min(vwgt[k]-in[k],excess[v]);
        if (cap > 0) {
          /* induce my new flow */
          excess[v] -= cap;
          excess[k] += cap;
          flow[j] += cap;
          flow[radj[j]] -= cap;
          in[k] += cap;
          /* mark this vertex as active */
          if (!active[k]) {
            active[k] = 1;
            q[nq] = k;
            nq = __inc_cycle(nq,nvtxs);
          }
        }
      } else if (height[k] < minh) {
        minh = height[k];
      }
      if (excess[v] <= 0) {
        /* if we have no more excess, exit */
        break;
      }
    }
    if (minh < nvtxs && excess[v] > 0) {
      /* perform a relabel operation */
      height[v] = minh+1;
      /* put this vertex back in teh queue */
      active[v] = 1;
      q[nq] = v;
      nq = __inc_cycle(nq,nvtxs);
    }
  }

  maxflow = excess[dst];

  dl_free(height);
  dl_free(excess);
  dl_free(active);
  dl_free(q);
  dl_free(radj);

  return maxflow;
}





#endif


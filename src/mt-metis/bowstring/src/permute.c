/**
 * @file permute.c
 * @brief Permutation functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-07-18
 */




#ifndef PERMUTE_C
#define PERMUTE_C




#include "permute.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void permute_graph(
    const vtx_t nvtxs, 
    adj_t * const xadj, 
    vtx_t * const adjncy,
    wgt_t * const vwgt, 
    wgt_t * const adjwgt)
{
  /* generate my permutation array */
  vtx_t * perm = vtx_alloc(nvtxs);

  vtx_incset(perm,0,1,nvtxs);
  dl_init_rand();
  perm = vtx_shuffle(perm,nvtxs);

  reorder_graph(nvtxs,xadj,adjncy,vwgt,adjwgt,perm);
  dl_free(perm);
}


void reorder_graph(
    const vtx_t nvtxs, 
    adj_t * const xadj, 
    vtx_t * const adjncy,
    wgt_t * const vwgt, 
    wgt_t * const adjwgt, 
    const vtx_t * const perm) 
{
  vtx_t i, k, pi;
  adj_t j;
  const size_t nedges = xadj[nvtxs];

  /* allocate my permuted junk */
  adj_t * rxadj = adj_alloc(nvtxs+1);
  vtx_t * radjncy = vtx_alloc(nedges);

  wgt_t * rvwgt = NULL, * radjwgt = NULL; 
  if (vwgt) {
    rvwgt = wgt_alloc(nvtxs);
  }
  if (adjwgt) {
    radjwgt = wgt_alloc(nedges);
  }

  vtx_t * rename = vtx_alloc(nvtxs);
  
  /* permute the xadj and the vwgt */
  rxadj[0] = 0;
  for (i=0;i<nvtxs;++i) {
    pi = perm[i];
    rename[pi] = i;
    rxadj[i+1] = (xadj[pi+1]-xadj[pi]) + rxadj[i];
    if (vwgt) {
      rvwgt[i] = vwgt[pi];
    }
  }

  /* permute the adjncy and adjwgt */
  k = 0;
  for (i=0;i<nvtxs;++i) {
    pi = perm[i];
    for (j=xadj[pi];j<xadj[pi+1];++j) {
      radjncy[k] = rename[adjncy[j]];
      if (adjwgt) {
        radjwgt[k] = adjwgt[j];
      }
      ++k;
    }
  }

  dl_free(rename);

  adj_copy(xadj,rxadj,nvtxs+1);
  vtx_copy(adjncy,radjncy,nedges);
  dl_free(rxadj);
  dl_free(radjncy);
  
  if (vwgt) {
    wgt_copy(vwgt,rvwgt,nvtxs);
    dl_free(rvwgt);
  }

  if (adjwgt) {
    wgt_copy(adjwgt,radjwgt,nedges);
    dl_free(radjwgt);
  }
}




#endif


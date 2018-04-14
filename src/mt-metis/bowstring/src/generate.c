/**
 * @file generate.c
 * @brief Functions for generating graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-10-19
 */




#ifndef BOWSTRING_GENERATE_C
#define BOWSTRING_GENERATE_C




#include "generate.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void generate_complete_graph(
    const vtx_t nvtxs, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  vtx_t i, k;
  adj_t j;

  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;

  const adj_t nedges = nvtxs*(nvtxs-1);

  xadj = adj_alloc(nvtxs+1);
  adjncy = vtx_alloc(nedges);

  /* populate the edges */
  j=0;
  xadj[0] = 0;
  for (i=0;i<nvtxs;++i) {
    for (k=0;k<i;++k) {
      adjncy[j++] = k;
    }
    for (k=i+1;k<nvtxs;++k) {
      adjncy[j++] = k;
    }
    xadj[i+1] = j;
  }

  if (r_xadj) {
    *r_xadj = xadj;
  } else {
    dl_free(xadj);
  }

  if (r_adjncy) {
    *r_adjncy = adjncy;
  } else {
    dl_free(adjncy);
  }

  if (r_vwgt) {
    *r_vwgt = wgt_init_alloc(1.0,nvtxs);
  }

  if (r_adjwgt) {
    *r_adjwgt = wgt_init_alloc(1.0,nedges);
  }
}


void generate_grid_graph(
    const vtx_t nvx, 
    const vtx_t nvy, 
    const vtx_t nvz, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  vtx_t x,y,z,i;
  adj_t j;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;

  const vtx_t nvtxs = nvx*nvy*nvz;
  const adj_t nedges = ((nvx-1)*(nvy*nvz) + (nvy-1)*(nvx*nvz) + 
      (nvz-1)*(nvx*nvy)) * 2;

  xadj = adj_calloc(nvtxs+1);
  adjncy = vtx_alloc(nedges);

  xadj[0] = j = 0;
  for (x=0;x<nvx;++x) {
    for (y=0;y<nvy;++y) {
      for (z=0;z<nvz;++z) {
        i = (x*nvy*nvz) + (y*nvz) + z;
        if (x>0) {
          /* -x */
          adjncy[j++] = i - (nvy*nvz);
        }
        if (x<nvx-1) {
          /* +x */
          adjncy[j++] = i + (nvy*nvz);
        }
        if (y>0) {
          /* -y */
          adjncy[j++] = i - nvz;
        }
        if (y<nvy-1) {
          /* +y */
          adjncy[j++] = i + nvz;
        }
        if (z>0) {
          /* -z */
          adjncy[j++] = i - 1;
        }
        if (z<nvz-1) {
          /* +z */
          adjncy[j++] = i + 1;
        }
        xadj[i+1] = j;
      }
    }
  }

  if (r_xadj) {
    *r_xadj = xadj;
  } else {
    dl_free(xadj);
  }

  if (r_adjncy) {
    *r_adjncy = adjncy;
  } else {
    dl_free(adjncy);
  }

  if (r_vwgt) {
    *r_vwgt = wgt_init_alloc(1.0,nvtxs);
  }

  if (r_adjwgt) {
    *r_adjwgt = wgt_init_alloc(1.0,nedges);
  }
}




#endif

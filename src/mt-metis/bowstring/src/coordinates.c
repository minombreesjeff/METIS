/**
 * @file coordinates.c
 * @brief Functions for assigning coordinates
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2014-01-26
 */




#ifndef BOWSTRING_COORDINATES_C
#define BOWSTRING_COORDINATES_C




#include "coordinates.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int bfs_coordinates(const vtx_t nvtxs, const adj_t * const xadj, 
    const vtx_t * const adjncy, const wgt_t * const adjwgt, const size_t ndim, 
    coord_t ** const coords)
{
  vtx_t start, i, k, sq, nq;
  adj_t j;
  size_t dim;

  int * visited;
  vtx_t * q;

  visited = int_calloc(nvtxs);
  q = vtx_alloc(nvtxs);

  for (dim=0;dim<ndim;++dim) {
    sq = nq = 0;
    start = vtx_rand(0,nvtxs);
    q[nq++] = start;
    visited[start] = (int)(dim+1);
    coords[dim][start] = 0.0;

    while (sq < nvtxs) {
      i = q[sq++];
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        if (visited[k] == (int)dim) {
          if (adjwgt) {
            coords[dim][k] = coords[dim][i]+(1.0/adjwgt[j]);
          } else {
            coords[dim][k] = coords[dim][i]+1.0;
          }
          coords[dim][k] += coord_rand(0,0.1);
          q[nq++] = k;
          visited[k] = (int)(dim+1);
        }
      }
    }
  }

  dl_free(q);
  dl_free(visited);

  return BOWSTRING_SUCCESS;
}




#endif


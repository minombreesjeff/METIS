/**
 * @file cut.c
 * @brief Functions for cutting graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2014-01-26
 */




#ifndef BOWSTRING_CUT_C
#define BOWSTRING_CUT_C




#include "cut.h"
#include "analyze.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static vlbl_t const NSPLIT = 2;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static int __split_graph(
    vtx_t const nvtxs, 
    adj_t const * xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    wgt_t const * const adjwgt, 
    vlbl_t const * const where, 
    coord_t const ** const coords, 
    size_t const ndim, 
    vtx_t * snvtxs, 
    adj_t ** sxadj, 
    vtx_t ** sadjncy, 
    wgt_t ** const svwgt, 
    wgt_t ** const sadjwgt, 
    coord_t *** const scoords, 
    vtx_t ** const lbl)
{
  size_t d;
  vtx_t i,k,v;
  adj_t j;
  vlbl_t w,x;
  
  adj_t snedges[NSPLIT];
  vtx_t * rename;

  rename = vtx_alloc(nvtxs);

  /* initialize counters */
  for (w=0;w<NSPLIT;++w) {
    snvtxs[w] = 0;
    snedges[w] = 0;
  }

  /* split coords/vwgt and count vertices */
  for (i=0;i<nvtxs;++i) {
    w = where[i];
    if (lbl) {
      lbl[w][snvtxs[w]] = i;
    }
    rename[i] = snvtxs[w];
    if (svwgt) {
      svwgt[w][snvtxs[w]] = vwgt[i];
    }
    snedges[w] += xadj[i+1] - xadj[i];
    ++snvtxs[w];
  }

  /* allocate edges */
  for (w=0;w<NSPLIT;++w) {
    scoords[w] = r_coord_sym_alloc(snvtxs[w],ndim);
    sxadj[w] = adj_alloc(snvtxs[w]+1);
    sxadj[w][0] = 0;
    sadjncy[w] = vtx_alloc(snedges[w]);
    if (adjwgt) {
      sadjwgt[w] = wgt_alloc(snedges[w]);
    }
    snedges[w] = 0;
  }

  /* populate edge arrays */
  for (i=0;i<nvtxs;++i) {
    w = where[i];
    v = rename[i];

    /* set coordinates */
    for (d=0;d<ndim;++d) {
      scoords[w][d][v] = coords[d][i];
    }

    /* avoid prefix sum */
    sxadj[w][v+1] = sxadj[w][v];   

    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      x = where[k];
      if (x == w) {
        sadjncy[w][sxadj[w][v+1]] = rename[k];
        if (adjwgt) {
          sadjwgt[w][sxadj[w][v+1]] = adjwgt[j];
        }
        ++sxadj[w][v+1];
      } else {
        /* throw away */
      }
    }
  }

  dl_free(rename);

  return BOWSTRING_SUCCESS;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void coordinate_bisection(
    vtx_t const nvtxs, 
    wgt_t const * const vwgt, 
    coord_t const * const x, 
    wgt_t const lwgt, 
    vlbl_t * where)
{
  vtx_t pivot,n,s,i,j;
  coord_t pv;
  wgt_t left, sum;
  coord_t * mx;
  wgt_t * mwgt = NULL;

  mx = coord_duplicate(x,nvtxs);
  if (vwgt) {
    mwgt = wgt_duplicate(vwgt,nvtxs);
  }
  n = nvtxs;
  s = 0;
  left = 0.0;

  while (1) {
    sum = 0.0;
    j = s; 
    pivot = vtx_rand(s,n);
    
    pv = mx[pivot];
    mx[pivot] = mx[n-1];
    for (i=s;i<n-1;++i) {
      if (mx[i] <= pv) {
        dl_swap(mx[j],mx[i]);
        if (mwgt) {
          dl_swap(mwgt[j],mwgt[i]);
        }
        if (mwgt) {
          sum += mwgt[i];
        } else {
          sum += 1;
        }
        ++j;
      }
    }
    mx[n-1] = mx[j];
    mx[j] = pv;
    if (left+sum == lwgt || s >= n-1) {
      break;
    } else if (left + sum < lwgt) {
      s = j;
      left += sum;
    } else {
      n = j;
    }
  }
  /* j is the median */
  pv = mx[j];

  dl_free(mx);
  if (mwgt) {
    dl_free(mwgt);
  }

  for (i=0;i<nvtxs;++i) {
    if (x[i] < pv) {
      where[i] = 1;
    } else {
      where[i] = 0;
    }
  }
}


void recursive_coordinate_bisection(
    vtx_t const nvtxs, 
    wgt_t const * const vwgt, 
    coord_t const ** coords, 
    size_t ndim, 
    wgt_t const * twgts, 
    vlbl_t const nparts, 
    vlbl_t * const where) 
{
  #ifdef XXX
  size_t d;
  vtx_t v;
  vlbl_t i, lnpart, rnpart;
  wgt_t mincut, cut, lwgt;
  vlbl_t * minwhere;

  /* graph parts */
  vtx_t snvtxs[2];
  adj_t * sxadj[2];
  vtx_t * sadjncy[2];
  wgt_t * svwgt[2];
  wgt_t * sadjwgt[2];
  vtx_t * slbl[2];
  coord_t ** scoords[2];

  minwhere = vlbl_alloc(nvtxs);

  lnpart = nparts/2;
  rnpart = nparts - lnpart;

  lwgt = 0;
  for (i=0;i<lnpart;++i) {
    lwgt += twgts[i];
  }

  mincut = -1;

  for (d=0;d<ndim;++d) {
    coordinate_bisection(nvtxs,vwgt,coords[d],lwgt,where);
    cut = calc_edgecut(nvtxs,xadj,adjncy,adjwgt,where);
    if (mincut == -1 || cut < mincut) {
      vlbl_copy(minwhere,where,nvtxs);
      mincut = cut;
    }
  }

  if (nparts > 2) {
    __split_graph(nvtxs,xadj,adjncy,vwgt,adjwgt,minwhere,coords,ndim,snvtxs, \
        sxadj,sadjncy,svwgt,sadjwgt,scoords,slbl);

    if (lnpart > 1) {
      recursive_coordinate_bisection(snvtxs[0],sxadj[0],sadjncy[0],svwgt[0], \
          sadjwgt[0],(coord_t const **)scoords[0],ndim,twgts,lnpart,minwhere);
      for (v=0;v<snvtxs[0];++v) {
        where[slbl[0][v]] = minwhere[v];
      }
    } else {
      for (v=0;v<snvtxs[1];++v) {
        where[slbl[1][v]] = 0;
      }
    }

    if (rnpart > 1) {
      recursive_coordinate_bisection(snvtxs[1],sxadj[1],sadjncy[1],svwgt[1], \
          sadjwgt[1],(coord_t const **)scoords[1],ndim,twgts+lnpart,rnpart, \
          minwhere);
      for (v=0;v<snvtxs[0];++v) {
        where[slbl[0][v]] = minwhere[v]+lnpart;
      }
    } else {
      for (v=0;v<snvtxs[1];++v) {
        where[slbl[1][v]] = lnpart;
      }
    }

    /* free subgraphs */
    for (i=0;i<2;++i) {
      dl_free(sxadj[i]);
      dl_free(sadjncy[i]);
      dl_free(svwgt[i]);
      dl_free(sadjwgt[i]);
      dl_free(slbl[i]);
      for (d=0;d<ndim;++d) {
        dl_free(scoords[i][d]);
      }
      dl_free(scoords[i]);
    }
  }

  vlbl_copy(where,minwhere,nvtxs);
  dl_free(minwhere);
  #endif
}




#endif

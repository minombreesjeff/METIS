/**
 * @file contract.c
 * @brief Functions for performing contraction
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2015-01-21
 */




#ifndef MTMETIS_CONTACT_C
#define MTMETIS_CONTACT_C




#include "contract.h"
#include "check.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct edge_t {
  vtx_t dst;
  wgt_t wgt;
} edge_t;




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


#define DEF_MASK_SIZE (0x2000)
static uint32_t const MASK_SIZE = DEF_MASK_SIZE;
static uint32_t const MASK = DEF_MASK_SIZE-1;
static vtx_t const MASK_MAX_DEG = DEF_MASK_SIZE >> 3;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLSORT_PREFIX edge
#define DLSORT_TYPE_T edge_t
#define DLSORT_COMPARE(a,b) ((a).dst < (b).dst)
#define DLSORT_STATIC
#include "dlsort_headers.h"
#undef DLSORT_STATIC
#undef DLSORT_COMPARE
#undef DLSORT_TYPE_T
#undef DLSORT_PREFIX






/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


/**
 * @brief Reverse the bits of a key for a given mask. Use for hashing into a
 * secondary hash table to prevent collisions.
 *
 * @param n The key to reverse.
 * @param mask The mask of the bits to be reversed.
 *
 * @return The reversed key.
 */
static inline vtx_t __reverse(
    vtx_t const n, 
    vtx_t const mask)
{
  vtx_t r = vtx_reversebits(n);
  int const mb = vtx_downlog2(mask);
  int const vs = sizeof(vtx_t)*8;
  if (vs >= 2*mb) {
    return (r >> (vs - (2*mb))) & mask;
  } else {
    return r >> (vs - mb);
  }
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Adjust the coarse vertex map given new graph distribution paramters. 
 *
 * @param cmap The coarse vertex map.
 * @param mynvtxs The number of vertices this thread owns.
 * @param olddist The old graph distribution parameters.
 * @param newdist The new graph distribution parameters.
 */
static void __adjust_cmap(
    vtx_t * const cmap, 
    vtx_t const mynvtxs, 
    graphdist_t const olddist, 
    graphdist_t const newdist)
{
  vtx_t i,k;
  tid_t o;

  for (i=0;i<mynvtxs;++i) {
    if (cmap[i] >= olddist.offset) { /* remote vertex */
      k = gvtx_to_lvtx(cmap[i],olddist);
      o = gvtx_to_tid(cmap[i],olddist);
      cmap[i] = lvtx_to_gvtx(k,o,newdist);
    }
  }
}


/**
 * @brief Perform contraction using a dense vector.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param mycnvtxs The number of coarse vertices owned by this thread.
 * @param gmatch The global match array.
 * @param fcmap The first fine vertex for each coarse vertex.
 */
static void __par_contract_DENSE(
    ctrl_t * const ctrl, 
    graph_t * const graph, 
    vtx_t const mycnvtxs, 
    vtx_t const * const * const gmatch, 
    vtx_t const * const fcmap)
{
  adj_t cnedges, l, maxdeg, j, i;
  tid_t o, t;
  vtx_t v, c, cg, k;
  wgt_t ewgt;
  adj_t * table;
  graph_t * cgraph;
  graphdist_t dist;

  tid_t const myid = dlthread_get_id(ctrl->comm);
  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);

  /* make accessing my old graph easy */
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const * const gxadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const gvwgt = (wgt_t const * const *)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const * const *)graph->adjwgt;

  vtx_t const * const * const gcmap = (vtx_t const **)graph->cmap;

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.contraction));
  }

  cgraph = par_graph_setup_coarse(graph,mycnvtxs);

  dist = cgraph->dist;

  __adjust_cmap(graph->cmap[myid],mynvtxs,graph->dist,dist);

  /* count possible edges */
  cnedges = 0;
  maxdeg = 0;
  for (c=0;c<mycnvtxs;++c) {
    v = fcmap[c];
    o = myid;
    l = 0;
    do {
      l += gxadj[o][v+1] - gxadj[o][v];
      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(o == myid && v == fcmap[c]));
    dl_storemax(maxdeg,l);
    cnedges += l;
  }

  adj_t * const mycxadj = cgraph->xadj[myid];
  vtx_t * const mycadjncy = cgraph->adjncy[myid] = vtx_alloc(cnedges);
  wgt_t * const mycvwgt = cgraph->vwgt[myid];
  wgt_t * const mycadjwgt = cgraph->adjwgt[myid] = wgt_alloc(cnedges);

  table = NULL;

  table = adj_init_alloc(NULL_ADJ,graph->gnvtxs);

  cnedges = 0;
  mycxadj[0] = 0;

  dlthread_barrier(ctrl->comm);

  /* set max degree for the coarse graph */
  for (c=0;c<mycnvtxs;++c) {
    cg = lvtx_to_gvtx(c,myid,dist);
    /* initialize the coarse vertex */
    mycvwgt[c] = 0;

    v = fcmap[c];
    o = myid;
    DL_ASSERT_EQUALS(myid,gvtx_to_tid(lvtx_to_gvtx(v,myid,graph->dist),
        graph->dist),"%"PF_TID_T);
    do { 
      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[o][v],dist),"%"PF_VTX_T);

      /* transfer over vertex stuff from v and u */
      mycvwgt[c] += gvwgt[o][v];

      for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
        k = gadjncy[o][j];
        if (k < graph->mynvtxs[o]) {
          t = o;
        } else {
          t = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        k = gcmap[t][k];
        if (gvtx_to_tid(k,dist) == myid) {
          k = gvtx_to_lvtx(k,dist);
        }
        if (k == c || k == cg) {
          /* internal edge */
        } else {
          /* external edge */
          i = table[k];
          ewgt = graph->uniformadjwgt ? 1 : gadjwgt[o][j];
          if (i == NULL_ADJ) {
            /* new edge */
            mycadjncy[cnedges] = k;
            mycadjwgt[cnedges] = ewgt;
            table[k] = cnedges++; 
          } else {
            /* duplicate edge */
            mycadjwgt[i] += ewgt;
          }
        }
      }

      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(myid == o && v == fcmap[c]));

    /* clear the table */
    for (j = cnedges;j > mycxadj[c];) {
      --j;
      k = mycadjncy[j];
      table[k] = NULL_ADJ;
    }

    mycxadj[c+1] = cnedges;
  }

  dl_free(table);

  cgraph->mynedges[myid] = cnedges;

  //graph_readjust_memory(cgraph,adjsize);

  dlthread_barrier(ctrl->comm);
  if (myid == 0) {
    cgraph->nedges = adj_sum(cgraph->mynedges,nthreads);
  }

  par_graph_setup_twgts(cgraph);
  
  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.contraction));
  }

  DL_ASSERT(check_graph(cgraph) == 1, "Bad graph generated in coarsening\n");
}


/**
 * @brief Perform contraction using a single hash table, performing a linear
 * scan on the edge list in the case of a collsion.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param mycnvtxs The number of coarse vertices owned by this thread.
 * @param gmatch The global match array.
 * @param fcmap The first fine vertex for each coarse vertex.
 */
static void __par_contract_CLS(
    ctrl_t * const ctrl, 
    graph_t * const graph, 
    vtx_t const mycnvtxs, 
    vtx_t const * const * const gmatch, 
    vtx_t const * const fcmap)
{
  adj_t cnedges, l, maxdeg, j, i, jj, start;
  tid_t o, t;
  vtx_t v, c, cg, k;
  wgt_t ewgt;
  graph_t * cgraph;
  graphdist_t dist;
  offset_t * htable;

  tid_t const myid = dlthread_get_id(ctrl->comm);
  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);

  /* make accessing my old graph easy */
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const * const gxadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const gvwgt = (wgt_t const * const *)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const * const *)graph->adjwgt;

  vtx_t const * const * const gcmap = (vtx_t const **)graph->cmap;

  /* count possible edges */
  cnedges = 0;
  maxdeg = 0;
  for (c=0;c<mycnvtxs;++c) {
    v = fcmap[c];
    o = myid;
    l = 0;
    do {
      l += gxadj[o][v+1] - gxadj[o][v];
      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(o == myid && v == fcmap[c]));
    dl_storemax(maxdeg,l);
    cnedges += l;
  }

  if (maxdeg > MASK_MAX_DEG) {
    __par_contract_DENSE(ctrl,graph,mycnvtxs,gmatch,fcmap);
    return;
  }

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.contraction));
  }

  cgraph = par_graph_setup_coarse(graph,mycnvtxs);

  dist = cgraph->dist;

  __adjust_cmap(graph->cmap[myid],mynvtxs,graph->dist,dist);

  adj_t * const mycxadj = cgraph->xadj[myid];
  vtx_t * const mycadjncy = cgraph->adjncy[myid] = vtx_alloc(cnedges);
  wgt_t * const mycvwgt = cgraph->vwgt[myid];
  wgt_t * const mycadjwgt = cgraph->adjwgt[myid] = wgt_alloc(cnedges);

  htable = offset_init_alloc(NULL_OFFSET,MASK_SIZE);

  cnedges = 0;
  mycxadj[0] = 0;

  dlthread_barrier(ctrl->comm);

  /* set max degree for the coarse graph */
  for (c=0;c<mycnvtxs;++c) {
    cg = lvtx_to_gvtx(c,myid,dist);
    /* initialize the coarse vertex */
    mycvwgt[c] = 0;

    v = fcmap[c];
    o = myid;
    DL_ASSERT_EQUALS(myid,gvtx_to_tid(lvtx_to_gvtx(v,myid,graph->dist),
        graph->dist),"%"PF_TID_T);
    start = cnedges;
    do {
      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[o][v],dist),"%"PF_VTX_T);

      /* transfer over vertex stuff from v and u */
      mycvwgt[c] += graph->uniformvwgt ? 1 : gvwgt[o][v];

      for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
        k = gadjncy[o][j];
        if (k < graph->mynvtxs[o]) {
          t = o;
        } else {
          t = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        k = gcmap[t][k];
        if (gvtx_to_tid(k,dist) == myid) {
          k = gvtx_to_lvtx(k,dist);
        }
        if (k == c || k == cg) {
          /* internal edge */
        } else {
          /* external edge */
          l = k&MASK;
          i = htable[l];
          ewgt = graph->uniformadjwgt ? 1 : gadjwgt[o][j];
          if (i == NULL_OFFSET) {
            /* new edge */
            mycadjncy[cnedges] = k;
            mycadjwgt[cnedges] = ewgt;
            htable[l] = (offset_t)(cnedges - start); 
            ++cnedges;
          } else {
            i += start;
            /* search for existing edge */
            for (jj=i;jj<cnedges;++jj) {
              if (mycadjncy[jj] == k) {
                mycadjwgt[jj] += ewgt;
                break;
              }
            }
            if (jj == cnedges) {
              /* we didn't find the edge, so add it */
              mycadjncy[cnedges] = k;
              mycadjwgt[cnedges++] = ewgt;
            }
          }
        }
      }

      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(myid == o && v == fcmap[c]));

    /* clear the htable */
    for (j = cnedges;j > mycxadj[c];) {
      --j;
      k = mycadjncy[j];
      l = (k&MASK);
      htable[l] = NULL_OFFSET;
    }

    mycxadj[c+1] = cnedges;
  }

  dl_free(htable);

  cgraph->mynedges[myid] = cnedges;

  //graph_readjust_memory(cgraph,adjsize);

  dlthread_barrier(ctrl->comm);
  if (myid == 0) {
    cgraph->nedges = adj_sum(cgraph->mynedges,nthreads);
  }

  par_graph_setup_twgts(cgraph);
  
  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.contraction));
  }

  DL_ASSERT(check_graph(cgraph) == 1, "Bad graph generated in coarsening\n");
}


/**
 * @brief Perform contraction by sorting and merging lists.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param mycnvtxs The number of coarse vertices owned by this thread.
 * @param gmatch The global match array.
 * @param fcmap The first fine vertex for each coarse vertex.
 */
static void __par_contract_SORT(
    ctrl_t * const ctrl, 
    graph_t * const graph, 
    vtx_t const mycnvtxs, 
    vtx_t const * const * const gmatch, 
    vtx_t const * const fcmap)
{
  adj_t cnedges, maxdeg, j, l;
  tid_t o, t;
  vtx_t v, c, cg, k, nlst;
  graph_t * cgraph;
  graphdist_t dist;
  edge_t * lst;

  tid_t const myid = dlthread_get_id(ctrl->comm);
  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);

  /* make accessing my old graph easy */
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const * const gxadj = (adj_t const * const *)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const * const *)graph->adjncy;
  wgt_t const * const * const gvwgt = (wgt_t const * const *)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const * const *)graph->adjwgt;

  vtx_t const * const * const gcmap = (vtx_t const **)graph->cmap;

  if (myid == 0) {
    dl_start_timer(&(ctrl->timers.contraction));
  }

  cgraph = par_graph_setup_coarse(graph,mycnvtxs);

  dist = cgraph->dist;

  __adjust_cmap(graph->cmap[myid],mynvtxs,graph->dist,dist);

  /* count possible edges */
  cnedges = 0;
  maxdeg = 0;
  for (c=0;c<mycnvtxs;++c) {
    v = fcmap[c];
    o = myid;
    l = 0;
    do {
      l += gxadj[o][v+1] - gxadj[o][v];
      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(o == myid && v == fcmap[c]));
    dl_storemax(maxdeg,l);
    cnedges += l;
  }

  adj_t * const mycxadj = cgraph->xadj[myid];
  vtx_t * const mycadjncy = cgraph->adjncy[myid] = vtx_alloc(cnedges);
  wgt_t * const mycvwgt = cgraph->vwgt[myid];
  wgt_t * const mycadjwgt = cgraph->adjwgt[myid] = wgt_alloc(cnedges);

  lst = malloc(sizeof(edge_t)*maxdeg);

  cnedges = 0;
  mycxadj[0] = 0;

  dlthread_barrier(ctrl->comm);

  /* set max degree for the coarse graph */
  for (c=0;c<mycnvtxs;++c) {
    cg = lvtx_to_gvtx(c,myid,dist);
    /* initialize the coarse vertex */
    mycvwgt[c] = 0;

    v = fcmap[c];
    o = myid;
    nlst = 0;
    DL_ASSERT_EQUALS(myid,gvtx_to_tid(lvtx_to_gvtx(v,myid,graph->dist),
        graph->dist),"%"PF_TID_T);
    do { 
      DL_ASSERT_EQUALS(c,gvtx_to_lvtx(gcmap[o][v],dist),"%"PF_VTX_T);

      /* transfer over vertex stuff from v and u */
      mycvwgt[c] += gvwgt[o][v];

      /* add coarse edges it list */
      for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
        k = gadjncy[o][j];
        if (k < graph->mynvtxs[o]) {
          t = o;
        } else {
          t = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        k = gcmap[t][k];
        if (gvtx_to_tid(k,dist) == myid) {
          k = gvtx_to_lvtx(k,dist);
        }
        if (k == c || k == cg) {
          /* internal edge -- ignore */
        } else {
          lst[nlst].dst = k;
          lst[nlst].wgt = gadjwgt[o][j];
          ++nlst;
        }
      }

      v = gmatch[o][v];
      if (v >= graph->mynvtxs[o]) {
        o = gvtx_to_tid(v,graph->dist);
        v = gvtx_to_lvtx(v,graph->dist);
      }
    } while (!(myid == o && v == fcmap[c]));

    if (nlst > 0) {
      /* sort and process edges */
      edge_quicksort(lst,nlst);

      /* add first edge */
      --nlst;
      mycadjncy[cnedges] = lst[nlst].dst;
      mycadjwgt[cnedges] = lst[nlst].wgt;
      ++cnedges;

      /* process the rest */
      while (nlst > 0) {
        --nlst;
        if (mycadjncy[cnedges-1] == lst[nlst].dst) {
          /* combine edges */
          mycadjwgt[cnedges-1] += lst[nlst].wgt;
        } else {
          /* add new edge */
          mycadjncy[cnedges] = lst[nlst].dst;
          mycadjwgt[cnedges] = lst[nlst].wgt;
          ++cnedges;
        }
      }
    }

    mycxadj[c+1] = cnedges;
  }

  dl_free(lst);

  cgraph->mynedges[myid] = cnedges;

  //graph_readjust_memory(cgraph,adjsize);

  dlthread_barrier(ctrl->comm);
  if (myid == 0) {
    cgraph->nedges = adj_sum(cgraph->mynedges,nthreads);
  }

  par_graph_setup_twgts(cgraph);
  
  if (myid == 0) {
    dl_stop_timer(&(ctrl->timers.contraction));
  }

  DL_ASSERT(check_graph(cgraph) == 1, "Bad graph generated in coarsening\n");
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void par_contract_graph(
    ctrl_t * const ctrl, 
    graph_t * const graph, 
    vtx_t const mycnvtxs, 
    vtx_t const * const * const gmatch, 
    vtx_t const * const fcmap)
{
  switch (ctrl->contype) {
    case MTMETIS_CONTYPE_CLS:
      __par_contract_CLS(ctrl,graph,mycnvtxs,gmatch,fcmap);
      break;
    case MTMETIS_CONTYPE_DENSE:
      __par_contract_DENSE(ctrl,graph,mycnvtxs,gmatch,fcmap);
      break;
    case MTMETIS_CONTYPE_SORT:
      __par_contract_SORT(ctrl,graph,mycnvtxs,gmatch,fcmap);
      break;
    default:
      dl_error("Unknown contraction type '%d'\n",ctrl->contype);
  }
}




#endif

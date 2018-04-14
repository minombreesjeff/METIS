/**
 * @file coarsen.c
 * @brief Coarsening functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */




#ifndef MTMETIS_COARSEN_C
#define MTMETIS_COARSEN_C




#include "coarsen.h"
#include "check.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLSORTKV_PREFIX vv
#define DLSORTKV_KEY_T vtx_t
#define DLSORTKV_VAL_T vtx_t
#define DLSORTKV_STATIC
#include "dlsortkv_headers.h"
#undef DLSORTKV_STATIC
#undef DLSORTKV_VAL_T
#undef DLSORTKV_KEY_T
#undef DLSORTKV_PREFIX




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


/**
* @brief Collapse an without a hash table.
*
* @param k The endpoint vertex. 
* @param ewgt The edge weight.
* @param cadjncy The coarse adjacency list.
* @param cadjwgt The coarse edge weight.
* @param htable The hashtable.
* @param nedges The number of coarse edges.
*
* @return The new number of coarse edges.
*/
static inline adj_t __collapse_edge_nomask(
    vtx_t const k,
    wgt_t const ewgt,
    vtx_t * const cadjncy,
    wgt_t * const cadjwgt,
    adj_t * const htable,
    adj_t nedges)
{
  adj_t m;

  if ((m = htable[k]) == NULL_ADJ) {
    cadjncy[nedges] = k;
    cadjwgt[nedges] = ewgt;
    htable[k] = nedges++;
  } else {
    cadjwgt[m] += ewgt;
  }

  return nedges;
}


/**
 * @brief Collapse an edge using a hash table
 *
 * @param start The start of the edges for this vertex
 * @param k The endpoint of the current edge.
 * @param ewgt The edge weight.
 * @param mask The hashtable mask.
 * @param cadjncy The coarse adjacency list.
 * @param cadjwgt The coarse edge weights.
 * @param htable The hash table for collapsing edges.
 * @param nedges The current number of coarse edges.
 *
 * @return The new number of coarse edges.
 */
static inline adj_t __collapse_edge_mask(
    adj_t const start,
    vtx_t const k,
    wgt_t const ewgt,
    vtx_t const mask,
    vtx_t * const cadjncy,
    wgt_t * const cadjwgt,
    adj_t * const htable,
    adj_t nedges)
{
  vtx_t kk;
  adj_t jj, m;

  kk = k&mask;
  if ((m = htable[kk]) == NULL_ADJ) {
    cadjncy[nedges] = k;
    cadjwgt[nedges] = ewgt;
    htable[kk]      = nedges++;
  } else if (cadjncy[m] == k) {
    cadjwgt[m] += ewgt;
  } else {
    for (jj=start; jj<nedges; ++jj) {
      if (cadjncy[jj] == k) {
        cadjwgt[jj] += ewgt;
        break;
      }
    }
    if (jj == nedges) {
      cadjncy[nedges]   = k;
      cadjwgt[nedges++] = ewgt;
    }
  }
  return nedges;
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
* @brief Create a coarse graph given a matching using hash table to identify
*   edges to be merged
*
* @param ctrl The control structure.
* @param graph The fine graph.
* @param cnvtxs The number of coarse vertices in the coarse graph.
* @param match The matchings of vertices (match[match[v]] = v).
* @param fcmap The mapping of the coarse vertex number to the lowest number
*   fine vertex in the coarse vertex.
*/
static void __coarsen_contract_graph(
    ctrl_t * const ctrl, 
    graph_t * const graph, 
    vtx_t const cnvtxs, 
    vtx_t const * const * const gmatch, 
    vtx_t const * const fcmap)
{
  vtx_t mask;
  vtx_t i, k, pv, pi, v, u, lv, lu, lvtx, nbrid, deg, maxdeg, c;
  tid_t utid;
  adj_t j, istart, iend, nedges, adjsize;
  adj_t * htable;
  vtx_t * cadjncy;
  wgt_t * cadjwgt;

  tid_t const myid = omp_get_thread_num();
  tid_t const nthreads = omp_get_num_threads();

  vtx_t const * const * const gcmap = (vtx_t const **)graph->cmap;

  /* fine graph pointers */
  vtx_t const mynvtxs = graph->mynvtxs[myid];

  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;
  wgt_t const * const * const gvwgt = (wgt_t const **)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const **)graph->adjwgt;

  graphdist_t const dist = graph->dist;

  /* fine graph parts */
  adj_t const * xadj = gxadj[myid];
  wgt_t const * vwgt = gvwgt[myid];
  vtx_t const * adjncy = gadjncy[myid];
  wgt_t const * adjwgt = gadjwgt[myid];

  /* implicit barrier */
  graph_t * const cgraph = graph_setup_coarse(graph, cnvtxs);

  /* coarse graph pointers */
  adj_t ** const gcxadj = cgraph->xadj;
  vtx_t ** const gcadjncy = cgraph->adjncy;
  wgt_t ** const gcvwgt = cgraph->vwgt;
  wgt_t ** const gcadjwgt = cgraph->adjwgt;

  /* coarse graph parts */
  adj_t * const cxadj = gcxadj[myid];
  wgt_t * const cvwgt = gcvwgt[myid];

  graphdist_t const cdist = cgraph->dist;

  DL_ASSERT(check_graph(graph),"Invalid graph");
  DL_ASSERT_EQUALS(cnvtxs,cgraph->mynvtxs[myid],"%"PF_VTX_T);

  #pragma omp master
  {
    dl_start_timer(&(ctrl->timers.contraction));
  }

  dprintf("[%"PF_TID_T"] Contracting %"PF_VTX_T" vertices to %"PF_VTX_T \
      " vertices\n",myid,mynvtxs,cnvtxs);
  
  maxdeg = adjsize = deg = 0;
  /* initialize the mmap */
  for (pi=0;pi<cnvtxs;++pi) {
    pv = lvtx_to_gvtx(pi,myid,cdist);
    v = fcmap[pi];
    lv = gvtx_to_lvtx(v,dist);
    DL_ASSERT(myid == gvtx_to_tid(v,cdist),"Vertex not owned by " \
        "processing thread");
    u = gmatch[myid][lv];
    if (v != u) {
      utid = gvtx_to_tid(u,dist);
      lu = gvtx_to_lvtx(u,dist);
      DL_ASSERT(gmatch[utid][lu] == v,"Bad match vector");
      DL_ASSERT(gcmap[utid][lu] == pv,"Bad cmap");
      deg = xadj[lv+1] - xadj[lv] + gxadj[utid][lu+1] - gxadj[utid][lu];
    } else {
      deg = xadj[lv+1] - xadj[lv];
    }
    if (deg > maxdeg) {
      maxdeg = deg;
    }
    adjsize += deg;
  }

  mask = 0x1000;

  /* adaptive */
  while (maxdeg > (mask >> 3) || graph->nedges/graph->nvtxs > mask/20) {
    if (cgraph->nvtxs < mask*2) {
      mask = 0;
      break;
    } else {
      mask <<= 1;
    }
  }

  cadjncy = gcadjncy[myid] = vtx_alloc(adjsize);
  cadjwgt = gcadjwgt[myid] = wgt_alloc(adjsize);

  maxdeg = nedges = 0; 
  if (mask) {
    dprintf("[%"PF_TID_T"] Using mask of %"PF_VTX_T"\n",myid,mask);
    htable = vtx_init_alloc(NULL_ADJ,mask); 
    --mask; /* make 1000 -> 0FFF */
    for (pi=0; pi<cnvtxs; ++pi) {
      pv = lvtx_to_gvtx(pi,myid,dist);

      cxadj[pi] = nedges;
      v = fcmap[pi];

      i = gvtx_to_lvtx(v,dist);

      DL_ASSERT_EQUALS(myid,gvtx_to_tid(v,dist),"%"PF_TID_T);
      DL_ASSERT_EQUALS(lvtx_to_gvtx(i,myid,dist),v,"%"PF_VTX_T);
      DL_ASSERT_EQUALS(gcmap[myid][i],pv,"%"PF_VTX_T);
      DL_ASSERT(i<mynvtxs,"Invalid local vertex number %" \
          PF_VTX_T"/%"PF_VTX_T,i,mynvtxs);

      istart = xadj[i];
      iend = xadj[i+1];
      for (j=istart; j<iend; ++j) {
        k = adjncy[j];
        if (k < mynvtxs) {
          nbrid = myid;
          lvtx = k;
        } else {
          nbrid = gvtx_to_tid(k,dist);
          lvtx = gvtx_to_lvtx(k,dist);
        }
        c = gcmap[nbrid][lvtx];
        if (c == pv) {
          continue;
        } else if (gvtx_to_tid(c,dist) == myid) {
          c = gvtx_to_lvtx(c,dist);
        }
        DL_ASSERT(c < cgraph->gnvtxs,"Invalid global vertex number: %" \
            PF_VTX_T"/%"PF_VTX_T,c,cgraph->gnvtxs);

        nedges = __collapse_edge_mask(cxadj[pi],c,adjwgt[j],mask,cadjncy, \
            cadjwgt,htable,nedges);
      }
      cvwgt[pi] = vwgt[i];

      if ((u = gmatch[myid][i]) != v) {
        lu = gvtx_to_lvtx(u,dist);
        utid = gvtx_to_tid(u,dist);
        DL_ASSERT(gcmap[utid][lu] == gcmap[myid][i],"Bad cmap");
        DL_ASSERT(gmatch[utid][lu] == v,"Bad match");
        /* combine u and v stuff */
        cvwgt[pi] += gvwgt[utid][lu];
        istart = gxadj[utid][lu];
        iend = gxadj[utid][lu+1];
        for (j=istart; j<iend; ++j) {
          k = gadjncy[utid][j];
          if (k < graph->mynvtxs[utid]) {
            nbrid = utid;
            lvtx = k;
          } else {
            nbrid = gvtx_to_tid(k,dist);
            lvtx = gvtx_to_lvtx(k,dist);
          }
          c = gcmap[nbrid][lvtx];
          if (c == pv) {
            continue;
          } else if (gvtx_to_tid(c,dist) == myid) {
            c = gvtx_to_lvtx(c,dist);
          }
          DL_ASSERT(c < cgraph->gnvtxs,"Invalid global vertex number: %" \
              PF_VTX_T"/%"PF_VTX_T,c,cgraph->gnvtxs);

          nedges = __collapse_edge_mask(cxadj[pi],c,gadjwgt[utid][j],mask, \
              cadjncy,cadjwgt,htable,nedges);
        }
      }

      /* Zero out the htable */
      for (j=cxadj[pi]; j<nedges; j++) {
        k = cadjncy[j];
        htable[k&mask] = NULL_ADJ;  
      }

      cxadj[pi+1] = nedges;

      if (cxadj[pi+1] - cxadj[pi] > maxdeg) {
         maxdeg = cxadj[pi+1] - cxadj[pi];
      }
    }
    dl_free(htable);
  } else {
    dprintf("[%"PF_TID_T"] Using no-mask with htable of size %"PF_VTX_T"\n", \
        myid,cgraph->gnvtxs);
    htable = vtx_init_alloc(NULL_ADJ,cgraph->gnvtxs); 
    for (pi=0;pi<cnvtxs;++pi) {
      pv = lvtx_to_gvtx(pi,myid,cdist);

      cxadj[pi] = nedges;
      v = fcmap[pi];

      i = gvtx_to_lvtx(v,dist);

      istart = xadj[i];
      iend = xadj[i+1];
      for (j=istart; j<iend; ++j) {
        k = adjncy[j];
        if (k < mynvtxs) {
          nbrid = myid;
          lvtx = k;
        } else {
          nbrid = gvtx_to_tid(k,dist);
          lvtx = gvtx_to_lvtx(k,dist);
        }
        c = gcmap[nbrid][lvtx];
        if (c == pv) {
          continue;
        } else if (gvtx_to_tid(c,dist) == myid) {
          c = gvtx_to_lvtx(c,dist);
        }
        DL_ASSERT(k < cgraph->gnvtxs,"Invalid global vertex number: %" \
            PF_VTX_T"/%"PF_VTX_T,c,cgraph->gnvtxs);

        nedges = __collapse_edge_nomask(c,adjwgt[j],cadjncy,cadjwgt,htable, \
          nedges);
      }
      cvwgt[pi] = vwgt[i];

      if ((u = gmatch[myid][i]) != v) {
        lu = gvtx_to_lvtx(u,dist);
        utid = gvtx_to_tid(u,dist);
        DL_ASSERT(gcmap[utid][lu] == gcmap[myid][i],"Bad cmap");
        DL_ASSERT(gmatch[utid][lu] == v,"Bad match");
        /* combine u and v stuff */
        cvwgt[pi] += gvwgt[utid][lu];
        istart = gxadj[utid][lu];
        iend = gxadj[utid][lu+1];
        for (j=istart; j<iend; ++j) {
          k = gadjncy[utid][j];
          if (k < graph->mynvtxs[utid]) {
            nbrid = utid;
            lvtx = k;
          } else {
            nbrid = gvtx_to_tid(k,dist);
            lvtx = gvtx_to_lvtx(k,dist);
          }
          c = gcmap[nbrid][lvtx];
          if (c == pv) {
            continue;
          } else if (gvtx_to_tid(c,dist) == myid) {
            c = gvtx_to_lvtx(c,dist);
          }

          nedges = __collapse_edge_nomask(c,gadjwgt[utid][j],cadjncy,cadjwgt, \
              htable,nedges);
        }
      }

      /* Zero out the htable and find the heaveist edge */
      for (j=cxadj[pi]; j<nedges; j++) {
        k = cadjncy[j];
        htable[k] = NULL_ADJ;
      }

      cxadj[pi+1] = nedges;
      
      if (cxadj[pi+1] - cxadj[pi] > maxdeg) {
         maxdeg = cxadj[pi+1] - cxadj[pi];
      }
    }
    dl_free(htable);
  }

  cgraph->mynedges[myid] = nedges;

  graph_readjust_memory(cgraph,adjsize);

  #pragma omp barrier
  #pragma omp master
  {
    cgraph->nedges = adj_sum(cgraph->mynedges,nthreads);
  }

  graph_setup_twgts(cgraph);
  #pragma omp barrier
  #pragma omp master
  {
    dl_stop_timer(&(ctrl->timers.contraction));
  }

  DL_ASSERT(check_graph(cgraph),"Bad graph");
}


static vtx_t __coarsen_match_RM(
    ctrl_t * const ctrl, 
    graph_t const * const graph,
    vtx_t * const * const gmatch, 
    vtx_t * const fcmap) 
{
  vtx_t i, pi, k, maxidx, last_unmatched, seed,nbrid,lvtx,gvtx,cnvtxs;
  adj_t j;
  vtx_t * match, * perm, * cmap;

  tid_t const myid = omp_get_thread_num();

  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  wgt_t const * const * const gvwgt = (wgt_t const **)graph->vwgt;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;
  vtx_t ** const gcmap = graph->cmap;

  wgt_t const maxvwgt = ctrl->maxvwgt;

  /* local graph pointers */
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = gxadj[myid];
  wgt_t const * const vwgt = gvwgt[myid];
  vtx_t const * const adjncy = gadjncy[myid];

  DL_ASSERT_EQUALS(graph->dist.nthreads,(tid_t)omp_get_num_threads(), \
      "%"PF_TID_T);

  #pragma omp master
  {
    dl_start_timer(&(ctrl->timers.matching));
  }

  perm = vtx_alloc(mynvtxs);

  k = 0;
  seed = ctrl->seed + myid;

  vtx_incset(perm,0,1,mynvtxs);
  vtx_pseudo_shuffle_r(perm,mynvtxs/8,mynvtxs,&seed);

  last_unmatched=0; /* this is private but ... */

  match = gmatch[myid];

  for (pi=0; pi<mynvtxs;++pi) {
    /* request my matches */
    i = perm[pi];

    if (match[i] == NULL_VTX) { /* Unmatched */
      gvtx = lvtx_to_gvtx(i,myid,graph->dist);
      maxidx = gvtx;
      
      if (vwgt[i] < maxvwgt) {
        /* Deal with island vertices. Match locally */
        if (xadj[i+1] == xadj[i]) { 
          last_unmatched = dl_max(pi, last_unmatched)+1;
          for (; last_unmatched<mynvtxs; last_unmatched++) {
            k = perm[last_unmatched];
            if (match[k] == NULL_VTX) {
              maxidx = lvtx_to_gvtx(k,myid,graph->dist);
              break;
            }
          }
        } else {
          /* Find a heavy-edge matching, subject to maxvwgt constraints */
          for (j=xadj[i]; j<xadj[i+1]; ++j) {
            k = adjncy[j];
            if (k < mynvtxs) {
              lvtx = k;
              nbrid = myid;
            } else {
              nbrid = gvtx_to_tid(k,graph->dist);
              lvtx = gvtx_to_lvtx(k,graph->dist);
            }
            if (vwgt[i]+gvwgt[nbrid][lvtx] <= maxvwgt && 
                (gmatch[nbrid][lvtx] == NULL_VTX)) {
              if (nbrid == myid) {
                /* I own it, lets go with it */
                maxidx = k;
                break;
              } else if (maxidx == gvtx) {
                maxidx = k;
              }
            }
          }
        }
      }
      if (maxidx < mynvtxs) {
        nbrid = myid;
        lvtx = maxidx;
        maxidx = lvtx_to_gvtx(maxidx,myid,graph->dist);
      } else {
        nbrid = gvtx_to_tid(maxidx,graph->dist);
        lvtx = gvtx_to_lvtx(maxidx,graph->dist);
      }
      /* handle unmatched vertices later -- and orphaned vertices */
      if (gvtx < maxidx) {
        match[i] = maxidx;
        gmatch[nbrid][lvtx] = gvtx;
      } else if (gvtx > maxidx) {
        gmatch[nbrid][lvtx] = gvtx;
        match[i] = maxidx;
      }
    }
  } /* outer match loop */
  #pragma omp barrier

  cmap = gcmap[myid] = perm;  /* re-use perm */

  cnvtxs = 0;
  for (i=0;i<mynvtxs;++i) {
    gvtx = lvtx_to_gvtx(i,myid,graph->dist);
    DL_ASSERT(gvtx_to_lvtx(gvtx,graph->dist) == i,"Mismatch local-global " \
        "vertex id");
    if (match[i] == NULL_VTX) {
      /* match any unmatched vertices with themselves */
      match[i] = gvtx;
    } else {
      nbrid = gvtx_to_tid(match[i],graph->dist);
      lvtx = gvtx_to_lvtx(match[i],graph->dist);
      if (gmatch[nbrid][lvtx] != gvtx) {
        /* match vertices in broken matches with themselves */
        match[i] = gvtx;
      }
    }

    if (MY_CVTX(gvtx,match[i],xadj[i+1]-xadj[i],
          gxadj[nbrid][lvtx+1]-gxadj[nbrid][lvtx])) {
      /* use global cmap[i] id */
      cmap[i] = lvtx_to_gvtx(cnvtxs,myid,graph->dist);
      fcmap[cnvtxs++] = gvtx;
    }
  }

  #pragma omp barrier
  /* tries to avoid false sharing */
  for (i=mynvtxs;i>0;) {
    --i;
    gvtx = lvtx_to_gvtx(i,myid,graph->dist);
    nbrid = gvtx_to_tid(match[i],graph->dist);
    lvtx = gvtx_to_lvtx(match[i],graph->dist);
    if (!MY_CVTX(gvtx,match[i],xadj[i+1]-xadj[i],
          gxadj[nbrid][lvtx+1]-gxadj[nbrid][lvtx])) { 
      cmap[i] = gcmap[nbrid][lvtx];
    }
  }

  #pragma omp barrier
  #pragma omp master
  {
    dl_stop_timer(&(ctrl->timers.matching));
  }

  return cnvtxs;
}


/**
* @brief Match the vertices in a graph using the SHEM coarsening scheme in
*   parallel and create a coarse graph and store it in graph->coarse
*
* @param ctrl THe control structure to use
* @param graph The fine graph to coarsen
* @param nthreads The number of vertices to use
*
* @return The number of vertices in the coarse graph 
*/
static vtx_t __coarsen_match_SHEM(
    ctrl_t * const ctrl, 
    graph_t const * const graph,
    vtx_t * const * const gmatch, 
    vtx_t * const fcmap) 
{
  unsigned int seed;
  vtx_t cnvtxs, i, pi, k, maxidx, maxwgt, last_unmatched, \
      lvtx, gvtx, lcl;
  wgt_t mywgt;
  tid_t nbrid;
  adj_t j, avgdegree;
  vtx_t * perm, * tperm, * cmap, * degrees;

  tid_t const myid = omp_get_thread_num();

  vtx_t ** const gcmap = graph->cmap;

  wgt_t const maxvwgt  = ctrl->maxvwgt;

  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;
  wgt_t const * const * const gvwgt = (wgt_t const **)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const **)graph->adjwgt;

  /* thread local graph pointers */
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = gxadj[myid];
  vtx_t const * const adjncy = gadjncy[myid];
  wgt_t const * const vwgt = gvwgt[myid];
  wgt_t const * const adjwgt = gadjwgt[myid];
  vtx_t * const match = gmatch[myid];

  #pragma omp master
  {
    dl_start_timer(&(ctrl->timers.matching));
  }

  /* matching vectors */
  perm = vtx_alloc(mynvtxs);
  tperm = vtx_alloc(mynvtxs);
  degrees = vtx_alloc(mynvtxs);

  k = 0;
  cnvtxs = 0;

  /* calculate the degree of each vertex, truncating to the average */
  avgdegree = (0.7*(graph->mynedges[myid]/mynvtxs))+1;
  for (i=0;i<mynvtxs;++i) {
    j = xadj[i+1] - xadj[i];
    degrees[i] = (j > avgdegree ? avgdegree : j);
  }

  /* create a pre-permutation array */
  vtx_incset(tperm,0,1,mynvtxs);

  /* shuffle permutation array and degrees the same */
  seed = ctrl->seed + myid;
  vtx_pseudo_shuffle_r(tperm,mynvtxs/8,mynvtxs,&seed);
  seed = ctrl->seed + myid;
  vtx_pseudo_shuffle_r(degrees,mynvtxs/8,mynvtxs,&seed);

  DL_ASSERT_EQUALS(degrees[0], \
      dl_min(avgdegree,xadj[tperm[0]+1]-xadj[tperm[0]]),"%"PF_ADJ_T);

  /* create permutation */
  vv_countingsort_kv(degrees,tperm,0,avgdegree,mynvtxs,perm,NULL);

  DL_ASSERT(mynvtxs < 2 || xadj[perm[0]+1] - xadj[perm[0]] <= \
      xadj[perm[mynvtxs-1]+1] - xadj[perm[mynvtxs-1]],"Sorting failed\n");

  /* free scratch space */
  dl_free(tperm);
  dl_free(degrees);

  last_unmatched=0; 

  for (pi=0; pi<mynvtxs;++pi) {
    /* request my matches */
    i = perm[pi];
      
    if (match[i] == NULL_VTX) {  /* Unmatched */
      mywgt = vwgt[i];
      maxwgt = 0;
      gvtx = lvtx_to_gvtx(i,myid,graph->dist);
      maxidx = gvtx;
      lcl = 0;

      if (mywgt < maxvwgt) {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i+1] == xadj[i]) { 
          last_unmatched = dl_max(pi, last_unmatched)+1;
          for (; last_unmatched<mynvtxs; last_unmatched++) {
            k = perm[last_unmatched];
            if (match[k] == NULL_VTX) {
              maxidx = k;
              break;
            }
          }
        } else {
          /* Find a heavy-edge matching, subject to maxvwgt constraints */
          for (j=xadj[i]; j<xadj[i+1]; ++j) {
            k = adjncy[j];
            if (k < mynvtxs) {
              lvtx = k;
              nbrid = myid;
              lcl = adjwgt[j] * LWGT;
            } else {
              nbrid = gvtx_to_tid(k,graph->dist);
              lvtx = gvtx_to_lvtx(k,graph->dist);
              lcl = adjwgt[j];
            }
            if (maxwgt < lcl && mywgt+gvwgt[nbrid][lvtx] <= maxvwgt && \
                gmatch[nbrid][lvtx] == NULL_VTX) {
              maxidx = k;
              maxwgt = lcl;
            }
          }
        }
      }
      /* handle unmatched vertices later */
      if (maxidx < mynvtxs) {
        nbrid = myid;
        lvtx = maxidx;
        maxidx = lvtx_to_gvtx(maxidx,myid,graph->dist);
      } else {
        nbrid = gvtx_to_tid(maxidx,graph->dist);
        lvtx = gvtx_to_lvtx(maxidx,graph->dist);
      }
      if (gvtx < maxidx) {
        match[i] = maxidx;
        gmatch[nbrid][lvtx] = gvtx;
      } else if (gvtx > maxidx) {
        gmatch[nbrid][lvtx] = gvtx;
        match[i] = maxidx;
      }
    }
  } /* outer match loop */
  #pragma omp barrier

  cmap = gcmap[myid] = perm;

  cnvtxs = 0;
  /* match any unmatched vertices with themselves */
  for (i=0;i<mynvtxs;++i) {
    gvtx = lvtx_to_gvtx(i,myid,graph->dist);
    DL_ASSERT(gvtx_to_lvtx(gvtx,graph->dist) == i,"Mismatch local-global " \
        "vertex id");
    if (match[i] == NULL_VTX) {
      match[i] = gvtx;
    } else {
      nbrid = gvtx_to_tid(match[i],graph->dist);
      lvtx = gvtx_to_lvtx(match[i],graph->dist);
      if (gmatch[nbrid][lvtx] != gvtx) {
        match[i] = gvtx;
      }
    }
    if (MY_CVTX(gvtx,match[i],xadj[i+1]-xadj[i],
          gxadj[nbrid][lvtx+1]-gxadj[nbrid][lvtx])) {
      /* use global cmap[i] id */
      cmap[i] = lvtx_to_gvtx(cnvtxs,myid,graph->dist);
      fcmap[cnvtxs++] = gvtx;
    }
  }

  #pragma omp barrier
  /* tries to avoid false sharing */
  for (i=mynvtxs;i>0;) {
    --i;
    gvtx = lvtx_to_gvtx(i,myid,graph->dist);
    nbrid = gvtx_to_tid(match[i],graph->dist);
    lvtx = gvtx_to_lvtx(match[i],graph->dist);
    if (!MY_CVTX(gvtx,match[i],xadj[i+1]-xadj[i],
          gxadj[nbrid][lvtx+1]-gxadj[nbrid][lvtx])) { 
      cmap[i] = gcmap[nbrid][lvtx];
    }
  }

  #pragma omp barrier
  #pragma omp master
  {
    dl_stop_timer(&(ctrl->timers.matching));
  }

  return cnvtxs;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


static int * __cg_eqwgts;
static vtx_t ** __cg_match;
graph_t * coarsen_graph(
    ctrl_t * ctrl,
    graph_t * graph)
{
  vtx_t i, cnvtxs;
  wgt_t base;
  int neqewgts;

  tid_t const nthreads = omp_get_num_threads();
  tid_t const myid = omp_get_thread_num();

  vtx_t * const fcmap = vtx_alloc(graph->mynvtxs[myid]);
  vtx_t * const match = vtx_init_alloc(NULL_VTX,graph->mynvtxs[myid]);

  /* ward off compiler warnings */
  cnvtxs = 0;

  #pragma omp master 
  {
    dl_start_timer(&ctrl->timers.coarsening);
  }

  DL_ASSERT(check_graph(graph),"Invalid graph");

  neqewgts=0;

  #pragma omp master
  {
    __cg_match = r_vtx_alloc(nthreads);
    __cg_eqwgts = int_alloc(nthreads);
  } 

  if (ctrl->ctype == MTMETIS_CTYPE_SHEM) {
    if (graph->mynvtxs[myid] > 0) {
      base = graph->adjwgt[myid][0]; 
      for (i=1;i<graph->xadj[myid][graph->mynvtxs[myid]];++i) {
        if (graph->adjwgt[myid][i] != base) {
          neqewgts = 1;
          break;
        }
      }
    }
  }
  #pragma omp barrier
  __cg_eqwgts[myid] = neqewgts;
  __cg_match[myid] = match;
  #pragma omp barrier
  neqewgts = int_sum(__cg_eqwgts,nthreads);

  /* set the maximum allowed coarsest vertex weight */
  ctrl->maxvwgt = 1.5*graph->tvwgt/ctrl->coarsen_to;

  do {
    par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH,"Graph{%zu} has %" \
        PF_VTX_T" vertices, %"PF_ADJ_T" edges, and %"PF_WGT_T" exposed edge " \
        "weight.\n",graph->level,graph->nvtxs,graph->nedges,graph->tadjwgt);

    /* allocate memory for cmap, if it has not already been done due to
       multiple cuts */
    #pragma omp master
    {
      if (graph->cmap == NULL) {
        /* threads need to allocate their own chunk inside the matching 
         * functions */
        graph->cmap = r_vtx_alloc(nthreads);
      }
    }
    #pragma omp barrier

    /* coarsening scheme selection used to go here */
    switch(ctrl->ctype) {
      case MTMETIS_CTYPE_RM:
        cnvtxs = __coarsen_match_RM(ctrl,graph,__cg_match,fcmap);
        break;
      case MTMETIS_CTYPE_SHEM:
        if (!neqewgts) {
          cnvtxs = __coarsen_match_RM(ctrl,graph,__cg_match,fcmap);
        } else {
          cnvtxs = __coarsen_match_SHEM(ctrl,graph,__cg_match,fcmap);
        }
        break;
      default:
        dl_error("Unknown ctype: %d\n",ctrl->ctype);
        break;
    }

    __coarsen_contract_graph(ctrl,graph,cnvtxs,(vtx_t const **)__cg_match, \
        fcmap);

    vtx_set(match,NULL_VTX,cnvtxs);

    graph = graph->coarser;

    neqewgts = 1;
  } while (graph->nvtxs > ctrl->coarsen_to && \
           graph->nvtxs < PAR_COARSEN_FRACTION*graph->finer->nvtxs && \
           graph->nedges > graph->nvtxs/2);

  par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH,"Graph{%zu} has %" \
      PF_VTX_T" vertices, %"PF_ADJ_T" edges, and %"PF_WGT_T" exposed edge " \
      "weight.\n",graph->level,graph->nvtxs,graph->nedges,graph->tadjwgt);

  dl_free(fcmap);
  dl_free(match);

  #pragma omp master
  {
    dl_free(__cg_eqwgts);
    dl_free(__cg_match);
  }

  DL_ASSERT(check_graph(graph),"Invalid graph");

  #pragma omp master
  {
    dl_stop_timer(&ctrl->timers.coarsening);
  }

  return graph;
}




#endif

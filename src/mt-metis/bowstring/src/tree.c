/**
 * @file tree.c
 * @brief Functions for generating trees in a graph
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-08-06
 */




#ifndef BOWSTRING_TREE_C
#define BOWSTRING_TREE_C




#include "tree.h"
#include "graph.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const MIN_PSEUDO_SHUFFLE = 64;
static adj_t const NULL_PARENT_EDGE = (adj_t)-1;
static vtx_t const NULL_PARENT = (vtx_t)-1;




/******************************************************************************
* PRIVATE TYPES ***************************************************************
******************************************************************************/


#define DLISET_PREFIX vtx
#define DLISET_TYPE_T vtx_t
#define DLISET_STATIC
#include "dliset_headers.h"
#undef DLISET_STATIC
#undef DLISET_TYPE_T
#undef DLISET_PREFIX


#define DLDJSET_PREFIX vtx
#define DLDJSET_TYPE_T vtx_t
#define DLDJSET_STATIC
#include "dldjset_headers.h"
#undef DLDJSET_STATIC
#undef DLDJSET_TYPE_T
#undef DLDJSET_PREFIX


#define DLPQ_KEY_T adj_t
#define DLPQ_VAL_T wgt_t
#define DLPQ_PREFIX aw
#define DLPQ_MIN 
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_MIN
#undef DLPQ_PREFIX
#undef DLPQ_KEY_T
#undef DLPQ_VAL_T




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __kruskal_mst(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    int * const adjmask)
{
  vtx_t i,k;
  adj_t j;
  adj_t * radj = NULL;
  vtx_t * alist = NULL;
  vtx_djset_t * set = NULL;
  aw_pq_t * q = NULL;

  adj_t const nedges = xadj[nvtxs];

  DL_ASSERT(check_graph(nvtxs,xadj,adjncy,adjwgt),"Bad graph passed to " \
      "__kruskal_mst()\n");

  set = vtx_djset_create(0,nvtxs);

  /* vertex lookup arrays */
  alist = vtx_alloc(nedges);

  /* build edge lookup */
  radj = adj_alloc(nedges);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  /* sort edges */
  q = aw_pq_create(0,nedges);
  for (i=0;i<nvtxs;++i) {
    vtx_set(alist+xadj[i],i,xadj[i+1]-xadj[i]);
    for (j=xadj[i];j<xadj[i+1];++j) {
      aw_pq_push(adjwgt[j],j,q);
    }
  }

  /* set all edges to 0 in the adjmask */
  int_set(adjmask,0,xadj[nvtxs]);

  while (q->size > 0) {
    j = aw_pq_pop(q);
    i = alist[j];
    k = adjncy[j];
    if (vtx_djset_find(i,set) != vtx_djset_find(k,set)) {
      adjmask[j] = adjmask[radj[j]] = 1;
      vtx_djset_join(i,k,set);
    }
  }

  aw_pq_free(q);
  vtx_djset_free(set);
  dl_free(alist);
  dl_free(radj);
}


static void __kruskal_rst(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    int * const adjmask)
{
  vtx_t i,k;
  adj_t j, lnedges, l;
  adj_t * perm = NULL, * aadj = NULL, * badj = NULL, * radj = NULL;
  vtx_t * alist = NULL, * blist = NULL;
  vtx_djset_t * set = NULL;

  adj_t const nedges = xadj[nvtxs];

  /* setup permutation */
  perm = adj_alloc(nedges/2);
  adj_incset(perm,0,1,nedges/2);
  if (((size_t)nedges)/2 < MIN_PSEUDO_SHUFFLE) {
    adj_shuffle(perm,nedges/2);
  } else {
    adj_pseudo_shuffle(perm,nedges/16,nedges/2);
  }

  radj = adj_alloc(nedges);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  /* set all edges to 0 in the adjmask */
  int_set(adjmask,0,xadj[nvtxs]);

  /* set up adjancy arrays */
  lnedges = 0;
  alist = vtx_alloc(nedges/2);
  aadj = adj_alloc(nedges/2);
  blist = vtx_alloc(nedges/2);
  badj = adj_alloc(nedges/2);
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (i < k) {
        l = perm[lnedges++];
        alist[l] = i;
        blist[l] = k;
        /* add the forward edge */
        aadj[l] = j;
        badj[l] = radj[j];
      }
    }
  }
  dl_free(perm);
  dl_free(radj);

  set = vtx_djset_create(0,nvtxs);
  
  while (lnedges > 0) {
    l = --lnedges;
    i = alist[l];
    k = blist[l];
    if (vtx_djset_find(i,set) != vtx_djset_find(k,set)) {
      adjmask[aadj[l]] = adjmask[badj[l]] = 1;
      vtx_djset_join(i,k,set);
    }
  }

  vtx_djset_free(set);
  dl_free(alist);
  dl_free(blist);
  dl_free(aadj);
  dl_free(badj);
}



/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void build_dfs_tree(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    vtx_t const start, 
    vtx_t * const label, 
    vtx_t * const rank, 
    adj_t * const parent, 
    int * const adjmask)
{
  vtx_t i, k, n, size, cnvtxs;
  adj_t j, r;
  vtx_t * stack;
  vtx_iset_t * vtxs;
  adj_t * nbrs;
  vtx_t * nnbrs;

  adj_t const nedges = xadj[nvtxs];

  vtxs = vtx_iset_create(0,nvtxs);
  nbrs = adj_alloc(nedges);
  nnbrs = vtx_alloc(nvtxs);
  vtx_iset_populate(vtxs);
  stack = vtx_alloc(nvtxs);
  cnvtxs = size = 0;

  /* prep my neighbor lists */
  for (i=0;i<nvtxs;++i) {
    nnbrs[i] = xadj[i+1]-xadj[i];
  }
  for (j=0;j<nedges;++j) {
    nbrs[j] = j;
  }

  /* add start to the tree */
  stack[size++] = start;
  if (label) {
    label[cnvtxs++] = start;
  } else {
    ++cnvtxs;
  }
  
  vtx_iset_remove(start,vtxs);

  /* select my tree */
  i = start;
  if (parent) {
    parent[i] = (adj_t)-1;
  }
  if (rank) {
    rank[i] = 0;
  }

  if (adjmask) {
    /* set all edges to 0 in the adjmask */
    int_set(adjmask,0,xadj[nvtxs]);
  }

  while (cnvtxs < nvtxs) {
    i = stack[size-1];
    /* find an unvisited neighbor */
    k = i;
    while (nnbrs[i] > 0) {
      r = adj_rand(0,nnbrs[i]);
      j = nbrs[xadj[i]+r];
      --nnbrs[i];
      nbrs[xadj[i]+r] = nbrs[xadj[i]+nnbrs[i]];
      n = adjncy[j];
      if (vtx_iset_contains(n,vtxs)) {
        vtx_iset_remove(n,vtxs);
        k = n;
        break;
      }
    } 
    if (k != i) {
      DL_ASSERT(k == adjncy[j],"Mismatch k and j\n");
      if (adjmask) {
        adjmask[j] = 1;
      }
      vtx_iset_remove(k,vtxs);
      /* find the reverse direction edge */
      if (adjmask || parent) {
        for (j=xadj[k];j<xadj[k+1];++j) {
          if (adjncy[j] == i) {
            if (adjmask) {
              adjmask[j] = 1;
            }
            if (parent) {
              parent[k] = j;
            }
            break;
          }
        }
      }
      DL_ASSERT(!parent || adjncy[parent[k]] == i,"Bad parent assignment\n");
      i = k;
      if (rank) {
        rank[i] = size;
      }
      stack[size++] = i;
      if (label) {
        label[cnvtxs++] = i;
      } else {
        ++cnvtxs;
      }
    } else {
      DL_ASSERT(nnbrs[i] == 0,"Moved past "PF_VTX_T" while it still has "
          PF_VTX_T" neighbors\n",i,nnbrs[i]);
      if (size > 1) {
        /* backtrack */
        --size;
      } else {
        DL_ASSERT(stack[0] == i,"Jumping from non-root vertex "PF_VTX_T"\n",i);

        /* jump to a random vertex */
        i = vtx_iset_get(vtx_rand(0,vtxs->size),vtxs);
        vtx_iset_remove(i,vtxs);
        if (parent) {
          parent[i] = (adj_t)-1;
        }
        if (rank) {
          rank[i] = 0;
        }
        stack[size++] = i;
        if (label) {
          label[cnvtxs++] = i;
        } else {
          ++cnvtxs;
        }
      }
    }
  } 
  DL_ASSERT(vtxs->size == 0,"Stopped building tree with vertices left\n");

  dl_free(nbrs);
  dl_free(nnbrs);
  vtx_iset_free(vtxs);
  dl_free(stack);
}


void build_bfs_tree(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    vtx_t const start, 
    vtx_t * const label, 
    vtx_t * const rank, 
    adj_t * const eparent, 
    int * const adjmask)
{
  vtx_t i, k, nq, sq, p;
  adj_t j;
  vtx_t * q;

  int * visited = int_calloc(nvtxs);
  vtx_t * parent = NULL;

  if (!label) {
    q = vtx_alloc(nvtxs);
  } else {
    q = label;
  }
  
  p = NULL_PARENT;
  sq = nq = 0;

  visited[start] = 1;
  q[nq++] = start;

  if (adjmask || eparent) {
    parent = vtx_alloc(nvtxs);
    parent[start] = p;
  }

  if (rank) {
    rank[start] = 0;
  }

  if (adjmask) {
    /* set all edges to 0 in the adjmask */
    int_set(adjmask,0,xadj[nvtxs]);
  }

  while (sq < nvtxs) {
    i = q[sq++];
    if (parent) {
      p = parent[i];
    }
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (!visited[k]) {
        if (parent) {
          parent[k] = i;
        }
        q[nq++] = k;
        visited[k] = 1;
        if (rank) {
          rank[k] = rank[i]+1;
        }
        if (adjmask) {
          adjmask[j] = 1;
        }
      } else if (parent && p == k) {
        if (eparent) {
          eparent[i] = j;
        }
        if (adjmask) {
          adjmask[j] = 1;
        }
      }
    }
  }

  if (parent) {
    dl_free(parent);
  }
  if (!label) {
    dl_free(q);
  }

  dl_free(visited);
}


void build_mst_tree(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    int * const adjmask)
{
  if (adjwgt) {
    __kruskal_mst(nvtxs,xadj,adjncy,adjwgt,adjmask); 
  } else {
    __kruskal_rst(nvtxs,xadj,adjncy,adjmask);
  }
}


void build_rst_tree(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    int * const adjmask)
{
  __kruskal_rst(nvtxs,xadj,adjncy,adjmask);
}




/* un-obscure */
#undef vtx_iset_t 
#undef vtx_iset_create 
#undef vtx_iset_get 
#undef vtx_iset_add 
#undef vtx_iset_remove 
#undef vtx_iset_free 
#undef vtx_djset_t
#undef vtx_djset_create 
#undef vtx_djset_find
#undef vtx_djset_join
#undef vtx_djset_free 



#endif

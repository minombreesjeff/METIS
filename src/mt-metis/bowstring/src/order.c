/**
 * @file order.c
 * @brief Permutation functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-07-18
 */




#ifndef BOWSTRING_ORDER_C
#define BOWSTRING_ORDER_C




#include "order.h"
#include "graph.h"
#include "tree.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLPQ_PREFIX vv
#define DLPQ_KEY_T vtx_t
#define DLPQ_VAL_T vtx_t
#define DLPQ_MIN
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_MIN
#undef DLPQ_VAL_T
#undef DLPQ_KEY_T
#undef DLPQ_PREFIX




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Create a random permutation vector.
 *
 * @param nvtxs The number of vertices to permute.
 * @param perm The permutation vector (output).
 */
static void __perm_RANDOM(
    vtx_t const nvtxs,
    vtx_t * const perm)
{
  vtx_incset(perm,0,1,nvtxs);
  vtx_shuffle(perm,nvtxs);
}


/**
 * @brief Create a post-order DFS permutation of the graph/tree.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacecny list pointer.
 * @param adjncy The adjacency list.
 * @param perm The permutation vector (output).
 */
static void __perm_POST(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    vtx_t * const perm)
{
  vtx_t nparent, i, n, p, k;
  adj_t j;
  int * visited;
  vtx_t * parent;
  adj_t * mark;

  parent = vtx_alloc(nvtxs);
  mark = adj_duplicate(xadj,nvtxs);
  visited = int_init_alloc(0,nvtxs);

  /* post-order traversal */
  n = 0;
  nparent = 0;
  i = 0;
  visited[i] = 1;
  parent[nparent++] = i;
  while (nparent > 0) {
    i = parent[nparent-1];
    if (nparent > 1) {
      p = parent[nparent-2];
    } else {
      p = NULL_VTX;
    }
    for (j=mark[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (p != k && !visited[k]) {
        mark[i] = j+1;
        visited[k] = 1;
        parent[nparent++] = k;
        goto NEXT;
      }
    }
    /* visit */
    perm[n++] = i;
    --nparent;
    NEXT:;
  }

  DL_ASSERT_EQUALS(n,nvtxs,PF_VTX_T);

  dl_free(visited);
  dl_free(mark);
  dl_free(parent);
}


/**
 * @brief Create a pre-order DFS permutation of the graph/tree.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacecny list pointer.
 * @param adjncy The adjacency list.
 * @param perm The permutation vector (output).
 */
static void __perm_DFS(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    vtx_t * const perm)
{
  build_dfs_tree(nvtxs,xadj,adjncy,vtx_rand(0,nvtxs),perm,NULL,NULL,NULL);
}


/**
 * @brief Create a BFS permutation of the graph/tree.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacecny list pointer.
 * @param adjncy The adjacency list.
 * @param perm The permutation vector (output).
 */
static void __perm_BFS(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    vtx_t * const perm)
{
  build_bfs_tree(nvtxs,xadj,adjncy,vtx_rand(0,nvtxs),perm,NULL,NULL,NULL);
}


/**
 * @brief Create a cuthillmckee ordering of the graph.
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param perm The permutation vector (output).
 */
static void __perm_RCM(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    vtx_t * const perm)
{
  vtx_t i,k,d,nordered,sr;
  adj_t j;
  vtx_t * deg;
  vv_pq_t * q, * rem;

  q = vv_pq_create(0,nvtxs);
  rem = vv_pq_create(0,nvtxs);
  /* offset pointer */
  deg = vtx_alloc(nvtxs);

  /* find my lowest degree vertex */
  for (i=0;i<nvtxs;++i) {
    d = 0;
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j]; 
      ++d;
    }
    deg[i] = d;
    vv_pq_push(d,i,rem);
  }

  sr = nordered = 0;

  /* loop through connected components */
  while (rem->size > 0) {
    i = vv_pq_pop(rem);
    perm[nordered++] = i;

    /* perform bfs */
    while (sr < nordered) {
      i = perm[sr++];
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j]; 
        if (vv_pq_contains(k,rem)) {
          /* local non-zero */
          vv_pq_remove(k,rem);
          vv_pq_push(deg[k],k,q);
        }
      }
      /* add rows/vertices in ascending order of local degree */
      while (q->size > 0) {
        k = vv_pq_pop(q);
        perm[nordered++] = k;
      }
    }
  }

  /* reverse the CM ordering */
  for (i=0;i<nvtxs/2;++i) {
    k = perm[i];
    perm[i] = perm[nvtxs-i-1];
    perm[nvtxs-i-1] = k;
  }

  vv_pq_free(q);
  vv_pq_free(rem);

  dl_free(deg);
}



/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void order_graph(
    vtx_t const nvtxs, 
    adj_t * const xadj, 
    vtx_t * const adjncy,
    wgt_t * const vwgt, 
    wgt_t * const adjwgt, 
    vtx_t const * const perm) 
{
  vtx_t i, k, pi;
  adj_t j;
  size_t const nedges = xadj[nvtxs];

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


void order_permutation(
    int const ordering,
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy,
    wgt_t const * const vwgt,
    wgt_t const * const adjwgt,
    vtx_t * const perm)
{
  DL_ASSERT(check_graph(nvtxs,xadj,adjncy,adjwgt),"Bad graph before " \
      "permutation\n");

  switch (ordering) {
    case BOWSTRING_ORDER_BFS:
      __perm_BFS(nvtxs,xadj,adjncy,perm);
      break;
    case BOWSTRING_ORDER_DFS:
      __perm_DFS(nvtxs,xadj,adjncy,perm);
      break;
    case BOWSTRING_ORDER_RANDOM:
      __perm_RANDOM(nvtxs,perm);
      break;
    case BOWSTRING_ORDER_POST:
      __perm_POST(nvtxs,xadj,adjncy,perm);
      break;
    case BOWSTRING_ORDER_RCM:
      __perm_RCM(nvtxs,xadj,adjncy,perm);
      break;
    default:
      dl_error("Unknown/Unimplemented ordering '%d'\n",ordering);
  }

  DL_ASSERT(check_graph(nvtxs,xadj,adjncy,adjwgt),"Bad graph after " \
      "permutation\n");
}








#endif


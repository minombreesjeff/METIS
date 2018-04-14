/**
 * @file dynamic.c
 * @brief Functions for dynamic graph operations.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Dominique LaSalle
 * @version 1
 * @date 2014-10-22
 */




#ifndef BOWSTRING_DYNAMIC_C
#define BOWSTRING_DYNAMIC_C




#include <bowstring_dynamic.h>
#include "base.h"
#include "graph.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static vtx_t const MIN_EDGES = 64ULL/sizeof(vtx_t);




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLMEM_PREFIX dynnode
#define DLMEM_TYPE_T dynnode_t
#define DLMEM_STATIC 1
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Determine the size of the current (or prospective) edge allocation.
 *
 * @param n The number of edges.
 */
static inline vtx_t __edge_size(
    vtx_t n)
{
  return dl_max(MIN_EDGES,vtx_uppow2(n));
}



/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


dyngraph_t * bowstring_dg_create(
    vtx_t const max_nvtxs)
{
  vtx_t i;
  dyngraph_t * graph;

  graph = (dyngraph_t*)malloc(sizeof(dyngraph_t));

  graph->nvtxs = 0;
  graph->nedges = 0;
  graph->max_nvtxs = max_nvtxs;

  /* initialize htable */
  graph->htable = vtx_init_alloc(NULL_VTX,max_nvtxs);

  /* initialize vertices */
  graph->nodes = dynnode_alloc(max_nvtxs); 
  for (i=0;i<max_nvtxs;++i) {
    graph->nodes[i].nadj = NULL_VTX;
  }

  DL_ASSERT(bowstring_dg_check(graph),"Bad graph after creation");

  return graph;
}


void bowstring_dg_addvtx(
    vtx_t const v,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt,
    vtx_t const nadj,
    dyngraph_t * const graph)
{
  vtx_t rnadj;

  ++graph->nvtxs;
  graph->nodes[v].nadj = nadj;

  if (nadj != NULL_VTX) {
    rnadj = __edge_size(nadj);

    /* allocate my edge arrays */
    graph->nodes[v].adj = vtx_alloc(rnadj);
    graph->nodes[v].wgt = wgt_alloc(rnadj);

    if (nadj > 0) {
      /* populate my edge arrays */
      vtx_copy(graph->nodes[v].adj,adjncy,nadj);
      wgt_copy(graph->nodes[v].wgt,adjwgt,nadj);
    }
  }
}


void bowstring_dg_addedge(
    vtx_t const v,
    vtx_t const u,
    wgt_t const w,
    dyngraph_t * const graph)
{
  vtx_t vnadj, unadj;

  DL_ASSERT(bowstring_dg_check(graph),"Bad graph before adding edge");

  vnadj = graph->nodes[v].nadj;
  unadj = graph->nodes[u].nadj;

  if (vnadj == NULL_VTX) {
    bowstring_dg_addvtx(v,NULL,NULL,0,graph);
    vnadj = 0;
  } else if (vnadj == __edge_size(vnadj)) {
    graph->nodes[v].adj = vtx_realloc(graph->nodes[v].adj,vnadj+1);
    graph->nodes[v].wgt = wgt_realloc(graph->nodes[v].wgt,vnadj+1);
  }

  if (unadj == NULL_VTX) {
    bowstring_dg_addvtx(u,NULL,NULL,0,graph);
    unadj = 0;
  } else if (unadj == __edge_size(unadj)) {
    graph->nodes[u].adj = vtx_realloc(graph->nodes[u].adj,unadj+1);
    graph->nodes[u].wgt = wgt_realloc(graph->nodes[u].wgt,unadj+1);
  }

  DL_ASSERT_EQUALS(bowstring_dg_connected(v,u,graph),(wgt_t)(-1.0),PF_WGT_T);

  graph->nodes[v].adj[vnadj] = u;
  graph->nodes[v].wgt[vnadj] = w;
  ++graph->nodes[v].nadj;

  graph->nodes[u].adj[unadj] = v;
  graph->nodes[u].wgt[unadj] = w;
  ++graph->nodes[u].nadj;

  ++graph->nedges;

  DL_ASSERT(bowstring_dg_check(graph),"Bad graph after adding edge "PF_VTX_T \
      ":"PF_VTX_T" of weight "PF_WGT_T,v,u,w);
}


void bowstring_dg_delvtx(
    vtx_t const v,
    dyngraph_t * const graph)
{
  vtx_t i, j, k;
  vtx_t * adj;

  adj = graph->nodes[v].adj;

  dl_free(graph->nodes[v].wgt);

  /* remove edges */
  for (i=0;i<graph->nodes[v].nadj;++i) {
    k = adj[i];
    for (j=0;j<graph->nodes[k].nadj;++j) {
      if (graph->nodes[k].adj[j] == v) {
        --graph->nodes[k].nadj;
        graph->nodes[k].adj[j] = graph->nodes[k].adj[graph->nodes[k].nadj];
        graph->nodes[k].wgt[j] = graph->nodes[k].wgt[graph->nodes[k].nadj];
        break;
      }
    }
  }

  graph->nedges -= graph->nodes[v].nadj;

  /* mark the vertex as not in the graph */
  graph->nodes[v].nadj = NULL_VTX;

  dl_free(adj); 

  --graph->nvtxs;

  DL_ASSERT(bowstring_dg_check(graph),"Bad graph after deleting vertex " \
      PF_VTX_T,v);
}


void bowstring_dg_deledge(
    vtx_t const v,
    vtx_t const u,
    dyngraph_t * const graph)
{
  vtx_t i, k;

  /* remove from v */
  for (i=0;i<graph->nodes[v].nadj;++i) {
    k = graph->nodes[v].adj[i];
    if (k == u) {
      --graph->nodes[v].nadj;
      graph->nodes[v].adj[i] = graph->nodes[v].adj[graph->nodes[v].nadj]; 
      graph->nodes[v].wgt[i] = graph->nodes[v].wgt[graph->nodes[v].nadj]; 
      break;
    }
  }
  DL_ASSERT_EQUALS(k,u,PF_VTX_T);
  
  /* remove from u */
  for (i=0;i<graph->nodes[u].nadj;++i) {
    k = graph->nodes[u].adj[i];
    if (k == v) {
      --graph->nodes[u].nadj;
      graph->nodes[u].adj[i] = graph->nodes[u].adj[graph->nodes[u].nadj]; 
      graph->nodes[u].wgt[i] = graph->nodes[u].wgt[graph->nodes[u].nadj]; 
      break;
    }
  }
  DL_ASSERT_EQUALS(k,v,PF_VTX_T);

  --graph->nedges;

  DL_ASSERT(bowstring_dg_check(graph),"Bad graph after deleting edge " \
      PF_VTX_T":"PF_VTX_T,v,u);
}


void bowstring_dg_updateedge(
    vtx_t const v,
    vtx_t const u,
    wgt_t const w,
    dyngraph_t * const graph)
{
  vtx_t i, k;

  /* remove from v */
  for (i=0;i<graph->nodes[v].nadj;++i) {
    k = graph->nodes[v].adj[i];
    if (k == u) {
      graph->nodes[v].wgt[i] = w; 
      break;
    }
  }
  
  /* remove from u */
  for (i=0;i<graph->nodes[u].nadj;++i) {
    k = graph->nodes[u].adj[i];
    if (k == v) {
      graph->nodes[u].wgt[i] = w; 
      break;
    }
  }

  DL_ASSERT(bowstring_dg_check(graph),"Bad graph after updating edge");
}


wgt_t bowstring_dg_connected(
    vtx_t const v,
    vtx_t const u,
    dyngraph_t const * const graph)
{
  vtx_t i, vnadj, unadj;
  wgt_t w;

  vnadj = graph->nodes[v].nadj;
  unadj = graph->nodes[u].nadj;

  w = -1;

  if (vnadj != NULL_VTX && unadj != NULL_VTX) {
    if (vnadj < unadj) {
      for (i=0;i<vnadj;++i) {
        if (graph->nodes[v].adj[i] == u) {
          w = graph->nodes[v].wgt[i];
          break;
        }
      }
    } else {
      for (i=0;i<unadj;++i) {
        if (graph->nodes[u].adj[i] == v) {
          w = graph->nodes[u].wgt[i];
          break;
        }
      }
    }
  }

  return w;
}


vtx_t bowstring_dg_contract(
    int const method,
    vtx_t const * const vtxs,
    vtx_t const nvtxs,
    dyngraph_t * const graph)
{
  vtx_t i, j, k, m, v, nadj, min, l, lnadj;
  vtx_t * htable, * adj;
  wgt_t * wgt;
  dynnode_t node;

  /* we assume the htable is full of NULL_VTX */
  htable = graph->htable;

  /* decide how many edges I might have */
  nadj = 0;
  for (i=0;i<nvtxs;++i) {
    v = vtxs[i];
    nadj += graph->nodes[v].nadj; 
  }

  min = vtx_min_value(vtxs,nvtxs);

  /* the size of my new arrays */
  nadj = __edge_size(dl_min(graph->nvtxs-nvtxs,nadj));

  adj = vtx_alloc(nadj);
  wgt = wgt_alloc(nadj);

  nadj = 0;
  switch (method) {
    case BOWSTRING_CONTRACT_SUM:
      for (i=0;i<nvtxs;++i) {
        v = vtxs[i];
        node = graph->nodes[v];
        for (j=0;j<node.nadj;++j) {
          k = node.adj[j];
          if ((m = htable[k]) == NULL_VTX) {
            /* add edge */
            adj[nadj] = k;
            wgt[nadj] = node.wgt[j];
            htable[k] = nadj++;
          } else {
            /* update edge */
            wgt[m] += node.wgt[j];
          }
          /* delete the reverse edge */
          for (l=0;l<graph->nodes[k].nadj;++l) {
            if (v == graph->nodes[k].adj[l]) {
              graph->nodes[k].adj[l] = \
                  graph->nodes[k].adj[--graph->nodes[k].nadj];
              break;
            }
          }
        }
      }
      break;
    case BOWSTRING_CONTRACT_MAX:
      for (i=0;i<nvtxs;++i) {
        v = vtxs[i];
        node = graph->nodes[v];
        for (j=0;j<node.nadj;++j) {
          k = node.adj[j];
          if ((m = htable[k]) == NULL_VTX) {
            /* add edge */
            adj[nadj] = k;
            wgt[nadj] = node.wgt[j];
            htable[k] = nadj++;
          } else if (wgt[m] < node.wgt[j]) {
            /* update edge */
            wgt[m] = node.wgt[j];
          }
          /* delete the reverse edge */
          for (l=0;l<graph->nodes[k].nadj;++l) {
            if (v == graph->nodes[k].adj[l]) {
              graph->nodes[k].adj[l] = \
                  graph->nodes[k].adj[--graph->nodes[k].nadj];
              break;
            }
          }
        }
      }
      break;
    case BOWSTRING_CONTRACT_MIN:
      for (i=0;i<nvtxs;++i) {
        v = vtxs[i];
        node = graph->nodes[v];
        for (j=0;j<node.nadj;++j) {
          k = node.adj[j];
          if ((m = htable[k]) == NULL_VTX) {
            /* add edge */
            adj[nadj] = k;
            wgt[nadj] = node.wgt[j];
            htable[k] = nadj++;
          } else if (wgt[m] > node.wgt[j]) {
            /* update edge */
            wgt[m] = node.wgt[j];
          }
          /* delete the reverse edge */
          for (l=0;l<graph->nodes[k].nadj;++l) {
            if (v == graph->nodes[k].adj[l]) {
              graph->nodes[k].adj[l] = \
                  graph->nodes[k].adj[--graph->nodes[k].nadj];
              break;
            }
          }
        }
      }
      break;
    default:
      dl_error("Unknown contraction method '%d'\n",method);
  }

  for (i=0;i<nadj;++i) {
    k = adj[j];

    /* add edges to new node -- we know we have room because we deleted */
    lnadj = graph->nodes[k].nadj++;
    graph->nodes[k].adj[lnadj] = min;
    graph->nodes[k].wgt[lnadj] = wgt[j];

    /* clear htable */
    htable[k] = NULL_VTX;
  }

  /* delete contracted nodes */
  for (i=0;i<nvtxs;++i) {
    v = vtxs[i]; 
    graph->nodes[v].nadj = NULL_VTX;
    dl_free(graph->nodes[v].adj);
    dl_free(graph->nodes[v].wgt);
  }

  /* save new node */
  graph->nodes[min].nadj = nadj;
  graph->nodes[min].adj = adj; 
  graph->nodes[min].wgt = wgt; 

  graph->nvtxs -= nvtxs-1;

  DL_ASSERT(bowstring_dg_check(graph),"Bad graph after contracting vertices");

  return min;
}


void bowstring_dg_tostatic(
    bowstring_dyngraph_t const * const graph,
    bowstring_adj_t ** const r_xadj,
    bowstring_vtx_t ** const r_adjncy,
    bowstring_wgt_t ** const r_adjwgt,
    bowstring_vtx_t ** const r_alias)
{
  vtx_t i, k;
  adj_t nedges, j;
  vtx_t * adjncy, * alias, * rename;
  adj_t * xadj;
  wgt_t * adjwgt;
  bowstring_dynnode_t node;
  
  xadj = adj_alloc(graph->nvtxs+1);
  alias = vtx_alloc(graph->nvtxs);
  rename = vtx_alloc(graph->max_nvtxs);

  /* setup xadj and count edges */
  xadj[0] = 0;
  k = 0;
  for (i=0;i<graph->max_nvtxs;++i) {
    if (graph->nodes[i].nadj != NULL_VTX) {
      xadj[k+1] = xadj[k] + graph->nodes[i].nadj;
      rename[i] = k;
      alias[k++] = i;
    }
  }
  nedges = xadj[graph->nvtxs];

  DL_ASSERT_EQUALS(nedges,2*graph->nedges,PF_ADJ_T);

  dprintf("Converting dynamic graph of "PF_VTX_T"/"PF_VTX_T" vertices and " \
      PF_ADJ_T" edges\n",graph->nvtxs,graph->max_nvtxs,nedges);

  if (nedges > 0) {
    /* allocate edges */
    adjncy = vtx_alloc(nedges);
    adjwgt = wgt_alloc(nedges);

    /* populate edges */
    for (i=0;i<graph->nvtxs;++i) {
      node = graph->nodes[alias[i]];
      DL_ASSERT(node.nadj != NULL_VTX,"Non present node being copied");
      for (j=0;j<node.nadj;++j) {
        adjncy[xadj[i]+j] = rename[node.adj[j]];
        adjwgt[xadj[i]+j] = node.wgt[j];
      }
    }
  } else {
    adjncy = NULL;
    adjwgt = NULL;
  }

  DL_ASSERT(check_graph(graph->nvtxs,xadj,adjncy,adjwgt),"Bad static graph " \
      "after converting from dynamic graph");

  dl_free(rename);

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
  if (r_adjwgt) {
    *r_adjwgt = adjwgt;
  } else {
    dl_free(adjwgt);
  }
  if (r_alias) {
    *r_alias = alias;
  } else {
    dl_free(alias);
  }
}


bowstring_dyngraph_t * bowstring_dg_todynamic(
    bowstring_vtx_t const nvtxs,
    bowstring_adj_t const * const xadj,
    bowstring_vtx_t const * const adjncy,
    bowstring_wgt_t const * const adjwgt)
{
  vtx_t i;
  wgt_t * ones;

  bowstring_dyngraph_t * graph;

  graph = bowstring_dg_create(nvtxs);

  if (adjwgt) {
    for (i=0;i<nvtxs;++i) {
      bowstring_dg_addvtx(i,adjncy+xadj[i],adjwgt+xadj[i],xadj[i+1]-xadj[i], \
          graph);
    }
  } else {
    ones = wgt_init_alloc(1,nvtxs);
    for (i=0;i<nvtxs;++i) {
      bowstring_dg_addvtx(i,adjncy+xadj[i],ones,xadj[i+1]-xadj[i],graph);
    }
    dl_free(ones);
  }
  graph->nedges = xadj[nvtxs]/2;

  DL_ASSERT(bowstring_dg_check(graph),"Bad dynamic graph converted from " \
      "static");

  return graph;
}


int bowstring_dg_check(
    bowstring_dyngraph_t const * const graph)
{
  vtx_t i, k, j, l, n;
  adj_t m;
  wgt_t wi, wk;

  n = 0;
  m = 0;
  for (i=0;i<graph->max_nvtxs;++i) {
    if (graph->nodes[i].nadj != NULL_VTX) {
      ++n;
      for (j=0;j<graph->nodes[i].nadj;++j) {
        k = graph->nodes[i].adj[j];
        if (k > graph->max_nvtxs) {
          eprintf("Invalid connection to vertex "PF_VTX_T"/"PF_VTX_T" from " \
              "vertex "PF_VTX_T"\n",k,graph->max_nvtxs,i);
          return 0;
        }
        if (graph->nodes[k].nadj == NULL_VTX) {
          eprintf("Invalid connection to missing vertex "PF_VTX_T"/"PF_VTX_T \
              " from vertex "PF_VTX_T"\n",k,graph->max_nvtxs,i);
          return 0;
        }
        wi = graph->nodes[i].wgt[j];
        for (l=0;l<graph->nodes[k].nadj;++l) {
          if (graph->nodes[k].adj[l] == i) {
            wk = graph->nodes[k].wgt[l];
            if (!dl_near_equal(wi,wk)) {
              eprintf("Unbalanced edge weight "PF_WGT_T"("PF_VTX_T"/"PF_VTX_T \
                  "):"PF_WGT_T"("PF_VTX_T"/"PF_VTX_T") for the edge from " \
                  "vertex "PF_VTX_T" to "PF_VTX_T"\n",wi,j, \
                  graph->nodes[i].nadj,wk,l,graph->nodes[k].nadj,i,k);
              return 0;
            }
            goto NEXT_EDGE;
          }
        }
        eprintf("Cound not find reverse edge going from "PF_VTX_T" to " \
            PF_VTX_T"\n",i,k);
        return 0;
        NEXT_EDGE:
        ++m;
      }
    }
  }
  m /= 2;

  if (n != graph->nvtxs) {
    eprintf("Incorrect number of vertices found "PF_VTX_T"/"PF_VTX_T"\n", \
        n,graph->nvtxs);
    return 0;
  }
  if (m != graph->nedges) {
    eprintf("Incorrect number of edges found "PF_ADJ_T"/"PF_ADJ_T"\n",m, \
        graph->nedges);
    return 0;
  }

  return 1;
}


vtx_t bowstring_dg_ncomp(
    bowstring_dyngraph_t const * const graph)
{
  vtx_t start, nq, sq, l, i, j, k;
  vtx_t * label, * q;

  q = vtx_alloc(graph->nvtxs);
  label = vtx_init_alloc(NULL_VTX,graph->max_nvtxs);

  start = 0;
  sq = 0;
  nq = 0;
  l = 0;
  while (start<graph->max_nvtxs) {
    if (graph->nodes[start].nadj != NULL_VTX && label[start] == NULL_VTX) {
      q[nq++] = start;
      label[start] = l;
      while (sq < nq) {
        i = q[sq++];
        for (j=0;j<graph->nodes[i].nadj;++j) {
          k = graph->nodes[i].adj[j];
          if (label[k] == NULL_VTX) {
            label[k] = l;
            q[nq++] = k;
          }
        }
      }
      ++l;
    }
    ++start;
  }
  dl_free(q);
  dl_free(label);

  return l;
}


void bowstring_dg_free(
    bowstring_dyngraph_t * graph)
{
  vtx_t i;

  for (i=0;i<graph->max_nvtxs;++i) {
    if (graph->nodes[i].nadj != NULL_VTX) {
      dl_free(graph->nodes[i].adj);
      dl_free(graph->nodes[i].wgt);
    }
  }

  dl_free(graph->nodes);
  dl_free(graph->htable);
  dl_free(graph);
}



#endif

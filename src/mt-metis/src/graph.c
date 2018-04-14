/**
* @file graph.c
* @brief routine for handling dgraphs
* @author Dominique LaSalle <lasalle@cs.umn.edu>
* @copyright 2014, Regents of the University of Minnesota
* @version 1
* @date 2012-06-12
*/




#ifndef MTMETIS_GRAPH_C
#define MTMETIS_GRAPH_C




#include <bowstring.h>
#include "graph.h"
#include "check.h"




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Free the part of the graph that thread 'myid' owns.
 *
 * @param graph The graph.
 * @param myid The thread id.
 */
static void __graph_free_part(
    graph_t * const graph,
    tid_t const myid)
{
  /* free graph structure */
  if (graph->free_xadj) {
    dl_free(graph->xadj[myid]);
  }
  if (graph->free_vwgt) { 
    dl_free(graph->vwgt[myid]);
  }
  if (graph->free_adjncy) {
    dl_free(graph->adjncy[myid]);
  }
  if (graph->free_adjwgt) {
    dl_free(graph->adjwgt[myid]);
  }

  /* free auxillery things */
  if (graph->cmap) {
    dl_free(graph->cmap[myid]);
  }
  if (graph->label) {
    dl_free(graph->label[myid]);
  }
}


/**
 * @brief Configure the distribution structure.
 *
 * @param maxnvtxs The maximum number of vertices owned by a thread.
 * @param nthreads The number of threads.
 * @param dist The distribution structure to configure.
 */
static void __graph_dist(
    vtx_t const maxnvtxs, 
    tid_t const nthreads,
    graphdist_t * const dist) 
{
  dist->nthreads = nthreads;
  dist->offset = vtx_uppow2(maxnvtxs+1);
  dist->mask = dist->offset - 1;
  dist->shift = vtx_downlog2(dist->offset);

  DL_ASSERT(maxnvtxs < (vtx_t)(1<<dist->shift),"Shift of %d for %"PF_VTX_T
      " vertices\n",dist->shift,maxnvtxs);
  DL_ASSERT_EQUALS(dist->offset,(vtx_t)(1<<dist->shift),"%"PF_VTX_T);
}


/**
 * @brief Distribute the vertices of a graph in continigous blocks. 
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param nthreads The number of threads to distribute the vertices over.
 * @param mynvtxs The number of vertices owned by each thread.
 * @param mynedges The number of edges owned by each thread.
 * @param lvtx The local vertex numbers.
 * @param owner The array assigning vertices to threads.
 */
static void __distribute_block(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    tid_t const nthreads,
    vtx_t * const mynvtxs,
    adj_t * const mynedges,
    vtx_t * const lvtx,
    tid_t * const owner)
{
  vtx_t v;
  tid_t myid;
  adj_t j, nedgesleft;

  adj_t const nedges = xadj[nvtxs];
  adj_t const avgnedges = nedges / nthreads;

  /* distribute vertices based on edge count */
  nedgesleft = nedges;
  myid = 0;
  for (v =0;v<nvtxs;++v) {
    if ((nthreads-myid-1)*avgnedges > nedgesleft) {
      ++myid;
    }
    owner[v] = myid;
    j = xadj[v+1] - xadj[v]; 
    mynedges[myid] += j;
    nedgesleft -= j;
  }
  for (v =0;v<nvtxs;++v) {
    myid = owner[v];
    lvtx[v] = mynvtxs[myid]++;
  }
}


/**
 * @brief Distribute the vertices of a graph cyclicly. 
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param nthreads The number of threads to distribute the vertices over.
 * @param mynvtxs The number of vertices owned by each thread.
 * @param mynedges The number of edges owned by each thread.
 * @param lvtx The local vertex numbers.
 * @param owner The array assigning vertices to threads.
 */
static void __distribute_cyclic(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    tid_t const nthreads,
    vtx_t * const mynvtxs,
    adj_t * const mynedges,
    vtx_t * const lvtx,
    tid_t * const owner)
{
  vtx_t i, v;
  tid_t myid;
  adj_t j, nedgesleft;
  vtx_t * perm;

  adj_t const nedges = xadj[nvtxs];
  adj_t const avgnedges = nedges / nthreads;

  perm = vtx_alloc(nvtxs);

  vtx_cyclicperm(perm,nthreads,nvtxs);

  /* distribute vertices based on edge count */
  nedgesleft = nedges;
  myid = 0;
  for (i =0;i<nvtxs;++i) {
    if ((nthreads-myid-1)*avgnedges > nedgesleft) {
      ++myid;
    }
    v = perm[i];
    owner[v] = myid;
    j = xadj[v+1] - xadj[v]; 
    mynedges[myid] += j;
    nedgesleft -= j;
  }
  for (i =0;i<nvtxs;++i) {
    myid = owner[i];
    lvtx[i] = mynvtxs[myid]++;
  }

  dl_free(perm);
}


/**
 * @brief Distribute the vertices of a graph block-cyclicly. 
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param nthreads The number of threads to distribute the vertices over.
 * @param mynvtxs The number of vertices owned by each thread.
 * @param mynedges The number of edges owned by each thread.
 * @param lvtx The local vertex numbers.
 * @param owner The array assigning vertices to threads.
 * @param block This size of a block.
 */
static void __distribute_blockcyclic(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    tid_t const nthreads,
    vtx_t * const mynvtxs,
    adj_t * const mynedges,
    vtx_t * const lvtx,
    tid_t * const owner,
    vtx_t block)
{
  vtx_t i, v;
  tid_t myid;
  adj_t j, nedgesleft;
  vtx_t * perm;

  adj_t const nedges = xadj[nvtxs];
  adj_t const avgnedges = nedges / nthreads;

  perm = vtx_alloc(nvtxs);

  /* create cyclic permutation */
  if (nthreads * block > nvtxs) {
    /* adjust the block if it will be imbalanced */
    block = dl_max(1,nvtxs/nthreads);
  }
  vtx_blockcyclicperm(perm,nthreads,block,nvtxs);

  /* distribute vertices based on edge count */
  nedgesleft = nedges;
  myid = 0;
  for (i =0;i<nvtxs;++i) {
    if ((nthreads-myid-1)*avgnedges > nedgesleft) {
      ++myid;
    }
    v = perm[i];
    owner[v] = myid;
    j = xadj[v+1] - xadj[v]; 
    mynedges[myid] += j;
    nedgesleft -= j;
  }
  for (i =0;i<nvtxs;++i) {
    myid = owner[i];
    lvtx[i] = mynvtxs[myid]++;
  }

  dl_free(perm);
}



/**
* @brief Initializes the graph
*
* @param graph The graph to initialize
* @param nthreads The number of threads to use
*/
static void __graph_init(
    graph_t * const graph,
    tid_t const nthreads) 
{
  DL_ASSERT(nthreads > 0,"__graph_init() called with 0 threads");
  #pragma omp master
  {
    graph->level = 0;

    /* graph size constants */
    graph->nvtxs = 0;
    graph->nedges = 0;
    graph->mincut = 0;
    graph->minvol = 0;
    graph->dist.nthreads = nthreads;

    /* memory for the graph structure */
    graph->mynvtxs = vtx_alloc(nthreads);
    graph->mynedges = adj_alloc(nthreads);
    graph->xadj = r_adj_alloc(nthreads);
    graph->vwgt = r_wgt_alloc(nthreads);
    graph->adjncy = r_vtx_alloc(nthreads);
    graph->adjwgt = r_wgt_alloc(nthreads);
    graph->label = NULL;
    graph->cmap = NULL;
    graph->where = NULL;
    graph->pwgts = NULL;
    graph->tvwgt = 0;
    graph->invtvwgt = 0;

    /* by default these are set to true, but the can be explicitly changed 
     * afterwards */
    graph->free_xadj = 1;
    graph->free_vwgt = 1;
    graph->free_adjncy = 1;
    graph->free_adjwgt = 1;

    /* linked-list structure */
    graph->coarser = NULL;
    graph->finer = NULL;
  }
  #pragma omp barrier
}




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


static graph_t * __gc_graph;
graph_t * graph_create(
    tid_t const nthreads)
{
  if (omp_in_parallel()) {
    #pragma omp master
    {
      __gc_graph = (graph_t*)calloc(1,sizeof(graph_t));
    }
    #pragma omp barrier
  } else {
    __gc_graph = (graph_t*)calloc(1,sizeof(graph_t));
  }
  __graph_init(__gc_graph,nthreads);

  return __gc_graph; 
}


static wgt_t * __gs_vwgt;
graph_t * graph_setup(
    vtx_t const nvtxs, 
    adj_t * const xadj, 
    vtx_t * const adjncy, 
    wgt_t * const adjwgt, 
    wgt_t * const vwgt) 
{
  graph_t * graph;

  tid_t const nthreads = omp_get_num_threads(); 
  tid_t const myid = omp_get_thread_num();

  graph = graph_create(nthreads);

  graph->mynvtxs[myid] = nvtxs;
  graph->mynedges[myid] = xadj[nvtxs];
  graph->xadj[myid] = xadj;
  graph->adjncy[myid] = adjncy;
  if (adjwgt) {
    graph->adjwgt[myid] = adjwgt;
  } else {
    graph->adjwgt[myid] = wgt_init_alloc(1,xadj[nvtxs]);
    graph->free_adjwgt = 1;
  }
  if (vwgt) {
    graph->vwgt[myid] = vwgt;
  } else {
    graph->vwgt[myid] = wgt_init_alloc(1,nvtxs);
    graph->free_vwgt = 1;
  }

  #pragma omp barrier
  #pragma omp master
  {
    graph->nvtxs = vtx_sum(graph->mynvtxs,nthreads);
    graph->nedges = adj_sum(graph->mynedges,nthreads);
    __graph_dist(vtx_max_value(graph->mynvtxs,nthreads),nthreads, \
        &(graph->dist));
    __gs_vwgt = wgt_alloc(nthreads);
  }
  #pragma omp barrier

  if (vwgt) {
    __gs_vwgt[myid] = wgt_sum(vwgt,nvtxs);
  }

  graph_setup_twgts(graph);

  DL_ASSERT(check_graph(graph),"Bad graph");

  return graph;
}


graph_t * graph_distribute(
    int distribution,
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    wgt_t const * const vwgt,
    tid_t const nthreads)
{
  vtx_t i,k,v,deg,mynvtxs;
  adj_t j,l;
  tid_t myid;
  vtx_t * dmynvtxs;
  adj_t * dmynedges;
  adj_t ** dxadj;
  vtx_t ** dadjncy, ** dlabel;
  wgt_t ** dadjwgt = NULL, ** dvwgt = NULL;
  graphdist_t dist;
  graph_t * graph;

  DL_ASSERT(nvtxs>0, "setup_graph() called with nvtxs = %"PF_VTX_T"\n",nvtxs);

  vtx_t * lvtx = vtx_alloc(nvtxs);
  tid_t * owner = tid_alloc(nvtxs);

  graph = graph_create(nthreads);

  /* set arrays from graph*/
  dmynvtxs = graph->mynvtxs;
  dmynedges = graph->mynedges;
  dxadj = graph->xadj;
  dadjncy = graph->adjncy;
  dadjwgt = graph->adjwgt;
  dvwgt = graph->vwgt;

  /* labels must be explicityl allocated */
  dlabel = graph->label = r_vtx_alloc(nthreads);

  /* zero out vertices and edges */
  vtx_set(dmynvtxs,0,nthreads);
  adj_set(dmynedges,0,nthreads);

  switch(distribution) {
    case MTMETIS_DISTRIBUTION_BLOCK:
      __distribute_block(nvtxs,xadj,nthreads,dmynvtxs,dmynedges, \
          lvtx,owner);
      break;
    case MTMETIS_DISTRIBUTION_CYCLIC:
      __distribute_cyclic(nvtxs,xadj,nthreads,dmynvtxs,dmynedges, \
          lvtx,owner);
      break;
    case MTMETIS_DISTRIBUTION_BLOCKCYCLIC:
      __distribute_blockcyclic(nvtxs,xadj,nthreads,dmynvtxs, \
          dmynedges,lvtx,owner,4096);
      break;
    default:
      dl_error("Unknown distribution '%d'\n",distribution);
  }

  __graph_dist(vtx_max_value(dmynvtxs,nthreads),nthreads,&dist);

  /* allocate arrays */
  for (myid =0;myid<nthreads;++myid) {
    mynvtxs = dmynvtxs[myid];
    dxadj[myid] = adj_alloc(mynvtxs+1);
    dxadj[myid][0] = 0;
    dlabel[myid] = vtx_alloc(mynvtxs);
    dadjncy[myid] = vtx_alloc(dmynedges[myid]);
    dvwgt[myid] = wgt_alloc(mynvtxs);
    dadjwgt[myid] = wgt_alloc(dmynedges[myid]);
    /* zero counts for insertion later */
    dmynvtxs[myid] = 0;
    dmynedges[myid] = 0;
  }

  /* set xadj and iadjwgt */
  for (v =0;v<nvtxs;++v) { 
    myid = owner[v];
    i = dmynvtxs[myid]++; 
    dlabel[myid][i] = v;
    DL_ASSERT_EQUALS(i,lvtx[v],"%"PF_VTX_T);
    deg = xadj[v+1] - xadj[v];
    dxadj[myid][i+1] = dxadj[myid][i] + deg;
    l = dmynedges[myid];
    if (adjwgt) {
      for (j =xadj[v];j<xadj[v+1];++j) {
        k = adjncy[j];
        if (owner[k] == myid) {
          dadjncy[myid][l] = lvtx[k];
        } else {
          dadjncy[myid][l] = lvtx_to_gvtx(lvtx[k],owner[k],dist);
        }
        dadjwgt[myid][l++] = adjwgt[j]; 
      }
    } else {
      for (j =xadj[v];j<xadj[v+1];++j) {
        k = adjncy[j];
        if (owner[k] == myid) {
          dadjncy[myid][l] = lvtx[k];
        } else {
          dadjncy[myid][l] = lvtx_to_gvtx(lvtx[k],owner[k],dist);
        }
        dadjwgt[myid][l++] = 1;
      }
    }
    DL_ASSERT_EQUALS(dxadj[myid][i+1],l,"%"PF_ADJ_T);
    dmynedges[myid] = l;
    if (vwgt) {
      dvwgt[myid][i] = vwgt[v];
    } else {
      dvwgt[myid][i] = 1;
    }
  }

  dl_free(owner);
  dl_free(lvtx);

  /* setup the graph */
  graph->gnvtxs = lvtx_to_gvtx(vtx_max_value(dmynvtxs,nthreads),nthreads-1, \
      dist);
  graph->nvtxs = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->dist = dist;

  /* setup tvwgt */
  graph->tvwgt = 0;
  for (myid=0;myid<nthreads;++myid) {
    graph->tvwgt += wgt_sum(graph->vwgt[myid],graph->mynvtxs[myid]);
  }
  graph->invtvwgt = 1.0/(graph->tvwgt > 0 ? graph->tvwgt : 1);

  /* setup tadjwgt */
  graph->tadjwgt = 0;
  for (myid=0;myid<nthreads;++myid) {
    graph->tadjwgt += wgt_sum(graph->adjwgt[myid],dmynedges[myid]);
  }

  /* set free configuration */
  graph->free_xadj = 1;
  graph->free_adjncy = 1;
  graph->free_adjwgt = 1;
  graph->free_vwgt = 1;

  DL_ASSERT(check_graph(graph),"Bad graph");

  return graph;
}


static adj_t * __gg_xadj;
static vtx_t * __gg_adjncy;
static wgt_t * __gg_vwgt;
static wgt_t * __gg_adjwgt;
static vtx_t * __gg_prefix;
static adj_t * __gg_nedges;
void graph_gather(
  graph_t const * const graph,
  adj_t ** const r_xadj,
  vtx_t ** const r_adjncy,
  wgt_t ** const r_vwgt,
  wgt_t ** const r_adjwgt,
  vtx_t * const r_voff)
{
  vtx_t i, k, voff, v;
  adj_t j, eoff;
  tid_t t;
  adj_t * gxadj;
  vtx_t * gadjncy;
  wgt_t * gvwgt, * gadjwgt;

  tid_t const myid = omp_get_thread_num();
  tid_t const nthreads = omp_get_num_threads();

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const mynedges = graph->mynedges[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const vwgt = graph->vwgt[myid];
  wgt_t const * const adjwgt = graph->adjwgt[myid];

  DL_ASSERT_EQUALS(mynedges,xadj[mynvtxs],"%"PF_ADJ_T);
  DL_ASSERT_EQUALS(adj_sum(graph->mynedges,nthreads),graph->nedges, \
      "%"PF_ADJ_T);

  #pragma omp master
  {
    __gg_xadj = adj_alloc(graph->nvtxs+1);
    __gg_adjncy = vtx_alloc(graph->nedges);
    __gg_vwgt = wgt_alloc(graph->nvtxs);
    __gg_adjwgt = wgt_alloc(graph->nedges);
    __gg_prefix = vtx_alloc(nthreads+1);
    __gg_nedges = adj_alloc(nthreads+1);
    /* create a prefix sum for placing vertices */
    __gg_prefix[0] = 0;
    for (t=0;t<nthreads;++t) {
      __gg_prefix[t+1] = __gg_prefix[t] + graph->mynvtxs[t];
    }
    /* create a prefix sum for placing edges */
    __gg_nedges[0] = 0;
    for (t=0;t<nthreads;++t) {
      __gg_nedges[t+1] = __gg_nedges[t] + graph->mynedges[t];
    }
    /* cap the xadj array */
    __gg_xadj[graph->nvtxs] = graph->nedges;
  }
  #pragma omp barrier
  voff = __gg_prefix[myid];
  eoff = __gg_nedges[myid];
  gxadj = __gg_xadj + voff;
  gadjncy = __gg_adjncy + eoff;
  gvwgt = __gg_vwgt + voff;
  gadjwgt = __gg_adjwgt + eoff;

  /* vertex ids are purely based on thread offsets */
  eoff = gxadj[0] = __gg_nedges[myid];
  for (i=1;i<mynvtxs;++i) {
    gxadj[i] = xadj[i] + eoff; 
  }

  /* insert edges into graph */
  for (i=0;i<mynvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        t = myid;
        v = k;
      } else {
        t = gvtx_to_tid(k,graph->dist);
        v = gvtx_to_lvtx(k,graph->dist);
      }
      gadjncy[j] = v + __gg_prefix[t];
    }
  }
  wgt_copy(gadjwgt,adjwgt,mynedges);

  /* copy vwgt */
  wgt_copy(gvwgt,vwgt,mynvtxs);

  #pragma omp barrier
  #pragma omp master
  {
    dl_free(__gg_nedges);
    dl_free(__gg_prefix);
  }

  /* assign pointers */
  *r_xadj = __gg_xadj;
  *r_adjncy = __gg_adjncy;
  *r_vwgt = __gg_vwgt;
  *r_adjwgt = __gg_adjwgt;
  *r_voff = voff;

  DL_ASSERT(bowstring_check_graph(graph->nvtxs,__gg_xadj,__gg_adjncy, \
        (bowstring_wgt_t*)__gg_adjwgt),"Bad graph gathered");
}


graph_t * graph_setup_coarse(
    graph_t * const graph, 
    vtx_t const cnvtxs)
{
  vtx_t mynvtxs;
  graph_t * cgraph;

  tid_t const myid = omp_get_thread_num();
  tid_t const nthreads = omp_get_num_threads();

  cgraph = graph_create(nthreads);

  #pragma omp master
  {
    graph->coarser = cgraph;
    cgraph->finer = graph;

    cgraph->level = graph->level + 1;

    DL_ASSERT_EQUALS(nthreads,graph->dist.nthreads,"%"PF_TID_T);

    cgraph->dist = graph->dist;

    cgraph->tvwgt = graph->tvwgt;
    cgraph->invtvwgt = graph->invtvwgt;
  }
  #pragma omp barrier

  cgraph = graph->coarser;

  DL_ASSERT(cgraph != NULL,"cgraph is NULL");

  cgraph->mynvtxs[myid] = cnvtxs;

  /* Allocate memory for the coarser graph */
  mynvtxs = cnvtxs;

  cgraph->xadj[myid] = adj_alloc(mynvtxs+1);
  if (mynvtxs > 0) {
    cgraph->vwgt[myid] = wgt_alloc(mynvtxs);
  } else {
    cgraph->xadj[myid][0] = 0;
    cgraph->vwgt[myid] = NULL;
  }

  cgraph->adjncy[myid] = NULL;
  cgraph->adjwgt[myid] = NULL;

  mynvtxs = vtx_omp_sumreduce(mynvtxs);

  #pragma omp master
  {
    cgraph->gnvtxs = graph->gnvtxs;
    cgraph->nvtxs = mynvtxs;
    DL_ASSERT(cgraph->gnvtxs >= cgraph->nvtxs,"Bad gnvtxs");
  }
  #pragma omp barrier

  return cgraph;
}


void graph_setup_twgts(
    graph_t * const graph)
{
  vtx_t vsum,asum;

  tid_t const myid = omp_get_thread_num();

  vsum = wgt_sum(graph->vwgt[myid],graph->mynvtxs[myid]);
  asum = wgt_sum(graph->adjwgt[myid],graph->mynedges[myid]);
  vsum = wgt_omp_sumreduce(vsum);
  asum = wgt_omp_sumreduce(asum);

  #pragma omp master 
  {
    graph->tvwgt = vsum;
    graph->tadjwgt = asum;
    graph->invtvwgt = 1.0/(graph->tvwgt > 0 ? graph->tvwgt : 1);
  }
  #pragma omp barrier
}


void graph_alloc_partmemory(
    ctrl_t * const ctrl,
    graph_t * const graph)
{
  tid_t const nthreads = graph->dist.nthreads;
  tid_t const myid = omp_get_thread_num();

  DL_ASSERT_EQUALS((tid_t)omp_get_num_threads(),graph->dist.nthreads, \
      "%"PF_TID_T);

  #pragma omp master
  {
    /* memory for the partition/refinement structure */
    graph->where = r_pid_alloc(nthreads);
    graph->pwgts = wgt_alloc(ctrl->nparts);
  }
  #pragma omp barrier
  graph->where[myid] = pid_alloc(graph->mynvtxs[myid]);
}


void graph_free(
    graph_t * graph)
{
  tid_t myid;

  if (omp_in_parallel()) {
    myid = omp_get_thread_num();
    __graph_free_part(graph,myid);
  } else {
    for (myid=0;myid<graph->dist.nthreads;++myid) {
      __graph_free_part(graph,myid);
    }
  }

  /* free partition/refinement structure */
  graph_free_rdata(graph);

  #pragma omp barrier
  #pragma omp master
  {
    dl_free(graph->xadj);
    dl_free(graph->adjncy);
    dl_free(graph->vwgt);
    dl_free(graph->adjwgt);
    dl_free(graph->mynvtxs);
    dl_free(graph->mynedges);

    if (graph->cmap) {
      dl_free(graph->cmap);
    }
    if (graph->label) {
      dl_free(graph->label);
    }

    dl_free(graph);
  }
}

void graph_free_rdata(
    graph_t * graph)
{
  tid_t myid;
  
  if (omp_in_parallel()) {
    myid = omp_get_thread_num();
    if (graph->where) {
      dl_free(graph->where[myid]);
    }
    if (graph->rename) {
      dl_free(graph->rename[myid]);
    }
  } else {
    for (myid=0;myid<graph->dist.nthreads;++myid) {
      if (graph->where) {
        dl_free(graph->where[myid]);
      }
      if (graph->rename) {
        dl_free(graph->rename[myid]);
      }
    }
  }

  #pragma omp barrier
  #pragma omp master
  {
    /* free partition/refinement structure */
    if (graph->pwgts) {
      dl_free(graph->pwgts);
      graph->pwgts = NULL;
    }
    if (graph->where) {
      dl_free(graph->where);
      graph->where = NULL;
    }
    if (graph->rename) {
      dl_free(graph->rename);
      graph->rename = NULL;
    }
  }
}


double graph_imbalance(
    graph_t const * const graph,
    pid_t const nparts,
    real_t const * const pijbm)
{
  vtx_t k;
  double max, cur;

  DL_ASSERT_EQUALS(wgt_sum(graph->pwgts,nparts),graph->tvwgt,"%"PF_WGT_T);

  max = 0;

  if (nparts > 256) {
    #pragma omp for schedule(static,128)
    for (k =0;k<nparts;++k) {
      cur = graph->pwgts[k]*pijbm[k];
      if (cur > max) {
        max = cur;
      }
    }
    max = double_omp_maxreduce_value(max);
  } else {
    for (k =0;k<nparts;++k) {
      cur = graph->pwgts[k]*pijbm[k];
      if (cur > max) {
        max = cur;
      }
    }
  }

  return max;
}


double graph_imbalance_diff(
    graph_t const * const graph,
    pid_t const nparts,
    real_t const * const pijbm,
    real_t const ubfactor)
{
  vtx_t k;
  double max, cur;

  DL_ASSERT_EQUALS(wgt_sum(graph->pwgts,nparts),graph->tvwgt,"%"PF_WGT_T);

  max = 0;

  if (nparts > 256) {
    #pragma omp for schedule(static,128)
    for (k =0;k<nparts;++k) {
      cur = graph->pwgts[k]*pijbm[k]-ubfactor;
      if (cur > max) {
        max = cur;
      }
    }
    max = double_omp_maxreduce_value(max);
  } else {
    for (k =0;k<nparts;++k) {
      cur = graph->pwgts[k]*pijbm[k]-ubfactor;
      if (cur > max) {
        max = cur;
      }
    }
  }

  return max;
}


wgt_t graph_cut(
    ctrl_t const * const ctrl,
    graph_t const * const graph,
    pid_t const * const * const where)
{
  vtx_t i, k, lvtx, nbrid;
  adj_t j;
  wgt_t cut;

  tid_t const myid = omp_get_thread_num();

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const adjwgt = graph->adjwgt[myid];
  pid_t const * const mywhere = where[myid];

  DL_ASSERT_EQUALS((int)graph->dist.nthreads,omp_get_num_threads(),"%d");

  cut = 0;

  if (graph->adjwgt == NULL) {
    for (i =0; i<mynvtxs; ++i) {
      for (j =xadj[i]; j<xadj[i+1]; ++j) {
        k = adjncy[j]; 
        if (k < mynvtxs) {
          lvtx = k;
          nbrid = myid;
        } else {
          nbrid = gvtx_to_tid(adjncy[j],graph->dist);
          lvtx = gvtx_to_lvtx(adjncy[j],graph->dist);
        }
        if (mywhere[i] != where[nbrid][lvtx]) {
          ++cut;
        }
      }
    }
  } else {
    for (i =0; i<mynvtxs; ++i) {
      for (j =xadj[i]; j<xadj[i+1]; ++j) {
        k = adjncy[j];
        if (k < mynvtxs) {
          lvtx = k;
          nbrid = myid;
        } else {
          nbrid = gvtx_to_tid(adjncy[j],graph->dist);
          lvtx = gvtx_to_lvtx(adjncy[j],graph->dist);
        }
        if (mywhere[i] != where[nbrid][lvtx]) {
          cut += adjwgt[j];
        }
      }
    }
  }
  cut = wgt_omp_sumreduce(cut);

  return cut/2;
}

int graph_isbalanced(
    ctrl_t const * const ctrl, 
    graph_t const * const graph, 
    real_t const ffactor)
{
  return (graph_imbalance_diff(graph,ctrl->nparts,ctrl->pijbm, \
      ctrl->ubfactor) <= ffactor);
}


void graph_readjust_memory(
    graph_t * const graph,
    adj_t adjsize)
{
  int const myid = omp_get_thread_num();
  adj_t const nedges = graph->xadj[myid][graph->mynvtxs[myid]];

  if (adjsize > 4096 && adjsize * 0.75 > nedges) {
    graph->adjncy[myid] = adj_realloc(graph->adjncy[myid],nedges);
    graph->adjwgt[myid] = wgt_realloc(graph->adjwgt[myid],nedges);
  }
}




#endif

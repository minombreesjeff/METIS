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
* PRIVATE TYPES ***************************************************************
******************************************************************************/


typedef struct gau_ptrs_t {
  adj_t * xadj;
  vtx_t * adjncy;
  wgt_t * vwgt;
  wgt_t * adjwgt;
  vtx_t * prefix;
  adj_t * nedges;
} gau_ptrs_t;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLPQ_PREFIX vv
#define DLPQ_KEY_T vtx_t
#define DLPQ_VAL_T vtx_t
#define DLPQ_MIN 1
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_MIN
#undef DLPQ_KEY_T
#undef DLPQ_VAL_T
#undef DLPQ_PREFIX




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


/* if we can get a 1 percent increase in our balance constraint, it is probably
 * worth it */
static double const MIN_ISLAND_WEIGHT = 0.01;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Allocate the structures for a unified graph (from a distributed one).
 *
 * @param graph The distributed graph.
 * @param r_uxadj A reference to the unified adjacency list pointer.
 * @param r_uadjncy A reference to the unified adjacecny indexes.
 * @param r_uvwgt A reference to the unified vertex weights.
 * @param r_uadjwgt A reference to the unified edges weights.
 * @param r_uprefix A reference to the prefix sum of vertices (per thread).
 * @param r_unedges A reference to the prefix sum of edges (per thread).
 */
static void __par_graph_alloc_unified(
    graph_t const * const graph,
    adj_t ** const r_uxadj,
    vtx_t ** const r_uadjncy,
    wgt_t ** const r_uvwgt,
    wgt_t ** const r_uadjwgt,
    vtx_t ** const r_uprefix,
    adj_t ** const r_unedges)
{
  tid_t t;
  gau_ptrs_t * ptr;

  tid_t const myid = dlthread_get_id(graph->comm);
  tid_t const nthreads = dlthread_get_nthreads(graph->comm);

  ptr = dlthread_get_shmem(sizeof(gau_ptrs_t),graph->comm);
  if (myid == 0) {
    ptr->xadj = adj_alloc(graph->nvtxs+1); 
    ptr->adjncy = vtx_alloc(graph->nedges); 
    ptr->vwgt = wgt_alloc(graph->nvtxs); 
    ptr->adjwgt = wgt_alloc(graph->nedges);
    ptr->prefix = vtx_alloc(nthreads+1);
    ptr->nedges = adj_alloc(nthreads+1);

    /* create a prefix sum for placing vertices */
    ptr->prefix[0] = 0;
    for (t=0;t<nthreads;++t) {
      ptr->prefix[t+1] = ptr->prefix[t] + graph->mynvtxs[t];
    }
    /* create a prefix sptr->m for placing edges */
    ptr->nedges[0] = 0;
    for (t=0;t<nthreads;++t) {
      ptr->nedges[t+1] = ptr->nedges[t] + graph->mynedges[t];
    }
    /* cap the xadj array */
    ptr->xadj[graph->nvtxs] = graph->nedges;
  }
  dlthread_barrier(graph->comm);

  *r_uxadj = ptr->xadj; 
  *r_uadjncy = ptr->adjncy; 
  *r_uvwgt = ptr->vwgt; 
  *r_uadjwgt = ptr->adjwgt; 
  *r_uprefix = ptr->prefix; 
  *r_unedges = ptr->nedges; 

  dlthread_barrier(graph->comm);
  dlthread_free_shmem(ptr,graph->comm);
}


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
  if (graph->label) {
    dl_free(graph->label[myid]);
  }
  if (graph->group) {
    dl_free(graph->group[myid]);
  }
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
  graph->level = 0;

  /* graph size constants */
  graph->nvtxs = 0;
  graph->nedges = 0;
  graph->mincut = 0;
  graph->minvol = 0;
  graph->dist.nthreads = nthreads;

  /* pre partitioning info */
  graph->ngroup = 0;
  graph->group = NULL;

  /* memory for the graph structure */
  if (nthreads > 0) {
    graph->mynvtxs = vtx_alloc(nthreads);
    graph->mynedges = adj_alloc(nthreads);
    graph->xadj = r_adj_alloc(nthreads);
    graph->vwgt = r_wgt_alloc(nthreads);
    graph->adjncy = r_vtx_alloc(nthreads);
    graph->adjwgt = r_wgt_alloc(nthreads);
  }
  graph->label = NULL;
  graph->rename = NULL;
  graph->cmap = NULL;
  graph->nislands = NULL;
  graph->uniformvwgt = 0;
  graph->uniformadjwgt = 0;
  graph->tvwgt = 0;
  graph->invtvwgt = 0;

  /* nullify partition data */
  graph->where = NULL;
  graph->pwgts = NULL;
  graph->vsinfo = NULL;
  graph->esinfo = NULL;
  graph->kwinfo = NULL;

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


static void __cuthillmckee(
    vtx_t const mynvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    pid_t * const perm)
{
  vtx_t i,k,d,nordered,sr;
  adj_t j;
  vtx_t * deg;
  vv_pq_t * q, * rem;

  q = vv_pq_create(0,mynvtxs);
  rem = vv_pq_create(0,mynvtxs);

  /* offset pointer */
  deg = vtx_alloc(mynvtxs);

  /* find my lowest degree vertex */
  for (i=0;i<mynvtxs;++i) {
    d = 0;
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j]; 
      if (k >= mynvtxs) {
        ++d;
      }
    }
    deg[i] = d;
    vv_pq_push(d,i,rem);
  }

  sr = nordered = mynvtxs;

  /* loop through connected components */
  while (rem->size > 0) {
    i = vv_pq_pop(rem);
    perm[nordered++] = i;

    /* perform bfs */
    while (sr < nordered) {
      i = perm[sr++];
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j]; 
        if (k >= mynvtxs && vv_pq_contains(k,rem)) {
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

  vv_pq_free(q);
  vv_pq_free(rem);
  /* un-offset */
  dl_free(deg);
}





/******************************************************************************
* PRIVATE PARALLEL FUNCTIONS **************************************************
******************************************************************************/


static void __par_graph_init(
    graph_t * const graph,
    dlthread_comm_t comm)
{
  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  if (myid == 0) {
    __graph_init(graph,nthreads);
    graph->comm = comm;
  }
  dlthread_barrier(comm);
}


static graph_t * __par_graph_bndgraph(
    ctrl_t const * const ctrl,
    graph_t const * const graph,
    int const * const * const include)
{
  vtx_t i, v, k, ninc, lvtx, maxnvtxs;
  adj_t j, l, con[2], incon[2];
  pid_t p;
  tid_t nbrid;
  adj_t nedges;
  wgt_t hwgt[2];
  adj_t * myxadj;
  vtx_t * myadjncy, * rename, ** label;
  wgt_t * myvwgt, * myadjwgt;
  graphdist_t dist;
  graph_t * bndgraph;

  tid_t const myid = dlthread_get_id(ctrl->comm);
  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const vwgt = graph->vwgt[myid];
  wgt_t const * const adjwgt = graph->adjwgt[myid];
  pid_t const * const * const gwhere = (pid_t const **)graph->where;

  /* construct boundary graph --
   *   the general strategy for constructing this graph is that the two super
   *   nodes will go first [0] and [1], followed [0]'s neighbors, and then
   *   [1]'s neighbors, then followed copying the interior points from the
   *   original graph. */

  rename = vtx_alloc(mynvtxs+2);

  label = dlthread_get_shmem(sizeof(vtx_t*)*nthreads,ctrl->comm);
  label[myid] = vtx_alloc(mynvtxs);

  /* count number of vertices each thread has in the boundary graph */
  ninc = 2;
  hwgt[0] = hwgt[1] = 0;
  for (i=0;i<mynvtxs;++i) {
    if (include[myid][i]) {
      label[myid][i] = ninc;
      rename[ninc++] = i;
    } else {
      hwgt[gwhere[myid][i]] += vwgt[i];
    }
  }

  myxadj = adj_init_alloc(0,ninc+1);
  myvwgt = wgt_alloc(ninc);
  myvwgt[0] = hwgt[0];
  myvwgt[1] = hwgt[1];

  nedges = 0;

  /* we'll add these as we insert the other end of edges */
  for (i=2;i<ninc;++i) {
    v = rename[i];
    myvwgt[i] = vwgt[v];
    con[0] = con[1] = NULL_ADJ;
    for (j=xadj[v];j<xadj[v+1];++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        lvtx = k;
        nbrid = myid;
      } else {
        lvtx = gvtx_to_lvtx(k,graph->dist);
        nbrid = gvtx_to_tid(k,graph->dist);
      }
      if (include[nbrid][lvtx]) { 
        ++myxadj[i+1];
      } else {
        p = gwhere[nbrid][lvtx];
        if (con[p] == NULL_ADJ) {
          con[p] = 1;
          ++myxadj[i+1];
          /* increase the size of the super-vertex */
          ++myxadj[p+1];
        }
      }
    }
  }

  /* implicit barrier neccessary to ensure label is populated */
  maxnvtxs = vtx_dlthread_maxreduce_value(ninc,ctrl->comm);

  graph_calc_dist(maxnvtxs,nthreads,&dist);

  nedges = adj_prefixsum_exc(myxadj+1,ninc);

  myadjncy = vtx_alloc(nedges);
  myadjwgt = wgt_alloc(nedges);

  for (i=2;i<ninc;++i) {
    v = rename[i];
    con[0] = con[1] = NULL_ADJ;
    l = myxadj[i+1];
    for (j=xadj[v];j<xadj[v+1];++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        lvtx = k;
        nbrid = myid;
      } else {
        lvtx = gvtx_to_lvtx(k,graph->dist);
        nbrid = gvtx_to_tid(k,graph->dist);
      }
      if (include[nbrid][lvtx]) { 
        if (nbrid == myid) {
          myadjncy[l] = label[nbrid][lvtx];
        } else {
          myadjncy[l] = lvtx_to_gvtx(label[nbrid][lvtx],nbrid,dist);
        }
        myadjwgt[l] = adjwgt[j];
        ++l;
      } else {
        p = gwhere[nbrid][lvtx];
        if (con[p] == NULL_ADJ) {
          myadjncy[l] = p; 
          myadjwgt[l] = adjwgt[j];
          con[p] = l++;
          /* reverse edge -- supernode */
          incon[p] = myxadj[p+1]++;
          myadjncy[incon[p]] = i;
          myadjwgt[incon[p]] = adjwgt[j];
        } else {
          myadjwgt[con[p]] += adjwgt[j];
          myadjwgt[incon[p]] += adjwgt[j];
        }
        DL_ASSERT_EQUALS(myadjwgt[con[p]],myadjwgt[incon[p]],"%"PF_WGT_T);
      }
    }
    myxadj[i+1] = l;
  }
  DL_ASSERT_EQUALS(myxadj[ninc],nedges,"%"PF_ADJ_T);
  dlthread_barrier(ctrl->comm);

  dl_free(label[myid]);

  dlthread_free_shmem(label,ctrl->comm);

  bndgraph = par_graph_setup(ninc,myxadj,myadjncy,myvwgt,myadjwgt,ctrl->comm);

  if (myid == 0) {
    bndgraph->label = r_vtx_alloc(nthreads);
  }
  dlthread_barrier(ctrl->comm);

  bndgraph->label[myid] = rename;

  dlthread_barrier(ctrl->comm);

  return bndgraph;
}


static tid_t __par_graph_extract_halves_nogroup(
    graph_t * const graph,
    pid_t const * const * const where,
    graph_t ** const halves)
{
  vtx_t i, k, hnvtxs, deg, g, v, u;
  adj_t j, l;
  pid_t w;
  tid_t o;
  tid_t hmyid, hnthreads, mygroup;
  dlthread_comm_t hcomm;
  vtx_t ** vprefix;
  vtx_t ** vsuffix;
  vtx_t nvtxs[2], hmynvtxs[2];
  adj_t nedges[2];

  tid_t const myid = dlthread_get_id(graph->comm);
  tid_t const nthreads = dlthread_get_nthreads(graph->comm);

  tid_t const fnthreads = nthreads / 2;

  /* create my vertex and edge counts */
  nvtxs[0] = nvtxs[1] = 0;
  nedges[0] = nedges[1] = 0;
  for (i=0;i<graph->mynvtxs[myid];++i) {
    w = where[myid][i];
    if (w < 2) {
      ++nvtxs[w];
    }
  }

  /* determine group assignements */
  if (myid < fnthreads) {
    mygroup = 0;
    hnthreads = fnthreads;
  } else {
    mygroup = 1;
    hnthreads = nthreads - fnthreads;
  }

  vprefix = dlthread_get_shmem(sizeof(vtx_t*)*4,graph->comm);
  vsuffix = vprefix+2;

  /* split the thread groups */
  hcomm = dlthread_comm_split(mygroup,2,graph->comm);

  DL_ASSERT_EQUALS((size_t)hnthreads,dlthread_get_nthreads(hcomm),"%zu");

  hmyid = dlthread_get_id(hcomm);

  halves[mygroup] = par_graph_create(hcomm);

  if (hmyid == 0) {
    /* lead threads only */
    vprefix[mygroup] = vtx_alloc(nthreads+1);
    halves[mygroup]->label = r_vtx_alloc(hnthreads);
  }

  if (myid == 0) {
    graph->rename = r_vtx_alloc(nthreads);
  }

  dlthread_barrier(graph->comm);

  /* vertices */
  vprefix[0][myid] = nvtxs[0]; 
  vprefix[1][myid] = nvtxs[1]; 

  /* rename vectors */
  graph->rename[myid] = vtx_alloc(graph->mynvtxs[myid]);

  dlthread_barrier(graph->comm);
  if (hmyid == 0) {
    /* create prefixsums of the vertices and edges */
    vprefix[mygroup][nthreads] = 0;
    vtx_prefixsum_exc(vprefix[mygroup],nthreads+1);

    hnvtxs = vprefix[mygroup][nthreads];

    /* create prefixsums for actual insertion into split graphs */
    vsuffix[mygroup] = vtx_alloc(hnthreads);
    for (o=0;o<hnthreads;++o) {
      vsuffix[mygroup][o] = vtx_chunkstart(o,hnthreads,hnvtxs);
      halves[mygroup]->mynvtxs[o] = vtx_chunksize(o,hnthreads,hnvtxs);
    }

    graph_calc_dist(vtx_max_value(halves[mygroup]->mynvtxs,hnthreads),hnthreads, \
        &(halves[mygroup]->dist));
  }
  dlthread_barrier(graph->comm);

  /* allocate vwgt and xadj */
  hnvtxs = halves[mygroup]->mynvtxs[hmyid];
  halves[mygroup]->xadj[hmyid] = adj_alloc(hnvtxs+1);
  halves[mygroup]->vwgt[hmyid] = wgt_alloc(hnvtxs);
  halves[mygroup]->label[hmyid] = vtx_alloc(hnvtxs);

  dlthread_barrier(graph->comm);
  /* insert vertex information into graphs */
  hmynvtxs[0] = hmynvtxs[1] = 0;
  for (i=0;i<graph->mynvtxs[myid];++i) {
    w = where[myid][i];
    if (w < 2) {
      u = hmynvtxs[w]++;
      hnvtxs = vprefix[w][nthreads];
      hnthreads = (w == 0 ? fnthreads : nthreads - fnthreads);
      
      /* get total vertex number */
      g = u + vprefix[w][myid];

      DL_ASSERT(g < hnvtxs,"Got vertex number of %"PF_VTX_T"/%"PF_VTX_T" " \
          "from %"PF_VTX_T" and %"PF_VTX_T"\n",g,hnvtxs,u, \
          vprefix[w][myid]);

      /* get new local vertex number */
      hmyid = vtx_chunkid(g,hnthreads,hnvtxs); 

      DL_ASSERT(hmyid < hnthreads,"Got chunk id of %"PF_TID_T"/%"PF_TID_T" " \
          "from %"PF_VTX_T", %"PF_TID_T", %"PF_VTX_T"\n",hmyid,hnthreads, \
          g,hnthreads,hnvtxs);

      v = g - vsuffix[w][hmyid];

      /* set rename */
      graph->rename[myid][i] = lvtx_to_gvtx(v,hmyid,halves[w]->dist);

      /* set alias as global vertex ID */
      halves[w]->label[hmyid][v] = lvtx_to_gvtx(i,myid,graph->dist);

      /* copy vertex weight */
      halves[w]->vwgt[hmyid][v] = graph->vwgt[myid][i];

      /* copy xadj info */
      deg = 0;
      for (j=graph->xadj[myid][i];j<graph->xadj[myid][i+1];++j) {
        k = graph->adjncy[myid][j];
        if (k < graph->mynvtxs[myid]) {
          o = myid;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        if (where[o][k] == w) {
          ++deg;
        }
      }
      halves[w]->xadj[hmyid][v] = deg;
    }
  }

  dlthread_barrier(graph->comm);

  /* fix respective xadj's */
  if (myid < fnthreads) {
    mygroup = 0;
    hmyid = myid; 
    hnthreads = fnthreads;
  } else {
    mygroup = 1;
    hmyid = myid - fnthreads;
    hnthreads = nthreads - fnthreads;
  }

  hnvtxs = halves[mygroup]->mynvtxs[hmyid];
  halves[mygroup]->xadj[hmyid][hnvtxs] = 0;
  adj_prefixsum_exc(halves[mygroup]->xadj[hmyid],hnvtxs+1);

  halves[mygroup]->adjncy[hmyid] = \
      vtx_alloc(halves[mygroup]->xadj[hmyid][hnvtxs]);
  halves[mygroup]->adjwgt[hmyid] = \
      wgt_alloc(halves[mygroup]->xadj[hmyid][hnvtxs]);

  halves[mygroup]->mynedges[hmyid] = halves[mygroup]->xadj[hmyid][hnvtxs];

  dlthread_barrier(graph->comm);
  /* insert edge information into graphs */
  hmynvtxs[0] = hmynvtxs[1] = 0;
  for (i=0;i<graph->mynvtxs[myid];++i) {
    w = where[myid][i];
    if (w < 2) {
      u = hmynvtxs[w]++;
      hnvtxs = vprefix[w][nthreads];
      hnthreads = (w == 0 ? fnthreads : nthreads - fnthreads);
      
      /* get total vertex number */
      g = u + vprefix[w][myid];

      /* get new local vertex number */
      hmyid = vtx_chunkid(g,hnthreads,hnvtxs); 
      v = g - vsuffix[w][hmyid];

      l = halves[w]->xadj[hmyid][v];
      for (j=graph->xadj[myid][i];j<graph->xadj[myid][i+1];++j) {
        k = graph->adjncy[myid][j];
        if (k < graph->mynvtxs[myid]) {
          o = myid;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        if (where[o][k] == w) {
          /* calcuate new endpoint */
          u = graph->rename[o][k];
          o = gvtx_to_tid(u,halves[w]->dist);
          if (o == hmyid) {
            k = gvtx_to_lvtx(u,halves[w]->dist);
          } else {
            k = u;
          }

          /* process edge */
          halves[w]->adjncy[hmyid][l] = k;
          halves[w]->adjwgt[hmyid][l] = graph->adjwgt[myid][j];
          ++l;
        }
      }
    }
  }
  dlthread_barrier(graph->comm);

  if (myid < fnthreads) {
    hmyid = myid; 
    mygroup = 0;
    hnthreads = fnthreads;
  } else {
    hmyid = myid - fnthreads;
    mygroup = 1;
    hnthreads = nthreads - fnthreads;
  }

  if (hmyid == 0) {
    /* lead threads only */
    halves[mygroup]->nvtxs = vtx_sum(halves[mygroup]->mynvtxs,hnthreads);
    halves[mygroup]->gnvtxs = max_gvtx(halves[mygroup]); 
    halves[mygroup]->nedges = adj_sum(halves[mygroup]->mynedges,hnthreads);

    dl_free(vprefix[mygroup]);
    dl_free(vsuffix[mygroup]);
  }

  par_graph_setup_twgts(halves[mygroup]);

  /* implicit barrier */
  dlthread_free_shmem(vprefix,graph->comm); /* vsuffix is part of this */

  return mygroup;
}


static tid_t __par_graph_extract_halves_group(
    graph_t * const graph,
    pid_t const * const * const gwhere,
    graph_t ** const halves)
{
  vtx_t i, k, hnvtxs, g, v, u, lvtx, maxgvtx;
  adj_t j, hnedges;
  pid_t d, other;
  tid_t hmyid, hnthreads, mygroup, o, t, nbrid;
  dlthread_comm_t hcomm;
  graph_t * mygraph;
  vtx_t * gvtxs;
  vtx_t * mydist[2];

  tid_t const myid = dlthread_get_id(graph->comm);
  tid_t const nthreads = dlthread_get_nthreads(graph->comm);

  tid_t const fnthreads = nthreads / 2;

  vtx_t const * const gmynvtxs = graph->mynvtxs;
  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;
  wgt_t const * const * const gvwgt = (wgt_t const **)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const **)graph->adjwgt;

  pid_t const ngroup = graph->ngroup;
  pid_t const * const * const group = (pid_t const **)graph->group;

  /* don't bother with weights on the base graph */
  int const do_wgt = graph->level > 0; /* this is a bug! */

  /* Create a global array of vertices by distribution id here ***************/

  /* allocate my rename vector first to use natural vectors */
  if (myid == 0) {
    graph->rename = r_vtx_alloc(nthreads);
  }

  gvtxs = dlthread_get_shmem(graph->nvtxs*sizeof(vtx_t),graph->comm);

  mydist[0] = vtx_init_alloc(0,ngroup*2+1);
  mydist[1] = mydist[0] + ngroup;

  /* create my vertex and edge counts */
  for (i=0;i<gmynvtxs[myid];++i) {
    other = gwhere[myid][i];
    if (other < 2) {
      /* count my vertices per dist */
      ++mydist[other][group[myid][i]];
    }
  }

  /* each thread owns disjoint distributions, so there really is a faster way
   * to do this */
  vtx_dlthread_sumareduce(mydist[0],2*ngroup,graph->comm);

  /* create a prefix sum of each mydist array to index into a global array of
   * vertices */
  vtx_prefixsum_exc(mydist[0],(2*ngroup)+1);

  /* insert my vertices global id numbers into the global array */
  for (i=0;i<gmynvtxs[myid];++i) {
    g = lvtx_to_gvtx(i,myid,graph->dist);
    other = gwhere[myid][i];
    if (other < 2) {
      v = mydist[other][group[myid][i]]++;
      gvtxs[v] = g;
    }
  }

  /* re-synchronize the mydist vector */
  vtx_dlthread_minareduce(mydist[0],2*ngroup,graph->comm);


  /* work on splitting the comm and graph here *******************************/

  /* determine group assignements */
  if (myid < fnthreads) {
    mygroup = 0;
    hnthreads = fnthreads;
  } else {
    mygroup = 1;
    hnthreads = nthreads - fnthreads;
  }

  /* split the thread groups */
  hcomm = dlthread_comm_split(mygroup,2,graph->comm);

  DL_ASSERT_EQUALS((size_t)hnthreads,dlthread_get_nthreads(hcomm),"%zu");

  hmyid = dlthread_get_id(hcomm);


  /* start constructing each half here ***************************************/

  mygraph = halves[mygroup] = par_graph_create(hcomm);

  if (hmyid == 0) {
    /* lead threads only */
    mygraph->label = r_vtx_alloc(hnthreads);
    mygraph->group = r_pid_alloc(hnthreads);
    mygraph->ngroup = graph->ngroup;
  }

  /* rename vectors */
  graph->rename[myid] = vtx_alloc(graph->mynvtxs[myid]);

  dlthread_barrier(graph->comm);

  hnvtxs = 0;
  hnedges = 0;

  /* each thread is now responsible for building its own part of the new 
   * graph -- first count my edges and vertices */
  for (d=hmyid;d<ngroup;d+=hnthreads) {
    /* go through each distribution and copy it */
    for (i=mydist[mygroup][d];i<mydist[mygroup][d+1];++i) {
      g = gvtxs[i];
      v = gvtx_to_lvtx(g,graph->dist);
      o = gvtx_to_tid(g,graph->dist);

      ++hnvtxs;
      for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
        k = gadjncy[o][j];
        if (k < gmynvtxs[o]) {
          lvtx = k;
          nbrid = o;
        } else {
          lvtx = gvtx_to_lvtx(k,graph->dist);
          nbrid = gvtx_to_tid(k,graph->dist);
        }
        other = gwhere[nbrid][lvtx];
        if (other == mygroup) {
          ++hnedges;
        }
      }
    }
  }

  maxgvtx = vtx_dlthread_sumreduce(hnvtxs,hcomm);

  if (hmyid == 0) {
    graph_calc_dist(maxgvtx,hnthreads,&(mygraph->dist));
  }

  /* allocate vwgt and xadj */
  mygraph->mynvtxs[hmyid] = hnvtxs;
  mygraph->mynedges[hmyid] = hnedges;
  mygraph->xadj[hmyid] = adj_alloc(hnvtxs+1);
  mygraph->vwgt[hmyid] = wgt_alloc(hnvtxs);
  mygraph->label[hmyid] = vtx_alloc(hnvtxs);
  mygraph->adjncy[hmyid] = vtx_alloc(hnedges);
  mygraph->adjwgt[hmyid] = wgt_alloc(hnedges);
  mygraph->group[hmyid] = pid_alloc(hnvtxs);

  /* sync so we can write to rename */
  dlthread_barrier(graph->comm);

  hnvtxs = 0;
  for (d=hmyid;d<ngroup;d+=hnthreads) {
    for (i=mydist[mygroup][d];i<mydist[mygroup][d+1];++i) {
      g = gvtxs[i];
      v = gvtx_to_lvtx(g,graph->dist);
      o = gvtx_to_tid(g,graph->dist);

      /* set rename */
      graph->rename[o][v] = lvtx_to_gvtx(hnvtxs,hmyid,mygraph->dist);
      ++hnvtxs;
    }
  }

  /* sync so we can read from rename */
  dlthread_barrier(graph->comm);

  /* populate edges */
  hnvtxs = 0;
  hnedges = 0;

  mygraph->xadj[hmyid][0] = 0;

  if (do_wgt) {
    for (d=hmyid;d<ngroup;d+=hnthreads) {
      for (i=mydist[mygroup][d];i<mydist[mygroup][d+1];++i) {
        g = gvtxs[i];
        v = gvtx_to_lvtx(g,graph->dist);
        o = gvtx_to_tid(g,graph->dist);

        mygraph->vwgt[hmyid][hnvtxs] = gvwgt[o][v];
        mygraph->label[hmyid][hnvtxs] = g;
        mygraph->group[hmyid][hnvtxs] = group[o][v];

        for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
          k = gadjncy[o][j];
          if (k < gmynvtxs[o]) {
            lvtx = k;
            nbrid = o;
          } else {
            lvtx = gvtx_to_lvtx(k,graph->dist);
            nbrid = gvtx_to_tid(k,graph->dist);
          }
          other = gwhere[nbrid][lvtx];
          if (other == mygroup) {
            /* calcuate new endpoint */
            u = graph->rename[nbrid][lvtx];
            t = gvtx_to_tid(u,mygraph->dist);
            if (t == hmyid) {
              k = gvtx_to_lvtx(u,mygraph->dist);
            } else {
              k = u;
            }
            
            mygraph->adjncy[hmyid][hnedges] = k;
            mygraph->adjwgt[hmyid][hnedges] = gadjwgt[o][j];
            ++hnedges;
          }
        }

        ++hnvtxs;
        mygraph->xadj[hmyid][hnvtxs] = hnedges;
      }
    }
  } else {
    for (d=hmyid;d<ngroup;d+=hnthreads) {
      for (i=mydist[mygroup][d];i<mydist[mygroup][d+1];++i) {
        g = gvtxs[i];
        v = gvtx_to_lvtx(g,graph->dist);
        o = gvtx_to_tid(g,graph->dist);

        mygraph->label[hmyid][hnvtxs] = g;
        mygraph->group[hmyid][hnvtxs] = group[o][v];

        for (j=gxadj[o][v];j<gxadj[o][v+1];++j) {
          k = gadjncy[o][j];
          if (k < gmynvtxs[o]) {
            lvtx = k;
            nbrid = o;
          } else {
            lvtx = gvtx_to_lvtx(k,graph->dist);
            nbrid = gvtx_to_tid(k,graph->dist);
          }
          other = gwhere[nbrid][lvtx];
          if (other == mygroup) {
            /* calcuate new endpoint */
            u = graph->rename[nbrid][lvtx];
            t = gvtx_to_tid(u,mygraph->dist);
            if (t == hmyid) {
              k = gvtx_to_lvtx(u,mygraph->dist);
            } else {
              k = u;
            }
            
            mygraph->adjncy[hmyid][hnedges] = k;
            ++hnedges;
          }
        }

        ++hnvtxs;
        mygraph->xadj[hmyid][hnvtxs] = hnedges;
      }
    }
    wgt_set(mygraph->vwgt[hmyid],1,hnvtxs);
    wgt_set(mygraph->adjwgt[hmyid],1,hnedges);
  }

  if (hmyid == 0) {
    /* lead threads only */
    mygraph->nvtxs = vtx_sum(mygraph->mynvtxs,hnthreads);
    mygraph->gnvtxs = max_gvtx(mygraph); 
    mygraph->nedges = adj_sum(mygraph->mynedges,hnthreads);
  }

  par_graph_setup_twgts(mygraph);

  dl_free(mydist[0]);
  dlthread_free_shmem(gvtxs,graph->comm);

  DL_ASSERT(check_graph(mygraph),"Bad graph extracted");

  return mygroup;
}




/******************************************************************************
* PUBLIC SERIAL FUNCTIONS *****************************************************
******************************************************************************/


graph_t * graph_create(
    tid_t const nthreads)
{
  graph_t * graph;

  graph = (graph_t*)calloc(1,sizeof(graph_t));

  __graph_init(graph,nthreads);

  return graph; 
}


graph_t * graph_setup(
    vtx_t * const nvtxs, 
    adj_t ** const xadj, 
    vtx_t ** const adjncy, 
    wgt_t ** const adjwgt, 
    wgt_t ** const vwgt,
    tid_t const nthreads) 
{
  tid_t myid;
  graph_t * graph;

  graph = graph_create(0);

  graph->mynvtxs = nvtxs;
  graph->xadj = xadj;
  graph->adjncy = adjncy;
  if (adjwgt) {
    graph->adjwgt = adjwgt;
  } else {
    graph->adjwgt = r_wgt_alloc(nthreads);
    for (myid=0;myid<nthreads;++myid) {
      graph->adjwgt[myid] = wgt_init_alloc(1,xadj[myid][nvtxs[myid]]);
    }
    graph->uniformadjwgt = 1;
    graph->free_adjwgt = 1;
  }
  if (vwgt) {
    graph->vwgt = vwgt;
  } else {
    graph->vwgt = r_wgt_alloc(nthreads);
    for (myid=0;myid<nthreads;++myid) {
      graph->vwgt[myid] = wgt_init_alloc(1,nvtxs[myid]);
    }
    graph->uniformvwgt = 1;
    graph->free_vwgt = 1;
  }

  for (myid=0;myid<nthreads;++myid) {
    graph->mynedges[myid] = xadj[myid][nvtxs[myid]];
  }

  graph->nvtxs = vtx_sum(graph->mynvtxs,nthreads);
  graph->nedges = adj_sum(graph->mynedges,nthreads);

  graph_calc_dist(vtx_max_value(graph->mynvtxs,nthreads),nthreads,&(graph->dist));

  graph_setup_twgts(graph);

  DL_ASSERT(check_graph(graph),"Bad graph");

  return graph;
}


graph_t * graph_distribute(
    int const distribution,
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt,
    wgt_t const * const adjwgt, 
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

  DL_ASSERT(nvtxs>0, "distribute_graph() called with nvtxs = %"PF_VTX_T"\n", \
      nvtxs);

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

  graph_calc_dist(vtx_max_value(dmynvtxs,nthreads),nthreads,&dist);

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
  if (vwgt) {
    for (myid=0;myid<nthreads;++myid) {
      graph->tvwgt += wgt_sum(graph->vwgt[myid],graph->mynvtxs[myid]);
    }
  } else {
    graph->tvwgt = graph->nvtxs;
    graph->uniformvwgt = 1;
  }
  graph->invtvwgt = 1.0/(graph->tvwgt > 0 ? graph->tvwgt : 1);

  /* setup tadjwgt */
  graph->tadjwgt = 0;
  if (adjwgt) {
    for (myid=0;myid<nthreads;++myid) {
      graph->tadjwgt += wgt_sum(graph->adjwgt[myid],dmynedges[myid]);
    }
  } else {
    graph->tadjwgt = graph->nedges;
    graph->uniformadjwgt = 1;
  }

  /* set free configuration */
  graph->free_xadj = 1;
  graph->free_adjncy = 1;
  graph->free_adjwgt = 1;
  graph->free_vwgt = 1;

  DL_ASSERT(check_graph(graph),"Bad graph");

  return graph;
}


void graph_gather(
  graph_t const * const graph,
  adj_t ** const r_xadj,
  vtx_t ** const r_adjncy,
  wgt_t ** const r_vwgt,
  wgt_t ** const r_adjwgt,
  vtx_t ** const r_voff)
{
  vtx_t i, k, v, g;
  adj_t j, l;
  tid_t t, myid;
  adj_t * gxadj;
  vtx_t * gadjncy, * prefix;
  wgt_t * gvwgt, * gadjwgt;

  tid_t const nthreads = graph->dist.nthreads;

  vtx_t const * const mynvtxs = graph->mynvtxs;
  adj_t const * const mynedges = graph->mynedges;
  adj_t const * const * const xadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const adjncy = (vtx_t const **)graph->adjncy;
  wgt_t const * const * const vwgt = (wgt_t const **)graph->vwgt;
  wgt_t const * const * const adjwgt = (wgt_t const **)graph->adjwgt;

  int const do_vwgt = !graph->uniformvwgt;
  int const do_adjwgt = !graph->uniformadjwgt;

  /* DL_ASSERT_EQUALS(mynedges,xadj[mynvtxs],"%"PF_ADJ_T); */
  DL_ASSERT_EQUALS(adj_sum(graph->mynedges,nthreads),graph->nedges, \
      "%"PF_ADJ_T);

  gxadj = adj_alloc(graph->nvtxs+1);
  gadjncy = vtx_alloc(graph->nedges);
  gvwgt = wgt_alloc(graph->nvtxs);
  gadjwgt = wgt_alloc(graph->nedges);

  /* create the vertex offsets */
  prefix = vtx_duplicate(mynvtxs,nthreads);
  vtx_prefixsum_exc(prefix,nthreads);

  /* vertex ids are purely based on thread offsets */
  gxadj[0] = 0;
  g = 0;
  l = 0;
  for (myid=0;myid<nthreads;++myid) {
    if (do_vwgt) {
      wgt_copy(gvwgt+g,vwgt[myid],mynvtxs[myid]);
    } else {
      wgt_set(gvwgt+g,1,mynvtxs[myid]);
    }
    if (do_adjwgt) {
      wgt_copy(gadjwgt+l,adjwgt[myid],mynedges[myid]);
    } else {
      wgt_set(gadjwgt+l,1,mynedges[myid]);
    }
    /* insert edges into graph */
    for (i=0;i<mynvtxs[myid];++i) {
      for (j=xadj[myid][i];j<xadj[myid][i+1];++j) {
        k = adjncy[myid][j];
        if (k < mynvtxs[myid]) {
          t = myid;
          v = k;
        } else {
          t = gvtx_to_tid(k,graph->dist);
          v = gvtx_to_lvtx(k,graph->dist);
        }
        gadjncy[l++] = v + prefix[t];
      }
      gxadj[++g] = l;
    }
  }

  /* assign pointers */
  *r_xadj = gxadj;
  *r_adjncy = gadjncy;
  *r_vwgt = gvwgt;
  *r_adjwgt = gadjwgt;

  if (r_voff) {
    *r_voff = prefix;
  } else {
    dl_free(prefix);
  }

  DL_ASSERT(bowstring_check_graph(graph->nvtxs,gxadj,gadjncy, \
        (bowstring_wgt_t*)gadjwgt),"Bad graph gathered");
}


graph_t * graph_setup_coarse(
    graph_t * const graph, 
    vtx_t * const cnvtxs)
{
  graph_t * cgraph;
  tid_t myid;

  tid_t const nthreads = graph->dist.nthreads;

  cgraph = graph_create(nthreads);

  graph->coarser = cgraph;
  cgraph->finer = graph;

  cgraph->level = graph->level + 1;

  DL_ASSERT_EQUALS(nthreads,graph->dist.nthreads,"%"PF_TID_T);

  graph_calc_dist(vtx_max_value(cgraph->mynvtxs,nthreads),nthreads, \
      &cgraph->dist);

  cgraph->tvwgt = graph->tvwgt;
  cgraph->invtvwgt = graph->invtvwgt;

  DL_ASSERT(cgraph != NULL,"cgraph is NULL");

  for (myid=0;myid<nthreads;++myid) {
    cgraph->mynvtxs[myid] = cnvtxs[myid];

    cgraph->xadj[myid] = adj_alloc(cnvtxs[myid]+1);
    if (cgraph->mynvtxs[myid] > 0) {
      cgraph->vwgt[myid] = wgt_alloc(cnvtxs[myid]);
    } else {
      cgraph->xadj[myid][0] = 0;
      cgraph->vwgt[myid] = NULL;
    }

    cgraph->adjncy[myid] = NULL;
    cgraph->adjwgt[myid] = NULL;
  }

  cgraph->gnvtxs = cgraph->gnvtxs;
  cgraph->nvtxs = vtx_sum(cgraph->mynvtxs,nthreads);

  DL_ASSERT(cgraph->gnvtxs >= cgraph->nvtxs,"Bad gnvtxs");

  return cgraph;
}


void graph_setup_twgts(
    graph_t * const graph)
{
  vtx_t i;
  adj_t j;
  twgt_t vsum,asum;
  tid_t myid;

  tid_t const nthreads = graph->dist.nthreads;

  if (graph->uniformvwgt) {
    vsum = graph->nvtxs;
  } else {
    vsum = 0;
    for (myid=0;myid<nthreads;++myid) {
      for (i=0;i<graph->mynvtxs[myid];++i) {
        vsum += graph->vwgt[myid][i];
      }
    }
  }

  if (graph->uniformadjwgt) {
    asum = graph->nedges;
  } else {
    asum = 0;
    for (myid=0;myid<nthreads;++myid) {
      for (j=0;j<graph->mynedges[myid];++j) {
        asum += graph->adjwgt[myid][j];
      }
    }
  }

  graph->tvwgt = vsum;
  graph->tadjwgt = asum;
  graph->invtvwgt = 1.0/(graph->tvwgt > 0 ? graph->tvwgt : 1);
}


void graph_alloc_partmemory(
    ctrl_t * const ctrl,
    graph_t * const graph)
{
  tid_t myid;

  tid_t const nthreads = graph->dist.nthreads;

  /* memory for the partition/refinement structure */
  graph->where = r_pid_alloc(nthreads);
  graph->pwgts = wgt_alloc(ctrl->nparts);

  for (myid=0;myid<nthreads;++myid) {
    graph->where[myid] = pid_alloc(graph->mynvtxs[myid]);
  }
}


void graph_free(
    graph_t * graph)
{
  tid_t myid;

  /* free partition/refinement structure */
  graph_free_rdata(graph);

  /* free individual arrays */
  for (myid=0;myid<graph->dist.nthreads;++myid) {
    __graph_free_part(graph,myid);
  }

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
  if (graph->group) {
    dl_free(graph->group);
  }

  dl_free(graph);
}


void graph_free_rdata(
    graph_t * graph)
{
  tid_t myid;
  
  for (myid=0;myid<graph->dist.nthreads;++myid) {
    if (graph->where) {
      dl_free(graph->where[myid]);
    }
    if (graph->rename) {
      dl_free(graph->rename[myid]);
    }
    if (graph->cmap) {
      dl_free(graph->cmap[myid]);
    }
  }

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
  if (graph->cmap) {
    dl_free(graph->cmap);
    graph->cmap = NULL;
  }
  if (graph->vsinfo) {
    vsinfo_free(graph);
  }
  if (graph->esinfo) {
    esinfo_free(graph);
  }
  if (graph->kwinfo) {
    par_kwinfo_free(graph);
  }
}


double graph_imbalance(
    graph_t const * const graph,
    pid_t const nparts,
    real_t const * const pijbm)
{
  vtx_t k;
  double max, cur;

  DL_ASSERT_EQUALS(wgt_lsum(graph->pwgts,nparts),graph->tvwgt,"%"PF_TWGT_T);

  max = 0;

  for (k=0;k<nparts;++k) {
    cur = graph->pwgts[k]*pijbm[k];
    if (cur > max) {
      max = cur;
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

  DL_ASSERT_EQUALS(wgt_lsum(graph->pwgts,nparts),graph->tvwgt,"%"PF_TWGT_T);

  max = 0;

  for (k =0;k<nparts;++k) {
    cur = graph->pwgts[k]*pijbm[k]-ubfactor;
    if (cur > max) {
      max = cur;
    }
  }

  return max;
}


wgt_t graph_cut(
    graph_t const * const graph,
    pid_t const * const * const gwhere)
{
  vtx_t i, k, lvtx, nbrid;
  adj_t j;
  wgt_t cut;
  tid_t myid;

  tid_t const nthreads = dlthread_get_nthreads(graph->comm);

  vtx_t const * const gmynvtxs = graph->mynvtxs;
  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;
  wgt_t const * const * const gadjwgt = (wgt_t const **)graph->adjwgt;

  DL_ASSERT_EQUALS((int)graph->dist.nthreads,(int)nthreads,"%d");

  cut = 0;
  if (graph->uniformadjwgt || graph->adjwgt == NULL) {
    for (myid=0;myid<nthreads;++myid) {
      for (i=0;i<gmynvtxs[myid];++i) {
        for (j=gxadj[myid][i];j<gxadj[myid][i+1];++j) {
          k = gadjncy[myid][j]; 
          if (k < gmynvtxs[myid]) {
            lvtx = k;
            nbrid = myid;
          } else {
            nbrid = gvtx_to_tid(gadjncy[myid][j],graph->dist);
            lvtx = gvtx_to_lvtx(gadjncy[myid][j],graph->dist);
          }
          if (gwhere[myid][i] != gwhere[nbrid][lvtx]) {
            ++cut;
          }
        }
      }
    }
  } else {
    for (myid=0;myid<nthreads;++myid) {
      for (i=0;i<gmynvtxs[myid];++i) {
        for (j=gxadj[myid][i];j<gxadj[myid][i+1];++j) {
          k = gadjncy[myid][j]; 
          if (k < gmynvtxs[myid]) {
            lvtx = k;
            nbrid = myid;
          } else {
            nbrid = gvtx_to_tid(gadjncy[myid][j],graph->dist);
            lvtx = gvtx_to_lvtx(gadjncy[myid][j],graph->dist);
          }
          if (gwhere[myid][i] != gwhere[nbrid][lvtx]) {
            cut += gadjwgt[myid][j];
          }
        }
      }
    }
  }

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
  int const myid = dlthread_get_id(graph->comm);
  adj_t const nedges = graph->xadj[myid][graph->mynvtxs[myid]];

  if (adjsize > 4096 && adjsize * 0.75 > nedges) {
    graph->adjncy[myid] = vtx_realloc(graph->adjncy[myid],nedges);
    graph->adjwgt[myid] = wgt_realloc(graph->adjwgt[myid],nedges);
  }
}


size_t graph_size(
    graph_t const * const graph)
{
  size_t nbytes;

  tid_t const nthreads = graph->dist.nthreads;
  vtx_t const nvtxs = graph->nvtxs;
  adj_t const nedges = graph->nedges;

  nbytes = sizeof(graph_t);

  if (graph->group) {
    nbytes += (sizeof(pid_t)*nvtxs) + (sizeof(pid_t*)*nthreads);
  }

  if (graph->mynvtxs) {
    nbytes += sizeof(vtx_t)*nthreads;
  }

  if (graph->mynedges) {
    nbytes += sizeof(adj_t)*nthreads;
  }

  if (graph->xadj) {
    nbytes += (sizeof(adj_t)*(nvtxs+1)) + (sizeof(adj_t*)*nthreads);
  }

  if (graph->vwgt) {
    nbytes += (sizeof(wgt_t)*nvtxs) + (sizeof(wgt_t*)*nthreads);
  }

  if (graph->adjncy) {
    nbytes += (sizeof(vtx_t)*nedges) + (sizeof(vtx_t*)*nthreads);
  }

  if (graph->adjwgt) {
    nbytes += (sizeof(wgt_t)*nedges) + (sizeof(wgt_t*)*nthreads);
  }

  if (graph->cmap) {
    nbytes += (sizeof(vtx_t)*nvtxs) + (sizeof(vtx_t*)*nthreads);
  }

  /* where and pwgts will be ignored for now as they should be moved to a
   * different structure */ 

  if (graph->nislands) {
    nbytes += sizeof(vtx_t)*nthreads;
  }

  if (graph->rename) {
    nbytes += (sizeof(vtx_t)*nvtxs) + (sizeof(vtx_t*)*nthreads);
  }

  if (graph->label) {
    nbytes += (sizeof(vtx_t)*nvtxs) + (sizeof(vtx_t*)*nthreads);
  }

  return nbytes;
}


void ser_graph_extract_halves(
    graph_t * const graph,
    pid_t const * const * const gwhere,
    graph_t ** const halves)
{
  vtx_t i, k, u;
  adj_t j, l;
  pid_t w, side;
  vtx_t hmynvtxs[2];
  adj_t hmynedges[2];

  vtx_t const mynvtxs = graph->mynvtxs[0];
  adj_t const * const xadj = graph->xadj[0];
  vtx_t const * const adjncy = graph->adjncy[0];
  wgt_t const * const vwgt = graph->vwgt[0];
  wgt_t const * const adjwgt = graph->adjwgt[0];

  pid_t const * const where = gwhere[0];

  int const do_vwgt = !graph->uniformvwgt;
  int const do_adjwgt = !graph->uniformadjwgt;

  DL_ASSERT_EQUALS(graph->dist.nthreads,1,"%"PF_TID_T);

  for (side=0;side<2;++side) {
    halves[side] = graph_create(1);
    halves[side]->comm = graph->comm;
    halves[side]->label = r_vtx_alloc(1);
  }

  /* rename vectors */
  graph->rename = r_vtx_alloc(1);
  graph->rename[0] = vtx_alloc(graph->mynvtxs[0]);

  /* count vertices and edges */
  hmynvtxs[0] = hmynvtxs[1] = 0;
  hmynedges[0] = hmynedges[1] = 0;
  for (i=0;i<mynvtxs;++i) {
    w = where[i];
    if (w < 2) {
      ++hmynvtxs[w];
      /* over estimates */
      hmynedges[w] += xadj[i+1] - xadj[i];
    }
  }

  /* allocate vertex arrays */
  halves[0]->xadj[0] = adj_alloc(hmynvtxs[0]+2);
  halves[1]->xadj[0] = adj_alloc(hmynvtxs[1]+2);
  halves[0]->vwgt[0] = wgt_alloc(hmynvtxs[0]);
  halves[1]->vwgt[0] = wgt_alloc(hmynvtxs[1]);
  halves[0]->label[0] = vtx_alloc(hmynvtxs[0]);
  halves[1]->label[0] = vtx_alloc(hmynvtxs[1]);

  /* allocate edge arrays */
  halves[0]->adjncy[0] = vtx_alloc(hmynedges[0]);
  halves[1]->adjncy[0] = vtx_alloc(hmynedges[1]);
  halves[0]->adjwgt[0] = wgt_alloc(hmynedges[0]);
  halves[1]->adjwgt[0] = wgt_alloc(hmynedges[1]);

  /* insert vertices and edges into graphs */
  hmynvtxs[0] = hmynvtxs[1] = 0;
  halves[0]->xadj[0][0] = halves[1]->xadj[0][0] = 0;
  for (i=0;i<mynvtxs;++i) {
    w = where[i];
    if (w < 2) {
      u = hmynvtxs[w]++;
      if (do_vwgt) {
        halves[w]->vwgt[0][u] = vwgt[u];
      }

      /* aliases */
      halves[w]->label[0][u] = lvtx_to_gvtx(i,0,graph->dist);
      graph->rename[0][i] = u;

      /* handl xadj */
      l = halves[w]->xadj[0][u];
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        if (where[k] == w) {
          halves[w]->adjncy[0][l] = k;
          if (do_adjwgt) {
            halves[w]->adjwgt[0][l] = adjwgt[j];
          }
          ++l;
        }
      }
      halves[w]->xadj[0][u+1] = l;
    }
  }
  hmynedges[0] = halves[0]->xadj[0][hmynvtxs[0]];
  hmynedges[1] = halves[1]->xadj[0][hmynvtxs[1]];

  if (!do_vwgt) {
    wgt_set(halves[0]->vwgt[0],1,hmynvtxs[0]);
    wgt_set(halves[1]->vwgt[0],1,hmynvtxs[1]);
    halves[0]->uniformvwgt = 1;
    halves[1]->uniformvwgt = 1;
  }

  if (!do_adjwgt) {
    wgt_set(halves[0]->adjwgt[0],1,hmynedges[0]);
    wgt_set(halves[1]->adjwgt[0],1,hmynedges[1]);
    halves[0]->uniformadjwgt = 1;
    halves[1]->uniformadjwgt = 1;
  }

  /* rename vertices */
  for (side=0;side<2;++side) {
    for (j=0;j<halves[side]->xadj[0][hmynvtxs[side]];++j) {
      k = halves[side]->adjncy[0][j];
      halves[side]->adjncy[0][j] = graph->rename[0][k];
    }
  }

  for (side=0;side<2;++side) {
    halves[side]->nvtxs = halves[side]->mynvtxs[0] = hmynvtxs[side];
    halves[side]->nedges = halves[side]->mynedges[0] = hmynedges[side];

    graph_calc_dist(hmynvtxs[side],1,&halves[side]->dist);

    halves[side]->gnvtxs = max_gvtx(halves[side]);

    graph_setup_twgts(halves[side]);

    DL_ASSERT_EQUALS(check_graph(halves[side]),1,"%d");
  }
}


void graph_calc_dist(
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




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


graph_t * par_graph_create(
    dlthread_comm_t const comm)
{
  graph_t * graph;

  graph = dlthread_get_shmem(sizeof(graph_t),comm);

  __par_graph_init(graph,comm);

  return graph; 
}


graph_t * par_graph_setup(
    vtx_t const nvtxs, 
    adj_t * const xadj, 
    vtx_t * const adjncy, 
    wgt_t * const vwgt,
    wgt_t * const adjwgt, 
    dlthread_comm_t const comm) 
{
  wgt_t asum, vsum;
  graph_t * graph;

  tid_t const myid = dlthread_get_id(comm);
  tid_t const nthreads = dlthread_get_nthreads(comm);

  graph = par_graph_create(comm);


  graph->mynvtxs[myid] = nvtxs;
  graph->mynedges[myid] = xadj[nvtxs];
  graph->xadj[myid] = xadj;
  graph->adjncy[myid] = adjncy;
  if (adjwgt) {
    graph->adjwgt[myid] = adjwgt;
    asum = wgt_sum(graph->adjwgt[myid],graph->mynedges[myid]);
    asum = wgt_dlthread_sumreduce(asum,graph->comm);
    graph->free_adjwgt = 0;
  } else {
    asum = xadj[nvtxs];
    graph->adjwgt[myid] = wgt_init_alloc(1,xadj[nvtxs]);
    graph->uniformadjwgt = 1;
  }
  if (vwgt) {
    graph->vwgt[myid] = vwgt;
    vsum = wgt_sum(graph->vwgt[myid],graph->mynvtxs[myid]);
    vsum = wgt_dlthread_sumreduce(vsum,graph->comm);
    graph->free_vwgt = 0;
  } else {
    vsum = nvtxs;
    graph->vwgt[myid] = wgt_init_alloc(1,nvtxs);
    graph->uniformvwgt = 1;
  }


  dlthread_barrier(comm);
  if (myid == 0) {
    graph->nvtxs = vtx_sum(graph->mynvtxs,nthreads);
    graph->nedges = adj_sum(graph->mynedges,nthreads);
    graph_calc_dist(vtx_max_value(graph->mynvtxs,nthreads),nthreads, \
        &(graph->dist));

    graph->gnvtxs = lvtx_to_gvtx(vtx_max_value(graph->mynvtxs,nthreads), \
        nthreads-1,graph->dist);

    if (graph->free_adjwgt) {
      /* we have all 1's for edge weight */
      graph->tadjwgt = graph->nedges;
    }
    if (graph->free_vwgt) {
      /* we have all 1's for vertex wegiht */
      graph->tvwgt = graph->nvtxs;
    }
    graph->invtvwgt = 1.0/(graph->tvwgt > 0 ? graph->tvwgt : 1);
  }
  dlthread_barrier(comm);

  par_graph_setup_twgts(graph);

  DL_ASSERT(check_graph(graph),"Bad graph");

  return graph;
}


void par_graph_setup_twgts(
    graph_t * const graph)
{
  vtx_t i;
  adj_t j;
  twgt_t vsum, asum;

  tid_t const myid = dlthread_get_id(graph->comm);

  if (graph->uniformvwgt) {
    vsum = graph->nvtxs;
  } else {
    vsum = 0;
    for (i=0;i<graph->mynvtxs[myid];++i) {
      vsum += graph->vwgt[myid][i];
    }
    vsum = twgt_dlthread_sumreduce(vsum,graph->comm);
  }
  if (graph->uniformadjwgt) {
    asum = graph->nedges;
  } else {
    asum = 0;
    for (j=0;j<graph->mynedges[myid];++j) {
      asum += graph->adjwgt[myid][j];
    }
    asum = twgt_dlthread_sumreduce(asum,graph->comm);
  }

  if (myid == 0) {
    graph->tvwgt = vsum;
    graph->tadjwgt = asum;
    graph->invtvwgt = 1.0/(graph->tvwgt > 0 ? graph->tvwgt : 1);
  }
  dlthread_barrier(graph->comm);
}


void par_graph_free(
    graph_t * graph)
{
  tid_t const myid = dlthread_get_id(graph->comm);

  __graph_free_part(graph,myid);

  /* free partition/refinement structure */
  par_graph_free_rdata(graph);

  dlthread_barrier(graph->comm);
  if (myid == 0) {
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
    if (graph->group) {
      dl_free(graph->group);
    }

    dl_free(graph);
  }
}


void par_graph_free_rdata(
    graph_t * graph)
{
  tid_t const myid = dlthread_get_id(graph->comm);

  if (graph->where) {
    dl_free(graph->where[myid]);
  }
  if (graph->rename) {
    dl_free(graph->rename[myid]);
  }
  if (graph->cmap) {
    dl_free(graph->cmap[myid]);
  }

  if (graph->vsinfo) {
    par_vsinfo_free(graph);
  }
  if (graph->esinfo) {
    par_esinfo_free(graph);
  }
  if (graph->kwinfo) {
    par_kwinfo_free(graph);
  }

  dlthread_barrier(graph->comm);

  if (myid == 0) {
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
    if (graph->cmap) {
      dl_free(graph->cmap);
      graph->cmap = NULL;
    }
  }
}


void par_graph_gather(
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

  /* unified graph parts */
  adj_t * uxadj;
  vtx_t * uadjncy;
  wgt_t * uvwgt;
  wgt_t * uadjwgt;
  vtx_t * uprefix;
  adj_t * unedges;

  tid_t const myid = dlthread_get_id(graph->comm);

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const mynedges = graph->mynedges[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const vwgt = graph->vwgt[myid];
  wgt_t const * const adjwgt = graph->adjwgt[myid];

  int const do_vwgt = !graph->uniformvwgt;
  int const do_adjwgt = !graph->uniformadjwgt;

  DL_ASSERT_EQUALS(mynedges,xadj[mynvtxs],"%"PF_ADJ_T);
  DL_ASSERT_EQUALS(adj_sum(graph->mynedges,graph->dist.nthreads), \
      graph->nedges,"%"PF_ADJ_T);

  __par_graph_alloc_unified(graph,&uxadj,&uadjncy,&uvwgt,&uadjwgt,&uprefix, \
      &unedges);

  voff = uprefix[myid];
  eoff = unedges[myid];
  gxadj = uxadj + voff;
  gadjncy = uadjncy + eoff;
  gvwgt = uvwgt + voff;
  gadjwgt = uadjwgt + eoff;

  /* vertex ids are purely based on thread offsets */
  eoff = gxadj[0] = unedges[myid];
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
      gadjncy[j] = v + uprefix[t];
    }
  }

  /* propagate weights */
  if (do_vwgt) {
    wgt_copy(gvwgt,vwgt,mynvtxs);
  } else {
    wgt_set(gvwgt,1,mynvtxs);
  }
  if (do_adjwgt) {
    wgt_copy(gadjwgt,adjwgt,mynedges);
  } else {
    wgt_set(gadjwgt,1,mynedges);
  }

  dlthread_barrier(graph->comm);
  if (myid == 0) {
    dl_free(unedges);
    dl_free(uprefix);
  }

  /* assign pointers */
  *r_xadj = uxadj;
  *r_adjncy = uadjncy;
  *r_vwgt = uvwgt;
  *r_adjwgt = uadjwgt;
  *r_voff = voff;

  DL_ASSERT(bowstring_check_graph(graph->nvtxs,uxadj,uadjncy, \
        (bowstring_wgt_t*)uadjwgt),"Bad graph gathered");
}


void par_graph_shuffle(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t const * const * const gwhere,
    int const wgts)
{
  vtx_t v, g, i, k, l, lvtx, smynvtxs, maxnvtxs;
  adj_t j, smynedges;
  tid_t nbrid, t, o;
  pid_t me;
  graphdist_t dist;
  pid_t * group = NULL;
  vtx_t * adjncy, * label, * vold;
  adj_t * myeprefix, * xadj;
  wgt_t * vwgt, * adjwgt;
  vtx_t ** vprefix, ** vlist, ** grename, ** myvprefix;
  adj_t ** eprefix;

  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);
  tid_t const myid = dlthread_get_id(ctrl->comm);

  vtx_t const * const gmynvtxs = graph->mynvtxs;
  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy; 
  wgt_t const * const * const gvwgt = (wgt_t const **)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const **)graph->adjwgt;
  vtx_t const * const * const glabel = (vtx_t const **)graph->label;

  /* should change this to intelligently think about the weights  --
  int const do_vwgt = !graph->uniformvwgt;
  int const do_adjwgt = !graph->uniformadjwgt; */

  vprefix = dlthread_get_shmem((sizeof(vtx_t*)*nthreads) + \
      (sizeof(adj_t*)*nthreads) + (sizeof(vtx_t*)*nthreads) + \
      (sizeof(vtx_t*)*nthreads) + (sizeof(vtx_t*)*nthreads),ctrl->comm);
  eprefix = (adj_t**)(vprefix+nthreads);
  vlist = (vtx_t**)(eprefix+nthreads);
  grename = (vtx_t**)(vlist+nthreads);
  myvprefix = (vtx_t**)(grename+nthreads);

  vprefix[myid] = vtx_init_alloc(0,nthreads+1);
  eprefix[myid] = adj_init_alloc(0,nthreads);
  myvprefix[myid] = vtx_alloc(nthreads+1);

  vprefix[myid][nthreads] = 0;

  vlist[myid] = vtx_alloc(gmynvtxs[myid]);
  grename[myid] = vtx_alloc(gmynvtxs[myid]);

  /* We first go through and count how many vertices and edges each thread is 
   * going to own. */
  for (i=0;i<gmynvtxs[myid];++i) {
    me = gwhere[myid][i]; 
    ++vprefix[myid][me];
    eprefix[myid][me] += gxadj[myid][i+1] - gxadj[myid][i];
  }

  myeprefix = adj_alloc(nthreads+1);
  vold = vtx_alloc(nthreads);

  dlthread_barrier(ctrl->comm);

  /* each thread then gathers its totals */
  for (o=0;o<nthreads;++o) {
    t = (o + myid) % nthreads;
    myvprefix[myid][t] = vprefix[t][myid];
    myeprefix[t] = eprefix[t][myid];
  }
  myvprefix[myid][nthreads] = 0;
  myeprefix[nthreads] = 0;

  vtx_prefixsum_exc(myvprefix[myid],nthreads+1);
  adj_prefixsum_exc(myeprefix,nthreads+1);

  smynvtxs = myvprefix[myid][nthreads];
  smynedges = myeprefix[nthreads];

  /* implicit barrier */
  maxnvtxs = vtx_dlthread_maxreduce_value(smynvtxs,ctrl->comm);
  graph_calc_dist(maxnvtxs,nthreads,&dist);

  /* thread creates an incoming and outgoing prefixsum */
  vtx_prefixsum_exc(vprefix[myid],nthreads+1);

  /* copy our initial offsets */
  vtx_copy(vold,vprefix[myid],nthreads);

  /* create outgoing vlist and rename vectors */
  for (i=0;i<gmynvtxs[myid];++i) {
    me = gwhere[myid][i]; 
    l = myvprefix[me][myid] + (vprefix[myid][me] - vold[me]);
    DL_ASSERT(l<myvprefix[me][nthreads],"Bad local vertex number: %"PF_VTX_T \
        "/%"PF_VTX_T"\n",l,myvprefix[me][nthreads]);
    /* the global vertex is determined by the myvprefix of the destination
     * thread */
    grename[myid][i] = lvtx_to_gvtx(l,me,dist);
    vlist[myid][vprefix[myid][me]++] = i;
  }

  /* de-shift vprefix */
  for (t=nthreads;t>0;--t) {
    vprefix[myid][t] = vprefix[myid][t-1];
  }
  vprefix[myid][0] = 0;

  dlthread_barrier(ctrl->comm);

  /* allocate arrays for my new graph */
  xadj = adj_alloc(smynvtxs+1);
  adjncy = vtx_alloc(smynedges);
  label = vtx_alloc(smynvtxs);

  /* copy group information */
  if (graph->group) {
    group = pid_alloc(smynvtxs);
    smynvtxs = 0;
    for (t=0;t<nthreads;++t) {
      for (i=vprefix[t][myid];i<vprefix[t][myid+1];++i) {
        v = vlist[t][i];
        group[smynvtxs] = graph->group[t][v];
        ++smynvtxs;
      }
    }
  }

  /* build my graph using the prefix sum arrays */
  if (wgts) {
    vwgt = wgt_alloc(smynvtxs);
    adjwgt = wgt_alloc(smynedges);
    smynvtxs = 0;
    smynedges = 0;
    for (t=0;t<nthreads;++t) {
      DL_ASSERT_EQUALS(smynvtxs,myvprefix[myid][t],"%"PF_VTX_T);
      for (i=vprefix[t][myid];i<vprefix[t][myid+1];++i) {
        v = vlist[t][i];
        label[smynvtxs] = glabel[t][v];
        vwgt[smynvtxs] = gvwgt[t][v];
        xadj[smynvtxs] = smynedges;
        for (j=gxadj[t][v];j<gxadj[t][v+1];++j) {
          k = gadjncy[t][j];
          if (k < gmynvtxs[t]) {
            lvtx = k;
            nbrid = t;
          } else {
            lvtx = gvtx_to_lvtx(k,graph->dist);
            nbrid = gvtx_to_tid(k,graph->dist);
          }
          /* the rename array contains the global vertex number */
          g = grename[nbrid][lvtx];
          if (gvtx_to_tid(g,dist) == myid) {
            adjncy[smynedges] = gvtx_to_lvtx(g,dist);
          } else {
            adjncy[smynedges] = g;
          }
          adjwgt[smynedges] = gadjwgt[nbrid][lvtx];
          ++smynedges;
        }
        ++smynvtxs;
      }
    }
  } else {
    vwgt = NULL;
    adjwgt = NULL;
    smynvtxs = 0;
    smynedges = 0;
    for (t=0;t<nthreads;++t) {
      DL_ASSERT_EQUALS(smynvtxs,myvprefix[myid][t],"%"PF_VTX_T);
      for (i=vprefix[t][myid];i<vprefix[t][myid+1];++i) {
        v = vlist[t][i];
        label[smynvtxs] = glabel[t][v];
        xadj[smynvtxs] = smynedges;
        for (j=gxadj[t][v];j<gxadj[t][v+1];++j) {
          k = gadjncy[t][j];
          if (k < gmynvtxs[t]) {
            lvtx = k;
            nbrid = t;
          } else {
            lvtx = gvtx_to_lvtx(k,graph->dist);
            nbrid = gvtx_to_tid(k,graph->dist);
          }
          /* the rename array contains the global vertex number */
          g = grename[nbrid][lvtx];
          if (gvtx_to_tid(g,dist) == myid) {
            adjncy[smynedges] = gvtx_to_lvtx(g,dist);
          } else {
            adjncy[smynedges] = g;
          }
          ++smynedges;
        }
        ++smynvtxs;
      }
    }
  }
  xadj[smynvtxs] = smynedges;

  dlthread_barrier(ctrl->comm);

  /* free intermediate data */
  dl_free(vold);
  dl_free(myeprefix);
  dl_free(myvprefix[myid]);
  dl_free(grename[myid]);
  dl_free(vprefix[myid]);
  dl_free(eprefix[myid]);
  dl_free(vlist[myid]);

  /* free the old graph components */
  dl_free(graph->xadj[myid]);
  dl_free(graph->adjncy[myid]);
  dl_free(graph->vwgt[myid]);
  dl_free(graph->adjwgt[myid]);
  dl_free(graph->label[myid]);

  /* set the new components */
  graph->mynvtxs[myid] = smynvtxs;
  graph->mynedges[myid] = smynedges;
  graph->xadj[myid] = xadj;
  graph->adjncy[myid] = adjncy;
  if (vwgt) {
    graph->vwgt[myid] = vwgt;
  } else {
    graph->vwgt[myid] = wgt_init_alloc(1,smynvtxs);
  }
  if (adjwgt) {
    graph->adjwgt[myid] = adjwgt;
  } else {
    graph->adjwgt[myid] = wgt_init_alloc(1,smynedges);
  }
  graph->label[myid] = label;
  graph->dist = dist;

  if (graph->group) {
    dl_free(graph->group[myid]);
    graph->group[myid] = group;
  }

  /* implicit barrier */
  dlthread_free_shmem(vprefix,ctrl->comm);

  DL_ASSERT(check_graph(graph),"Invalid graph after shuffle");
}


graph_t * par_graph_setup_coarse(
    graph_t * const graph, 
    vtx_t const cnvtxs)
{
  vtx_t mynvtxs;
  graph_t * cgraph;

  tid_t const myid = dlthread_get_id(graph->comm);
  tid_t const nthreads = dlthread_get_nthreads(graph->comm);

  cgraph = par_graph_create(graph->comm);

  if (myid == 0) {
    graph->coarser = cgraph;
    cgraph->finer = graph;

    cgraph->level = graph->level + 1;

    cgraph->tvwgt = graph->tvwgt;
    cgraph->invtvwgt = graph->invtvwgt;
  }
  dlthread_barrier(graph->comm);

  cgraph = graph->coarser;

  DL_ASSERT(cgraph != NULL,"cgraph is NULL");

  cgraph->mynvtxs[myid] = cnvtxs;

  /* Allocate memory for the coarser graph */
  mynvtxs = cnvtxs;

  cgraph->xadj[myid] = adj_alloc(mynvtxs+1);
  cgraph->vwgt[myid] = wgt_alloc(mynvtxs);

  cgraph->adjncy[myid] = NULL;
  cgraph->adjwgt[myid] = NULL;

  mynvtxs = vtx_dlthread_sumreduce(mynvtxs,graph->comm);

  if (myid == 0) {
    graph_calc_dist(vtx_max_value(cgraph->mynvtxs,nthreads),nthreads, \
        &cgraph->dist);

    cgraph->gnvtxs = max_gvtx(cgraph);
    cgraph->nvtxs = mynvtxs;
    DL_ASSERT(cgraph->gnvtxs >= cgraph->nvtxs,"Bad gnvtxs of %"PF_VTX_T"/%" \
        PF_VTX_T,cgraph->gnvtxs,cgraph->nvtxs);
    cgraph->comm = graph->comm;
  }
  dlthread_barrier(graph->comm);

  DL_ASSERT(cgraph->finer != NULL,"Failed to set cgraph->finer");
  DL_ASSERT(graph->coarser != NULL,"Failed to set graph->coarser");

  return cgraph;
}


void par_graph_alloc_partmemory(
    ctrl_t * const ctrl,
    graph_t * const graph)
{
  tid_t const nthreads = dlthread_get_nthreads(graph->comm);
  tid_t const myid = dlthread_get_id(graph->comm);

  DL_ASSERT_EQUALS(nthreads,graph->dist.nthreads,"%"PF_TID_T);

  if (myid == 0) {
    /* memory for the partition/refinement structure */
    graph->where = r_pid_alloc(nthreads);
    graph->pwgts = wgt_alloc(ctrl->nparts);
  }
  dlthread_barrier(graph->comm);

  graph->where[myid] = pid_alloc(graph->mynvtxs[myid]);
}


tid_t par_graph_extract_halves(
    graph_t * const graph,
    pid_t const * const * const gwhere,
    graph_t ** const halves)
{
  tid_t mygroup;

  if (graph->group) {
    mygroup = __par_graph_extract_halves_group(graph,gwhere,halves);
  } else {
    mygroup = __par_graph_extract_halves_nogroup(graph,gwhere,halves);
  }

  return mygroup;
}


adj_t * par_graph_build_radj(
    graph_t const * const graph)
{
  vtx_t i, k, kk, lvtx, olvtx;
  tid_t nbrid, onbrid;
  adj_t j, jj;
  adj_t * radj;

  tid_t const myid = dlthread_get_id(graph->comm);

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const nedges = graph->xadj[myid][mynvtxs];

  vtx_t const * const gmynvtxs = graph->mynvtxs;
  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;

  radj = adj_alloc(nedges);

  /* populate my radj */
  for (i=0;i<mynvtxs;++i) {
    for (j=gxadj[myid][i];j<gxadj[myid][i+1];++j) {
      k = gadjncy[myid][j];
      if (k < mynvtxs) {
        lvtx = k;
        nbrid = myid;
      } else {
        lvtx = gvtx_to_lvtx(k,graph->dist);
        nbrid = gvtx_to_tid(k,graph->dist);
      }
      for (jj=gxadj[nbrid][lvtx];jj<gxadj[nbrid][lvtx+1];++jj) {
        kk = gadjncy[nbrid][jj];
        if (kk < gmynvtxs[nbrid]) {
          olvtx = kk;
          onbrid = nbrid;
        } else {
          olvtx = gvtx_to_lvtx(kk,graph->dist);
          onbrid = gvtx_to_tid(kk,graph->dist);
        }
        if (onbrid == myid && olvtx == i) {
          radj[j] = jj;
          break;
        }
      }
    }
  }

  return radj;
}


void par_graph_intext_vtx(
    graph_t const * const graph,
    vtx_t * const r_nint,
    vtx_t * const r_next)
{
  vtx_t i, k;
  adj_t j;
  vtx_t next, nint;
  
  tid_t const myid = dlthread_get_id(graph->comm);

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];

  next = 0;
  nint = 0;

  for (i=0;i<mynvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k >= mynvtxs) {
        ++next;
        break;
      }
    }
    if (j == xadj[i+1]) {
      ++nint;
    }
  }

  if (r_nint) {
    *r_nint = nint;
  }
  if (r_next) {
    *r_next = next;
  }
}


wgt_t par_graph_cut(
    graph_t const * const graph,
    pid_t const * const * const where)
{
  vtx_t i, k, lvtx, nbrid;
  adj_t j;
  wgt_t cut;

  tid_t const myid = dlthread_get_id(graph->comm);

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const adjwgt = graph->adjwgt[myid];
  pid_t const * const mywhere = where[myid];

  DL_ASSERT_EQUALS((int)graph->dist.nthreads, \
      (int)dlthread_get_nthreads(graph->comm),"%d");

  cut = 0;

  if (graph->uniformadjwgt || graph->adjwgt == NULL) {
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
  cut = wgt_dlthread_sumreduce(cut,graph->comm);

  return cut/2;
}


void par_graph_removeislands(
    ctrl_t * const ctrl,
    graph_t * const graph)
{
  vtx_t i, k, nislands, nvtxs, lvtx, mynvtxs;
  adj_t j;
  pid_t p;
  tid_t nbrid;
  wgt_t iwgt;
  adj_t * xadj; 
  wgt_t * vwgt, * ivwgt;
  vtx_t * rename, * adjncy;
  vtx_t ** grename, ** glabel;

  tid_t const myid = dlthread_get_id(ctrl->comm);
  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);

  int const do_vwgt = !graph->uniformvwgt;

  mynvtxs = graph->mynvtxs[myid];
  xadj = graph->xadj[myid];
  adjncy = graph->adjncy[myid];
  vwgt = graph->vwgt[myid];

  /* see if there are islands */
  iwgt = 0;
  for (i=0;i<mynvtxs;++i) {
    if (xadj[i+1] == xadj[i]) {
      if (do_vwgt) {
        iwgt += vwgt[i]; 
      } else {
        ++iwgt;
      }
    }
  }

  /* implicit barrier */
  iwgt = wgt_dlthread_sumreduce(iwgt,ctrl->comm);

  if (iwgt < MIN_ISLAND_WEIGHT * graph->tvwgt) {
    /* not worth it */
    par_dprintf("Not removing islands: %0.03lf%%\n", \
        100.0*iwgt/(double)graph->tvwgt);
    return;
  }

  par_dprintf("Removing islands: %0.03lf%%\n", \
      100.0*iwgt/(double)graph->tvwgt);

  grename = dlthread_get_shmem(nthreads*sizeof(vtx_t*),ctrl->comm);
  grename[myid] = rename = vtx_alloc(mynvtxs);

  glabel = dlthread_get_shmem(nthreads*sizeof(vtx_t*),ctrl->comm);
  glabel[myid] = vtx_alloc(mynvtxs);

  if (do_vwgt) {
    ivwgt = wgt_alloc(mynvtxs);
  } else {
    ivwgt = NULL;
  }

  /* remove islands */
  iwgt = 0;
  nvtxs = 0;
  nislands = 0;
  for (i=0;i<mynvtxs;++i) {
    if (xadj[i+1] == xadj[i]) {
      ++nislands;
      iwgt += vwgt[i]; 
      /* starts at mynvtxs-1 */
      rename[i] = mynvtxs-nislands;
      if (graph->label) {
        glabel[myid][mynvtxs-nislands] = graph->label[myid][i];
      } else {
        glabel[myid][mynvtxs-nislands] = i;
      }
      if (do_vwgt) {
        ivwgt[mynvtxs-nislands] = vwgt[i];
      }
    } else {
      if (do_vwgt) {
        vwgt[nvtxs] = vwgt[i];
      }
      xadj[nvtxs+1] = xadj[i+1];
      rename[i] = nvtxs;
      if (graph->label) {
        glabel[myid][nvtxs] = graph->label[myid][i];
      } else {
        glabel[myid][nvtxs] = i;
      }
      ++nvtxs;
    }
  }

  dlthread_barrier(ctrl->comm);

  /* adjust edges */
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        adjncy[j] = rename[k];
      } else {
        lvtx = gvtx_to_lvtx(k,graph->dist);
        nbrid = gvtx_to_tid(k,graph->dist);
        adjncy[j] = lvtx_to_gvtx(grename[nbrid][lvtx],nbrid,graph->dist);  
      }
    }
  }

  /* insert islands at the end */
  if (do_vwgt) {
    wgt_copy(vwgt+nvtxs,ivwgt+nvtxs,nislands);
    dl_free(ivwgt);
  }

  if (graph->label) {
    dl_free(graph->label[myid]);
    dlthread_barrier(ctrl->comm);
  }

  /* adjust ctrl */
  if (myid == 0) {
    ctrl->ubfactor *= graph->tvwgt/(real_t)(graph->tvwgt - iwgt);

    graph->tvwgt -= iwgt;
    graph->invtvwgt = 1.0 / graph->tvwgt;

    if (ctrl->pijbm) {
      for (p=0;p<ctrl->nparts;++p) {
        ctrl->pijbm[p] = graph->invtvwgt / ctrl->tpwgts[p];
      }
    }

    if (graph->label) {
      dl_free(graph->label);
    }
    graph->label = glabel;
    graph->nislands = vtx_alloc(nthreads);
  }
  dlthread_barrier(ctrl->comm);

  /* adjust numbers */
  graph->mynvtxs[myid] = nvtxs;
  graph->nislands[myid] = nislands;

  /* implicit barrier */
  dl_free(rename);
  dlthread_free_shmem(grename,ctrl->comm);

  if (myid == 0) {
    graph->nvtxs = vtx_sum(graph->mynvtxs,nthreads);
  }

  par_dprintf("Removed %"PF_VTX_T" islands for new balance constraint of %" \
      PF_REAL_T"\n",nislands,ctrl->ubfactor);
}


void par_graph_restoreislands(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const * const gwhere)
{
  vtx_t i, mynvtxs;
  pid_t p, nparts;
  wgt_t iwgt, excess, twgt, uwgt;
  wgt_t * lpwgts;
  double * fpwgts;

  tid_t const myid = dlthread_get_id(ctrl->comm);
  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);

  wgt_t const * const vwgt = graph->vwgt[myid];

  int const do_vwgt = !graph->uniformvwgt;

  switch (ctrl->ptype) {
    case MTMETIS_PTYPE_ND:
    case MTMETIS_PTYPE_VSEP:
      nparts = 2;
      break;
    case MTMETIS_PTYPE_ESEP:
    case MTMETIS_PTYPE_RB:
    case MTMETIS_PTYPE_KWAY:
    default:
      nparts = ctrl->nparts;
  }

  mynvtxs = graph->mynvtxs[myid];
  lpwgts = wgt_init_alloc(0,nparts);

  /* restore */
  if (do_vwgt) {
    iwgt = wgt_sum(vwgt+mynvtxs,graph->nislands[myid]);
    twgt = wgt_dlthread_sumreduce(iwgt,ctrl->comm);
  } else {
    iwgt = graph->nislands[myid];
    twgt = vtx_sum(graph->nislands,nthreads);
  }
  graph->mynvtxs[myid] += graph->nislands[myid];

  if (myid == 0) {
    ctrl->ubfactor *= graph->tvwgt/(real_t)(graph->tvwgt + twgt);

    graph->tvwgt += twgt;
    graph->invtvwgt = 1.0 / graph->tvwgt;

    if (ctrl->pijbm) {
      for (p=0;p<ctrl->nparts;++p) {
        ctrl->pijbm[p] = graph->invtvwgt / ctrl->tpwgts[p];
      }
    }
  }

  /* calculate fillable pwgts */
  fpwgts = dlthread_get_shmem(nparts*sizeof(double),ctrl->comm);
  if (myid == 0) {
    excess = 0;
    for (p=0;p<nparts;++p) {
      uwgt = (ctrl->tpwgts[p]*graph->tvwgt)-graph->pwgts[p];
      if (uwgt < 0) {
        uwgt = 0;
      }
      fpwgts[p] = uwgt;
      excess += uwgt;
    }
    for (p=0;p<nparts;++p) {
      fpwgts[p] /= excess;
    }
  }
  dlthread_barrier(ctrl->comm);

  /* assign my island vertices to partitions */
  p = myid % nparts;
  for (i=0;i<graph->nislands[myid];++i) {
    while ((lpwgts[p]/(double)iwgt) >= fpwgts[p]) {
      p = (p+1)%nparts;
    }
    gwhere[myid][i+mynvtxs] = p;
    if (do_vwgt) {
      lpwgts[p] += vwgt[i+mynvtxs];
    } else {
      lpwgts[p] += 1;
    }
  }

  wgt_dlthread_sumareduce(lpwgts,nparts,ctrl->comm);
  
  if (myid == 0) {
    graph->nvtxs = vtx_sum(graph->mynvtxs,nthreads);
    for (p=0;p<nparts;++p) {
      graph->pwgts[p] += lpwgts[p];
    }
    dl_free(graph->nislands);
  }

  dl_free(lpwgts);

  /* free unused stuff */
  dlthread_free_shmem(fpwgts,ctrl->comm);
}


void par_graph_extract_parts(
    graph_t * const graph,
    pid_t const * const * const gwhere,
    pid_t const nparts,
    graph_t ** const parts)
{
  /* This funciton needs to handle cases where nparts < nthreads and 
   * nparts > nthreads. The latter will imply that it can be executed serially
   * without issue. We will assume that in hte case that nthreads > nparts,
   * nthreads is divisible by nparts. */
  vtx_t i, k, pnvtxs, deg, g, v, u;
  adj_t j, l;
  pid_t w;
  tid_t o, pmyid, pid, pnthreads;
  dlthread_comm_t hcomm;
  vtx_t ** vprefix;
  vtx_t ** vsuffix;
  vtx_t * nvtxs, * pmynvtxs;

  tid_t const myid = dlthread_get_id(graph->comm);
  tid_t const nthreads = dlthread_get_nthreads(graph->comm);

  /* create my vertex and edge counts */
  nvtxs = vtx_init_alloc(0,nparts);
  for (i=0;i<graph->mynvtxs[myid];++i) {
    w = gwhere[myid][i];
    if (w < nparts) { /* allow for the excluding of separators */
      ++nvtxs[w];
    }
  }

  /* allocate space for prefix and suffix arrays for all parts */
  vprefix = dlthread_get_shmem(sizeof(vtx_t*)*2*nparts,graph->comm);
  vsuffix = vprefix+nparts;

  for (pid=myid%nparts;pid<nparts;pid+=nthreads) { /* everythread */

    /* some make this evenly distribute threads based on partition size */
    pnthreads = (nthreads / nparts) + (nthreads % nparts > 0 ? 1 : 0);

    /* create communicators for my new graphs */
    if (pnthreads > 1) {
      /* each thread will only work on its partition */
      hcomm = dlthread_comm_split(pid,nparts,graph->comm);
    } else {
      /* each thread may work on more than one partition */
      hcomm = DLTHREAD_COMM_SINGLE;
    }

    /* all threads with their respective part(s) */
    parts[pid] = par_graph_create(hcomm);
  }

  for (pid=myid;pid<nparts;pid+=nthreads) { /* lead threads only */
    vprefix[pid] = vtx_alloc(nthreads+1);
    pnthreads = dlthread_get_nthreads(parts[pid]->comm);
    parts[pid]->label = r_vtx_alloc(pnthreads);
  }

  /* setup rename */
  if (myid == 0) {
    graph->rename = r_vtx_alloc(nthreads);
  }
  dlthread_barrier(graph->comm);
  graph->rename[myid] = vtx_alloc(graph->mynvtxs[myid]);

  /* assign prefixes */
  for (pid=0;pid<nparts;++pid) {
    vprefix[pid][myid] = nvtxs[pid]; 
  }

  dlthread_barrier(graph->comm);
  for (pid=myid;pid<nparts;pid+=nthreads) { /* lead threads only */
    pnthreads = dlthread_get_nthreads(parts[pid]->comm);

    /* create prefixsums of the vertices and edges */
    vprefix[pid][nthreads] = 0;
    vtx_prefixsum_exc(vprefix[pid],nthreads+1);

    pnvtxs = vprefix[pid][nthreads];

    /* create prefixsums for actual insertion into split graphs */
    vsuffix[pid] = vtx_alloc(pnthreads);
    for (o=0;o<pnthreads;++o) {
      vsuffix[pid][o] = vtx_chunkstart(o,pnthreads,pnvtxs);
      parts[pid]->mynvtxs[o] = vtx_chunksize(o,pnthreads,pnvtxs);
    }

    graph_calc_dist(vtx_max_value(parts[pid]->mynvtxs,pnthreads),pnthreads, \
        &(parts[pid]->dist));
  }
  dlthread_barrier(graph->comm);

  /* allocate vwgt and xadj */
  for (pid=myid%nparts;pid<nparts;pid+=nthreads) {
    pmyid = dlthread_get_id(parts[pid]->comm);
    pnvtxs = parts[pid]->mynvtxs[pmyid];
    parts[pid]->xadj[pmyid] = adj_alloc(pnvtxs+1);
    parts[pid]->vwgt[pmyid] = wgt_alloc(pnvtxs);
    parts[pid]->label[pmyid] = vtx_alloc(pnvtxs);
  }

  dlthread_barrier(graph->comm);
  /* insert vertex information into graphs */
  pmynvtxs = vtx_init_alloc(0,nparts);
  for (i=0;i<graph->mynvtxs[myid];++i) {
    w = gwhere[myid][i];
    if (w < nparts) {
      u = pmynvtxs[w]++;

      DL_ASSERT_EQUALS(parts[w]->dist.nthreads, \
          (tid_t)dlthread_get_nthreads(parts[w]->comm),"%"PF_TID_T);

      pnvtxs = vprefix[w][nthreads];
      
      /* get total vertex number */
      g = u + vprefix[w][myid];

      DL_ASSERT(g < pnvtxs,"Got vertex number of %"PF_VTX_T"/%"PF_VTX_T" " \
          "from %"PF_VTX_T" and %"PF_VTX_T"\n",g,pnvtxs,u, \
          vprefix[w][myid]);

      /* get new local vertex number */
      pmyid = vtx_chunkid(g,parts[w]->dist.nthreads,pnvtxs); 

      DL_ASSERT(pmyid < parts[w]->dist.nthreads,"Got chunk id of %"PF_TID_T \
          "/%"PF_TID_T" from %"PF_VTX_T", %"PF_TID_T", %"PF_VTX_T"\n",pmyid, \
          parts[w]->dist.nthreads,g,pnthreads,pnvtxs);

      v = g - vsuffix[w][pmyid];

      /* set rename */
      graph->rename[myid][i] = lvtx_to_gvtx(v,pmyid,parts[w]->dist);

      /* set alias as global vertex ID */
      parts[w]->label[pmyid][v] = lvtx_to_gvtx(i,myid,graph->dist);

      /* copy vertex weight */
      parts[w]->vwgt[pmyid][v] = graph->vwgt[myid][i];

      /* copy xadj info */
      deg = 0;
      for (j=graph->xadj[myid][i];j<graph->xadj[myid][i+1];++j) {
        k = graph->adjncy[myid][j];
        if (k < graph->mynvtxs[myid]) {
          o = myid;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        if (gwhere[o][k] == w) {
          ++deg;
        }
      }
      parts[w]->xadj[pmyid][v] = deg;
    }
  }

  dlthread_barrier(graph->comm);

  /* fix respective xadj's */
  for (pid=myid%nparts;pid<nparts;pid+=nthreads) { /* every thread */
    pmyid = dlthread_get_id(parts[pid]->comm);
    pnvtxs = parts[pid]->mynvtxs[pmyid];
    parts[pid]->xadj[pmyid][pnvtxs] = 0;
    adj_prefixsum_exc(parts[pid]->xadj[pmyid],pnvtxs+1);

    parts[pid]->adjncy[pmyid] = \
        vtx_alloc(parts[pid]->xadj[pmyid][pnvtxs]);
    parts[pid]->adjwgt[pmyid] = \
        wgt_alloc(parts[pid]->xadj[pmyid][pnvtxs]);

    parts[pid]->mynedges[pmyid] = parts[pid]->xadj[pmyid][pnvtxs];
  }

  dlthread_barrier(graph->comm);
  /* insert edge information into graphs */
  vtx_set(pmynvtxs,0,nparts);
  for (i=0;i<graph->mynvtxs[myid];++i) {
    w = gwhere[myid][i];
    if (w < nparts) {
      u = pmynvtxs[w]++;
      pnvtxs = vprefix[w][nthreads];
      
      /* get total vertex number */
      g = u + vprefix[w][myid];

      /* get new local vertex number */
      pmyid = vtx_chunkid(g,parts[w]->dist.nthreads,pnvtxs); 
      v = g - vsuffix[w][pmyid];

      l = parts[w]->xadj[pmyid][v];
      for (j=graph->xadj[myid][i];j<graph->xadj[myid][i+1];++j) {
        k = graph->adjncy[myid][j];
        if (k < graph->mynvtxs[myid]) {
          o = myid;
        } else {
          o = gvtx_to_tid(k,graph->dist);
          k = gvtx_to_lvtx(k,graph->dist);
        }
        if (gwhere[o][k] == w) {
          /* calcuate new endpoint */
          u = graph->rename[o][k];
          o = gvtx_to_tid(u,parts[w]->dist);
          if (o == pmyid) {
            k = gvtx_to_lvtx(u,parts[w]->dist);
          } else {
            k = u;
          }

          /* process edge */
          parts[w]->adjncy[pmyid][l] = k;
          parts[w]->adjwgt[pmyid][l] = graph->adjwgt[myid][j];
          ++l;
        }
      }
    }
  }
  dl_free(nvtxs);
  dl_free(pmynvtxs);
  dlthread_barrier(graph->comm);

  for (pid=myid;pid<nparts;pid+=nthreads) { /* lead threads only */
    parts[pid]->nvtxs = vtx_sum(parts[pid]->mynvtxs,parts[pid]->dist.nthreads);
    parts[pid]->gnvtxs = max_gvtx(parts[pid]); 
    parts[pid]->nedges = adj_sum(parts[pid]->mynedges, \
        parts[pid]->dist.nthreads);

    dl_free(vprefix[pid]);
    dl_free(vsuffix[pid]);
  }

  for (pid=myid%nparts;pid<nparts;pid+=nthreads) {
    par_graph_setup_twgts(parts[pid]);
  }

  /* implicit barrier */
  dlthread_free_shmem(vprefix,graph->comm); /* vsuffix is part of this */
}



#endif

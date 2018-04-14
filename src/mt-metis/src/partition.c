/**
 * @file partition.c
 * @brief Partitioning functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */




#ifndef MTMETIS_PARTITION_C
#define MTMETIS_PARTITION_C




#include "partition.h"
#include "coarsen.h"
#include "initpart.h"
#include "uncoarsen.h"




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Calcuate teh edgecut of a partitioning, serially.
 *
 * @param graph The graph.
 * @param where The partition ids.
 *
 * @return  The total weight of cut edges.
 */
static wgt_t __partition_calc_cut(
    graph_t const * const graph,
    pid_t const * const * const where)
{
  vtx_t i,k,l,mynvtxs;
  adj_t j;
  wgt_t cut;
  pid_t me, other;
  tid_t myid, o;

  tid_t const nthreads = graph->dist.nthreads;

  cut = 0;
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    for (i=0;i<mynvtxs;++i) {
      me = where[myid][i];
      for (j=graph->xadj[myid][i];j<graph->xadj[myid][i+1];++j) {
        k = graph->adjncy[myid][j];
        if (k > mynvtxs) {
          l = gvtx_to_lvtx(k,graph->dist);
          o = gvtx_to_tid(k,graph->dist);
        } else {
          l = k;
          o = myid;
        }
        other = where[o][l];
        if (other != me) {
          cut += graph->adjwgt[myid][j];
        }
      }
    }
  }

  return cut/2;
}


/**
 * @brief Calculate the communication volume of partitioning (serially).
 *
 * @param graph The partitioned graph.
 * @param where The partition IDs.
 * @param nparts The number of partitions.
 *
 * @return The communication volume. 
 */
static vtx_t __partition_calc_comvol(
    graph_t const * const graph,
    pid_t const * const * const where,
    pid_t const nparts)
{
  vtx_t vol, i, k, l, mynvtxs, g;
  adj_t j;
  pid_t me, other;
  tid_t o, myid;
  vtx_t * marker;

  tid_t const nthreads = graph->dist.nthreads;

  marker = vtx_init_alloc(NULL_VTX,nparts);

  vol = 0;
  for (myid=0;myid<nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    for (i=0;i<mynvtxs;++i) {
      me = where[myid][i];
      g = lvtx_to_gvtx(i,myid,graph->dist);
      marker[me] = g; 
      for (j=graph->xadj[myid][i];j<graph->xadj[myid][i+1];++j) {
        k = graph->adjncy[myid][j];
        if (k > mynvtxs) {
          l = gvtx_to_lvtx(k,graph->dist);
          o = gvtx_to_tid(k,graph->dist);
        } else {
          l = k;
          o = myid;
        }
        other = where[o][l];
        if (marker[other] != g) {
          marker[other] = g;
          ++vol;
        }
      }
    }
  }

  dl_free(marker);

  return vol;
}


/**
 * @brief Perform a multilevel kway partitioning in parallel (to be called by
 * each thread in a parallel region).
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param where The allocated where vector.
 *
 * @return Total weight of cut edges.
 */
static wgt_t __partition_mlevel_kway(
    ctrl_t * const ctrl,
    graph_t * const graph,
    pid_t * const * const where)
{
  vtx_t i;
  wgt_t curobj, bestobj;
  double curbal, bestbal;
  graph_t * cgraph;

  tid_t const myid = omp_get_thread_num();

  curobj = 0;
  bestobj = 0;
  curbal=0.0;
  bestbal=0.0;

  for (i=0;i<ctrl->ncuts;++i) {
    cgraph = coarsen_graph(ctrl,graph);

    initpart_kway(ctrl,cgraph);

    DL_ASSERT_EQUALS(graph_cut(ctrl,cgraph, \
        (const pid_t **)cgraph->where),cgraph->mincut,"%"PF_WGT_T); 

    uncoarsen_kway(ctrl, graph, cgraph);

    curobj = graph->mincut;

    curbal = graph_imbalance(graph,ctrl->nparts,ctrl->pijbm);

    if (i == 0  \
        || (curbal <= 0.0005 && bestobj > curobj) \
        || (bestbal > 0.0005 && curbal < bestbal)) {
      pid_copy(where[myid],graph->where[myid],graph->mynvtxs[myid]);
      bestobj = curobj;
      bestbal = curbal;
    }

    graph_free_rdata(graph);

    if (bestobj == 0)
      break;
  }

  return bestobj;
}





/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void partition_kway(
    ctrl_t * const ctrl, 
    graph_t * const graph, 
    pid_t ** const where)
{
  pid_t i;

  if (!ctrl->pijbm) {
    #pragma omp master
    {
      ctrl->pijbm = real_alloc(ctrl->nparts);
    }
  }
  #pragma omp barrier

  /* set up multipliers for making balance computations easier */
  #pragma omp for 
  for (i=0;i<ctrl->nparts;++i) {
    ctrl->pijbm[i] = graph->invtvwgt / ctrl->tpwgts[i];
  }

  graph->mincut = __partition_mlevel_kway(ctrl,graph,where); 
}


void partition_print_info(
    ctrl_t const * ctrl,
    graph_t const * graph,
    pid_t const * const * where)
{
  vtx_t i, k, mynvtxs;
  pid_t nparts;
  tid_t myid;
  wgt_t tvwgt, mvwgt;
  wgt_t * kpwgts;
  real_t * tpwgts;
  double unbalance;
  const wgt_t *vwgt; 
  const pid_t *mywhere;

  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  dl_print_footer('*');
  printf(" size of vtx_t: %zu, adj_t: %zu, wgt_t: %zu, pid_t: %zu, tid_t: " \
      "%zu, real_t: %zu\n",8*sizeof(vtx_t), 8*sizeof(adj_t), 8*sizeof(wgt_t), \
      8*sizeof(pid_t), 8*sizeof(tid_t), 8*sizeof(real_t));
  printf("\n");
  dl_print_header("Graph Information",'-');
  printf("#Vertices: %"PF_VTX_T", #Edges: %"PF_ADJ_T", #Parts: %"PF_PID_T"\n", 
    graph->nvtxs, graph->nedges/2, ctrl->nparts);

  printf("\n");
  printf("\n");
  dl_print_header("Direct k-way Partitioning",'-');


  /* Compute objective-related infomration */
  printf(" - Edgecut: %"PF_ADJ_T", communication volume: %"PF_VTX_T".\n\n", \
    __partition_calc_cut(graph,where),__partition_calc_comvol(graph,where, \
      nparts));


  /* Compute constraint-related information */
  kpwgts = wgt_init_alloc(0,nparts);

  for (myid=0;myid<graph->dist.nthreads;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    vwgt = graph->vwgt[myid];
    mywhere = where[myid];
    for (i=0; i<mynvtxs; ++i) {
      kpwgts[mywhere[i]] += vwgt[i];
    }
  }

  /* Report on balance */
  printf(" - Balance:\n");
  tvwgt = wgt_sum(kpwgts,nparts);
  k = 0;
  unbalance = 1.0*kpwgts[k]/(tpwgts[k]*tvwgt);
  for (i=1;i<nparts;++i) {
    if (unbalance < 1.0*kpwgts[i]/(tpwgts[i]*tvwgt)) {
      unbalance = 1.0*kpwgts[i]/(tpwgts[i]*tvwgt);
      k = i;
    }
  }
  mvwgt = 0;
  for (myid=0;myid<graph->dist.nthreads;++myid) {
    if ((mynvtxs = graph->mynvtxs[myid]) > 0) {
      vwgt = graph->vwgt[myid];
      k = wgt_max_index(vwgt,mynvtxs);
      if (vwgt[k] > mvwgt) {
        mvwgt = vwgt[k];
      }
    }
  }

  printf("     constraint #0:  %5.3lf out of %5.3lf\n", \
      unbalance, 1.0*nparts*mvwgt/ (1.0*tvwgt));
  printf("\n");

  for (k=0, unbalance=kpwgts[k]/(tpwgts[k]*tvwgt), i=1; i<nparts; i++) {
    if (unbalance < kpwgts[i]/(tpwgts[i]*tvwgt)) {
      unbalance = kpwgts[i]/(tpwgts[i]*tvwgt);
      k = i;
    }
  }

  printf(" - Most overweight partition:\n");
  printf("     pid: %"PF_PID_T", actual: %"PF_WGT_T", desired: %"PF_WGT_T \
         ", ratio: %.2lf\n",k,kpwgts[k],(wgt_t)(tvwgt*tpwgts[k]),unbalance);
  printf("\n");

  dl_free(kpwgts);


  #ifdef XXX
  /* Compute subdomain adjacency information */
  pptr = imalloc(nparts+1, "ComputePartitionInfo: pptr");
  pind = imalloc(nvtxs, "ComputePartitionInfo: pind");
  pdom = imalloc(nparts, "ComputePartitionInfo: pdom");

  iarray2csr(nvtxs, nparts, (idx_t*)where, pptr, pind);

  maxndom = nparts+1;
  minndom = 0;
  for (tndom=0, pid=0; pid<nparts; pid++) {
    iset(nparts, 0, pdom);
    for (ii=pptr[pid]; ii<pptr[pid+1]; ii++) {
      i = pind[ii];
      for (j=xadj[i].start; j<xadj[i].end; j++)
        pdom[where[adjncy[j]]] += adjwgt[j];
    }
    pdom[pid] = 0;
    for (ndom=0, i=0; i<nparts; i++)
      ndom += (pdom[i] > 0 ? 1 : 0);
    tndom += ndom;
    if (pid == 0 || maxndom < ndom)
      maxndom = ndom;
    if (pid == 0 || minndom > ndom)
      minndom = ndom;
  }

  printf(" - Subdomain connectivity: max: %"PRIDX", min: %"PRIDX \
      ", avg: %.2"PRREAL"\n\n",maxndom, minndom, 1.0*tndom/nparts);
      
  gk_free((void **)&pptr, &pind, &pdom, LTERM);
  #endif


  dl_print_footer('*');
}            




#endif


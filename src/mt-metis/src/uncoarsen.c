/**
 * @file uncoarsen.c
 * @brief Uncoarsening functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */





#ifndef MTMETIS_UNCOARSEN_C
#define MTMETIS_UNCOARSEN_C




#include "uncoarsen.h"
#include "ucinfo.h"
#include "refine.h"
#include "check.h"




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static wgt_t ** __upk_gpwgts;
/**
 * @brief Compute the parameters for performing kway partitioning.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * 
 * @return The newly setup ucinfo.
 */
static ucinfo_t * __uncoarsen_partparams_kway(
    ctrl_t * const ctrl,
    graph_t * const graph)
{
  vtx_t other,me,i,k,lvtx,nbrid,na;
  adj_t j, l;
  wgt_t mincut;
  nbrinfo_t * nbrinfo; 
  nbrinfo_t * myrinfo;
  adjinfo_t * mynbrs;
  vtx_iset_t * bnd;
  ucinfo_t * ucinfo;
  wgt_t * mypwgts;

  pid_t const nparts = ctrl->nparts;

  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;
  wgt_t const * const * const gvwgt = (wgt_t const **)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const **)graph->adjwgt;

  pid_t const * const * const gwhere  = (pid_t const **)graph->where;

  tid_t const myid = omp_get_thread_num();
  tid_t const nthreads = omp_get_num_threads();

  wgt_t * const pwgts = graph->pwgts;

  vtx_t const mynvtxs = graph->mynvtxs[myid];

  adj_t const * const xadj = gxadj[myid];
  vtx_t const * const adjncy = gadjncy[myid];
  wgt_t const * const vwgt = gvwgt[myid];
  wgt_t const * const adjwgt = gadjwgt[myid];
  pid_t const * const where = gwhere[myid];

  DL_ASSERT(graph->pwgts != NULL,"Non-allocated pwgts");
  DL_ASSERT(graph->where != NULL,"Non-allocated where");

  #pragma omp master
  {
    __upk_gpwgts = r_wgt_alloc(nthreads);
  }
  #pragma omp barrier
  mypwgts = __upk_gpwgts[myid] = wgt_init_alloc(0,nparts);

  /* reset the size of neighbor infor array */
  ucinfo = ucinfo_create(ctrl,graph);
  nbrinfo = ucinfo->nbrinfo;
  ucinfo->nnbrpool = 0;
  bnd = ucinfo->bnd;

  mincut = 0;
  vtx_iset_clear(bnd);

  /* calculate partition weights */
  for (i=0; i<mynvtxs; ++i) {
    mypwgts[where[i]] += vwgt[i];
  }
  #pragma omp barrier
  #pragma omp for schedule(static) 
  for (i=0; i<nparts;++i) {
    k = 0;
    for (j=0; j<nthreads;++j) {
      k += __upk_gpwgts[j][i];
    }
    pwgts[i] = k;
  }

  /* clear neighbor info */ 
  memset(nbrinfo,0,sizeof(nbrinfo_t)*mynvtxs);

  /* calculate nbrinfo for vertices */
  for (i=0; i<mynvtxs; ++i) {
    myrinfo = nbrinfo+i;

    na = dl_min(nparts,xadj[i+1]-xadj[i]);
    myrinfo->nbrstart = ucinfo->nnbrpool;
    /* check to see if we need to expand the pool */
    if ((ucinfo->nnbrpool += na) > ucinfo->maxnnbrpool) { 
      ucinfo->maxnnbrpool *= NBRPOOL_EXP_RATE;
      ucinfo->nbrpool = adjinfo_realloc(ucinfo->nbrpool,ucinfo->maxnnbrpool);
    }

    mynbrs = ucinfo->nbrpool + myrinfo->nbrstart;
    me = where[i];

    for (j=xadj[i]; j<xadj[i+1]; ++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        lvtx = k;
        nbrid = myid;
      } else {
        nbrid = gvtx_to_tid(k,graph->dist);
        lvtx = gvtx_to_lvtx(k,graph->dist); 
      }
      if (me == gwhere[nbrid][lvtx]) {
        myrinfo->id += adjwgt[j];
      } else {
        myrinfo->ed += adjwgt[j];
      }
    }

    /* Time to compute the particular external degrees */
    if (myrinfo->ed > 0) {
      mincut += myrinfo->ed;
      for (j=xadj[i]; j<xadj[i+1]; ++j) {
        k = adjncy[j];
        if (k < mynvtxs) {
          lvtx = k;
          nbrid = myid;
        } else {
          nbrid = gvtx_to_tid(k,graph->dist);
          lvtx = gvtx_to_lvtx(k,graph->dist); 
        }
        other = gwhere[nbrid][lvtx];
        if (me != other) {
          for (l=0; l<myrinfo->nnbrs; l++) {
            if (mynbrs[l].pid == other) {
              mynbrs[l].ed += adjwgt[j];
              break;
            }
          }
          if (l == myrinfo->nnbrs) {
            mynbrs[l].pid = other;
            mynbrs[l].ed  = adjwgt[j];
            myrinfo->nnbrs++;
          }
        }
      }

      /* Only ed-id>=0 nodes are considered to be in the boundary */
      if (myrinfo->ed >= myrinfo->id) {
        vtx_iset_add(i,bnd);
      }
      DL_ASSERT(myrinfo->nnbrs > 0,"No neighbors.");
    } else if (myrinfo->id == 0) {
      vtx_iset_add(i,bnd);
    } else {
      myrinfo->nbrstart = NULL_ADJ;
      ucinfo->nnbrpool -= na;
      DL_ASSERT_EQUALS(myrinfo->nnbrs,0,"%"PF_ADJ_T);
    }
  } 

  mincut = wgt_omp_sumreduce(mincut);

  dl_free(__upk_gpwgts[myid]);

  #pragma omp barrier
  #pragma omp master
  {
    graph->mincut = mincut/2;

    dl_free(__upk_gpwgts);

  }

  /* the checks */
  DL_ASSERT(check_info(ucinfo,graph,(pid_t const **)gwhere),"Bad info");
  DL_ASSERT(check_bnd(ucinfo->bnd,graph),"Bad boundary");

  return ucinfo;
}


static ucinfo_t ** __upk_ucinfos;
/**
 * @brief Project a kway partitioning.
 *
 * @param ctrl The control structure containing runtime parameters.
 * @param graph The partitioned graph to project the partition from.
 */
static void __uncoarsen_project_kway(
    ctrl_t * const ctrl,
    graph_t * const graph,
    ucinfo_t * const ucinfo)
{
  vtx_t i, k, pi, lvtx, nbrid, ned, nid;
  adj_t j, l, istart, iend;
  pid_t me, other, na;
  wgt_t tid, ted;
  nbrinfo_t * myrinfo;
  adjinfo_t * mynbrs;
  vtx_t * id, * ed;
  pid_t * htable;
  vtx_iset_t * bnd;
  nbrinfo_t * nbrinfo;

  tid_t const myid = omp_get_thread_num();
  tid_t const nthreads = omp_get_num_threads();

  pid_t const nparts = ctrl->nparts;
  graph_t * const cgraph = graph->coarser;
  pid_t const * const * const gcwhere = (pid_t const **)cgraph->where;
  vtx_t * const * const gcmap = graph->cmap;
  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;
  wgt_t const * const * const gadjwgt = (wgt_t const **)graph->adjwgt;

  #pragma omp master
  {
    dl_start_timer(&(ctrl->timers.projection));
  }

  graph_alloc_partmemory(ctrl,graph);

  graph->mincut = cgraph->mincut;

  pid_t ** const gwhere = graph->where;
  pid_t * const where = gwhere[myid];

  vtx_t const mynvtxs = graph->mynvtxs[myid];

  vtx_t * const cmap = gcmap[myid];

  adj_t const * const xadj = gxadj[myid];
  vtx_t const * const adjncy = gadjncy[myid];
  wgt_t const * const adjwgt = gadjwgt[myid];

  #pragma omp master
  {
    __upk_ucinfos = (ucinfo_t**)malloc(sizeof(ucinfo_t*)*nthreads);
  }
  #pragma omp barrier
  __upk_ucinfos[myid] = ucinfo;
  #pragma omp barrier
  
  /* Compute the required info for refinement */
  ucinfo->nnbrpool = 0;

  ed = vtx_alloc(mynvtxs);
  id = vtx_alloc(mynvtxs);

  ned = nid = 0;
  for (i=0;i<mynvtxs;++i) {
    k = cmap[i];
    if (k < cgraph->mynvtxs[myid]) {
      lvtx = k;
      nbrid = myid;
    } else {
      lvtx = gvtx_to_lvtx(k,graph->dist);
      nbrid = gvtx_to_tid(k,graph->dist);
    }
    where[i] = gcwhere[nbrid][lvtx];
    if (__upk_ucinfos[nbrid]->nbrinfo[lvtx].ed > 0) {
      ed[ned++] = i;
    } else {
      id[nid++] = i;
    }
  }

  htable = pid_init_alloc(NULL_PID,nparts);


  #pragma omp master
  {
    wgt_copy(graph->pwgts,cgraph->pwgts,nparts);
  }
  #pragma omp barrier

  graph_free(graph->coarser);


  /* expanpd nbrinfo */
  dl_free(ucinfo->nbrinfo);
  nbrinfo = ucinfo->nbrinfo = nbrinfo_calloc(mynvtxs);

  /* expand boundary */
  vtx_iset_free(ucinfo->bnd);
  bnd = ucinfo->bnd = vtx_iset_create(0,mynvtxs);

  for (pi=0;pi<nid;++pi) {
    i = id[pi];
    istart = xadj[i];
    iend = xadj[i+1];

    myrinfo = nbrinfo+i;

    for (tid=0, j=istart; j<iend; j++) {
      tid += adjwgt[j];
    }
    myrinfo->id = tid;
    myrinfo->nbrstart = NULL_ADJ;
    if (tid == 0) {
      /* keep islands on the border */
      vtx_iset_add(i,bnd);
    }
    DL_ASSERT_EQUALS(myrinfo->ed,0,"%"PF_WGT_T);
    DL_ASSERT_EQUALS(myrinfo->nnbrs,0,"%"PF_PID_T);
  }
  for (pi=0;pi<ned;++pi) {
    i = ed[pi];

    istart = xadj[i];
    iend = xadj[i+1];

    myrinfo = nbrinfo+i;

    na = dl_min(nparts,xadj[i+1]-xadj[i]);
    myrinfo->nbrstart = ucinfo->nnbrpool;
    /* check to see if we need to expand the pool */
    if ((ucinfo->nnbrpool += na) > ucinfo->maxnnbrpool) { 
      ucinfo->maxnnbrpool *= NBRPOOL_EXP_RATE;
      ucinfo->nbrpool = adjinfo_realloc(ucinfo->nbrpool,ucinfo->maxnnbrpool);
    }

    mynbrs = ucinfo->nbrpool + myrinfo->nbrstart;

    me = where[i];
    tid = 0;
    ted = 0;
    for (j=istart; j<iend; ++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        lvtx = k;
        nbrid = myid;
      } else {
        lvtx = gvtx_to_lvtx(k,graph->dist);
        nbrid = gvtx_to_tid(k,graph->dist);
      }
      other = gwhere[nbrid][lvtx];
      if (me == other) {
        tid += adjwgt[j];
      } else {
        ted += adjwgt[j];
        if ((l = htable[other]) == NULL_PID) {
          htable[other] = myrinfo->nnbrs;
          mynbrs[myrinfo->nnbrs].pid = other;
          mynbrs[myrinfo->nnbrs++].ed = adjwgt[j];
        } else {
          mynbrs[l].ed += adjwgt[j];
        }
      }
    }
    myrinfo->id = tid;
    myrinfo->ed = ted;

    if (ted > 0) {
      if (ted >= tid) {
        vtx_iset_add(i,bnd);
      }
      for (j=0; j<myrinfo->nnbrs; ++j) {
        htable[mynbrs[j].pid] = NULL_ADJ;
      }
    } else if (tid == 0) {
      vtx_iset_add(i,bnd);
    }
    if (myrinfo->nnbrs == 0) {
      ucinfo->nnbrpool -= na;
      myrinfo->nbrstart = NULL_ADJ;
    }
  }

  dl_free(htable);
  dl_free(ed);
  dl_free(id);


  #pragma omp master
  {
    dl_free(__upk_ucinfos);
    graph->coarser = NULL;
    DL_ASSERT(check_info(ucinfo,graph,(pid_t const **)gwhere),"Bad info");
    DL_ASSERT(check_bnd(ucinfo->bnd,graph),"Bad boundary");
  }
  #pragma omp barrier
  #pragma omp master
  {
    dl_stop_timer(&(ctrl->timers.projection));
  }
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void uncoarsen_kway(
    ctrl_t * const ctrl,
    graph_t * const orggraph,
    graph_t * graph)
{
  vtx_t i;
  ucinfo_t * ucinfo;

  #pragma omp master
  {
    dl_start_timer(&ctrl->timers.uncoarsening);
  }

  /* Compute the parameters of the coarsest graph */
  ucinfo = __uncoarsen_partparams_kway(ctrl,graph);
  
  /* Refine each successively finer graph */
  for (i=0; ;i++) {
    (void)refine_kway(ctrl,graph,ctrl->nrefpass,5.0,ucinfo); 

    par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH,"Final partition on " \
        "graph %zu: %"PF_WGT_T" cut and %5.4lf balance\n",graph->level, \
        graph->mincut,graph_imbalance(graph,ctrl->nparts,ctrl->pijbm));


    /* if we finished uncoarsening, stop */
    if (graph == orggraph)
      break;

    graph = graph->finer;

    __uncoarsen_project_kway(ctrl,graph,ucinfo);
  }

  ucinfo_free(ucinfo);

  #pragma omp master
  {
    dl_stop_timer(&ctrl->timers.uncoarsening);
  }
}







#endif

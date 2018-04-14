/**
 * @file kwayrefine.c
 * @brief KWay refinement routines
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2015, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */





#ifndef MTMETIS_KWAYREFINE_C
#define MTMETIS_KWAYREFINE_C




#include "kwayrefine.h"
#include "check.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct move_t {
  pid_t to;
  vtx_t vtx;
} move_t;


typedef struct update_t {
  pid_t to;
  pid_t from;
  wgt_t ewgt;
  vtx_t nbr;
} update_t;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLPQ_PREFIX vw
#define DLPQ_KEY_T wgt_t
#define DLPQ_VAL_T vtx_t
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_VAL_T
#undef DLPQ_KEY_T
#undef DLPQ_PREFIX


#define DLLCB_PREFIX move
#define DLLCB_TYPE_T move_t
#define DLLCB_STATIC 1
#include "dllcb_headers.h"
#undef DLLCB_STATIC
#undef DLLCB_TYPE_T
#undef DLLCB_PREFIX


#define DLLCB_PREFIX update
#define DLLCB_TYPE_T update_t
#define DLLCB_STATIC 1
#include "dllcb_headers.h"
#undef DLLCB_STATIC
#undef DLLCB_TYPE_T
#undef DLLCB_PREFIX


#define DLPQ_PREFIX vt
#define DLPQ_KEY_T wgt_t
#define DLPQ_VAL_T vtx_t
#define DLPQ_USE_HT
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_USE_HT
#undef DLPQ_VAL_T
#undef DLPQ_KEY_T
#undef DLPQ_PREFIX


#define DLMSET_PREFIX vtx
#define DLMSET_TYPE_T vtx_t
#define DLMSET_STATIC
#include "dlmset_headers.h"
#undef DLMSET_STATIC
#undef DLMSET_TYPE_T
#undef DLMSET_PREFIX




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const MIN_HILL_SIZE = 3;
static vtx_t const MAX_HILLS = 64;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static inline pid_t __partner(
    pid_t const side,
    pid_t const offset,
    pid_t const nparts,
    pid_t const d)
{
  pid_t other;

  pid_t const cycle = offset*2;
  pid_t const block = ((nparts/cycle)+(nparts%cycle ? 1 : 0))*cycle;

  /* figure out our side */
  if ((side / offset) % 2 == d) {
    other = (side + offset) % block;
  } else {
    other = ((block + side) - offset) % block;
  }

  return other;
}


static inline int __right_side(
    int const dir,
    pid_t const to,
    pid_t const from)
{
  if (dir) {
    return to < from;
  } else {
    return to > from;
  }
}


/**
 * @brief Update a vertex incident the one being moved.
 *
 * @param ctrl The control structure.
 * @param k The vertex to update.
 * @param to The partition the vertex is being moved to.
 * @param from The partition the vertex is being moved from.
 * @param ewgt The edge connecting the moved vertex.
 * @param graph The graph.
 * @param bnd The boundary set of vertices.
 * @param queue The priority queue of vertices to move.
 *
 * @return The change in edgecut.
 */
static wgt_t __update_vertex(
    ctrl_t * ctrl, 
    vtx_t const k, 
    pid_t const to, 
    pid_t const from, 
    wgt_t const ewgt, 
    graph_t * graph, 
    kwinfo_t * const kwinfo,
    vw_pq_t * queue)
{
  vtx_t l;
  wgt_t oed;
  real_t rgain;
  vtx_iset_t * bnd;
  kwnbrinfo_t * myrinfo;
  adjinfo_t * mynbrs;

  tid_t const myid = dlthread_get_id(ctrl->comm);

  pid_t * const where = graph->where[myid];
  pid_t const nparts = ctrl->nparts;
  pid_t const me = where[k];

  int const greedy = ctrl->rtype == MTMETIS_RTYPE_GREEDY;

  bnd = kwinfo->bnd;

  /* create my workspace */
  myrinfo = kwinfo->nbrinfo+k;

  oed = myrinfo->ed;
  
  mynbrs = kwinfo_get_nbrs(kwinfo,k, \
      dl_min(nparts,graph->xadj[myid][k+1]-graph->xadj[myid][k]));

  if (me == to) {
    myrinfo->id += ewgt;
    myrinfo->ed -= ewgt;
  } else if (me == from) {
    myrinfo->id -= ewgt;
    myrinfo->ed += ewgt;
  }
  /* add it to the boundary if necessary */
  if (!vtx_iset_contains(k,bnd)) {
    if (is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
      vtx_iset_add(k,bnd);
    }
  } else if (!is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
    vtx_iset_remove(k,bnd);
  }

  /* update nbrs */
  for (l=0;l<myrinfo->nnbrs;++l) {
    if (mynbrs[l].pid == from) {
      if (mynbrs[l].ed == ewgt) {
        mynbrs[l] = mynbrs[--myrinfo->nnbrs];
      } else {
        mynbrs[l].ed -= ewgt;
      }
      break;
    }
  }
  for (l=0;l<myrinfo->nnbrs;++l) {
    if (mynbrs[l].pid == to) {
      mynbrs[l].ed += ewgt;
      break;
    }
  }
  if (to != me && l == myrinfo->nnbrs) {
    mynbrs[myrinfo->nnbrs].ed = ewgt;
    mynbrs[myrinfo->nnbrs++].pid = to;
  }

  if (queue) {
    rgain = ((myrinfo->nnbrs > 0 ? \
          ((real_t)myrinfo->ed) / sqrt(myrinfo->nnbrs) : 0.0) \
          - myrinfo->id);

    if ((me == to || me == from)) {
      if (vw_pq_contains(k,queue)) {
        if ((!greedy && myrinfo->ed > 0) || rgain >= 0) {
          vw_pq_update(rgain,k,queue);
        } else {
          vw_pq_remove(k,queue);
        }
      }
    }
  }

  return oed - myrinfo->ed;
}


static wgt_t __move_vertex(
    ctrl_t * const ctrl,
    graph_t * const graph,
    tid_t const myid,
    vtx_t const i,
    pid_t const to,
    kwinfo_t * const kwinfo,
    wgt_t * const pwgts,
    pid_t * const where,
    vw_pq_t * const q,
    update_combuffer_t * const combuffer)
{
  vtx_t k;
  adj_t j;
  wgt_t cut, ted, ewgt;
  tid_t nbrid;
  update_t up;
  adjinfo_t * mynbrs;

  pid_t const nparts = ctrl->nparts;
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = graph->xadj[myid];
  vtx_t const * const adjncy = graph->adjncy[myid];
  wgt_t const * const adjwgt = graph->adjwgt[myid];
  wgt_t const * const vwgt = graph->vwgt[myid];

  vtx_iset_t * const bnd = kwinfo->bnd;

  kwnbrinfo_t * const myrinfo = kwinfo->nbrinfo+i;
  pid_t const from = where[i];

  int const greedy = ctrl->rtype == MTMETIS_RTYPE_GREEDY;

  cut = 0;

  pwgts[to] += vwgt[i];
  pwgts[from] -= vwgt[i];
  where[i] = to;

  ted = myrinfo->ed;

  mynbrs = kwinfo_get_nbrs(kwinfo,i, \
      dl_min(nparts,graph->xadj[myid][i+1]-graph->xadj[myid][i]));

  for (k=0;k<myrinfo->nnbrs;++k) {
    if (mynbrs[k].pid == to) {
      break;
    }
  }
  if (k==myrinfo->nnbrs) {
    k = NULL_PID;
  }
  
  /* make the move */
  if (k != NULL_PID) {
    myrinfo->ed += myrinfo->id-mynbrs[k].ed;
    dl_swap(myrinfo->id,mynbrs[k].ed);
  } else if (myrinfo->id > 0) {
    myrinfo->ed += myrinfo->id;
    mynbrs[myrinfo->nnbrs].ed = myrinfo->id;
    k = myrinfo->nnbrs++;
    myrinfo->id = 0;
  }

  /* old minus new */
  cut += ted - myrinfo->ed;

  if (mynbrs[k].ed == 0) {
    mynbrs[k] = mynbrs[--myrinfo->nnbrs];
  } else {
    mynbrs[k].pid = from;
  }
 
  /* see if this vertex should be removed/added from the boundary */
  if (vtx_iset_contains(i,bnd)) {
    if (!is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
      vtx_iset_remove(i,bnd);
    }
  } else {
    if (is_bnd(myrinfo->id,myrinfo->ed,greedy)) {
      vtx_iset_add(i,bnd);
    }
  }

  /* update neighbors */
  for(j=xadj[i];j<xadj[i+1];++j) {
    k = adjncy[j];
    ewgt = adjwgt[j];
    if (k < mynvtxs) {
      /* I own it */
      cut += __update_vertex(ctrl,k,to,from,ewgt,graph,kwinfo,q);
    } else {
      /* notify my neighbor */
      nbrid = gvtx_to_tid(k,graph->dist);

      up.to = to;
      up.from = from;
      up.ewgt = ewgt;
      up.nbr = gvtx_to_lvtx(k,graph->dist);

      update_combuffer_add(nbrid,up,combuffer);
    }
  }

  return cut;
}


static inline void __par_sync_pwgts(
    tid_t const myid,
    pid_t const nparts,
    wgt_t * const gpwgts,
    wgt_t * const lpwgts,
    dlthread_comm_t const comm)
{
  pid_t p;

  /* turn local pwgts into deltas */
  for (p=0;p<nparts;++p) {
    lpwgts[p] -= gpwgts[p];
  }

  /* create global deltas */
  wgt_dlthread_sumareduce(lpwgts,nparts,comm);

  /* set local pwgts to be global pwgts */
  for (p=0;p<nparts;++p) {
    lpwgts[p] += gpwgts[p];
  }

  dlthread_barrier(comm);

  /* re-sync global pwgts */
  if (myid == 0) {
    for (p=0;p<nparts;++p) {
      gpwgts[p] = lpwgts[p];
    }
  }

  dlthread_barrier(comm);
}




/******************************************************************************
* REFINEMENT FUNCTIONS ********************************************************
******************************************************************************/


static vtx_t __par_kwayrefine_GREEDY(
    ctrl_t * const ctrl, 
    graph_t * const graph,
    size_t const niter, 
    kwinfo_t * const kwinfo)
{
  vtx_t c, i, k, nmoved;
  adj_t j;
  wgt_t gain, wgt, mycut, ewgt;
  pid_t to, from;
  size_t pass;
  real_t rgain;
  wgt_t * lpwgts;
  kwnbrinfo_t * myrinfo;
  adjinfo_t const * mynbrs;
  vtx_iset_t * bnd;
  vw_pq_t * q;
  wgt_t * minwgt, * maxwgt;
  update_t up;
  update_combuffer_t * combuffer;

  tid_t const myid = dlthread_get_id(ctrl->comm);

  /* Link the graph fields */
  vtx_t const mynvtxs = graph->mynvtxs[myid];
  wgt_t const * const vwgt = graph->vwgt[myid];

  pid_t ** const gwhere = graph->where;
  wgt_t * const pwgts = graph->pwgts;
  
  pid_t const nparts = ctrl->nparts;

  kwnbrinfo_t * const nbrinfo = kwinfo->nbrinfo;
  pid_t * const where = gwhere[myid];
  real_t const * const tpwgts = ctrl->tpwgts;

  combuffer = update_combuffer_create(graph->mynedges[myid],ctrl->comm);

  lpwgts = wgt_alloc(nparts);
  wgt_copy(lpwgts,pwgts,nparts);

  minwgt = wgt_alloc(nparts);
  maxwgt = wgt_alloc(nparts);


  bnd = kwinfo->bnd;

  /* setup max/min partition weights */
  for (i=0;i<nparts;++i) {
    maxwgt[i]  = ctrl->tpwgts[i]*graph->tvwgt*ctrl->ubfactor;
    minwgt[i]  = ctrl->tpwgts[i]*graph->tvwgt*(1.0/ctrl->ubfactor);
  }

  DL_ASSERT(check_kwinfo(kwinfo,graph,(pid_t const **)gwhere),"Bad kwinfo");
  DL_ASSERT(check_kwbnd(kwinfo->bnd,graph,1),"Bad boundary");

  q = vw_pq_create(0,mynvtxs); 

  nmoved = 0;
  for (pass=0; pass<niter; pass++) {
    mycut = 0;
    for (c=0;c<2;++c) {
      dlthread_barrier(ctrl->comm);

      /* fill up my queue with my vertices */
      vw_pq_clear(q);
      for (i=0;i<bnd->size;++i) {
        k = vtx_iset_get(i,bnd);
        DL_ASSERT(k < mynvtxs,"Invalid border vertex %"PF_VTX_T,k);
        if (nbrinfo[k].nnbrs > 0) {
          /* only insert vertices with external neighbors */
          rgain = (1.0*nbrinfo[k].ed/sqrt(nbrinfo[k].nnbrs)) - nbrinfo[k].id;
          vw_pq_push(rgain,k,q);
        }
      }
      /* make moves */

      do {
        /* perform updates */
        while (update_combuffer_next(&up,combuffer)) {
          k = up.nbr;
          ewgt = up.ewgt;
          to = up.to;
          from = up.from;

          mycut += __update_vertex(ctrl,k,to,from,ewgt,graph,kwinfo,q);
        }

        /* move a vertice */
        if (q->size > 0) {
          i = vw_pq_pop(q);

          myrinfo = kwinfo->nbrinfo+i;
          mynbrs = kwinfo_get_nbrs_ro(kwinfo,i, \
              dl_min(nparts,graph->xadj[myid][i+1]-graph->xadj[myid][i]));

          from = where[i];
          wgt = vwgt[i];

          if (myrinfo->id > 0 && lpwgts[from]-wgt < minwgt[from]) {
            continue;
          }

          /* find the first eligible partition */
          for (k=0;k<myrinfo->nnbrs; ++k) {
            to = mynbrs[k].pid;
            if (!__right_side(c,to,from)) {
              continue;
            }
            if (lpwgts[to]+wgt <= maxwgt[to]) {
              if (mynbrs[k].ed >= myrinfo->id) {
                break;
              }
            }
          }
          if (k == myrinfo->nnbrs) {
            /* if there aren't any eligable partitions, abort */
            continue;
          }

          /* see if there is a better one from the eligable one */
          for (j=k+1; j<myrinfo->nnbrs; ++j) {
            to = mynbrs[j].pid;
            if (!__right_side(c,to,from)) {
              continue;
            }
            if (mynbrs[j].ed >= mynbrs[k].ed) {
              gain = mynbrs[j].ed-myrinfo->id; 
              DL_ASSERT(gain >= 0, "Invalid gain of %"PF_WGT_T,gain);
              if ((gain > 0 && lpwgts[to]+wgt <= maxwgt[to]) \
                  || (mynbrs[j].ed == mynbrs[k].ed && \
                     tpwgts[mynbrs[k].pid]*lpwgts[to] < \
                     tpwgts[to]*lpwgts[mynbrs[k].pid])) {
                k = j;
              }
            }
          }
          to = mynbrs[k].pid;

          if (mynbrs[k].ed >= myrinfo->id) { 
            gain = mynbrs[k].ed-myrinfo->id;
            if (!(gain > 0 || (gain == 0 \
                      && (lpwgts[from] >= maxwgt[from]  \
                          || tpwgts[to]*lpwgts[from] > \
                          tpwgts[from]*(lpwgts[to]+wgt))))) {
              continue;
            }
          }

          if (lpwgts[to] + wgt > maxwgt[to] || 
              lpwgts[from] - wgt < minwgt[from]) {
            /* whatever you do, don't push the red button */
            continue;
          }

          /* make the move ***************************************************/
          ++nmoved;

          mycut += __move_vertex(ctrl,graph,myid,i,to,kwinfo,lpwgts, \
              where,q,combuffer);
        } 
      } while ((q->size > 0 && vw_pq_top(q) >= 0) || \
          !update_combuffer_finish(combuffer));

      DL_ASSERT_EQUALS(update_combuffer_next(NULL,combuffer),0,"%d");

      update_combuffer_clear(combuffer);

      /* update my partition weights */
      __par_sync_pwgts(myid,nparts,pwgts,lpwgts,ctrl->comm);

    } /* end directions */

    mycut = wgt_dlthread_sumreduce(mycut,ctrl->comm);

    par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH, \
        "Refinement pass %zu: %"PF_WGT_T" improvement\n",pass,mycut);

    if (mycut == 0) {
      break;
    }

    if (myid == 0) {
      graph->mincut -= (mycut/2);
    }
  } /* end passes */

  nmoved = vtx_dlthread_sumreduce(nmoved,ctrl->comm);

  vw_pq_free(q);

  dl_free(minwgt);
  dl_free(maxwgt);
  dl_free(lpwgts);

  DL_ASSERT(check_kwinfo(kwinfo,graph,(pid_t const **)gwhere),"Bad kwinfo");
  DL_ASSERT(check_kwbnd(kwinfo->bnd,graph,1),"Bad boundary");
  DL_ASSERT(graph->mincut >= 0,"Invalid mincut of %"PF_WGT_T, \
      graph->mincut);
  DL_ASSERT_EQUALS(graph_cut(graph,(pid_t const**)gwhere), \
      graph->mincut,"%"PF_WGT_T);

  /* implicit barrier */
  update_combuffer_free(combuffer);

  return nmoved;
}




/******************************************************************************
* PUBLIC PARALLEL FUNCTIONS ***************************************************
******************************************************************************/


vtx_t par_kwayrefine(
    ctrl_t * const ctrl,
    graph_t * const graph,
    kwinfo_t * const kwinfo)
{
  vtx_t nmoves; 

  nmoves = 0;

  switch (ctrl->rtype) {
    case MTMETIS_RTYPE_GREEDY:
      nmoves = __par_kwayrefine_GREEDY(ctrl,graph,ctrl->nrefpass,kwinfo);
      break;
    default:
      dl_error("Unsupported refinement type '%d' for K-Way partitions.", \
          ctrl->rtype);
  }

  par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH,"%zu) [%"PF_VTX_T" %" \
      PF_ADJ_T"] {%"PF_WGT_T" %"PF_VTX_T"}\n",graph->level,graph->nvtxs, \
      graph->nedges,graph->mincut,nmoves);

  return nmoves;
}




#endif

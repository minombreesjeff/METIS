/**
 * @file kwayfm.c
 * @brief KWay refinement routines
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */





#ifndef MTMETIS_REFINE_C
#define MTMETIS_REFINE_C




#include "refine.h"
#include "check.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct move_t {
  pid_t to;
  pid_t from;
  vtx_t mvid;
} move_t;


typedef struct update_t {
  pid_t to;
  pid_t from;
  wgt_t ewgt;
  adj_t adjid;
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


#define DLMEM_PREFIX move
#define DLMEM_TYPE_T move_t
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX update
#define DLMEM_TYPE_T update_t
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


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
static wgt_t __refine_update_vertex(
    ctrl_t * ctrl, 
    vtx_t const k, 
    pid_t const to, 
    pid_t const from, 
    wgt_t const ewgt, 
    graph_t * graph, 
    ucinfo_t * const ucinfo,
    vw_priority_queue_t * queue)
{
  vtx_t l;
  wgt_t oed;
  pid_t na;
  real_t rgain;
  vtx_iset_t * bnd;
  nbrinfo_t * myrinfo;

  tid_t const myid = omp_get_thread_num();

  pid_t * const where = graph->where[myid];
  pid_t const nparts = ctrl->nparts;
  pid_t const me = where[k];

  bnd = ucinfo->bnd;

  /* create my workspace */
  myrinfo = ucinfo->nbrinfo+k;

  oed = myrinfo->ed;
  
  if (myrinfo->nbrstart == NULL_ADJ) {
    na = dl_min(nparts,graph->xadj[myid][k+1]-graph->xadj[myid][k]);
    myrinfo->nbrstart = ucinfo->nnbrpool;
    /* check to see if we need to expand the pool */
    if ((ucinfo->nnbrpool += na) > ucinfo->maxnnbrpool) { 
      ucinfo->maxnnbrpool *= NBRPOOL_EXP_RATE;
      ucinfo->nbrpool = adjinfo_realloc(ucinfo->nbrpool, ucinfo->maxnnbrpool);
    }
  } else {
    na = -1;
  }

  adjinfo_t * mynbrs = ucinfo->nbrpool + myrinfo->nbrstart;

  if (me == to) {
    myrinfo->id += ewgt;
    myrinfo->ed -= ewgt;
  } else if (me == from) {
    myrinfo->id -= ewgt;
    myrinfo->ed += ewgt;
  }
  /* add it to the boundary if necessary */
  if (!vtx_iset_contains(k,bnd)) {
    if (myrinfo->ed >= myrinfo->id) {
      vtx_iset_add(k,bnd);
    }
  } else if (myrinfo->ed < myrinfo->id) {
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

  rgain = ((myrinfo->nnbrs > 0 ? \
        ((real_t)myrinfo->ed) / sqrt(myrinfo->nnbrs) : 0.0) \
        - myrinfo->id);

  if ((me == to || me == from)) {
    if (vw_maxpq_contains(k,queue)) {
      if (rgain >= 0) {
        vw_maxpq_update(rgain,k,queue);
      } else {
        vw_maxpq_delete(k,queue);
      }
    }
  }

  return oed - myrinfo->ed;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


static adj_t ** __rk_nnbrub; /* neighbor update buffer indexes */
static update_t *** __rk_nbrub; /* neighbor update buffer */
static wgt_t ** __rk_gpwgts;
static wgt_t ** __rk_poffs; 
static int ** __rk_done; /* signal when a thread is done */
vtx_t refine_kway(
    ctrl_t * const ctrl, 
    graph_t * const graph,
    size_t const niter, 
    real_t const ffactor,
    ucinfo_t * const ucinfo)
{
  int mydone;
  vtx_t v, c, i, iii, k, nop, nmoved, lvtx, nbrid, totalmoves;
  adj_t j, l, jjj;
  wgt_t gain, ted, vwgt, ewgt, nmincut, nsum;
  pid_t to, from;
  tid_t owner;
  size_t pass;
  real_t rgain, rwgt;
  wgt_t * mypwgts, * poff;
  nbrinfo_t * myrinfo;
  adjinfo_t * mynbrs;
  vtx_iset_t * bnd;
  vw_priority_queue_t * queue;
  wgt_t * minwgt, * maxwgt, * itpwgts;
  real_t * pwgtr;
  vtx_t * maxnub;
  adj_t * nupdate_buffer;
  update_t ** update_buffer;

  /* my pending moves */
  vtx_t nmove_buffer;
  move_t * move_buffer;

  update_t * otb; /* for later */

  tid_t const myid = omp_get_thread_num();
  tid_t const nthreads = omp_get_num_threads();

  /* Link the graph fields */
  adj_t const * const * const gxadj = (adj_t const **)graph->xadj;
  vtx_t const * const * const gadjncy = (vtx_t const **)graph->adjncy;
  wgt_t const * const * const gvwgt = (wgt_t const **)graph->vwgt;
  wgt_t const * const * const gadjwgt = (wgt_t const **)graph->adjwgt;

  pid_t ** const gwhere = graph->where;
  wgt_t * const pwgts = graph->pwgts;
  
  pid_t const nparts = ctrl->nparts;

  vtx_t const mynvtxs = graph->mynvtxs[myid];
  adj_t const * const xadj = gxadj[myid];
  vtx_t const * const adjncy = gadjncy[myid];
  wgt_t const * const adjwgt = gadjwgt[myid];

  nbrinfo_t * const nbrinfo = ucinfo->nbrinfo;
  pid_t * const where = gwhere[myid];

  #pragma omp master
  {
    dl_start_timer(&(ctrl->timers.refinement));
  }

  /* setup shared vectors */
  #pragma omp master
  {
    __rk_gpwgts = r_wgt_alloc(nthreads);
    __rk_poffs = r_wgt_alloc(nthreads);
    __rk_nbrub = (update_t***)malloc(sizeof(update_t**)*nthreads);
    __rk_nnbrub = r_adj_alloc(nthreads);
    __rk_done = r_int_alloc(nthreads);
  }
  #pragma omp barrier

  __rk_done[myid] = &mydone;

  poff = __rk_poffs[myid] = wgt_alloc(nparts);
  mypwgts = __rk_gpwgts[myid] = wgt_alloc(nparts);

  minwgt = wgt_alloc(nparts);
  maxwgt = wgt_alloc(nparts);
  itpwgts = wgt_alloc(nparts);
  pwgtr = real_alloc(nparts);

  maxnub = vtx_alloc(nthreads);
  nupdate_buffer = adj_alloc(nthreads);
  update_buffer = r_update_alloc(nthreads);
   
  totalmoves = 0;

  bnd = ucinfo->bnd;

  /* Setup the weight intervals of the various subdomains */

  for (i=0;i<nparts;++i) {
    itpwgts[i] = ctrl->tpwgts[i]*graph->tvwgt;
    maxwgt[i]  = ctrl->tpwgts[i]*graph->tvwgt*ctrl->ubfactor;
    minwgt[i]  = ctrl->tpwgts[i]*graph->tvwgt*(1.0/ctrl->ubfactor);
  }

  DL_ASSERT(check_info(ucinfo,graph,(pid_t const **)gwhere),"Bad ucinfo");
  DL_ASSERT(check_bnd(ucinfo->bnd,graph),"Bad boundary");

  queue = vw_priority_queue_create(0,mynvtxs); 

  /* my outgoing updates */
  for (i=0;i<nthreads;++i) {
    maxnub[i] = MOVEBUFFERSIZE;
    update_buffer[i] = update_alloc(maxnub[i]);
  }
  adj_set(nupdate_buffer,0,nthreads);


  nmove_buffer = 0;
  move_buffer = move_alloc(MOVEBUFFERSIZE);


  /* let other threads access my update buffer */
  __rk_nbrub[myid] = update_buffer;
  __rk_nnbrub[myid] = nupdate_buffer;

  /*=====================================================================*
  * The top-level refinement loop 
  *=====================================================================*/
  for (pass=0; pass<niter; pass++) {
    nmoved = nmincut = 0;
    for (c=0;c<2;++c) {
      mydone = 0;
      /* calculate pwgt ratios */
      for (i=0;i<nparts;++i) {
        pwgtr[i] = pow(pwgts[i] / (real_t)itpwgts[i],2);
      }
      #pragma omp barrier

      /* fill up my queue with my vertices */
      for (i=0;i<bnd->size;++i) {
        k = vtx_iset_get(i,bnd);
        DL_ASSERT(k < mynvtxs,"Invalid border vertex %"PF_VTX_T,k);
        rgain = ((nbrinfo[k].nnbrs > 0 ? 
             1.0*nbrinfo[k].ed/sqrt(nbrinfo[k].nnbrs) : 0.0) 
             - nbrinfo[k].id) * pwgtr[where[k]];
        vw_maxpq_push(rgain,k,queue);
      }
      /* make moves */

      for (iii=0;;++iii) {
        /* move a vertice */
        wgt_copy(mypwgts,pwgts,nparts);
        wgt_set(poff,0,nparts);
        #pragma omp barrier

        for (jjj=0;jjj<MOVEBUFFERSIZE;) {
          if (queue->size > 0) {
            i = vw_maxpq_pop(queue);
            myrinfo = ucinfo->nbrinfo+i;
            if (myrinfo->nbrstart == NULL_ADJ) {
              mynbrs = NULL;
            } else {
              mynbrs = ucinfo->nbrpool + myrinfo->nbrstart;
            }

            from = where[i];
            vwgt = gvwgt[myid][i];

            if (myrinfo->id > 0 && mypwgts[from]-vwgt < minwgt[from]) {
              continue;
            }

            for (k=0;k<myrinfo->nnbrs; ++k) {
              to = mynbrs[k].pid;
              if (!RIGHT_SIDE(c,to,from)) {
                continue;
              }
              if (mynbrs[k].ed >= myrinfo->id) {
                gain = mynbrs[k].ed-myrinfo->id; 
                if (mypwgts[to]+vwgt <= maxwgt[to]+ffactor*gain) {
                  break;
                }
              }
            }
            if (k == myrinfo->nnbrs) {
              continue;
            }

            for (j=k+1; j<myrinfo->nnbrs; ++j) {
              to = mynbrs[j].pid;
              if (!RIGHT_SIDE(c,to,from)) {
                continue;
              }
              if (mynbrs[j].ed >= mynbrs[k].ed) {
                gain = mynbrs[j].ed-myrinfo->id; 
                if ((gain > 0 && mypwgts[to]+vwgt <= maxwgt[to]+ffactor*gain) \
                    || (mynbrs[j].ed == mynbrs[k].ed && 
                       itpwgts[mynbrs[k].pid]*mypwgts[to] < 
                       itpwgts[to]*mypwgts[mynbrs[k].pid])) {
                  k = j;
                }
              }
            }

            to = mynbrs[k].pid;

            if (mynbrs[k].ed >= myrinfo->id) { 
              gain = mynbrs[k].ed-myrinfo->id;
              if (!(gain > 0 || (gain == 0 \
                        && (mypwgts[from] >= maxwgt[from]  \
                            || itpwgts[to]*mypwgts[from] > \
                            itpwgts[from]*(mypwgts[to]+vwgt))))) {
                continue;
              }
            }

            #ifdef XXX
            /* update my pwgts every once in a while */
            if ((jjj+myid)%3==0) {
              mypwgts[to] = pwgts[to];
              mypwgts[from] = pwgts[from];
              for (j=0;j<nthreads;++j) {
                mypwgts[to] += __rk_poffs[j][to];
                mypwgts[from] += __rk_poffs[j][from];
              }
            }
            #endif

            /****************************************************************
            * The move is good!                                            *
            ***************************************************************/
            if (mypwgts[to] + vwgt > maxwgt[to] || 
                mypwgts[from] - vwgt < minwgt[from]) {
              /* whatever you do, don't push the red button */
              continue;
            }
            mypwgts[to] += vwgt;
            mypwgts[from] -= vwgt;
            poff[to] += vwgt;
            poff[from] -= vwgt;
            where[i] = to;
            
            /* make the move */
            ted = myrinfo->ed;
            myrinfo->ed += myrinfo->id-mynbrs[k].ed;
            dl_swap(myrinfo->id,mynbrs[k].ed);

            /* old minus new */
            nmincut += ted - myrinfo->ed;

            if (mynbrs[k].ed == 0) {
              mynbrs[k] = mynbrs[--myrinfo->nnbrs];
            } else {
              mynbrs[k].pid = from;
            }
           
            if (!vtx_iset_contains(i,bnd)) {
              if (myrinfo->ed >= myrinfo->id) {
                vtx_iset_add(i,bnd);
              }
            } else if (myrinfo->ed < myrinfo->id) {
              vtx_iset_remove(i,bnd);
            }

            /* update vertices */
            v = lvtx_to_gvtx(i,myid,graph->dist);
            /* put it in my move buffer */
            move_buffer[nmove_buffer].mvid = v;
            move_buffer[nmove_buffer].to = to;
            move_buffer[nmove_buffer].from = from;
            for(j=xadj[i];j<xadj[i+1];++j) {
              k = adjncy[j];
              if (k < mynvtxs) {
                /* I own it */
                ewgt = adjwgt[j];
                nmincut += __refine_update_vertex(ctrl,k,to,from,ewgt, \
                    graph,ucinfo,queue);
              }
            }
            ++nmove_buffer;
            ++jjj;
          } else { 
            if (mydone == 0) {
              mydone = 1;
            }
            break;
          }
        }
        /* fix my mypwgts */
        for (i=0;i<nparts;++i) {
          mypwgts[i] = poff[i];
        }
        #pragma omp barrier
      
        /* figure out what we can move */
        #pragma omp for schedule(static)
        for (i=0;i<nparts;++i) {
          nsum = __rk_poffs[0][i];
          for (j=1;j<nthreads;++j) {
            nsum += __rk_poffs[j][i];
          }
          /* I should only scale weight that is in the wrong direction... */
          if (nsum > 0 && nsum+pwgts[i] > maxwgt[i]) {
            rwgt = (pwgts[i]+nsum - maxwgt[i])/(real_t)nsum;
            for(j=0;j<nthreads;++j) {
              __rk_gpwgts[j][i] = __rk_poffs[j][i]*rwgt;
            }
          } else if (nsum < 0 && pwgts[i]+nsum < minwgt[i]) {
            rwgt = (minwgt[i] -pwgts[i]-nsum)/(real_t)(-nsum);
            for(j=0;j<nthreads;++j) {
              __rk_gpwgts[j][i] = __rk_poffs[j][i]*rwgt;
            }
          } else {
            for (j=0;j<nthreads;++j) {
              __rk_gpwgts[j][i] = 0;
            }
          }
        }
        /* count my number of overweight partitions */
        nop = 0;
        for (i=0;i<nparts;++i) {
          if (mypwgts[i] != 0) {
            ++nop;
          }
        }

        /* undo moves */
        for (l=nmove_buffer;l>0&&nop>0;) {
          --l;
          to = move_buffer[l].to;
          from = move_buffer[l].from;
          v = move_buffer[l].mvid;
          i = gvtx_to_lvtx(v,graph->dist);
          if (mypwgts[to] > 0 || mypwgts[from] < 0) {
            /* undo it */
            move_buffer[l] = move_buffer[--nmove_buffer];
            if (mypwgts[to] != 0 && 
                (mypwgts[to] -= gvwgt[myid][i]) < 0) {
              mypwgts[to] = 0;
              --nop;
            }
            if (mypwgts[from] != 0 && 
                (mypwgts[from] += gvwgt[myid][i]) > 0) {
              mypwgts[from] = 0;
              --nop;
            }
            /* un-update the moved vertex */
            where[i] = from;
            myrinfo = nbrinfo+i;
            adjinfo_t * mynbrs = ucinfo->nbrpool + myrinfo->nbrstart;

            ted = myrinfo->ed;

            for (k=0;k<myrinfo->nnbrs;++k) {
              if (mynbrs[k].pid == from) {
                break;
              }
            }
            if (k < myrinfo->nnbrs) {
              myrinfo->ed += myrinfo->id-mynbrs[k].ed;
              dl_swap(myrinfo->id,mynbrs[k].ed);

              /* old minus new */
              nmincut += ted - myrinfo->ed;

              if (mynbrs[k].ed == 0) {
                mynbrs[k] = mynbrs[--myrinfo->nnbrs];
              } else {
                mynbrs[k].pid = to;
              }
            } else {
              /* I moved from a partition I wasn't connected to */
              myrinfo->ed += myrinfo->id;
              mynbrs[myrinfo->nnbrs].pid = to;
              mynbrs[myrinfo->nnbrs++].ed = myrinfo->id;
              nmincut -= myrinfo->id;
              myrinfo->id = 0;
            }
            if (!vtx_iset_contains(i,bnd)) {
              if (myrinfo->ed >= myrinfo->id) {
                vtx_iset_add(i,bnd);
              }
            } else if (myrinfo->ed < myrinfo->id) {
              vtx_iset_remove(i,bnd);
            }

            /* un-update the neighboring vertices I own */
            for(j=xadj[i];j<xadj[i+1];++j) {
              k = adjncy[j];
              if (k < mynvtxs) {
                /* I own it */
                ewgt = adjwgt[j];
                nmincut += __refine_update_vertex(ctrl,k,from,to,ewgt, \
                    graph,ucinfo,queue);
              } else {
                /* they don't know about it yet, and will never know */
              }
            }
          }
        }

        /* re-total the move weight and update remote neighbors */
        wgt_set(mypwgts,0,nparts);
        for (l=0;l<nmove_buffer;++l) {
          v = move_buffer[l].mvid;
          i = gvtx_to_lvtx(v,graph->dist);
          to = move_buffer[l].to;
          from = move_buffer[l].from;
          DL_ASSERT(where[i] == to,"Bad move in move buffer");
          mypwgts[to] += gvwgt[myid][i];
          mypwgts[from] -= gvwgt[myid][i];
          /* remote neighbor updates */
          for (j=xadj[i];j<xadj[i+1];++j) {
            k = adjncy[j];
            if (k >= mynvtxs) {
              nbrid = gvtx_to_tid(k,graph->dist);
              if (nupdate_buffer[nbrid] == maxnub[nbrid]) {
                maxnub[nbrid]*=2;
                update_buffer[nbrid] = update_realloc(update_buffer[nbrid], \
                    maxnub[nbrid]);
              }
              update_buffer[nbrid][nupdate_buffer[nbrid]].adjid = k; 
              update_buffer[nbrid][nupdate_buffer[nbrid]].ewgt = adjwgt[j];
              update_buffer[nbrid][nupdate_buffer[nbrid]].to = to;
              update_buffer[nbrid][nupdate_buffer[nbrid]++].from = from;
            }
          }
        }
        nmoved += nmove_buffer;
        nmove_buffer = 0;
        #pragma omp barrier

        /* update partition weights */
        #pragma omp for schedule (static) nowait
        for (i=0;i<nparts;++i) {
          nsum = __rk_gpwgts[0][i] + pwgts[i];
          for (j=1;j<nthreads;++j) {
            nsum += __rk_gpwgts[j][i];
          }
          pwgts[i] = nsum;
        }

        /* copy from the senders rather than to the recievers */
        nupdate_buffer[myid] = 0;
        for (owner=(myid+1)%nthreads;owner!=myid;owner=(owner+1)%nthreads) {
          if (__rk_nnbrub[owner][myid] > 0) { 
            otb = __rk_nbrub[owner][myid];
            for (j=0;j<__rk_nnbrub[owner][myid];++j) {
              ewgt = otb[j].ewgt;
              to = otb[j].to;
              from = otb[j].from;
              v = otb[j].adjid;
              lvtx = gvtx_to_lvtx(v,graph->dist);
              nmincut += __refine_update_vertex(ctrl,lvtx,to,from,ewgt,graph,
                  ucinfo,queue);
            }
            __rk_nnbrub[owner][myid] = 0;
            mydone = 0; /* make sure I know I'm not done */
          }
        }

        #pragma omp barrier
        for (i=0;i<nthreads;++i) {
          if (*(__rk_done[i]) != 1) {
            break;
          }
        }
        if (i == nthreads) {
          break;
        }
      }
      #pragma omp barrier
      DL_ASSERT_EQUALS(__rk_nnbrub[myid][myid],0,"%"PF_ADJ_T);
    }
    nmoved = vtx_omp_sumreduce(nmoved);
    nmincut = wgt_omp_sumreduce(nmincut);

    #pragma omp master 
    {
      /* update edgecut */
      graph->mincut -= (nmincut/2);

      /* keep track of total moves */
      totalmoves += nmoved;
    }

    DL_ASSERT(graph->mincut >= 0,"Invalid mincut of %"PF_WGT_T, \
        graph->mincut);
    DL_ASSERT(check_info(ucinfo,graph,(pid_t const **)gwhere),"Bad ucinfo");
    DL_ASSERT_EQUALS(graph_cut(ctrl,graph,(pid_t const**)gwhere), \
        graph->mincut,"%"PF_WGT_T);

    par_vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_HIGH, \
        "Refinement pass %zu: %"PF_VTX_T" moves, %"PF_WGT_T" cut " \
        "improvement (%"PF_WGT_T")\n",pass,nmoved,nmincut,graph->mincut);

    /* definitely want to sync before we break */
    #pragma omp barrier
    if (nmoved == 0 || nmincut == 0) {
      break;
    }
  } /* end passes */
  #pragma omp barrier

  vw_priority_queue_free(queue);
  for (i=0;i<nthreads;++i) {
    dl_free(update_buffer[i]);
  }

  dl_free(minwgt);
  dl_free(maxwgt);
  dl_free(itpwgts);
  dl_free(pwgtr);

  dl_free(maxnub);
  dl_free(nupdate_buffer);
  dl_free(update_buffer);
  dl_free(move_buffer);

  dl_free(__rk_gpwgts[myid]);
  dl_free(__rk_poffs[myid]);

  /* the reduction operation acts as a barrier */
  totalmoves = vtx_omp_sumreduce(totalmoves); /* actually redundant */

  DL_ASSERT(check_bnd(ucinfo->bnd,graph),"Bad boundary");

  #pragma omp master 
  {
    dl_free(__rk_gpwgts);
    dl_free(__rk_poffs);
    dl_free(__rk_nbrub);
    dl_free(__rk_nnbrub);
    dl_free(__rk_done);
  }
  #pragma omp barrier
  #pragma omp master
  {
    dl_stop_timer(&(ctrl->timers.refinement));
  }

  return totalmoves;
}



#endif

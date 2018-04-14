/**
 * @file kwayfm.c
 * @brief KWay refinement routines
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */


#include "includes.h"

/**
* @brief Parallel kway-refinement
*
* @param dctrl control strucutre
* @param graph the graph who's partition to refine
* @param niter number of iterations of refinement to perform
* @param ffactor fudge factor (allow positive gains at teh cost of balance)
* @param nthreads number of threads to use
*/
idx_t ParGreedy_KWayOptimize(dctrl_t * const dctrl, dgraph_t *const graph,
  const idx_t niter, const real_t ffactor, idx_t * const * const * const refcon) 
{
  INIT_PARALLEL();

  idx_t v,c,i, ii, iii, j,jj,jjj, k, l, m, pi,pass, gain, ted,
      tid,from, me, to, vwgt,ewgt,nop,nsum, mybnd, nnbrpool,
      seed,nmoved, nmincut,owner,oldnnbrs,lvtx,nbrid;

  real_t rgain,rwgt;
  move_t move;

  idx_t totalmoves = 0;
  ctrl_t * const ctrl = dctrl->ctrl;

  nnbrpool = dctrl->nnbrpool[myid];

  /* Link the graph fields */
  const idx_t nvtxs  = graph->nvtxs;
  const idx_t * const * const gvwgt = (const idx_t **)graph->vwgt;
  const idx_t * const * const gxadj = (const idx_t **)graph->xadj;
  const idx_t * const * const gadjncy = (const idx_t **)graph->adjncy;
  const idx_t * const * const gadjwgt = (const idx_t **)graph->adjwgt;

  idx_t ** const gwhere = graph->where;
  idx_t * const pwgts = graph->pwgts;
  
  const idx_t nparts = ctrl->nparts;

  idx_t * const buffer = dctrl->buffer1;
  idx_t * const partial = dctrl->buffer2;

  idx_t * const * const gpwgts = refcon[0];
  idx_t * const * const poffs = refcon[1];
  idx_t * const * const done = refcon[2];
  idx_t * const mydone = done[myid];

  /* Setup the weight intervals of the various subdomains */
  idx_t minwgt[nparts];
  idx_t maxwgt[nparts];
  idx_t itpwgts[nparts];

  for (i=0;i<nparts;++i) {
    itpwgts[i] = ctrl->tpwgts[i]*graph->tvwgt[0];
    maxwgt[i]  = ctrl->tpwgts[i]*graph->tvwgt[0]*ctrl->ubfactors[0];
    minwgt[i]  = ctrl->tpwgts[i]*graph->tvwgt[0]*(1.0/ctrl->ubfactors[0]);
  }

  real_t pwgtr[nparts];
  ASSERT(checkInfo(dctrl,graph,(const idx_t **)gwhere) == 1);
  ASSERT(checkBND(graph) == 1);

  idx_t nbnd = 0;

  startwctimer(dctrl->refineParPartTmr);

  seed = ctrl->seed + myid;

  const idx_t mynvtxs = graph->mynvtxs[myid];

  idx_t * const poff = poffs[myid];
  idx_t * const mypwgts = gpwgts[myid];
  idx_t * const where = gwhere[myid];

  const idx_t * const xadj = gxadj[myid];
  const idx_t * const adjncy = gadjncy[myid];
  const idx_t * const adjwgt = gadjwgt[myid];

  ckrinfo_t * const ckrinfo = graph->ckrinfo[myid];

  nbnd = graph->mynbnd[myid];
  idx_t * const bndind = graph->bndind[myid];
  idx_t * const bndptr = graph->bndptr[myid];

  ckrinfo_t * myrinfo;
  cnbr_t * mynbrs;

  rpq_t * const queue = rpqCreate(mynvtxs); 

  /* my outgoing updates */
  update_t * update_buffer[nthreads];
  idx_t maxnub[nthreads];
  for (i=0;i<nthreads;++i) {
    maxnub[i] = MOVEBUFFERSIZE;
    update_buffer[i] = (update_t*)
      gk_malloc(sizeof(update_t)*maxnub[i],"X");
  }
  idx_t nupdate_buffer[nthreads];
  iset(nthreads,0,nupdate_buffer);

  /* my pending moves */
  move_t * move_buffer = 
    (move_t*)gk_malloc(sizeof(move_t)*MOVEBUFFERSIZE,"X");
  idx_t nmove_buffer = 0;

  /* let other threads access my update buffer */
  dctrl->tptr_void[myid] = (void*)update_buffer;
  dctrl->tptr_idx[myid] = nupdate_buffer;

  /* give me friendly ways to access the other threads update_buffer's */
  idx_t ** nnbrub = dctrl->tptr_idx; /* neighbor update buffer indexes */
  update_t *** nbrub = (update_t***)dctrl->tptr_void;/* neighbor update buffer */
  update_t * otb; /* for later */

  /*=====================================================================*
  * The top-level refinement loop 
  *=====================================================================*/
  for (pass=0; pass<niter; pass++) {
    nmoved = nmincut = 0;
    for (c=0;c<2;++c) {
      *mydone = 0;
      /* calculate pwgt ratios */
      for (i=0;i<nparts;++i) {
        pwgtr[i] = pow(pwgts[i] / (real_t)itpwgts[i],2);
      }
      #pragma omp barrier
      startwctimer(dctrl->rpqTmr);

      /* fill up my queue with my vertices */
      for (i=0;i<nbnd;++i) {
        j = bndind[i];
        ASSERT(j >= 0 && j < mynvtxs);
        rgain = ((ckrinfo[j].nnbrs > 0 ? 
             1.0*ckrinfo[j].ed/sqrt(ckrinfo[j].nnbrs) : 0.0) 
             - ckrinfo[j].id) * pwgtr[where[j]];
        rpqInsert(queue,j,rgain);
      }
      /* make moves */
      stopwctimer(dctrl->rpqTmr);
      startwctimer(dctrl->emptyQueueTmr);

      for (iii=0;;++iii) {
        /* move a vertice */
        icopy(nparts,pwgts,mypwgts);
        iset(nparts,0,poff);
        #pragma omp barrier
        startwctimer(dctrl->checkMovTmr);

        for (jjj=0;jjj<MOVEBUFFERSIZE;) {
          if ((i = rpqGetTop(queue)) > -1) {
            myrinfo = graph->ckrinfo[myid]+i;
            if (myrinfo->inbr < 0) {
              mynbrs = NULL;
            } else {
              mynbrs = dctrl->nbrpool[myid] + myrinfo->inbr;
            }

            from = where[i];
            vwgt = gvwgt[myid][i];

            if (myrinfo->id > 0 && mypwgts[from]-vwgt < minwgt[from]) {
              continue;
            }

            for (k=myrinfo->nnbrs-1; k>=0; k--) {
              to = mynbrs[k].pid;
              if (!RIGHT_SIDE(c,to,from)) {
                continue;
              }
              gain = mynbrs[k].ed-myrinfo->id; 
              if (gain >= 0 && mypwgts[to]+vwgt <= maxwgt[to]+ffactor*gain)  
                break;
            }
            if (k < 0) {
              continue;
            }

            for (j=k-1; j>=0; j--) {
              to = mynbrs[j].pid;
              if (!RIGHT_SIDE(c,to,from)) {
                continue;
              }
              gain = mynbrs[j].ed-myrinfo->id; 
              if ((mynbrs[j].ed > mynbrs[k].ed && 
                    mypwgts[to]+vwgt <= maxwgt[to]+ffactor*gain) 
                  ||
                  (mynbrs[j].ed == mynbrs[k].ed && 
                   itpwgts[mynbrs[k].pid]*mypwgts[to] < 
                     itpwgts[to]*mypwgts[mynbrs[k].pid]))
                k = j;
            }

            to = mynbrs[k].pid;

            gain = mynbrs[k].ed-myrinfo->id;
            if (!(gain > 0 
                  || (gain == 0  
                      && (mypwgts[from] >= maxwgt[from] 
                          || itpwgts[to]*mypwgts[from] > 
                          itpwgts[from]*(mypwgts[to]+vwgt)
                         )
                     )
                 )
               ) {
              continue;
            }

            /* update my pwgts every once in a while */
            #ifdef XXX
            if ((jjj+myid)%3==0) {
              mypwgts[to] = pwgts[to];
              mypwgts[from] = pwgts[from];
              for (j=0;j<nthreads;++j) {
                mypwgts[to] += poffs[j][to];
                mypwgts[from] += poffs[j][from];
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
            tid = myrinfo->id;
            ted = myrinfo->ed;
            myrinfo->ed += myrinfo->id-mynbrs[k].ed;
            SWAP(myrinfo->id, mynbrs[k].ed, j);

            /* old minus new */
            nmincut += ted - myrinfo->ed;

            if (mynbrs[k].ed == 0) {
              mynbrs[k] = mynbrs[--myrinfo->nnbrs];
            } else {
              mynbrs[k].pid = from;
            }
           
            if (bndptr[i] == LISTEMPTY && myrinfo->ed - myrinfo->id >= 0) {
              BNDInsert(nbnd,bndind,bndptr,i);
            } else if (bndptr[i] != LISTEMPTY && myrinfo->ed-myrinfo->id < 0) {
              BNDDelete(nbnd,bndind,bndptr,i);
            }

            /* update vertices */
            v = LVTX_2_GVTX(i,myid,graph->dshift);
            /* put it in my move buffer */
            move_buffer[nmove_buffer].mvid = v;
            move_buffer[nmove_buffer].to = to;
            move_buffer[nmove_buffer].from = from;
            for(j=xadj[i];j<xadj[i+1];++j) {
              k = adjncy[j];
              ewgt = adjwgt[j];
              owner = GVTX_2_THRID(k,graph->dmask);
              if (owner != myid) {
                /* tell the other thread about later */
              } else { /* I own it */
                lvtx = GVTX_2_LVTX(k,graph->dshift);
                nmincut += updateVertex(dctrl,myid,lvtx,to,from,ewgt,graph,
                    &nbnd,&nnbrpool,queue,-1);
              }
            }
            ++nmove_buffer;
            ++jjj;
          } else { 
            if (*mydone == 0) {
              *mydone = 1;
            }
            break;
          }
        }
        /* fix my mypwgts */
        for (i=0;i<nparts;++i) {
          mypwgts[i] = poff[i];
        }
        #pragma omp barrier
        stopwctimer(dctrl->checkMovTmr);
        startwctimer(dctrl->makMovTmr);
        startwctimer(dctrl->pwgtsTmr);
      
        /* figure out what we can move */
        #pragma omp for schedule(static)
        for (i=0;i<nparts;++i) {
          nsum = poffs[0][i];
          for (j=1;j<nthreads;++j) {
            nsum += poffs[j][i];
          }
          /* I should only scale weight that is in the wrong direction... */
          if (nsum > 0 && nsum+pwgts[i] > maxwgt[i]) {
            rwgt = (pwgts[i]+nsum - maxwgt[i])/(real_t)nsum;
            for(j=0;j<nthreads;++j) {
              gpwgts[j][i] = poffs[j][i]*rwgt;
            }
          } else if (nsum < 0 && pwgts[i]+nsum < minwgt[i]) {
            rwgt = (minwgt[i] -pwgts[i]-nsum)/(real_t)(-nsum);
            for(j=0;j<nthreads;++j) {
              gpwgts[j][i] = poffs[j][i]*rwgt;
            }
          } else {
            for (j=0;j<nthreads;++j) {
              gpwgts[j][i] = 0;
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

        stopwctimer(dctrl->pwgtsTmr);
        startwctimer(dctrl->undoMovTmr);

        /* undo moves */
        for (l=nmove_buffer-1;l>=0&&nop>0;--l) {
          to = move_buffer[l].to;
          from = move_buffer[l].from;
          v = move_buffer[l].mvid;
          i = GVTX_2_LVTX(v,graph->dshift);
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
            myrinfo = ckrinfo+i;
            cnbr_t * mynbrs = dctrl->nbrpool[myid] + myrinfo->inbr;

            tid = myrinfo->id;
            ted = myrinfo->ed;

            for (k=0;k<myrinfo->nnbrs;++k) {
              if (mynbrs[k].pid == from) {
                break;
              }
            }
            if (k < myrinfo->nnbrs) {
              myrinfo->ed += myrinfo->id-mynbrs[k].ed;
              SWAP(myrinfo->id, mynbrs[k].ed, j);

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
            if (bndptr[i] == LISTEMPTY && myrinfo->ed - myrinfo->id >= 0) {
              BNDInsert(nbnd,bndind,bndptr,i);
            } else if (bndptr[i] != LISTEMPTY && myrinfo->ed-myrinfo->id < 0) {
              BNDDelete(nbnd,bndind,bndptr,i);
            }

            /* un-update the neighboring vertices I own */
            for(j=xadj[i];j<xadj[i+1];++j) {
              k = adjncy[j];
              ewgt = adjwgt[j];
              owner = GVTX_2_THRID(k,graph->dmask);
              if (owner != myid) {
                /* they don't know about it yet, and will never know */
              } else { /* I own it */
                lvtx = GVTX_2_LVTX(k,graph->dshift);
                nmincut += updateVertex(dctrl,myid,lvtx,from,to,ewgt,graph,
                    &nbnd,&nnbrpool,queue,1);
              }
            }
          }
        }

        stopwctimer(dctrl->undoMovTmr);
        startwctimer(dctrl->finMovTmr);

        /* re-total the move weight and update remote neighbors */
        iset(nparts,0,mypwgts);
        for (l=0;l<nmove_buffer;++l) {
          v = move_buffer[l].mvid;
          i= GVTX_2_LVTX(v,graph->dshift);
          to = move_buffer[l].to;
          from = move_buffer[l].from;
          ASSERT(where[i] == to);
          mypwgts[to] += gvwgt[myid][i];
          mypwgts[from] -= gvwgt[myid][i];
          /* remote neighbor updates */
          for (j=xadj[i];j<xadj[i+1];++j) {
            nbrid = GVTX_2_THRID(adjncy[j],graph->dmask);
            if (nbrid != myid) {
              if (nupdate_buffer[nbrid] == maxnub[nbrid]) {
                maxnub[nbrid]*=2;
                update_buffer[nbrid] = (update_t*)
                  gk_realloc(update_buffer[nbrid], 
                      sizeof(update_t)*maxnub[nbrid],"X");
              }
              update_buffer[nbrid][nupdate_buffer[nbrid]].adjid = adjncy[j]; 
              update_buffer[nbrid][nupdate_buffer[nbrid]].ewgt = adjwgt[j];
              update_buffer[nbrid][nupdate_buffer[nbrid]].to = to;
              update_buffer[nbrid][nupdate_buffer[nbrid]++].from = from;
            }
          }
        }
        nmoved += nmove_buffer;
        nmove_buffer = 0;
        #pragma omp barrier

        stopwctimer(dctrl->finMovTmr);
        stopwctimer(dctrl->makMovTmr);
        startwctimer(dctrl->updNeiTmr);

        /* update partition weights */
        #pragma omp for schedule (static) nowait
        for (i=0;i<nparts;++i) {
          nsum = gpwgts[0][i] + pwgts[i];
          for (j=1;j<nthreads;++j) {
            nsum += gpwgts[j][i];
          }
          pwgts[i] = nsum;
        }

        /* copy from the senders rather than to the recievers */
        nupdate_buffer[myid] = 0;
        for (owner=(myid+1)%nthreads;owner!=myid;) {
          ASSERT(owner != myid);
          if (nnbrub[owner][myid] > 0) { 
            otb = nbrub[owner][myid];
            for (j=0;j<nnbrub[owner][myid];++j) {
              ewgt = otb[j].ewgt;
              to = otb[j].to;
              from = otb[j].from;
              v = otb[j].adjid;
              nmincut += updateVertex(dctrl,myid,
                  GVTX_2_LVTX(v,graph->dshift),to,from,ewgt,graph,
                  &nbnd,&nnbrpool,queue,0);
            }
            nnbrub[owner][myid] = 0;
            *mydone = 0; /* make sure I know I'm not done */
          }
          owner = (owner+1)%nthreads;
        }

        stopwctimer(dctrl->updNeiTmr);

        #pragma omp barrier
        for (i=0;i<nthreads;++i) {
          if (*(done[i]) != 1) {
            break;
          }
        }
        if (i == nthreads) {
          break;
        }
      }
      #pragma omp barrier
      ASSERT(nnbrub[myid][myid] == 0);

      stopwctimer(dctrl->emptyQueueTmr);
    }
    dl_omp_sumreduce(myid,nmincut,buffer,nthreads);
    dl_omp_sumreduce(myid,nmoved,partial,nthreads);

    ASSERT(nnbrub[myid][myid] == 0);

    #pragma omp master 
    {
      graph->mincut -= (nmincut/2);
      ASSERT(graph->mincut >= 0);
      #ifdef DEBUG
      real_t bal = SerComputeLoadImbalance(dctrl,graph,nparts,
          dctrl->ctrl->pijbm);

      printf("PASS %"PRIDX": NMOVES = %"PRIDX" MINCUT = %"PRIDX"(%"PRIDX") BAL = %5.4f\n",pass,
          nmoved,nmincut, graph->mincut,bal);
      #endif
      ASSERT(checkInfo(dctrl,graph,(const idx_t **)gwhere) == 1);
      ASSERT(SerComputeCut(graph,(const idx_t**)gwhere) 
          == graph->mincut);
      /* keep track of total moves */
      totalmoves += nmoved;
    }

    /* definitely want to sync before we break */
    #pragma omp barrier
    if (nmoved == 0 || nmincut == 0) {
      break;
    }
  } /* end passes */
  #pragma omp barrier
  startwctimer(dctrl->bndTmr);
  stopwctimer(dctrl->bndTmr);

  dctrl->nnbrpool[myid] = nnbrpool;
  graph->mynbnd[myid] = nbnd;

  rpqDestroy(queue);
  for (i=0;i<nthreads;++i) {
    gk_free((void**)&update_buffer[i],LTERM);
  }

  gk_free((void**)&move_buffer,&update_buffer, LTERM);

  dl_omp_sumreduce(myid,nbnd,dctrl->buffer1,nthreads);
  dl_omp_sumreduce(myid,totalmoves,dctrl->buffer2,nthreads);

  stopwctimer(dctrl->refineParPartTmr);

  #pragma omp master 
  {
    graph->nbnd = nbnd;
    ASSERT(checkBND(graph) == 1);
  }
  #pragma omp barrier

  return totalmoves;
}

idx_t updateVertex(dctrl_t * dctrl, const idx_t myid,
    const idx_t k, const idx_t to, const idx_t from, const idx_t ewgt, 
    dgraph_t * graph, idx_t * rnbnd, idx_t * rnnbrpool, rpq_t * queue, 
    const idx_t done) 
{
  idx_t l,oldnnbrs,nbnd,oed,na,oid,nnbrpool;
  idx_t * where = graph->where[myid];

  const idx_t nparts = dctrl->ctrl->nparts;
  const idx_t me = where[k];

  ctrl_t * const ctrl = dctrl->ctrl;

  /* create my workspace */
  ckrinfo_t * myrinfo = graph->ckrinfo[myid]+k;

  nbnd = *rnbnd;
  nnbrpool = *rnnbrpool;
  idx_t * const bndptr = graph->bndptr[myid];
  idx_t * const bndind = graph->bndind[myid];

  oed = myrinfo->ed;
  oid = myrinfo->id;
  
  if (myrinfo->inbr < 0) {
    ASSERT(oldnnbrs == 0);
    na = gk_min(nparts,graph->xadj[myid][k+1]-graph->xadj[myid][k]);
    myrinfo->inbr = nnbrpool;
    /* check to see if we need to expand the pool */
    if ((nnbrpool += na) > dctrl->maxnnbrpool[myid]) { 
      dctrl->maxnnbrpool[myid] *= NBRPOOL_EXP_RATE;
      dctrl->nbrpool[myid] = (cnbr_t*)gk_realloc(dctrl->nbrpool[myid],
          sizeof(cnbr_t)*dctrl->maxnnbrpool[myid],"X");
    }
  } else {
    ASSERT(oldnnbrs >= 0);
    na = -1;
  }

  cnbr_t * mynbrs = dctrl->nbrpool[myid] + myrinfo->inbr;

  if (me == to) {
    myrinfo->id += ewgt;
    myrinfo->ed -= ewgt;
  } else if (me == from) {
    myrinfo->id -= ewgt;
    myrinfo->ed += ewgt;
  }
  /* add it to the boundary if necessary */
  if (bndptr[k] == LISTEMPTY && myrinfo->ed - myrinfo->id >= 0) {
    BNDInsert(nbnd,bndind,bndptr,k);
  } else if (bndptr[k] != LISTEMPTY && myrinfo->ed-myrinfo->id < 0) {
    BNDDelete(nbnd,bndind,bndptr,k);
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
  } else {
    ASSERT(na == -1);
  }

  real_t rgain = ((myrinfo->nnbrs > 0 ? 
        ((real_t)myrinfo->ed) / sqrt(myrinfo->nnbrs) : 0.0)
        - myrinfo->id);

  if ((me == to || me == from || oldnnbrs != myrinfo->nnbrs)) {
    if (queue->locator[k] != -1) {
      if (rgain >= 0) {
        rpqUpdate(queue,k,rgain);
      } else {
        rpqDelete(queue,k);
      }
    }
  }

  *rnbnd = nbnd;
  *rnnbrpool = nnbrpool;

  return oed - myrinfo->ed;
}

/**
 * @file kwayrefine.c
 * @brief KWay refinement routines
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */


#include "includes.h"

/*************************************************************************/
/*! This function is the entry point of cut-based refinement */
/*************************************************************************/
void ParRefineKWay(dctrl_t * dctrl, dgraph_t *orggraph, dgraph_t *graph)
{
  INIT_PARALLEL();

  idx_t i, moved = 0;

  ctrl_t * const ctrl = dctrl->ctrl;
  const idx_t nparts = ctrl->nparts;
  idx_t ** ptr[3];

  /* allocate the refcon */
  #pragma omp master
  {
    ptr[0] = pimalloc(nthreads,"X");
    ptr[1] = pimalloc(nthreads,"X");
    ptr[2] = pimalloc(nthreads,"X");
    dctrl->tptr_void[0] = (void*)ptr; 
  }
  #pragma omp barrier
  idx_t done = 0;
  idx_t ** const * const refcon = (idx_t ***)dctrl->tptr_void[0];
  refcon[0][myid] = imalloc(nparts,"X");
  refcon[1][myid] = imalloc(nparts,"X");
  refcon[2][myid] = &done;

  ASSERT(checkWhere(dctrl,graph,(const idx_t **)graph->where) == 1);

  /* Compute the parameters of the coarsest graph */
  ParComputeKWayPartitionParams(dctrl, graph,refcon[0]);
  
  /* Refine each successively finer graph */
  for (i=0; ;i++) {
    #ifdef DEBUG
    #pragma omp master
    {
      real_t bal = 0.0, curbal;
      idx_t j;
      for (j=0;j<nparts;++j) {
        if ((curbal = (nparts*graph->pwgts[j] / (real_t)graph->tvwgt[0] )) 
            > bal) {
          bal = curbal;
        }
      }
      printf("GRC %"PRIDX": NVTXS = %"PRIDX", NEDGES = %"PRIDX", MINCUT = %"PRIDX", BAL = %f\n",i,
          graph->nvtxs, graph->nedges, graph->mincut,bal);
    }
    #pragma omp barrier
    #endif

    startwctimer(ctrl->RefTmr);

    moved += ParGreedy_KWayOptimize(dctrl, graph, ctrl->niter, 5.0, 
        (idx_t *const * const *)refcon); 

    stopwctimer(ctrl->RefTmr);

    if (graph == orggraph)
      break;

    graph = graph->finer;

    startwctimer(ctrl->ProjectTmr);

    ParProjectKWayPartition(dctrl, graph);

    stopwctimer(ctrl->ProjectTmr);
  }

  gk_free((void**)&refcon[0][myid],&refcon[1][myid],LTERM);
  #pragma omp barrier
  #pragma omp master
  {
    gk_free((void**)&refcon[0],&refcon[1],&refcon[2],LTERM);
    printf("Total vertices moved = %"PRIDX"\n",moved);
  }
}


void ParComputeKWayPartitionParams(dctrl_t * const dctrl, 
    dgraph_t * const graph, idx_t * const * const gpwgts)
{
  INIT_PARALLEL();

  idx_t size[nthreads];

  ctrl_t * const ctrl = dctrl->ctrl;

  const idx_t nparts = ctrl->nparts;
  const idx_t nvtxs  = graph->nvtxs;

  const idx_t * const * const gxadj = (const idx_t **)graph->xadj;
  const idx_t * const * const gvwgt = (const idx_t **)graph->vwgt;
  const idx_t * const * const gadjncy = (const idx_t **)graph->adjncy;
  const idx_t * const * const gadjwgt = (const idx_t **)graph->adjwgt;

  const idx_t * const * const gwhere  = (const idx_t **)graph->where;

  idx_t * const * const gbndind = graph->bndind;
  idx_t * const * const gbndptr = graph->bndptr;

  idx_t * const buffer = dctrl->buffer1;
  idx_t * const partial = dctrl->buffer2;

  idx_t mincut = 0;
  idx_t nbnd = 0;

  idx_t * const pwgts = graph->pwgts;

  /* Compute the required info for refinement */
  ParCnbrpoolReset(dctrl);

  idx_t other,me,i,j,k,pi,v,lvtx,nbrid,na;

  idx_t * const mypwgts = iset(nparts,0,gpwgts[myid]);

  const idx_t mynvtxs = graph->mynvtxs[myid];

  const idx_t * const xadj = gxadj[myid];
  const idx_t * const vwgt = gvwgt[myid];
  const idx_t * const adjncy = gadjncy[myid];
  const idx_t * const adjwgt = gadjwgt[myid];
  const idx_t * const where = gwhere[myid];
  ckrinfo_t * const ckrinfo = graph->ckrinfo[myid];

  idx_t * const bndptr = gbndptr[myid];
  idx_t * const bndind = gbndind[myid];

  ckrinfo_t * myrinfo;
  cnbr_t * mynbrs;
  idx_t nnbrpool = dctrl->nnbrpool[myid];

  iset(mynvtxs,LISTEMPTY,bndptr);

  for (i=0; i<mynvtxs; ++i) {
    mypwgts[where[i]] += vwgt[i];
  }
  #pragma omp barrier
  #pragma omp for schedule(static) 
  for (i=0; i<nparts;++i) {
    k = 0;
    for (j=0; j<nthreads;++j) {
      k += gpwgts[j][i];
    }
    pwgts[i] = k;
  }
  memset(ckrinfo,0,sizeof(ckrinfo_t)*mynvtxs);

  for (i=0; i<mynvtxs; ++i) {
    v = LVTX_2_GVTX(i,myid,graph->dshift);

    myrinfo = ckrinfo+i;

    na = gk_min(nparts,xadj[i+1]-xadj[i]);
    myrinfo->inbr = nnbrpool;
    /* check to see if we need to expand the pool */
    if ((nnbrpool += na) > dctrl->maxnnbrpool[myid]) { 
      dctrl->maxnnbrpool[myid] *= NBRPOOL_EXP_RATE;
      dctrl->nbrpool[myid] = (cnbr_t*)gk_realloc(dctrl->nbrpool[myid],
          sizeof(cnbr_t)*dctrl->maxnnbrpool[myid],"X");
    }

    mynbrs = dctrl->nbrpool[myid] + myrinfo->inbr;
    me = where[i];

    for (j=xadj[i]; j<xadj[i+1]; ++j) {
      nbrid = GVTX_2_THRID(adjncy[j],graph->dmask);
      lvtx = GVTX_2_LVTX(adjncy[j],graph->dshift); 
      if (me == gwhere[nbrid][lvtx])
        myrinfo->id += adjwgt[j];
      else
        myrinfo->ed += adjwgt[j];
    }

    /* Time to compute the particular external degrees */
    if (myrinfo->ed > 0) {
      mincut += myrinfo->ed;
      for (j=xadj[i]; j<xadj[i+1]; ++j) {
        nbrid = GVTX_2_THRID(adjncy[j],graph->dmask);
        lvtx = GVTX_2_LVTX(adjncy[j],graph->dshift); 
        other = gwhere[nbrid][lvtx];
        if (me != other) {
          for (k=0; k<myrinfo->nnbrs; k++) {
            if (mynbrs[k].pid == other) {
              mynbrs[k].ed += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->nnbrs) {
            mynbrs[k].pid = other;
            mynbrs[k].ed  = adjwgt[j];
            myrinfo->nnbrs++;
          }
        }
      }

      /* Only ed-id>=0 nodes are considered to be in the boundary */
      if (myrinfo->ed-myrinfo->id >= 0) {
        BNDInsert(nbnd,bndind,bndptr, i);
      }
      ASSERT(myrinfo->nnbrs > 0);
    } else if (myrinfo->id == 0) {
      BNDInsert(nbnd,bndind,bndptr, i);
    } else {
      myrinfo->inbr = -1;
      nnbrpool -= na;
      ASSERT(myrinfo->nnbrs == 0);
    }
  } 
  graph->mynbnd[myid] = nbnd;
  dctrl->nnbrpool[myid] = nnbrpool;

  dl_omp_sumreduce(myid,mincut,dctrl->buffer1,nthreads);
  dl_omp_sumreduce(myid,nbnd,dctrl->buffer2,nthreads);

  #pragma omp master
  {
    graph->nbnd = nbnd;
    graph->mincut = mincut/2;

    /* the checks */
    ASSERT(checkInfo(dctrl,graph,(const idx_t **)gwhere) == 1);
    ASSERT(checkBND(graph) == 1);
  }
  #pragma omp barrier
}


/*************************************************************************/
/*! This function projects a partition, and at the same time computes the
 parameters for refinement. */
/*************************************************************************/
void ParProjectKWayPartition(dctrl_t * const dctrl, dgraph_t * const graph) 
{
  INIT_PARALLEL();

  ctrl_t * const ctrl = dctrl->ctrl;

  const idx_t nparts = ctrl->nparts;

  dgraph_t * const cgraph = graph->coarser;
  const idx_t *const * const gcwhere = (const idx_t **)cgraph->where;

  const idx_t nvtxs = graph->nvtxs;
  idx_t * const * const gcmap = graph->cmap;

  const idx_t * const * const gxadj = (const idx_t **)graph->xadj;
  const idx_t * const * const gadjncy = (const idx_t **)graph->adjncy;
  const idx_t * const * const gadjwgt = (const idx_t **)graph->adjwgt;

  ParAllocateKWayPartitionMemory(ctrl, graph);

  graph->mincut = cgraph->mincut;

  idx_t ** const gwhere = graph->where;
  
  idx_t * const buffer = dctrl->buffer1;
  idx_t * const partial = dctrl->buffer2;

  idx_t nbnd = 0;

  /* Compute the required info for refinement */
  ParCnbrpoolReset(dctrl);

  startwctimer(dctrl->projectParPartTmr);

  idx_t me,other,tid,ted,i,j,k,pi,lvtx,nbrid,na,istart,iend;

  const idx_t mynvtxs = graph->mynvtxs[myid];

  idx_t * const cmap = gcmap[myid];
  idx_t * const where = gwhere[myid];

  const idx_t * const xadj = gxadj[myid];
  const idx_t * const adjncy = gadjncy[myid];
  const idx_t * const adjwgt = gadjwgt[myid];

  idx_t * const bndind = graph->bndind[myid];
  idx_t * const bndptr = graph->bndptr[myid];

  ckrinfo_t * const ckrinfo = graph->ckrinfo[myid];

  ckrinfo_t * myrinfo;
  cnbr_t * mynbrs;

  startwctimer(dctrl->projectWhereTmr);

  idx_t ned, nid;
  idx_t * ed = imalloc(mynvtxs,"X");
  idx_t * id = imalloc(mynvtxs,"X");

  ned = nid = 0;
  for (i=0;i<mynvtxs;++i) {
    k = cmap[i];
    lvtx = GVTX_2_LVTX(k,graph->dshift);
    nbrid = GVTX_2_THRID(k,graph->dmask);
    where[i] = gcwhere[nbrid][lvtx];
    if (cgraph->ckrinfo[nbrid][lvtx].ed > 0) {
      ed[ned++] = i;
    } else {
      id[nid++] = i;
    }
  }

  idx_t nnbrpool = dctrl->nnbrpool[myid];
  idx_t * htable = ismalloc(nparts,-1,"X");
  iset(mynvtxs,LISTEMPTY,bndptr);
  memset(ckrinfo,0,sizeof(ckrinfo_t)*mynvtxs);

  #pragma omp master
  {
    icopy(nparts,cgraph->pwgts,graph->pwgts);
  }
  #pragma omp barrier

  stopwctimer(dctrl->projectWhereTmr);
  ParFreeGraph(&graph->coarser);
  startwctimer(dctrl->projectCKRTmr);

  startwctimer(dctrl->projectIdTmr);
  for (pi=0;pi<nid;++pi) {
    i = id[pi];
    istart = xadj[i];
    iend = xadj[i+1];

    myrinfo = ckrinfo+i;

    for (tid=0, j=istart; j<iend; j++) {
      tid += adjwgt[j];
    }
    myrinfo->id   = tid;
    myrinfo->inbr = -1;
    if (tid == 0) {
      /* keep islands on the border */
      BNDInsert(nbnd,bndind,bndptr,i);
    }
    ASSERT(myrinfo->ed == 0);
    ASSERT(myrinfo->nnbrs == 0);
  }
  stopwctimer(dctrl->projectIdTmr);
  startwctimer(dctrl->projectEdTmr);
  for (pi=0;pi<ned;++pi) {
    i = ed[pi];

    istart = xadj[i];
    iend = xadj[i+1];

    myrinfo = ckrinfo+i;

    na = gk_min(nparts,xadj[i+1]-xadj[i]);
    myrinfo->inbr = nnbrpool;
    /* check to see if we need to expand the pool */
    if ((nnbrpool += na) > dctrl->maxnnbrpool[myid]) { 
      dctrl->maxnnbrpool[myid] *= NBRPOOL_EXP_RATE;
      dctrl->nbrpool[myid] = (cnbr_t*) gk_realloc(dctrl->nbrpool[myid],
            sizeof(cnbr_t)*dctrl->maxnnbrpool[myid],"X");
    }

    mynbrs = dctrl->nbrpool[myid] + myrinfo->inbr;

    me = where[i];
    for (tid=0, ted=0, j=istart; j<iend; j++) {
      k = adjncy[j];
      lvtx = GVTX_2_LVTX(k,graph->dshift);
      nbrid = GVTX_2_THRID(k,graph->dmask);
      other = gwhere[nbrid][lvtx];
      if (me == other) {
        tid += adjwgt[j];
      } else {
        ted += adjwgt[j];
        if ((k = htable[other]) == -1) {
          htable[other]               = myrinfo->nnbrs;
          mynbrs[myrinfo->nnbrs].pid  = other;
          mynbrs[myrinfo->nnbrs++].ed = adjwgt[j];
        } else {
          mynbrs[k].ed += adjwgt[j];
        }
      }
    }
    myrinfo->id = tid;
    myrinfo->ed = ted;

    if (ted > 0) {
      if (ted-tid >= 0) {
        BNDInsert(nbnd,bndind,bndptr,i);
      }
      for (j=0; j<myrinfo->nnbrs; ++j) {
        htable[mynbrs[j].pid] = -1;
      }
    } else if (tid == 0) {
      BNDInsert(nbnd,bndind,bndptr,i);
    }
    if (myrinfo->nnbrs == 0) {
      nnbrpool -= na;
      myrinfo->inbr = -1;
    }
  }
  stopwctimer(dctrl->projectEdTmr);
  stopwctimer(dctrl->projectCKRTmr);

  dctrl->nnbrpool[myid] = nnbrpool;
  graph->mynbnd[myid] = nbnd;
  gk_free((void**)&htable,&ed,&id,LTERM);

  dl_omp_sumreduce(myid,nbnd,dctrl->buffer1,nthreads);


  stopwctimer(dctrl->projectParPartTmr);

  #pragma omp master
  {
    graph->nbnd = nbnd;
    graph->coarser = NULL;
    ASSERT(checkInfo(dctrl,graph,(const idx_t **)gwhere) == 1);
    ASSERT(checkBND(graph) == 1);
    ASSERT(checkWhere(dctrl,graph,(const idx_t**)gwhere) == 1);
  }
  #pragma omp barrier
}



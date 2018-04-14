/**
* @file mem.c
* @brief Supplimental memory functions from GKlib
* @author Dominique LaSalle <lasalle@cs.umn.edu>
* Copyright 2012, Regents of the University of Minnesota
* @version 1
* @date 2012-06-13
*/

#include "includes.h"

idx_t ** pimalloc(const size_t n, const char * msg)
{
  return (idx_t**)gk_malloc(sizeof(idx_t*)*n,(char*)msg);
}

/**
* @brief Initialize an array of move structures
*
* @param n The size of the array
* @param x The value to set move_t.to to
* @param ptr The pointer to teh array
*
* @return  Return the poitner to teh array
*/
move_t * mset(const idx_t n, idx_t x, move_t * ptr) {
  idx_t i;
  for(i=0;i<n;++i) {
    ptr[i].to = x;
  }
  return ptr;
}

/**
* @brief Allocate memory for making a KWayPartition
*
* @param ctrl The control structure
* @param graph The graph that will be partitioned
* @param nthreads The number of threads to use
*/
void ParAllocateKWayPartitionMemory(ctrl_t * const ctrl, 
    dgraph_t * const graph) 
{
  INIT_PARALLEL();

  #pragma omp master
  {
    graph->pwgts  = imalloc(ctrl->nparts, "X");
    graph->where  = pimalloc(nthreads,  "X");
    graph->mynbnd = imalloc(nthreads,"X");
    graph->bndptr = pimalloc(nthreads,  "X");
    graph->bndind = pimalloc(nthreads,  "X");
    graph->ckrinfo = (ckrinfo_t **)gk_malloc(nthreads*sizeof(ckrinfo_t*), "X");
  }
  #pragma omp barrier
  ASSERT(myid < nthreads);

  const idx_t mynvtxs = graph->mynvtxs[myid];
  graph->where[myid] = imalloc(mynvtxs,"X");
  graph->bndptr[myid] = imalloc(mynvtxs,"X");
  graph->bndind[myid] = imalloc(mynvtxs,"X");
  graph->ckrinfo[myid] = (ckrinfo_t *)gk_malloc(mynvtxs*sizeof(ckrinfo_t),"X");
  #pragma omp barrier
}

dctrl_t * AllocateDCtrl(const idx_t nvtxs, const idx_t nparts,
    const idx_t nthreads) 
{
  idx_t i,j;
  idx_t * options;
  real_t * tpwgts, *ubvec;
  ctrl_t * ctrl;
  dctrl_t * dctrl = (dctrl_t*)gk_malloc(sizeof(dctrl_t),"X");
  memset(dctrl,0,sizeof(dctrl_t));

  options = imalloc(METIS_NOPTIONS,"ParAllocateDCtrl: options");
  METIS_SetDefaultOptions(options);

  /* initialize tpwgts */
  tpwgts = rsmalloc(nparts, -1.0, "PartGraphKway: tpwgts");
  for (i=0;i<nparts;++i) {
    tpwgts[i] = 1.0 / nparts;
  }

  /* initialize ubvec */
  ubvec = rsmalloc(1,UBFACTOR,"PartGraphKway: ubvec");

  /* setup ctrl */
  dctrl->ctrl = ctrl = SetupCtrl(METIS_OP_KMETIS,options,1,nparts,tpwgts,
      ubvec);
  gk_free((void**)&options,&tpwgts,&ubvec,LTERM);
  if (!ctrl) {
    return NULL;
  }
  dctrl->freeCtrlParts=1;

  /* set various run parameters that depend on the graph */
  ctrl->seed = 1;
  ctrl->CoarsenTo = gk_min(PAR_COARSEN_FACTOR*nparts,
      gk_max(nvtxs/(20*gk_log2(nparts)),20*nparts));
  ctrl->nIparts = 5;

  /* communication buffers */
  dctrl->buffer1 = imalloc(falseShareArraySize(nthreads,idx_t),"X");
  dctrl->buffer2 = imalloc(falseShareArraySize(nthreads,idx_t),"X");
  dctrl->buffer3 = rmalloc(falseShareArraySize(nthreads,idx_t),"X");
  dctrl->tptr_idx = (idx_t **) gk_malloc(sizeof(idx_t*)*nthreads,"X");
  dctrl->tptr_real = (real_t **) gk_malloc(sizeof(real_t*)*nthreads,"X");
  dctrl->tptr_void = (void**) gk_malloc(sizeof(void*)*nthreads,"X");

  return dctrl;
}

void FreeDCtrl(dctrl_t ** pdctrl, const idx_t nthreads) {
  idx_t i;
  dctrl_t * dctrl = *pdctrl;
  gk_free((void**)&dctrl->buffer1,&dctrl->buffer2,&dctrl->tptr_idx,
      &dctrl->tptr_real,&dctrl->tptr_void,&dctrl->buffer3,LTERM);
  if (dctrl->freeCtrlParts) {
    FreeCtrl(&dctrl->ctrl);
  } else {
    FreeWorkSpace(dctrl->ctrl);
    gk_free((void**)&dctrl->ctrl,LTERM);
    dctrl->ctrl = NULL;
  }
  if (dctrl->nbrpool != NULL) {
    #pragma omp parallel num_threads(nthreads)
    {
      INIT_PARALLEL();
      gk_free((void**)&dctrl->nbrpool[myid],LTERM);
    }
    gk_free((void**)&dctrl->nbrpool,&dctrl->nnbrpool,
        &dctrl->maxnnbrpool,LTERM);
  }

  gk_free((void**)&dctrl,LTERM);

  *pdctrl = NULL;
}


void ParAllocateRefinementWorkSpace(dctrl_t * dctrl, const dgraph_t * graph)
{
  INIT_PARALLEL();
  #pragma omp master
  {
    dctrl->nbrpool = (cnbr_t**)gk_malloc(sizeof(cnbr_t*)*graph->ndist,"X");
    dctrl->maxnnbrpool = imalloc(nthreads,"X");
    dctrl->nnbrpool = imalloc(nthreads,"X");
  }
  #pragma omp barrier

  dctrl->maxnnbrpool[myid] = (4*graph->nedges) / nthreads;
  dctrl->nbrpool[myid] = (cnbr_t*)
    gk_malloc(sizeof(cnbr_t)*dctrl->maxnnbrpool[myid],"X");
  dctrl->nnbrpool[myid] = 0;
}


void ParCnbrpoolReset(dctrl_t * dctrl)
{
  INIT_PARALLEL();
  dctrl->nnbrpool[myid] = 0;
}

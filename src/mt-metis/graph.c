/**
* @file graph.c
* @brief routine for handling dgraphs
* @author Dominique LaSalle <lasalle@cs.umn.edu>
* @copyright 2013, Regents of the University of Minnesota
* @version 1
* @date 2012-06-12
*/

#include "includes.h"

/**
* @brief Create a graph
*
* @param nthreads Number of threads to use
*
* @return The created graph strucutre
*/
dgraph_t *SerCreateGraph(void)
{
  dgraph_t *graph;

  graph = (dgraph_t *)gk_malloc(sizeof(dgraph_t), "CreateGraph: graph");

  SerInitGraph(graph);

  return graph; 
}

/**
* @brief Initializes the graph
*
* @param graph The graph to initialize
* @param nthreads The number of threads to use
*/
void SerInitGraph(dgraph_t *graph) 
{
  memset((void *)graph, 0, sizeof(dgraph_t));

  /* graph size constants */
  graph->nvtxs     = -1;
  graph->nedges    = -1;
  graph->mincut    = -1;
  graph->minvol    = -1;
  graph->nbnd      = -1;
  graph->maxdeg    = -1;
  graph->ndist     = 0;

  /* memory for the graph structure */
  graph->xadj      = NULL;
  graph->vwgt      = NULL;
  graph->adjncy    = NULL;
  graph->adjwgt    = NULL;
  graph->label     = NULL;
  graph->cmap      = NULL;
  graph->tvwgt     = NULL;
  graph->invtvwgt  = NULL;

  /* by default these are set to true, but the can be explicitly changed afterwards */
  graph->free_xadj   = 1;
  graph->free_vwgt   = 1;
  graph->free_adjncy = 1;
  graph->free_adjwgt = 1;

  /* memory for the partition/refinement structure */
  graph->where     = NULL;
  graph->pwgts     = NULL;
  graph->id        = NULL;
  graph->ed        = NULL;
  graph->bndptr    = NULL;
  graph->bndind    = NULL;
  graph->ckrinfo   = NULL;

  /* linked-list structure */
  graph->coarser   = NULL;
  graph->finer     = NULL;
}

/**
* @brief Create a dgraph based on inputs
*
* @param dctrl The ctrl_t wrapper
* @param nvtxs Number of vertices
* @param xadj Adjacency list pointers
* @param adjncy Adjacency list
* @param vwgt Vertex weights 
* @param adjwgt Edge weights
* @param nthreads Number of threads to run with
*
* @return The Dgraph
*/
dgraph_t *ParSetupGraph(dgraph_t ** r_graph, idx_t * buffer, const idx_t nvtxs, 
    const idx_t * const oxadj, const idx_t * const oadjncy, 
    const idx_t * const oadjwgt, const idx_t * const ovwgt) 
{
  INIT_PARALLEL();

  dgraph_t * graph;

  #pragma omp master
  {
    /* allocate the graph and fill in the fields */
    *r_graph = graph = SerCreateGraph();

    graph->xadj = pimalloc(nthreads,"X");
    graph->vwgt = pimalloc(nthreads,"X");
    graph->adjncy = pimalloc(nthreads,"X");
    graph->adjwgt = pimalloc(nthreads,"X");

    graph->mynvtxs = imalloc(nthreads,"X"); 
    graph->mymaxdeg = imalloc(nthreads,"X");

    graph->gnvtxs = graph->nvtxs = nvtxs;
    graph->ndist = nthreads;
    if (nthreads > 1) {
      graph->dshift = (sizeof(idx_t)*8)-__builtin_clz(nthreads-1);
      graph->dsize = (1 << graph->dshift);
      graph->dmask = graph->dsize -1;
    } else {
      graph->dshift = 0;
      graph->dsize = 1;
      graph->dmask = 0x0;
    }
    *r_graph = graph;
  }
  #pragma omp barrier
  graph = *r_graph;

  idx_t ** const gxadj = graph->xadj;
  idx_t ** const gvwgt = graph->vwgt;
  idx_t ** const gadjncy = graph->adjncy;
  idx_t ** const gadjwgt = graph->adjwgt;

  idx_t i,j,jj,v, maxdeg, adjsize, cycle,cyclestart,cycleend;
  maxdeg = adjsize = 0;

  const idx_t cyclesize = nthreads*BLOCKSIZE;
  const idx_t ncycles = (nvtxs / cyclesize) +1;
  const idx_t lastblocksize = gk_min(
      gk_max((nvtxs%cyclesize) - (myid*BLOCKSIZE),0), BLOCKSIZE);
  const idx_t mynvtxs = graph->mynvtxs[myid] = 
    (nvtxs /cyclesize)*BLOCKSIZE + lastblocksize;

  /* we'll have an xadj of {0} if we don't own any vertices */
  idx_t * const xadj = gxadj[myid] = imalloc(mynvtxs+1,"X");
  if (mynvtxs > 0) {
    idx_t * const vwgt = gvwgt[myid] = imalloc(mynvtxs,"X");
    xadj[0] = j = 0;
    /* do xadj and vwgt */
    for (i=0;i<mynvtxs;++i) {
      /* v = LVTX_2_GVTX(i,myid,graph->dshift); */
      v = (i&BLOCKMASK)+(((LVTX_2_LBID(i)*nthreads)+myid)<<BLOCKSHIFT); 
      xadj[i+1] = (j += (oxadj[v+1]-oxadj[v])); 
      vwgt[i] = ovwgt[v];
      if (xadj[i+1] - xadj[i] > maxdeg) {
        maxdeg = xadj[i+1]-xadj[i];
      }
    }
    /* finish allocations */
    idx_t * const adjncy = gadjncy[myid] = imalloc(xadj[mynvtxs],"X"); 
    idx_t * const adjwgt = gadjwgt[myid] = imalloc(xadj[mynvtxs],"X"); 
    /* do adjncy and adjwgt */
    for (i=0;i<mynvtxs;++i) {
      v = (i&BLOCKMASK)+((LVTX_2_LBID(i)*nthreads+myid)<<BLOCKSHIFT); 
      jj = xadj[i];
      for (j=oxadj[v];j<oxadj[v+1];++j) {
        adjncy[jj] = oadjncy[j] + 
          ((oadjncy[j]/(nthreads*BLOCKSIZE))*(graph->dsize-nthreads)*BLOCKSIZE);  
        adjwgt[jj++] = oadjwgt[j]; 
      }
      ASSERT(jj == xadj[i+1]);
    }
  } else {
    xadj[0] = 0;
    gvwgt[myid] = NULL;
    gadjncy[myid] = NULL;
    gadjwgt[myid] = NULL;
  }
  graph->mymaxdeg[myid] = maxdeg;
  dl_omp_maxreduce(myid,maxdeg,buffer,nthreads);

  #pragma omp master
  {
    graph->maxdeg = maxdeg;
    graph->nedges = oxadj[nvtxs];

    graph->tvwgt    = imalloc(1, "SetupGraph: tvwgts");
    graph->invtvwgt = rmalloc(1, "SetupGraph: invtvwgts");

    graph->tvwgt[0]    = isum(nvtxs, (idx_t*)ovwgt, 1);
    graph->invtvwgt[0] = 1.0/(graph->tvwgt[0] > 0 ? graph->tvwgt[0] : 1);

  }

  #pragma omp barrier
  ASSERT(checkGraph(graph) == 1);

  return graph;
}


/**
* @brief Sets up the structure of a coarse graph
*
* @param graph The graph that has been matched
* @param cnvtxs The number of coarse vertices
* @param dovsize Do vertex size
* @param nthreads The number of threads to use
*
* @return The coarse graph structure
*/
dgraph_t * ParSetupCoarseGraph(idx_t * buffer, dgraph_t ** r_graph, 
    dgraph_t * const graph, const idx_t cnvtxs)
{
  INIT_PARALLEL();

  #pragma omp master
  {
    *r_graph = SerCreateGraph(); 

    (*r_graph)->finer = graph;
    graph->coarser = *r_graph;

    ASSERT(nthreads == graph->ndist);
    (*r_graph)->ndist = graph->ndist;
    (*r_graph)->dsize = graph->dsize;
    (*r_graph)->dshift = graph->dshift;
    (*r_graph)->dmask = graph->dmask;

    (*r_graph)->mymaxdeg = imalloc(nthreads,"X");
    (*r_graph)->mynvtxs = imalloc(nthreads,"X");

    (*r_graph)->xadj = pimalloc(nthreads,"X");
    (*r_graph)->vwgt = pimalloc(nthreads,"X");
    (*r_graph)->adjncy = pimalloc(nthreads,"X");
    (*r_graph)->adjwgt = pimalloc(nthreads,"X");

  }
  #pragma omp barrier
  dgraph_t * const cgraph = *r_graph;
  cgraph->mynvtxs[myid] = cnvtxs;

  idx_t ** const gxadj = cgraph->xadj;   
  idx_t ** const gvwgt = cgraph->vwgt;
  idx_t ** const gadjncy = cgraph->adjncy;
  idx_t ** const gadjwgt = cgraph->adjwgt;

  /* Allocate memory for the coarser graph */
  idx_t mynvtxs = cnvtxs;

  idx_t i,j,k,v,maxid;

  cgraph->xadj[myid]     = imalloc(mynvtxs+1,"X");
  if (mynvtxs > 0) {
    cgraph->vwgt[myid]     = imalloc(mynvtxs, "X");
  } else {
    cgraph->xadj[myid][0] = 0;
    cgraph->vwgt[myid] = NULL;
  }

  cgraph->adjncy[myid]   = NULL;
  cgraph->adjwgt[myid]   = NULL;
  maxid = LVTX_2_GVTX(mynvtxs-1,myid,graph->dshift)+1;

  dl_omp_maxreduce(myid,maxid,buffer,nthreads);
  dl_omp_sumreduce(myid,mynvtxs,buffer,nthreads);
  #pragma omp master
  {
    cgraph->gnvtxs = maxid;
    cgraph->nvtxs = mynvtxs;
    cgraph->tvwgt = imalloc(1, "X");
    cgraph->invtvwgt = rmalloc(1, "X");
    ASSERT(cgraph->gnvtxs >= cgraph->nvtxs);
  }


  #pragma omp barrier

  return cgraph;
}


/**
* @brief Setup derived tvwgt
*
* @param graph The graph to setup
* @param nthreads The number of threads to use
*/
void ParSetupGraph_tvwgt(idx_t * buffer, dgraph_t *const graph)
{
  INIT_PARALLEL();
  idx_t sum = 0;

  ASSERT(nthreads == graph->ndist);
  sum = isum(graph->mynvtxs[myid],graph->vwgt[myid],1);
  dl_omp_sumreduce(myid,sum,buffer,nthreads);

  #pragma omp master 
  {
    graph->tvwgt    = imalloc(1, "SetupGraph: tvwgts");
    graph->invtvwgt = rmalloc(1, "SetupGraph: invtvwgts");

    graph->tvwgt[0]    = sum;
    graph->invtvwgt[0] = 1.0/(graph->tvwgt[0] > 0 ? graph->tvwgt[0] : 1);
  }
}

/**
* @brief Sets up the graphs label info
*
* @param graph The graph to setup the label for
* @param nthreads The number of threads to use
*/
void ParSetupGraph_label(dgraph_t * const graph)
{
  idx_t i;

  #pragma omp master
  {
    if (graph->label == NULL) {
      graph->label = imalloc(graph->nvtxs, "SetupGraph_label: label");
    }
  }
  iincset(graph->nvtxs,0,graph->label);
  #pragma omp barrier
}

/**
* @brief Free a graph strucutre 
*
* @param r_graph The graph to free
* @param nthreads The number of threads to use
*/
void ParFreeGraph(dgraph_t ** r_graph) 
{
  INIT_PARALLEL();
  ASSERT(nthreads == (*r_graph)->ndist);

  idx_t i;
  dgraph_t * const graph = *r_graph;
  
  /* free graph structure */
  if (graph->free_xadj) {
    gk_free((void **)&graph->xadj[myid], LTERM);
  }
  if (graph->cmap != NULL) {
    gk_free((void**)&graph->cmap[myid],LTERM);
  }
  if (graph->free_vwgt) { 
    gk_free((void **)&graph->vwgt[myid], LTERM);
  }
  if (graph->free_adjncy) {
    gk_free((void **)&graph->adjncy[myid], LTERM);
  }
  if (graph->free_adjwgt) {
    gk_free((void **)&graph->adjwgt[myid], LTERM);
  }

  /* free partition/refinement structure */
  ParFreeRData(graph);

  #pragma omp barrier
  #pragma omp master
  {
    gk_free((void **)&graph->xadj,&graph->vwgt,&graph->adjncy,&graph->adjwgt,
        &graph->tvwgt, &graph->invtvwgt, &graph->label,&graph->mynvtxs,
        &graph->mymaxdeg, &graph->cmap,r_graph, LTERM);
  }
}

/**
* @brief Free refinement/partition arrays
*
* @param graph Structure containing arrays to be freed
* @param nthreads Number of threads to use
*/
void ParFreeRData(dgraph_t *graph) 
{
  INIT_PARALLEL();
  ASSERT(graph->ndist == nthreads);


  if (graph->where != NULL) {
    gk_free((void**)&graph->where[myid],LTERM);
  }
  if (graph->ckrinfo != NULL) {
    gk_free((void**)&graph->ckrinfo[myid],LTERM);
  }
  if (graph->id != NULL) {
    gk_free((void**)&graph->id[myid],LTERM);
  }
  if (graph->ed != NULL) {
    gk_free((void**)&graph->ed[myid],LTERM);
  }
  if (graph->bndind != NULL) {
    gk_free((void**)&graph->bndind[myid],LTERM);
  }
  if (graph->bndptr != NULL) {
    gk_free((void**)&graph->bndptr[myid],LTERM);
  }
  if (graph->ckrinfo != NULL) {
    gk_free((void**)&graph->ckrinfo[myid],LTERM);
  }
  if (graph->rename != NULL) {
    gk_free((void**)&graph->rename[myid],LTERM);
  }

  #pragma omp barrier
  #pragma omp master
  {
    /* free partition/refinement structure */
    gk_free((void **)&graph->where, &graph->pwgts, &graph->id, &graph->ed, 
        &graph->bndptr, &graph->bndind, &graph->mynbnd,&graph->ckrinfo,
        &graph->where,&graph->rename,LTERM);
  }
}


real_t ParComputeLoadImbalance(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t nparts, const real_t *pijbm)
{
  INIT_PARALLEL();

  idx_t i, *pwgts;
  real_t max, cur;

  max = 1.0;
  #pragma omp master
  {
    pwgts = graph->pwgts;

    for (i=0;i<nparts;++i) {
         cur = pwgts[i]*pijbm[i];
      if (cur > max) {
        max = cur;
      }
    }
  }
  dl_omp_maxreduce(myid,max,dctrl->buffer3,nthreads);

  return max;
}

real_t SerComputeLoadImbalance(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t nparts, const real_t *pijbm)
{
  idx_t i, *pwgts;
  real_t max, cur;

  max = 1.0;
  pwgts = graph->pwgts;

  for (i=0;i<nparts;++i) {
       cur = pwgts[i]*pijbm[i];
    if (cur > max) {
      max = cur;
    }
  }

  return max;
}

real_t ParComputeLoadImbalanceDiff(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t nparts, const real_t *pijbm, const real_t *ubvec)
{
  INIT_PARALLEL();

  idx_t i, j, k, *pwgts;
  real_t max, cur;
  
  max = -1.0;

  #pragma omp master
  {
    pwgts = graph->pwgts;
    for (k=0;k<nparts;++k) {
      cur = pwgts[k]*pijbm[k] - ubvec[0];
      if (cur > max) {
        max = cur;
      }
    }
  }
  dl_omp_maxreduce(myid,max,dctrl->buffer3,nthreads);

  return max;
}

real_t SerComputeLoadImbalanceDiff(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t nparts, const real_t *pijbm, const real_t *ubvec)
{
  idx_t i, j, k,  *pwgts;
  real_t max, cur;
  
  max = -1.0;

  pwgts = graph->pwgts;
  for (k=0;k<nparts;++k) {
    cur = pwgts[k]*pijbm[k] - ubvec[0];
    if (cur > max) {
      max = cur;
    }
  }

  return max;
}

/**
* @brief Compute the edge cut of a graph
*
* @param graph The graph to compute the edge cut of
* @param where The vector containing partition informatino for each vertex
* @param nthreads The number of threads to use 
*
* @return The total weight of the cut egdes 
*/
idx_t ParComputeCut(dctrl_t * dctrl,const dgraph_t *const graph, 
    const idx_t * const * const where)
{
  INIT_PARALLEL();

  idx_t cut = 0;
  ASSERT(graph->ndist == nthreads);

  idx_t i,j,gvtx, lvtx, nbrid;

  const idx_t mynvtxs = graph->mynvtxs[myid];
  const idx_t * const xadj = graph->xadj[myid];
  const idx_t * const adjncy = graph->adjncy[myid];
  const idx_t * const mywhere = where[myid];

  if (graph->adjwgt == NULL) {
    for (i=0; i<mynvtxs; ++i) {
      for (j=xadj[i]; j<xadj[i+1]; ++j) {
        nbrid = GVTX_2_THRID(adjncy[j],graph->dmask);
        lvtx = GVTX_2_LVTX(adjncy[j],graph->dshift);
        if (mywhere[i] != where[nbrid][lvtx]) {
          cut++;
        }
      }
    }
  } else {
    const idx_t * const adjwgt = graph->adjwgt[myid];
    for (i=0; i<mynvtxs; ++i) {
      for (j=xadj[i]; j<xadj[i+1]; ++j) {
        nbrid = GVTX_2_THRID(adjncy[j],graph->dmask);
        lvtx = GVTX_2_LVTX(adjncy[j],graph->dshift);
        if (mywhere[i] != where[nbrid][lvtx]) {
          cut+=adjwgt[j];
        }
      }
    }
  }
  dl_omp_sumreduce(myid,cut,dctrl->buffer1,nthreads);

  return cut/2;
}

idx_t SerComputeCut(const dgraph_t * const graph,
    const idx_t * const * const where)
{
  idx_t cut = 0;
  idx_t myid;
  const idx_t rnthreads = graph->ndist;
  for (myid=0;myid<graph->ndist;++myid) {
    idx_t i,j,gvtx, lvtx, nbrid;

    const idx_t mynvtxs = graph->mynvtxs[myid];
    const idx_t * const xadj = graph->xadj[myid];
    const idx_t * const adjncy = graph->adjncy[myid];
    const idx_t * const mywhere = where[myid];

    if (graph->adjwgt == NULL) {
      for (i=0; i<mynvtxs; ++i) {
        for (j=xadj[i]; j<xadj[i+1]; ++j) {
          nbrid = GVTX_2_THRID(adjncy[j],graph->dmask);
          lvtx = GVTX_2_LVTX(adjncy[j],graph->dshift);
          if (mywhere[i] != where[nbrid][lvtx]) {
            cut++;
          }
        }
      }
    } else {
      const idx_t * const adjwgt = graph->adjwgt[myid];
      for (i=0; i<mynvtxs; ++i) {
        for (j=xadj[i]; j<xadj[i+1]; ++j) {
          nbrid = GVTX_2_THRID(adjncy[j],graph->dmask);
          lvtx = GVTX_2_LVTX(adjncy[j],graph->dshift);
          if (mywhere[i] != where[nbrid][lvtx]) {
            cut+=adjwgt[j];
          }
        }
      }
    }
  }
  return cut/2;
}


/**
* @brief Computes the communication volume of a graph
*
* @param graph The graph to compute the comv of
* @param where The vector containing the partition informatin for each vertex
* @param nthreads The number of threads to use
*
* @return 
*/
idx_t ParComputeVolume(const dgraph_t *graph, const idx_t *const * const where)
{
  /* so helpful */
  return 0;
}

idx_t SerComputeVolume(const dgraph_t *graph, const idx_t *const * const where)
{
  /* so helpful */
  return 0;
}

/**
* @brief Check if the graph is within the balance constraints + ffactor
*
* @param ctrl Control structure
* @param graph The graph
* @param ffactor The fudge-factor
* @param nthreads The number of threads to use
*
* @return 0 if unbalanced, 1 if balanced
*/
idx_t ParIsBalanced(dctrl_t *dctrl, const dgraph_t *graph, const real_t ffactor)
{
  ctrl_t * ctrl = dctrl->ctrl;
  return (ParComputeLoadImbalanceDiff(dctrl,graph, ctrl->nparts, ctrl->pijbm, 
    ctrl->ubfactors) <= ffactor);
}

graph_t * ParConvertDGraph(dctrl_t * dctrl, dgraph_t * const ograph,
    graph_t ** r_graph)
{
  INIT_PARALLEL();
  
  ASSERT(checkGraph(ograph) == 1);
  const idx_t nvtxs = ograph->nvtxs;
  const idx_t gnvtxs = ograph->gnvtxs;
  const idx_t nedges = ograph->nedges;
  const idx_t ndist = ograph->ndist;

  #pragma omp master
  {
    void * ptr[4];
    ptr[0] = imalloc(nvtxs+1,"X");
    ptr[1] = imalloc(nvtxs,"X");
    ptr[2] = imalloc(nedges,"X");
    ptr[3] = imalloc(nedges,"X");
    dctrl->tptr_void[0] = ptr;
    ograph->rename = pimalloc(ndist,"X");
  }
  #pragma omp barrier
  void ** ptr = dctrl->tptr_void[0];
  idx_t * const xadj = (idx_t*)ptr[0];
  idx_t * const vwgt = (idx_t*)ptr[1];
  idx_t * const adjncy = (idx_t*)ptr[2];
  idx_t * const adjwgt = (idx_t*)ptr[3];

  idx_t ** const rename = ograph->rename;

  /* if this isn't called with the exact number of threads, do it serially */
  idx_t i,j,v,nbrid,u;

  const idx_t mynvtxs = ograph->mynvtxs[myid];

  const idx_t * const oxadj = ograph->xadj[myid];
  const idx_t * const ovwgt = ograph->vwgt[myid];
  const idx_t * const oadjncy = ograph->adjncy[myid];
  const idx_t * const oadjwgt = ograph->adjwgt[myid];

  rename[myid] = imalloc(mynvtxs,"X");

  v = mynvtxs;
  j = oxadj[mynvtxs];

  /* prefix sum for numbering */
  dl_omp_prefixsum(myid,v,dctrl->buffer1,nthreads);
  dl_omp_prefixsum(myid,j,dctrl->buffer2,nthreads);
  ASSERT(j<=nedges);
  ASSERT(v<=nvtxs);

  /* copy xadj and vwgt */
  if (mynvtxs >0) {
    icopy(mynvtxs,(idx_t*)ovwgt,vwgt+v);
    icopy(oxadj[mynvtxs],(idx_t*)oadjncy,adjncy+j);
    icopy(oxadj[mynvtxs],(idx_t*)oadjwgt,adjwgt+j);
    for (i=0;i<mynvtxs;++i) {
      xadj[v] = oxadj[i]+j; 
      rename[myid][i] = v++;
    }
    xadj[v] = oxadj[i]+j;
    ASSERT(v<=nvtxs);
  }
  #pragma omp barrier
  if (mynvtxs > 0) {
    for (i=j;i<xadj[v];++i) {
      u = GVTX_2_LVTX(adjncy[i],ograph->dshift);
      nbrid = GVTX_2_THRID(adjncy[i],ograph->dmask);
      adjncy[i] = rename[nbrid][u];
    }
  }

  #pragma omp barrier
  #pragma omp master
  {
    *r_graph = SetupGraph(dctrl->ctrl,nvtxs,1,xadj,adjncy,vwgt,
      NULL,adjwgt);
    (*r_graph)->free_adjncy = (*r_graph)->free_adjwgt = 1;
    (*r_graph)->free_xadj = (*r_graph)->free_vwgt =1;
    ASSERT((*r_graph)->nedges == ograph->nedges);
    ASSERT(ograph->nedges == (*r_graph)->xadj[(*r_graph)->nvtxs]);
    ASSERT((*r_graph)->nvtxs == ograph->nvtxs);
    (*r_graph)->label = imalloc(nvtxs,"X");

    if (!ograph->label) {
      iincset(nvtxs,0,(*r_graph)->label);
    } else {
      icopy(ograph->nvtxs,ograph->label,(*r_graph)->label);
    }
    ASSERT(CheckGraph((*r_graph),0,1) == 1);
  }
  #pragma omp barrier

  return *r_graph;
}

idx_t computeEEW(dgraph_t * const graph) 
{
  idx_t sum = 0;

  idx_t myid;
  for (myid=0;myid<graph->ndist;++myid) {
    idx_t i,j;

    for (i=0;i<graph->mynvtxs[myid];++i) {
      for (j=graph->xadj[myid][i];j<graph->xadj[myid][i+1];++j) {
        sum += graph->adjwgt[myid][j];
      }
    }
  }

  return sum;
}

void nParReAdjustMemory(const idx_t myid, dctrl_t *const dctrl, 
    dgraph_t * const graph, const idx_t adjsize) 
{
  const idx_t nedges = graph->xadj[myid][graph->mynvtxs[myid]];
  if (adjsize > 4096 && adjsize * 0.75 > nedges) {
    #ifndef NDEBUG
    printf("[%"PRIDX"] Shrinkng adjncy and adjwgt from %"PRIDX" to %"PRIDX"\n",myid,adjsize,
        nedges);
    #endif
    graph->adjncy[myid] = irealloc(graph->adjncy[myid],nedges,"X");
    graph->adjwgt[myid] = irealloc(graph->adjwgt[myid],nedges,"X");
  }
}

#ifdef XXX
void ParExtractGraphParts(dctrl_t * dctrl, dgraph_t * graph, 
    dgraph_t ** sgraphs, const idx_t nparts, const idx_t nthreads) 
{
  const idx_t nvtxs = graph->nvtxs;
  const idx_t nedges = graph->nedges;
  const idx_t * const where = graph->where;
  const xadj_t * const xadj = graph->xadj;
  const idx_t * const vwgt = graph->vwgt;
  const idx_t * const label = graph->label;
  const idx_t * const adjncy = graph->adjncy;
  const idx_t * const adjwgt = graph->adjwgt;

  idx_t gnvtxs[nparts], gnedges[nparts];
  idx_t psnvtxs[nparts][nthreads], psnedges[nparts][nthreads];
  idx_t * rename = imalloc(nvtxs,"X");

  #pragma omp parallel \
    shared(graph,sgraphs,dctrl,gnvtxs,gnedges,psnvtxs,psnedges,rename) \
    default(none) num_threads(nthreads)
  {
    INIT_PARALLEL(nvtxs);

    idx_t i,j, k,l,mypart,istart,iend;
    idx_t * myrename, * svwgt[nparts], * slabel[nparts], *sadjncy[nparts],
          *sadjwgt[nparts];
    const idx_t * mywhere;
    idx_t sadjsize[nparts], snvtxs[nparts], snedges[nparts];
    const xadj_t * myxadj;
    xadj_t * sxadj[nparts];
    xadj_t lxadj;

    /* zero out my counts */
    iset(nparts,0,snvtxs);
    iset(nparts,0,snedges);

    mywhere = where + mystart;
    myxadj = xadj + mystart;
    myrename = rename + mystart;

    /* determine the size of the parts of each graph we own */
    for (i=0;i<mysize;++i) {
      k = mywhere[i];
      ASSERT(k >= 0 && k < nparts);
      myrename[i] = snvtxs[k]++;
      snedges[k] += myxadj[i].end - myxadj[i].start;
    }
    icopy(nparts,snedges,sadjsize);

    /* let everyone know the totals */
    /* prefix sum my vertices -- assumes nparts >= nthreads */
    for (i=0;i<nparts;++i) {
      psnvtxs[i][myid] = snvtxs[i];
      psnedges[i][myid] = snedges[i];
    }
    #pragma omp barrier
    #pragma omp for
    for (i=0;i<nparts;++i) {
      /* add */
      for (j=1;j<rnthreads;++j) {
        psnvtxs[i][j] += psnvtxs[i][j-1];
        psnedges[i][j] += psnedges[i][j-1];
      }
      gnvtxs[i] = psnvtxs[i][rnthreads-1];
      gnedges[i] = psnedges[i][rnthreads-1];
      /* shift */
      for (j=rnthreads-1;j>0;--j) {
        psnvtxs[i][j] = psnvtxs[i][j-1];
        psnedges[i][j] = psnedges[i][j-1];
      }
      psnvtxs[i][0] = 0;
      psnedges[i][0] = 0;
    }

    ASSERT(isum(nparts,gnvtxs,1) == nvtxs);
    ASSERT(isum(nparts,gnedges,1) == nedges);

    /* create the split graphs */
    #pragma omp for
    for (i=0;i<nparts;++i) {
      /* create each graph with only 1 adj section */
      sgraphs[i] = ParSetupSplitGraph(graph,gnvtxs[i],gnedges[i],nthreads);
    }

    /* adjust my renaming */
    for (i=0;i<mysize;++i) {
      myrename[i] += psnvtxs[mywhere[i]][myid];
    }

    /* set my pointers based on the totals */
    for (i=0;i<nparts;++i) {
      svwgt[i] = sgraphs[i]->vwgt;
      slabel[i] = sgraphs[i]->label;
      sxadj[i] = sgraphs[i]->xadj;
      sadjncy[i] = sgraphs[i]->adjncy;
      sadjwgt[i] = sgraphs[i]->adjwgt;
      snvtxs[i] = psnvtxs[i][myid];
      snedges[i] = psnedges[i][myid];
    }

    /* main loop for splitting */
    for (i=mystart;i<myend;++i) {
      mypart = where[i];

      /* local xadj copy */
      sxadj[mypart][snvtxs[mypart]].start = snedges[mypart];

      /* copy over edges */
      istart = xadj[i].start;
      iend = xadj[i].end;
      if (graph->bndptr[i] == -1) { /* interior vertex */
        icopy(iend-istart,(idx_t*)adjncy+istart,
            sadjncy[mypart]+snedges[mypart]);
        icopy(iend-istart,(idx_t*)adjwgt+istart,
            sadjwgt[mypart]+snedges[mypart]);
        snedges[mypart] += iend - istart;
      } else { /* exterior vertex */
        for (j=istart;j<iend;++j) {
          k = adjncy[j];
          if (where[k] == mypart) { /* if its an internal edge */
            sadjncy[mypart][snedges[mypart]] = k;
            sadjwgt[mypart][snedges[mypart]++] = adjwgt[j];
          }
        }
      }

      /* copy over vertex weights */
      svwgt[mypart][snvtxs[mypart]] = vwgt[i];

      /* copy over vertex labels */
      slabel[mypart][snvtxs[mypart]] = label[i];

      /* finish the vertex */
      sxadj[mypart][snvtxs[mypart]++].end = snedges[mypart];
    }

    /* redump edges back to psnedges */
    for (i=0;i<nparts;++i) {
      psnedges[i][myid] = snedges[i] - psnedges[i][myid];
    }

    #pragma omp barrier

    /* update number of edges and tvwgts per graph */
    #pragma omp for
    for (i=0;i<nparts;++i) {
      sgraphs[i]->nedges = 0;
      for (j=0;j<nthreads;++j) {
        sgraphs[i]->nedges += psnedges[i][j];
      }
      ParSetupGraph_tvwgt(sgraphs[i],1);
    }

    /* rename edges */
    for (mypart=0;mypart<nparts;++mypart) {
      for (i=psnvtxs[mypart][myid];i<snvtxs[mypart];++i) {
        for (j=sxadj[mypart][i].start;j<sxadj[mypart][i].end;++j) {
          sadjncy[mypart][j] = rename[sadjncy[mypart][j]];
        }
      }
    }

    #ifndef NDEBUG
    #pragma omp barrier
    #pragma omp master
    {
      /* check nvtxs */
      idx_t ntvtxs, ntedges;
      ntedges = ntvtxs = 0;
      for (i=0;i<nparts;++i) {
        ntvtxs += sgraphs[i]->nvtxs;
        ntedges += sgraphs[i]->nedges;
      }
      ASSERT(ntvtxs == graph->nvtxs);
      ASSERT(ntedges <= graph->nedges);
    }
    #pragma omp for
    for (i=0;i<nparts;++i) {
      ASSERT(checkGraph(sgraphs[i]) == 1);
    }
    #endif
  }
  gk_free((void**)&rename,LTERM);
}


dgraph_t *ParSetupSplitGraph(dgraph_t *graph, const idx_t snvtxs, 
    const idx_t snedges, const idx_t nthreads)
{
  dgraph_t *sgraph;

  sgraph = ParCreateGraph(nthreads);

  sgraph->nvtxs  = snvtxs;
  sgraph->adjsize = sgraph->nedges = snedges;

  sgraph->maxdeg = -1;

  /* Allocate memory for the splitted graph */
  sgraph->xadj        = xmalloc(snvtxs, "SetupSplitGraph: xadj");
  sgraph->vwgt        = imalloc(sgraph->ncon*snvtxs, "SetupSplitGraph: vwgt");
  sgraph->adjncy      = imalloc(sgraph->adjsize,  "SetupSplitGraph: adjncy");
  sgraph->adjwgt      = imalloc(sgraph->adjsize,  "SetupSplitGraph: adjwgt");
  sgraph->label	      = imalloc(snvtxs,   "SetupSplitGraph: label");
  sgraph->tvwgt       = imalloc(sgraph->ncon, "SetupSplitGraph: tvwgt");
  sgraph->invtvwgt    = rmalloc(sgraph->ncon, "SetupSplitGraph: invtvwgt");

  if (graph->vsize)
    sgraph->vsize     = imalloc(snvtxs,   "SetupSplitGraph: vsize");

  return sgraph;
}
#endif


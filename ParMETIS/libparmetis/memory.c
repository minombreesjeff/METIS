/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * Started 2/24/96
 * George
 *
 * $Id: memory.c 10412 2011-06-25 23:12:57Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************/
/*! This function allocate various pools of memory */
/*************************************************************************/
void AllocateWSpace(ctrl_t *ctrl, graph_t *graph)
{
  ctrl->mcore = gk_mcoreCreate(10*graph->nvtxs*sizeof(idx_t));
}


/*************************************************************************/
/*! This function allocates refinement-specific memory for the workspace */
/*************************************************************************/
void AllocateRefinementWorkSpace(ctrl_t *ctrl, idx_t nbrpoolsize)
{
  ctrl->nbrpoolsize     = nbrpoolsize;
  ctrl->nbrpoolcpos     = 0;
  ctrl->nbrpoolreallocs = 0;

  ctrl->cnbrpool = (cnbr_t *)gk_malloc(ctrl->nbrpoolsize*sizeof(cnbr_t), 
                                  "AllocateRefinementWorkSpace: cnbrpool");

}


/*************************************************************************/
/*! This function de-allocate various pools of memory */
/**************************************************************************/
void FreeWSpace(ctrl_t *ctrl)
{
  gk_mcoreDestroy(&ctrl->mcore, 0);
  /*
  printf(" nbrpool statistics\n" 
         "        nbrpoolsize: %12zu   nbrpoolcpos: %12zu\n"
         "    nbrpoolreallocs: %12zu\n\n",
         ctrl->nbrpoolsize,  ctrl->nbrpoolcpos, ctrl->nbrpoolreallocs);
  */

  gk_free((void **)&ctrl->cnbrpool, LTERM);
  ctrl->nbrpoolsize = 0;
  ctrl->nbrpoolcpos = 0;

}


/*************************************************************************/
/*! This function de-allocates memory allocated for the control structures */
/*************************************************************************/
void FreeCtrl(ctrl_t *ctrl)
{
  FreeWSpace(ctrl);

  gkMPI_Comm_free(&(ctrl->gcomm));
}


/*************************************************************************/
/*! This function allocate space from the workspace/heap */
/*************************************************************************/
void *wspacemalloc(ctrl_t *ctrl, size_t nbytes)
{
  return gk_mcoreMalloc(ctrl->mcore, nbytes);
}


/*************************************************************************/
/*! This function allocate space from the core  */
/*************************************************************************/
idx_t *iwspacemalloc(ctrl_t *ctrl, size_t n)
{
  return (idx_t *)wspacemalloc(ctrl, n*sizeof(idx_t));
}

/*************************************************************************/
/*! This function resets the cnbrpool */
/*************************************************************************/
void cnbrpoolReset(ctrl_t *ctrl)
{
  ctrl->nbrpoolcpos = 0;
}


/*************************************************************************/
/*! This function gets the next free index from cnbrpool */
/*************************************************************************/
idx_t cnbrpoolGetNext(ctrl_t *ctrl, idx_t nnbrs)
{
  ctrl->nbrpoolcpos += nnbrs;

  if (ctrl->nbrpoolcpos > ctrl->nbrpoolsize) {
    ctrl->nbrpoolsize += gk_max(10*nnbrs, ctrl->nbrpoolsize/2);

    ctrl->cnbrpool = (cnbr_t *)gk_realloc(ctrl->cnbrpool,
                          ctrl->nbrpoolsize*sizeof(cnbr_t), "cnbrpoolGet: cnbrpool");
    ctrl->nbrpoolreallocs++;
  }

  return ctrl->nbrpoolcpos - nnbrs;
}

/*************************************************************************/
/*! This function allocate space from the core */
/*************************************************************************/
real_t *rwspacemalloc(ctrl_t *ctrl, size_t n)
{
  return (real_t *)wspacemalloc(ctrl, n*sizeof(real_t));
}


/*************************************************************************/
/*! This function allocate space from the core */
/*************************************************************************/
ikv_t *ikvwspacemalloc(ctrl_t *ctrl, size_t n)
{
  return (ikv_t *)wspacemalloc(ctrl, n*sizeof(ikv_t));
}


/*************************************************************************/
/*! This function allocate space from the core */
/*************************************************************************/
rkv_t *rkvwspacemalloc(ctrl_t *ctrl, size_t n)
{
  return (rkv_t *)wspacemalloc(ctrl, n*sizeof(rkv_t));
}


/*************************************************************************
* This function creates a Coarsegraph_t data structure and initializes
* the various fields
**************************************************************************/
graph_t *CreateGraph(void)
{
  graph_t *graph;

  graph = (graph_t *)gk_malloc(sizeof(graph_t), "CreateCoarseGraph: graph");

  InitGraph(graph);

  return graph;
}


/*************************************************************************
* This function creates a Coarsegraph_t data structure and initializes
* the various fields
**************************************************************************/
void InitGraph(graph_t *graph) 
{
  memset(graph, 0, sizeof(graph_t));

  graph->gnvtxs = graph->nvtxs = graph->nedges = graph->nsep = -1;
  graph->nnbrs = graph->nrecv = graph->nsend = graph->nlocal = -1;
  graph->xadj = graph->vwgt = graph->vsize = graph->adjncy = graph->adjwgt = NULL;
  graph->nvwgt = NULL;
  graph->vtxdist = NULL;
  graph->match = graph->cmap = NULL;
  graph->label = NULL;

  graph->peind = NULL;
  graph->sendptr = graph->sendind = graph->recvptr = graph->recvind = NULL;
  graph->imap = NULL;
  graph->pexadj = graph->peadjncy = graph->peadjloc = NULL;
  graph->lperm = NULL;

  graph->slens = graph->rlens = NULL;
  graph->rcand = NULL;

  graph->where = graph->home = graph->lpwgts = graph->gpwgts = NULL;
  graph->lnpwgts = graph->gnpwgts = NULL;
  graph->ckrinfo = NULL;

  graph->nrinfo  = NULL;
  graph->sepind  = NULL;

  graph->coarser = graph->finer = NULL;

}

/*************************************************************************/
/*! This function deallocates any memory stored in a graph */
/*************************************************************************/
void FreeGraph(graph_t *graph) 
{

  /* Graph structure fields */
  gk_free((void **)&graph->xadj, 
         (void **)&graph->vwgt,
         (void **)&graph->nvwgt,
         (void **)&graph->vsize,
         (void **)&graph->adjncy,
         (void **)&graph->adjwgt,
         (void **)&graph->vtxdist, 
         (void **)&graph->home, 
         LTERM);

  FreeNonGraphFields(graph); 

  gk_free((void **)&graph, LTERM);
}


/*************************************************************************/
/*! This function deallocates the non-graph structure fields of a graph
    data structure */
/*************************************************************************/
void FreeNonGraphFields(graph_t *graph) 
{

  gk_free(
      /* Coarsening fields */
      (void **)&graph->match, 
      (void **)&graph->cmap, 

      /* Initial partitioning fields */
      (void **)&graph->label, 

      /* Communication/Setup fields */
      (void **)&graph->peind, 
      (void **)&graph->sendptr, 
      (void **)&graph->sendind, 
      (void **)&graph->recvptr, 
      (void **)&graph->recvind, 
      (void **)&graph->imap,
      (void **)&graph->pexadj,
      (void **)&graph->peadjncy,
      (void **)&graph->peadjloc,
      (void **)&graph->lperm, 

      /* Projection fields */
      (void **)&graph->rlens,
      (void **)&graph->slens,
      (void **)&graph->rcand,

      /* Refinement fields */
      (void **)&graph->where, 
      (void **)&graph->lpwgts, 
      (void **)&graph->gpwgts, 
      (void **)&graph->lnpwgts, 
      (void **)&graph->gnpwgts, 
      (void **)&graph->ckrinfo, 
      (void **)&graph->nrinfo, 
      (void **)&graph->sepind,

      LTERM);
}


/*************************************************************************/
/*! This function deallocates the non-graph and non-setup structure fields 
    of a graph data structure */
/*************************************************************************/
void FreeNonGraphNonSetupFields(graph_t *graph) 
{

  gk_free(
      /* Coarsening fields */
      (void **)&graph->match, 
      (void **)&graph->cmap, 

      /* Initial partitioning fields */
      (void **)&graph->label, 

      /* Projection fields */
      (void **)&graph->rlens,
      (void **)&graph->slens,
      (void **)&graph->rcand,

      /* Refinement fields */
      (void **)&graph->where, 
      (void **)&graph->lpwgts, 
      (void **)&graph->gpwgts, 
      (void **)&graph->lnpwgts, 
      (void **)&graph->gnpwgts, 
      (void **)&graph->ckrinfo, 
      (void **)&graph->nrinfo, 
      (void **)&graph->sepind,

      LTERM);
}


/*************************************************************************/
/*! This function frees any memory allocated for storing the initial graph
    and performs the local to global (i.e., original numbering of the
    adjacency list)
*/
/*************************************************************************/
void FreeInitialGraphAndRemap(graph_t *graph, idx_t wgtflag, idx_t freevsize) 
{
  idx_t i, nedges;
  idx_t *adjncy, *imap;

  nedges = graph->nedges;
  adjncy = graph->adjncy;
  imap   = graph->imap;

  if (imap != NULL) {
    for (i=0; i<nedges; i++)
      adjncy[i] = imap[adjncy[i]];  /* Apply local to global transformation */
  }

  /* Free fields that are not related to the structure of the graph */
  FreeNonGraphFields(graph); 

  /* Free some derived graph-structure fields */
  gk_free((void **)&graph->nvwgt, &graph->home, LTERM);

  if (freevsize)
    gk_free((void **)&graph->vsize, LTERM);
  if ((wgtflag&2) == 0) 
    gk_free((void **)&graph->vwgt, LTERM);
  if ((wgtflag&1) == 0) 
    gk_free((void **)&graph->adjwgt, LTERM);

  gk_free((void **)&graph, LTERM);
}

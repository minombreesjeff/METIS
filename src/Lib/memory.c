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
 * $Id: memory.c,v 1.4 1997/12/22 21:28:38 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function allocates memory for the workspace
**************************************************************************/
void AllocateWorkSpace(CtrlType *ctrl, GraphType *graph, int nparts)
{
  if (ctrl->optype == OP_KMETIS)
    ctrl->wspace.edegrees = (EDegreeType *)GKmalloc(graph->nedges*sizeof(EDegreeType), "AllocateWorkSpace: edegrees");
  else
    ctrl->wspace.edegrees = (EDegreeType *)idxmalloc(graph->nedges, "AllocateWorkSpace: edegrees");

  ctrl->wspace.maxcore = 5*(graph->nvtxs+1) +                 /* Refinement vectors */
                         4*(nparts+1) +                       /* Partition weights etc */
                         2*graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* 2-way refinement */
                         2*(NEG_GAINSPAN+PLUS_GAINSPAN+1)*(sizeof(ListNodeType *)/sizeof(idxtype)) + /* 2-way refinement */
                         20  /* padding for 64 bit machines */
                         ;

  ctrl->wspace.core = idxmalloc(ctrl->wspace.maxcore, "AllocateWorkSpace: maxcore");
  ctrl->wspace.ccore = 0;
}


/*************************************************************************
* This function allocates memory for the workspace
**************************************************************************/
void FreeWorkSpace(CtrlType *ctrl, GraphType *graph)
{
  GKfree(&ctrl->wspace.edegrees, &ctrl->wspace.core, -1);

}

/*************************************************************************
* This function allocate space from the core 
**************************************************************************/
idxtype *idxwspacemalloc(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore += n;
  ASSERT(ctrl->wspace.ccore <= ctrl->wspace.maxcore);
  return ctrl->wspace.core + ctrl->wspace.ccore - n;
}

/*************************************************************************
* This function frees space from the core 
**************************************************************************/
void idxwspacefree(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore -= n;
  ASSERT(ctrl->wspace.ccore >= 0);
}



/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
GraphType *CreateGraph(void)
{
  GraphType *graph;

  graph = (GraphType *)GKmalloc(sizeof(GraphType), "CreateCoarseGraph: graph");

  InitGraph(graph);

  return graph;
}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
void InitGraph(GraphType *graph) 
{
  graph->gdata = graph->rdata = NULL;

  graph->nvtxs = graph->nedges = -1;

  graph->xadj = graph->vwgt = graph->adjncy = graph->adjwgt = NULL;
  graph->adjwgtsum = NULL;
  graph->label = NULL;
  graph->cmap = NULL;

  graph->where = graph->pwgts = NULL;
  graph->id = graph->ed = NULL;
  graph->bndptr = graph->bndind = NULL;
  graph->rinfo = NULL;
  graph->nrinfo = NULL;

  graph->coarser = graph->finer = NULL;

}

/*************************************************************************
* This function deallocates any memory stored in a graph
**************************************************************************/
void FreeGraph(GraphType *graph) 
{

  GKfree(&graph->gdata, &graph->rdata, -1);
  free(graph);
}


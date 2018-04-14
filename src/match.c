/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * match.c
 *
 * This file contains routines that perform the maximal graph matching
 *
 * Started 8/27/94
 * George
 *
 * $Id: match.c,v 1.5 1996/11/05 19:06:48 karypis Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Global variables 
**************************************************************************/
extern CtrlType *__Ctrl;		/* mlevelpart.c */


/*************************************************************************
* This function finds a maximal matching for a graph randomly.
**************************************************************************/
int RM_Match(CoarseGraphType *graph)
{
  int i, j, k, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int *cmap, *match;
  int idx;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomMatch: graph->cmap");
  perm = icoremalloc(graph->nvtxs, "RandomMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      idx = i;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        if (cmap[k] == UNMATCHED) {
          idx = k;
          break;
        }
      }

      cmap[i] = cmap[idx] = nvtxs++;  
      match[i] = idx;
      match[idx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}


/*************************************************************************
* This function implements the Sorted RM matching
**************************************************************************/
int SRM_Match(CoarseGraphType *graph, int flag)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int maxidx, maxvwgt;
  int *cmap, *match, *perm;
  KeyValueType *cand;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "SRM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "SRM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "SRM_Match: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  if (flag == DO_SORT) {
    cand = (KeyValueType *)icoremalloc(2*graph->nvtxs, "SHEM_Match: cand", 1);

    for (ii=graph->nvtxs-1; ii>=0; ii--) {
      i = perm[ii];
      cand[i].key = -graph->vtxs[i]->nedges;
      cand[i].val = i;
    }
    SortKeyValueNodesDec(cand, graph->nvtxs);

    for (ii=graph->nvtxs-1; ii>=0; ii--) 
      perm[ii] = cand[ii].val;

    icorefree(2*graph->nvtxs);
  }

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxidx = i;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          maxidx = k;
          break;
        }
      }

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and heaby edge heuristics
**************************************************************************/
int HEM_Match(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int maxewgt;		        /* The maximum edge-weight and the index */
  int maxidx;
  int *cmap, *match;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomHeavyEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomHeavyEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxidx = i;
      maxewgt = 0;

      for (j=cvtx->nedges-1; j>=0; j--) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && maxewgt < kwgt) {
          maxidx = k;
          maxewgt = kwgt;
        }
      }

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function implements the Sorted HEM matching scheme
**************************************************************************/
int SHEM_Match(CoarseGraphType *graph, int flag)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx, **vtxs;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int maxewgt;		        /* The maximum edge-weight and the index */
  int maxidx, maxvwgt;
  int *cmap, *match, *perm;
  KeyValueType *cand;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  vtxs = graph->vtxs;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "SHEM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "SHEM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "SHEM_Match: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  if (flag == DO_SORT) {
    cand = (KeyValueType *)icoremalloc(2*graph->nvtxs, "SHEM_Match: cand", 1);

    for (ii=graph->nvtxs-1; ii>=0; ii--) {
      i = perm[ii];
      cand[i].key = -vtxs[i]->nedges;
      cand[i].val = i;
    }
    SortKeyValueNodesDec(cand, graph->nvtxs);

    for (i=graph->nvtxs-1; i>=0; i--) 
      perm[i] = cand[i].val;

    icorefree(2*graph->nvtxs);
  }

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = vtxs[i];
      maxvwgt -= cvtx->vwgt;
      maxidx = i;
      maxewgt = 0;

      for (j=cvtx->nedges-1; j>=0; j--) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && vtxs[k]->vwgt < maxvwgt && maxewgt < kwgt) {
          maxidx = k;
          maxewgt = kwgt;
        }
      }

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;

      maxvwgt += cvtx->vwgt;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and light edge heuristics
**************************************************************************/
int LEM_Match(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int minewgt;		        /* The maximum edge-weight and the index */
  int minidx;
  int *cmap, *match;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomLightEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomLightEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      minewgt = 1000000;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED) {
          if (minewgt > kwgt) {
            minidx = k;
            minewgt = kwgt;
          }
        }
      }
      if (minewgt == 1000000)
        minidx = i;

      cmap[i] = cmap[minidx] = nvtxs++;  
      match[minidx] = i;
      match[i] = minidx;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and light edge heuristics
**************************************************************************/
int LEM_Match_W(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int minewgt;		        /* The maximum edge-weight and the index */
  int minidx, maxvwgt;
  int *cmap, *match;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomLightEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomLightEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      minewgt = 1000000;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          if (minewgt > kwgt) {
            minidx = k;
            minewgt = kwgt;
          }
        }
      }
      if (minewgt == 1000000)
        minidx = i;

      cmap[i] = cmap[minidx] = nvtxs++;  
      match[minidx] = i;
      match[i] = minidx;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}




/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and maximu clique edge heuristics
**************************************************************************/
int HCM_Match(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  float maxclique;	        /* The maximum edge-weight and the index */
  int maxidx;
  float tmp;
  int *cmap, *match;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomHeavyEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomHeavyEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxclique = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED) {
          tmp = ((float)(cvtx->cewgt + graph->vtxs[k]->cewgt + 2.0*kwgt)) /
                ((float)(cvtx->vwgt + graph->vtxs[k]->vwgt));
          if (maxclique < tmp) {
            maxidx = k;
            maxclique = tmp;
          }
        }
      }
      if (maxclique == -1)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and maximu clique edge heuristics
**************************************************************************/
int HCM_Match_W(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  float maxclique;	        /* The maximum edge-weight and the index */
  int maxidx, maxvwgt;
  float tmp;
  int *cmap, *match;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomHeavyEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomHeavyEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxclique = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          tmp = ((float)(cvtx->cewgt + graph->vtxs[k]->cewgt + 2.0*kwgt)) /
                ((float)(cvtx->vwgt + graph->vtxs[k]->vwgt));
          if (maxclique < tmp) {
            maxidx = k;
            maxclique = tmp;
          }
        }
      }
      if (maxclique == -1)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}




/*************************************************************************
* This function implements the fast version of Modified HEM
**************************************************************************/
int MHEM_Match(CoarseGraphType *graph)
{
  int i, ii, j, k, l, cewgt, kwgt, nedges;
  VertexType *cvtx;             /* The current vertex */
  EdgeType *edges;
  int nvtxs;                    /* The number of vertices in coersed graph */
  int *perm;                    /* Random permutation array */
  int maxewgt, maxaewgt;        /* The maximum edge-weight and the index */
  int maxidx;
  int *cmap, *match, *lookup;
  int *same, nsame;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "MHEM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "MHEM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "MHEM_Match: perm", 0);
  lookup = icoremalloc(graph->nvtxs, "MHEM_Match: lookup", 0);
  same = icoremalloc(1000, "MHEM_Match: same", 1);
  if (same == NULL) {
    icorefree(2*graph->nvtxs);
    return 0;
  }

  iset(graph->nvtxs, -1, lookup);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      nsame = maxewgt = maxaewgt = 0;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;

        if (cmap[k] == UNMATCHED) {
          if (maxewgt < kwgt) {
            maxidx = k;
            maxewgt = kwgt;
            same[0] = k;
            nsame = 1;
          }
          else { 
            if (maxewgt == kwgt) 
              same[nsame++] = k;
          }
        }
      }

      if (nsame > 2) {
        for (j=0; j<cvtx->nedges; j++)
          lookup[cvtx->edges[j].edge] = ii;

        for (j=0; j<nsame; j+=2) {
          k = same[j];
          edges = graph->vtxs[k]->edges;
          nedges = graph->vtxs[k]->nedges;

          cewgt = 0;
          for (l=0; l<nedges; l++) {
            if (lookup[edges[l].edge] == ii)
              cewgt += edges[l].ewgt; 
          }
          if (maxaewgt < cewgt) {
            maxidx = k;
            maxaewgt = cewgt;
          }
        }
      }

      if (maxewgt == 0)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(1000+2*graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function implements the fast version of Sorted Modified HEM
**************************************************************************/
int SMHEM_Match(CoarseGraphType *graph, int flag)
{
  int i, ii, j, k, l, cewgt, kwgt, nedges;
  VertexType *cvtx;             /* The current vertex */
  EdgeType *edges;
  int nvtxs;                    /* The number of vertices in coersed graph */
  int *perm;                    /* Random permutation array */
  int maxewgt, maxaewgt;        /* The maximum edge-weight and the index */
  int maxidx;
  int *cmap, *match, *lookup;
  int maxvwgt;
  int *same, nsame;
  KeyValueType *cand;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "SMHEM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "SMHEM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "SMHEM_Match: perm", 0);
  lookup = icoremalloc(graph->nvtxs, "SMHEM_Match: lookup", 0);
  same = icoremalloc(1000, "SMHEM_Match: same", 1);

  iset(graph->nvtxs, -1, lookup);
  RandomPermute(perm, graph->nvtxs, 1);

  if (flag == DO_SORT) {
    cand = (KeyValueType *)icoremalloc(2*graph->nvtxs, "SMHEM_Match: cand", 1);
    for (ii=0; ii<graph->nvtxs; ii++) {
      i = perm[ii];
      cand[i].key = -graph->vtxs[i]->nedges;
      cand[i].val = i;
    }
    SortKeyValueNodesDec(cand, graph->nvtxs);

    for (i=0; i<graph->nvtxs; i++)
      perm[i] = cand[i].val;

    icorefree(2*graph->nvtxs);
  }

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      nsame = maxewgt = maxaewgt = 0;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;

        if (cmap[k] == UNMATCHED) {
          if (maxewgt < kwgt && (maxvwgt > cvtx->vwgt+graph->vtxs[k]->vwgt)) {
            maxidx = k;
            maxewgt = kwgt;
            same[0] = k;
            nsame = 1;
          }
          else { 
            if (maxewgt == kwgt && (maxvwgt > cvtx->vwgt+graph->vtxs[k]->vwgt)) 
              same[nsame++] = k;
          }
        }
      }

      if (nsame > 2) {
        for (j=0; j<cvtx->nedges; j++)
          lookup[cvtx->edges[j].edge] = ii;

        for (j=0; j<nsame; j+=2) {
          k = same[j];
          edges = graph->vtxs[k]->edges;
          nedges = graph->vtxs[k]->nedges;

          cewgt = 0;
          for (l=0; l<nedges; l++) {
            if (lookup[edges[l].edge] == ii)
              cewgt += edges[l].ewgt; 
          }
          if (maxaewgt < cewgt) {
            maxidx = k;
            maxaewgt = cewgt;
          }
        }
      }

      if (maxewgt == 0)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(1000+2*graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
int CreateCoarseGraph(CoarseGraphType *graph, int nvtxs)
{
  int v, u, i, j, l, m, vtx, lastedge, nedges, status;
  CoarseGraphType *coarser;
  VertexType *oldvtx, *newvtx;
  EdgeType *oldedges, *edges;
  int *cmap, *match, *htable;
  
  coarser = CreateGraph();
  coarser->nvtxs = nvtxs;

  coarser->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*nvtxs, "CreateCoarseGraph: coarser->vtxs");
  coarser->allvtxs = (VertexType *)GKmalloc(sizeof(VertexType)*nvtxs, "CreateCoarseGraph: coarser->vtxs");
  for (i=0; i<nvtxs; i++)
    coarser->vtxs[i] = coarser->allvtxs+i;

  coarser->finer = graph;
  graph->coarser = coarser;
  cmap = graph->cmap;
  match = graph->match;

  htable = icoremalloc(nvtxs, "CreateCoarseGraph: htable", 0);
  iset(nvtxs, -1, htable);

  edges = GetEdgePool();
  lastedge = 0;

  for (v=0; v<graph->nvtxs; v++) {
    if (v > match[v])
      continue;  /* Already took care of this */

    vtx = cmap[v];
    newvtx = coarser->vtxs[vtx];

    newvtx->vwgt = graph->vtxs[v]->vwgt;
    newvtx->cewgt = graph->vtxs[v]->cewgt;
    newvtx->ewgtsum = graph->vtxs[v]->ewgtsum;
    newvtx->edges = edges + lastedge;
    oldedges = graph->vtxs[v]->edges;
    nedges = graph->vtxs[v]->nedges;

    m = lastedge;
    for (i=0; i<nedges; i++) {
      if ((j = cmap[oldedges[i].edge]) == vtx)
        continue;

      if ((l = htable[j]) == -1) { /* I'm seeing this vertex for the first time */
        htable[j] = m;
        edges[m].edge = j;
        edges[m].ewgt = oldedges[i].ewgt;
        m++;
      }
      else {
        edges[l].ewgt += oldedges[i].ewgt;
      }
    }

    u = match[v];
    if (v != u) { 
      newvtx->vwgt += graph->vtxs[u]->vwgt;
      newvtx->cewgt += graph->vtxs[u]->cewgt;
      newvtx->ewgtsum += graph->vtxs[u]->ewgtsum;
      oldedges = graph->vtxs[u]->edges;
      nedges = graph->vtxs[u]->nedges;

      for (i=0; i<nedges; i++) {
        if ((j = cmap[oldedges[i].edge]) == vtx) {
          /* Add/Remove the compacted edge weight */
          INC_DEC(newvtx->cewgt, newvtx->ewgtsum, 2*oldedges[i].ewgt);
          continue;
        }

        if ((l = htable[j]) != -1) {
          edges[l].ewgt += oldedges[i].ewgt;
        }
        else { /* I'm seeing this vertex for the first time */
          htable[j] = m;
          edges[m].edge = j;
          edges[m].ewgt = oldedges[i].ewgt;
          m++;
        }
      }
    }

    /* Reset htable */
    for (i=lastedge; i<m; i++)
      htable[edges[i].edge] = -1;

    newvtx->nedges = m-lastedge;
    lastedge = m;
  }

  icorefree(nvtxs);

  if ((status = SetEdgePool(lastedge)) == 0) {
    FreeGraph(coarser);
    graph->coarser = NULL;
    GKfree(graph->cmap, graph->match, -1);
    graph->cmap = graph->match = NULL;
  }
  else {
    coarser->nedges = lastedge;
  }

  return status;
}



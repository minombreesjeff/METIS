/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * graph.c
 *
 * This file contains functions that deal with setting up the graphs
 * for METIS.
 *
 * Started 7/25/97
 * George
 *
 * $Id: graph.c,v 1.1 1997/11/04 23:19:14 karypis Exp $
 *
 */

#include <metis.h>

/*************************************************************************
* This function sets up the graph from the user input
**************************************************************************/
void SetUpGraph(GraphType *graph, int OpType, int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int wgtflag)
{
  int i, j, sum, gsize;

  InitGraph(graph);

  graph->nvtxs = nvtxs;
  graph->nedges = xadj[nvtxs];
  graph->xadj = xadj;
  graph->adjncy = adjncy;

  gsize = 0; 
  if ((wgtflag&2) == 0)
    gsize += nvtxs;
  if ((wgtflag&1) == 0)
    gsize += graph->nedges;

  gsize += 2*nvtxs;

  graph->gdata = idxmalloc(gsize, "SetUpGraph: gdata");

  gsize = 0;
  if ((wgtflag&2) == 0) {
    vwgt = graph->vwgt = idxset(nvtxs, 1, graph->gdata);
    gsize += nvtxs;
  }
  else
    graph->vwgt = vwgt;

  if ((wgtflag&1) == 0) {
    adjwgt = graph->adjwgt = idxset(graph->nedges, 1, graph->gdata+gsize);
    gsize += graph->nedges;
  }
  else
    graph->adjwgt = adjwgt;


  graph->adjwgtsum = graph->gdata + gsize;
  gsize += nvtxs;

  for (i=0; i<nvtxs; i++) {
    sum = 0;
    for (j=xadj[i]; j<xadj[i+1]; j++)
      sum += adjwgt[j];
    graph->adjwgtsum[i] = sum;
  }

  graph->cmap = graph->gdata + gsize;
  gsize += nvtxs;

  if (OpType != OP_KMETIS) {
    graph->label = idxmalloc(nvtxs, "SetUpGraph: label");

    for (i=0; i<nvtxs; i++)
      graph->label[i] = i;
  }

  /* RandomizeGraph(graph);*/
}


/*************************************************************************
* This function randomly permutes the adjacency lists of a graph
**************************************************************************/
void RandomizeGraph(GraphType *graph)
{
  int i, j, k, l, tmp, nvtxs;
  idxtype *xadj, *adjncy, *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  for (i=0; i<nvtxs; i++) {
    l = xadj[i+1]-xadj[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = xadj[i] + RandomInRange(l);
      SWAP(adjncy[j], adjncy[k], tmp);
      SWAP(adjwgt[j], adjwgt[k], tmp);
    }
  }
}


/*************************************************************************
* This function checks whether or not partition pid is contigous
**************************************************************************/
int IsConnectedSubdomain(CtrlType *ctrl, GraphType *graph, int pid, int report)
{
  int i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
  idxtype *xadj, *adjncy, *where, *touched, *queue;
  idxtype *cptr;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue = idxmalloc(nvtxs, "IsConnected: queue");
  cptr = idxmalloc(nvtxs, "IsConnected: cptr");

  nleft = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] == pid) 
      nleft++;
  }

  for (i=0; i<nvtxs; i++) {
    if (where[i] == pid) 
      break;
  }

  touched[i] = 1;
  queue[0] = i;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (where[i] == pid && !touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (where[k] == pid && !touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  if (ncmps > 1 && report) {
    printf("The graph has %d connected components in partition %d:\t", ncmps, pid);
    for (i=0; i<ncmps; i++) {
      wgt = 0;
      for (j=cptr[i]; j<cptr[i+1]; j++)
        wgt += graph->vwgt[queue[j]];
      printf("[%5d %5d] ", cptr[i+1]-cptr[i], wgt);
      /*
      if (cptr[i+1]-cptr[i] == 1)
        printf("[%d %d] ", queue[cptr[i]], xadj[queue[cptr[i]]+1]-xadj[queue[cptr[i]]]);
      */
    }
    printf("\n");
  }

  GKfree(&touched, &queue, &cptr, -1);

  return (ncmps == 1 ? 1 : 0);
}


/*************************************************************************
* This function checks whether or not partition pid is contigous
**************************************************************************/
int IsConnected(CtrlType *ctrl, GraphType *graph, int report)
{
  int i, j, k, nvtxs, first, last;
  idxtype *xadj, *adjncy, *touched, *queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue = idxmalloc(nvtxs, "IsConnected: queue");

  touched[0] = 1;
  queue[0] = 0;
  first = 0; last = 1;

  while (first < last) {
    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }

  if (first != nvtxs && report)
    printf("The graph is not connected. It has %d disconnected vertices!\n", nvtxs-first);

  return (first == nvtxs ? 1 : 0);
}


/*************************************************************************
* This function checks whether or not partition pid is contigous
**************************************************************************/
int IsConnected2(GraphType *graph, int report)
{
  int i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
  idxtype *xadj, *adjncy, *where, *touched, *queue;
  idxtype *cptr;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: touched");
  queue = idxmalloc(nvtxs, "IsConnected: queue");
  cptr = idxmalloc(nvtxs, "IsConnected: cptr");

  nleft = nvtxs;
  touched[0] = 1;
  queue[0] = 0;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (!touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  if (ncmps > 1 && report) {
    printf("%d connected components:\t", ncmps);
    for (i=0; i<ncmps; i++) {
      if (cptr[i+1]-cptr[i] > 200)
        printf("[%5d] ", cptr[i+1]-cptr[i]);
    }
    printf("\n");
  }

  GKfree(&touched, &queue, &cptr, -1);

  return (ncmps == 1 ? 1 : 0);
}


/*************************************************************************
* This function returns the number of connected components in cptr,cind
* The separator of the graph is used to split it and then find its components.
**************************************************************************/
int FindComponents(CtrlType *ctrl, GraphType *graph, idxtype *cptr, idxtype *cind)
{
  int i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
  idxtype *xadj, *adjncy, *where, *touched, *queue;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  touched = idxsmalloc(nvtxs, 0, "IsConnected: queue");

  for (i=0; i<graph->nbnd; i++)
    touched[graph->bndind[i]] = 1;

  queue = cind;

  nleft = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2) 
      nleft++;
  }

  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2)
      break;
  }

  touched[i] = 1;
  queue[0] = i;
  first = 0; last = 1;

  cptr[0] = 0;  /* This actually points to queue */
  ncmps = 0;
  while (first != nleft) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      for (i=0; i<nvtxs; i++) {
        if (!touched[i])
          break;
      }
      queue[last++] = i;
      touched[i] = 1;
    }

    i = queue[first++];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (!touched[k]) {
        queue[last++] = k;
        touched[k] = 1;
      }
    }
  }
  cptr[++ncmps] = first;

  free(touched);

  return ncmps;
}



/*************************************************************************
* This function gathers some statistics for degree distribution 
**************************************************************************/
void AnalyzeGraph(GraphType *graph)
{
  int i, j, k, l, nvtxs, maxdegree, median;
  idxtype *xadj, *adjncy, *adjwgt, *vwgt;
  idxtype *dfreq;

  nvtxs = graph->nvtxs; 
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  vwgt = graph->vwgt;

  maxdegree = 0;
  for (i=0; i<nvtxs; i++) {
    if (xadj[i+1]-xadj[i] > maxdegree)
      maxdegree = xadj[i+1]-xadj[i];
  }

  dfreq = idxsmalloc(maxdegree+1, 0, "AnalyzeGraph: maxdegree");

  for (i=0; i<nvtxs; i++) 
    dfreq[xadj[i+1]-xadj[i]]++;
  for (k=j=i=0; i<=maxdegree; i++) {
    if (dfreq[i] > 0) {
      if (k <= nvtxs/2 && k+dfreq[i] > nvtxs/2)
        median = i;
      k += dfreq[i];
      j++;
      printf("%d:%d ", i, dfreq[i]);
      if (j%10 == 0)
        printf("\n");
    }
  }
  printf("\nMedian: %d, \tAverage: %d\n", median, graph->nedges/nvtxs);
  free(dfreq);
}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * serkmetis.c
 *
 * This file contains all the functions for serial kmetis
 *
 * Started 10/19/96
 * George
 *
 * $Id: serkmetis.c,v 1.3 1997/07/18 00:32:13 karypis Exp $
 *
 */

#include <par_kmetis.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
void Ser_KMetis(GraphType *graph, int nparts, float ubfactor)
{
  int i, j, k, maxvwgt;
  GraphType *cgraph;
  EdgeType *degrees;
  float balance;

  degrees = (EdgeType *)GKmalloc(sizeof(EdgeType)*graph->xadj[graph->nvtxs], "Ser_KMetis: degrees");

  maxvwgt = idxsum(graph->nvtxs, graph->vwgt)/(5*nparts);

  cgraph = Ser_Coarsen(graph, 10*nparts, maxvwgt);
  
  Ser_InitPartition(cgraph, nparts, ubfactor);

  Ser_Refine(graph, cgraph, nparts, ubfactor, degrees);

  free(degrees);

  balance = 1.0*nparts*graph->lpwgts[idxamax(nparts, graph->lpwgts)]/(1.0*idxsum(nparts, graph->lpwgts));

  /* printf("Cut: %5d, Balance: %7.4lf (%7.4lf)\n", graph->mincut, balance, ubfactor); */

  if (balance > ubfactor)
    graph->mincut = balance*graph->mincut;
}


/*************************************************************************
* Let the game begin
**************************************************************************/
void Ser_OMetis(GraphType *graph, float ubfactor)
{
  int i, j, k, maxvwgt;
  GraphType *cgraph;
  EdgeType *degrees;
  float balance;

  degrees = (EdgeType *)GKmalloc(sizeof(EdgeType)*graph->xadj[graph->nvtxs], "Ser_KMetis: degrees");

  maxvwgt = idxsum(graph->nvtxs, graph->vwgt) / 20;

  cgraph = Ser_Coarsen(graph, 50, maxvwgt);
  
  Ser_InitPartition(cgraph, 2, ubfactor);

  Ser_Refine(graph, cgraph, 2, ubfactor, degrees);

  GKfree(&degrees, &graph->lpwgts, -1);

  Ser_ConstructSeparator(graph);
  Ser_NodeComputePartitionParams(graph);
  Ser_NodeFM(graph, ubfactor, 2);

}



/*************************************************************************
* This function takes a graph and creates the link list of coarser graphs
**************************************************************************/
GraphType *Ser_Coarsen(GraphType *graph, int CoarsenTo, int maxvwgt)
{
  GraphType *cgraph;

  cgraph = graph;

  while (cgraph->nvtxs > CoarsenTo) {
    Ser_Match_HEM(cgraph, maxvwgt);

    cgraph = cgraph->coarser;

    /* Make sure that we exit if things do not coarsen!  */
    if (cgraph->finer->nvtxs*COARSEN_FRACTION < cgraph->nvtxs) {
      /* printf("\nSlow progress in coarsening. Aborted at %d nodes\n", cgraph->nvtxs); */
      break;
    }
  }

  return cgraph;
}


/*************************************************************************
* This function finds a matching
**************************************************************************/
void Ser_Match_HEM(GraphType *graph, int maxvwgt)
{
  int i, ii, j, nvtxs, cnvtxs, maxidx, maxwgt;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *match, *cmap, *perm;


  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  perm = idxmalloc(nvtxs, "Ser_Match_HEM: perm");
  RandomPermute(nvtxs, perm, 1);

  match = graph->match = idxsmalloc(nvtxs, UNMATCHED, "Ser_Match_HEM: match");
  cmap = graph->cmap = idxmalloc(nvtxs, "Ser_Match_HEM: cmap");

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;

      /* Find a heavy-edge matching */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (match[adjncy[j]] == UNMATCHED && maxwgt < adjwgt[j] && vwgt[i]+vwgt[adjncy[j]] <= maxvwgt) {
          maxwgt = adjwgt[j];
          maxidx = adjncy[j];
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  Ser_CreateCoarseGraph(graph, cnvtxs, perm);

  free(perm);
}




/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void Ser_CreateCoarseGraph(GraphType *graph, int cnvtxs, idxtype *perm)
{
  int i, j, k, nvtxs, nedges, cnedges, v, u;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *match, *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cadjncy, *cadjwgt;
  GraphType *cgraph;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  match = graph->match;
  cmap = graph->cmap;

  /* Initialize the coarser graph */
  cgraph = CreateGraph();
  cgraph->maxvwgt = graph->maxvwgt;
  cgraph->nvtxs = cnvtxs;

  cgraph->finer = graph;
  graph->coarser = cgraph;


  /*************************************************************
  * Create the coarser graph
  **************************************************************/
  /* Allocate memory for the coarser graph, and fire up coarsening */
  cxadj = cgraph->xadj = idxmalloc(cnvtxs+1, "Ser_CreateCoarserGraph: cxadj");
  cvwgt = cgraph->vwgt = idxmalloc(cnvtxs, "Ser_CreateCoarserGraph: cvwgt");
  cadjncy = cgraph->adjncy = idxmalloc(xadj[nvtxs], "Ser_CreateCoarserGraph: xadjncy");
  cadjwgt = cgraph->adjwgt = idxmalloc(xadj[nvtxs], "Ser_CreateCoarserGraph: xadjwgt");

  htable = idxsmalloc(cnvtxs, -1, "Ser_CreateCoarserGraph: htable");

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs)
      continue;

    cvwgt[cnvtxs] = vwgt[v];
    nedges = 0;

    for (j=xadj[v]; j<xadj[v+1]; j++) {
      k = cmap[adjncy[j]];
      if (k != cnvtxs) {  /* If this is not an internal edge */
        if (htable[k] == -1) {
          cadjncy[cnedges+nedges] = k;
          cadjwgt[cnedges+nedges] = adjwgt[j];
          htable[k] = nedges++;
        }
        else {
          cadjwgt[cnedges+htable[k]] += adjwgt[j];
        }
      }
    }

    u = match[v];
    if (v != u) { 
      cvwgt[cnvtxs] += vwgt[u];

      for (j=xadj[u]; j<xadj[u+1]; j++) {
        k = cmap[adjncy[j]];
        if (k != cnvtxs) {  /* If this is not an internal edge */
          if (htable[k] == -1) {
            cadjncy[cnedges+nedges] = k;
            cadjwgt[cnedges+nedges] = adjwgt[j];
            htable[k] = nedges++;
          }
          else {
            cadjwgt[cnedges+htable[k]] += adjwgt[j];
          }
        }
      }
    }

    cnedges += nedges;
    for (j=cxadj[cnvtxs]; j<cnedges; j++)
      htable[cadjncy[j]] = -1;  /* Zero out the htable */
    cxadj[++cnvtxs] = cnedges;
  }

  cgraph->nedges = cnedges;

  free(htable);

}






/*************************************************************************
* This function is the entry point of the initial partition algorithm
* that does recursive bissection.
**************************************************************************/
void Ser_InitPartition(GraphType *graph, int nparts, float ubfactor)
{
  int i, j, ndim, nvtxs, penum, firstvtx;
  idxtype *where, *label, *bisection;

  nvtxs = graph->nvtxs;
  graph->level = -2;

  where = idxmalloc(nvtxs, "Ser_InitPartition: where");
  label = graph->label = idxmalloc(nvtxs, "Ser_InitPartition: label");
  for (i=0; i<nvtxs; i++)
    label[i] = i;

  bisection = idxmalloc(nvtxs, "Ser_InitPartition: bisection");
  Ser_RecursiveBisection(graph, where, nparts, 1.01, 0, bisection);

  graph->where = where;

  free(bisection);
}




/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void Ser_RecursiveBisection(GraphType *graph, idxtype *part, int nparts, float ubfactor, int fpart, idxtype *bisection)
{
  int i, nvtxs, tvwgt, zeropwgt;
  GraphType lgraph, rgraph;

  nvtxs = graph->nvtxs;

  tvwgt = idxsum(nvtxs, graph->vwgt);
  zeropwgt = ((tvwgt+1)*(nparts>>1))/nparts;

  if (nparts%2 == 0 && nparts > 2 && nvtxs > 200)
    Ser_MlevelBisectGraph(graph, ubfactor, bisection);
  else
    Ser_BisectGraph(graph, zeropwgt, ubfactor, bisection);


  /* printf("ICut: %d\n", graph->mincut); */

  for (i=0; i<nvtxs; i++)
    part[graph->label[i]] = bisection[i] + fpart;

  if (nparts > 2)
    Ser_SplitGraphPart(graph, bisection, &lgraph, &rgraph);

  if (graph->level == -1)
    GKfree(&graph->xadj, &graph->vwgt, &graph->label, &graph->adjncy, &graph->adjwgt, -1); 

  if (nparts > 3) {
    Ser_RecursiveBisection(&rgraph, part, (nparts+1)/2, ubfactor, fpart+nparts/2, bisection);
    Ser_RecursiveBisection(&lgraph, part, nparts/2, ubfactor, fpart, bisection);
  }
  else if (nparts == 3) {
    Ser_RecursiveBisection(&rgraph, part, (nparts+1)/2, ubfactor, fpart+nparts/2, bisection);
    GKfree(&lgraph.xadj, &lgraph.vwgt, &lgraph.label, &lgraph.adjncy, &lgraph.adjwgt, -1);
  }

}


/*************************************************************************
* This function performs multilevel bisection
**************************************************************************/
void Ser_MlevelBisectGraph(GraphType *graph, float ubfactor, idxtype *part)
{
  Ser_KMetis(graph, 2, ubfactor);

  idxcopy(graph->nvtxs, graph->where, part);

  GKfree(&graph->where, &graph->match, &graph->cmap, &graph->rinfo, &graph->lpwgts, -1);
  graph->where = graph->match = graph->cmap = graph->lpwgts = NULL;
  graph->rinfo = NULL;
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void Ser_BisectGraph(GraphType *graph, int zeropwgt, float ubfactor, idxtype *part)
{
  int i, j, k, nvtxs, drain, nleft, first, last, tvwgt, pwgts[2], minpwgt[2], maxpwgt[2], from, bestcut, icut, mincut, me, pass, nbfs;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *queue, *touched, *gain, *bestpart;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  bestpart = idxmalloc(nvtxs, "BisectGraph: bestpart");
  queue = idxmalloc(nvtxs, "BisectGraph: queue");
  touched = idxmalloc(nvtxs, "BisectGraph: touched");

  tvwgt = idxsum(nvtxs, vwgt);
  maxpwgt[0] = ubfactor*zeropwgt;
  maxpwgt[1] = ubfactor*(tvwgt-zeropwgt);
  minpwgt[0] = (1.0/ubfactor)*zeropwgt;
  minpwgt[1] = (1.0/ubfactor)*(tvwgt-zeropwgt);

  for (nbfs=0; nbfs<NIPARTS; nbfs++) {
    idxset(nvtxs, 0, touched);

    pwgts[1] = tvwgt;
    pwgts[0] = 0;

    idxset(nvtxs, 1, part);

    queue[0] = RandomInRange(nvtxs);
    touched[queue[0]] = 1;
    first = 0; last = 1;
    nleft = nvtxs-1;
    drain = 0;

    /* Start the BFS from queue to get a partition */
    for (;;) {
      if (first == last) { /* Empty. Disconnected graph! */
        if (nleft == 0 || drain)
          break;

        k = RandomInRange(nleft);
        for (i=0; i<nvtxs; i++) {
          if (touched[i] == 0) {
            if (k == 0)
              break;
            else
              k--;
          }
        }

        queue[0] = i;
        touched[i] = 1;
        first = 0; last = 1;;
        nleft--;
      }

      i = queue[first++];
      if (pwgts[1]-vwgt[i] < minpwgt[1]) {
        drain = 1;
        continue;
      }

      part[i] = 0;
      INC_DEC(pwgts[0], pwgts[1], vwgt[i]);
      if (pwgts[1] <= maxpwgt[1])
        break;

      drain = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (touched[k] == 0) {
          queue[last++] = k;
          touched[k] = 1;
          nleft--;
        }
      }
    }

/*
    if (pwgts[0] < minpwgt[0] || pwgts[1] < minpwgt[1] || pwgts[0] > maxpwgt[0] || pwgts[1] > maxpwgt[1])
      printf("Unbalanced I-partitions! %d %d [%d %d] [%d %d]\n", pwgts[0], pwgts[1], minpwgt[0], maxpwgt[0], minpwgt[1], maxpwgt[1]);
*/

    /*************************************************************
    * Do some partition refinement 
    **************************************************************/
    graph->where = part;
    Ser_EdgeFM(graph, minpwgt, maxpwgt, NLGR_PASSES);

    if (nbfs == 0 || bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, part, bestpart);
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestpart, part);

  GKfree(&bestpart, &queue, &touched, -1);
}



/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void Ser_EdgeFM(GraphType *graph, int *minpwgt, int *maxpwgt, int npasses)
{
  int i, ii, j, k, nvtxs, pwgts[2], from, to, icut, mincut, pass, me;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *where, *gain, *perm;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  gain = idxsmalloc(nvtxs, 0, "Ser_EdgeFM: gain");
  perm = idxsmalloc(nvtxs, 0, "Ser_EdgeFM: perm");
  for (i=0; i<nvtxs; i++)
    perm[i] = i;

  /* Compute initial gains */
  mincut = pwgts[0] = pwgts[1] = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    pwgts[me] += vwgt[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[adjncy[j]] == me)
        gain[i] -= adjwgt[j];
      else {
        gain[i] += adjwgt[j];
        mincut += adjwgt[j];
      }
    }
  }
  mincut = mincut/2;

  /* printf("Initial Cut: %d\n", mincut); */

  for (pass=0; pass<npasses; pass++) { /* Do a number of passes */
    /* printf("Cut: %d [%d %d %d]\n", mincut, zeropwgt, pwgts[0], pwgts[1]); */
    icut = mincut;
    RandomPermute(nvtxs, perm, 0);

    for (ii=0; ii<nvtxs; ii++) {
      i = perm[ii];
      from = where[i];
      to = (from+1)%2;
      if (gain[i] >= 0 && pwgts[from]-vwgt[i] >= minpwgt[from] && pwgts[to]+vwgt[i] <= maxpwgt[to]) {
        mincut -= gain[i];
        where[i] = to;
        INC_DEC(pwgts[to], pwgts[from], vwgt[i]);

        gain[i] = -gain[i];

        /* Update gains of adjacent vertices */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (where[k] == from) 
            gain[k] += 2*adjwgt[j];
          else 
            gain[k] -= 2*adjwgt[j];
        }
      }
    }

    if (icut == mincut)
      break;
  }

/*
  if (pwgts[0] < minpwgt[0] || pwgts[1] < minpwgt[1] || pwgts[0] > maxpwgt[0] || pwgts[1] > maxpwgt[1])
    printf("Unbalanced partitions! %d %d [%d %d] [%d %d]\n", pwgts[0], pwgts[1], minpwgt[0], maxpwgt[0], minpwgt[1], maxpwgt[1]);
*/
  graph->mincut = mincut;

  GKfree(&gain, &perm, -1);
}


/*************************************************************************
* This function takes a graph and a bisection and splits it into two graphs.
**************************************************************************/
void Ser_SplitGraphPart(GraphType *graph, idxtype *part, GraphType *lgraph, GraphType *rgraph)
{
  int i, j, k, l, mypart, nvtxs, snvtxs[2], snedges[2];
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *label;
  idxtype *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *slabel[2];
  idxtype *rename;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  label = graph->label;

  rename = idxmalloc(nvtxs, "Ser_SplitGraphPart: rename");
  
  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  for (i=0; i<nvtxs; i++) {
    k = part[i];
    rename[i] = snvtxs[k]++;
    snedges[k] += xadj[i+1]-xadj[i];
  }

  /* printf("%d %d -> [%d %d] [%d %d]\n", nvtxs, graph->nedges, snvtxs[0], snedges[0], snvtxs[1], snedges[1]); */

  sxadj[0] = lgraph->xadj = idxmalloc(snvtxs[0]+1, "Ser_SplitGraphPart: lgraph->xadj");
  svwgt[0] = lgraph->vwgt = idxmalloc(snvtxs[0], "Ser_SplitGraphPart: lgraph->vwgt");
  slabel[0] = lgraph->label = idxmalloc(snvtxs[0], "Ser_SplitGraphPart: lgraph->label");
  sadjncy[0] = lgraph->adjncy = idxmalloc(snedges[0], "Ser_SplitGraphPart: lgraph->adjncy");
  sadjwgt[0] = lgraph->adjwgt = idxmalloc(snedges[0], "Ser_SplitGraphPart: lgraph->adjwgt");
  sxadj[1] = rgraph->xadj = idxmalloc(snvtxs[1]+1, "Ser_SplitGraphPart: rgraph->xadj");
  svwgt[1] = rgraph->vwgt = idxmalloc(snvtxs[1], "Ser_SplitGraphPart: rgraph->vwgt");
  slabel[1] = rgraph->label = idxmalloc(snvtxs[1], "Ser_SplitGraphPart: rgraph->label");
  sadjncy[1] = rgraph->adjncy = idxmalloc(snedges[1], "Ser_SplitGraphPart: rgraph->adjncy");
  sadjwgt[1] = rgraph->adjwgt = idxmalloc(snedges[1], "Ser_SplitGraphPart: rgraph->adjwgt");

  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  sxadj[0][0] = sxadj[1][0] = 0;
  for (i=0; i<nvtxs; i++) {
    mypart = part[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (part[k] == mypart) {
        sadjncy[mypart][snedges[mypart]] = rename[k]; 
        sadjwgt[mypart][snedges[mypart]++] = adjwgt[j]; 
      }
    }

    svwgt[mypart][snvtxs[mypart]] = vwgt[i];
    slabel[mypart][snvtxs[mypart]] = label[i];
    sxadj[mypart][++snvtxs[mypart]] = snedges[mypart];
  }

  lgraph->nvtxs = snvtxs[0];
  lgraph->nedges = snedges[0];
  lgraph->level = -1;
  rgraph->nvtxs = snvtxs[1];
  rgraph->nedges = snedges[1];
  rgraph->level = -1;


  /* printf("%d %d -> [%d %d] [%d %d]\n", nvtxs, graph->nedges, snvtxs[0], snedges[0], snvtxs[1], snedges[1]); */

  free(rename);
}





/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void Ser_Refine(GraphType *orggraph, GraphType *graph, int nparts, float ubfactor, EdgeType *degrees)
{
  int i;

  for (i=0; ; i++) {
    Ser_ComputePartitionParams(graph, nparts, degrees);
    Ser_KWayRefine(graph, nparts, ubfactor, 6); 

    if (graph == orggraph)
      break;

    graph = graph->finer;
    Ser_ProjectPartition(graph);
  }
}


/*************************************************************************
* This function projects a partition.
**************************************************************************/
void Ser_ProjectPartition(GraphType *graph)
{
  int i, j, nvtxs;
  idxtype *match, *cmap, *where, *cwhere;
  GraphType *cgraph;


  cgraph = graph->coarser;
  cwhere = cgraph->where;

  nvtxs = graph->nvtxs;
  match = graph->match;
  cmap = graph->cmap;
  where = graph->where = idxmalloc(nvtxs, "Ser_ProjectPartition: graph->where");

  for (i=0; i<nvtxs; i++) 
    where[i] = cwhere[cmap[i]];

  FreeGraph(graph->coarser);
  graph->coarser = NULL;

}



/*************************************************************************
* This function computes the initial id/ed 
**************************************************************************/
void Ser_ComputePartitionParams(GraphType *graph, int nparts, EdgeType *degrees)
{
  int i, j, k, l, nvtxs;
  idxtype *xadj, *adjncy, *adjwgt, *lpwgts;
  idxtype *where;
  RInfoType *rinfo, *myrinfo;
  EdgeType *edegrees;
  int me, other;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  rinfo = graph->rinfo = (RInfoType *)GKmalloc(sizeof(RInfoType)*nvtxs, "ComputePartitionParams: rinfo");
  lpwgts = graph->lpwgts = idxsmalloc(nparts, 0, "ComputePartitionParams: lpwgts");


  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  graph->mincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    myrinfo = rinfo+i;

    lpwgts[me] += graph->vwgt[i];

    myrinfo->degrees = degrees+xadj[i];
    myrinfo->ndegrees = myrinfo->id = myrinfo->ed = 0;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]])
        myrinfo->id += adjwgt[j];
      else
        myrinfo->ed += adjwgt[j];
    }


    if (myrinfo->ed > 0) {  /* Time to do some serious work */
      graph->mincut += myrinfo->ed;
      edegrees = myrinfo->degrees;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (edegrees[k].edge == other) {
              edegrees[k].ewgt += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            edegrees[k].edge = other;
            edegrees[k].ewgt = adjwgt[j];
            myrinfo->ndegrees++;
          }
        }
      }
    }
  }
  graph->mincut = graph->mincut/2;
}




/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Ser_KWayRefine(GraphType *graph, int nparts, float ubfactor, int npasses)
{
  int i, ii, iii, j, jj, k, l, pass, nvtxs, nmoves; 
  int from, me, other, vwgt, minpwgt, maxpwgt, badminpwgt, badmaxpwgt, oldcut, mincut;
  idxtype *xadj, *adjncy, *adjwgt;
  idxtype *where, *lpwgts, *perm;
  RInfoType *rinfo, *myrinfo;
  EdgeType *degrees;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  rinfo = graph->rinfo;
  lpwgts = graph->lpwgts;

  minpwgt = lpwgts[idxamin(nparts, lpwgts)];
  maxpwgt = lpwgts[idxamax(nparts, lpwgts)];
  badminpwgt = (1.0/ubfactor)*idxsum(nparts, lpwgts)/nparts;
  badmaxpwgt = ubfactor*idxsum(nparts, lpwgts)/nparts;

/*
  printf("K-way refinement [%5d, %5d], [%5d, %5d] Cut: %5d, %5d\n", minpwgt, maxpwgt, badminpwgt, badmaxpwgt, graph->mincut, graph->nvtxs);
*/

  perm = idxmalloc(nvtxs, "Ser_KWayRefine: perm");
  for (i=0; i<nvtxs; i++)
    perm[i] = i;

  for (pass=0; pass<npasses; pass++) {
    oldcut = graph->mincut;
    minpwgt = lpwgts[idxamin(nparts, lpwgts)];
    maxpwgt = lpwgts[idxamax(nparts, lpwgts)];

    RandomPermute(nvtxs, perm, 0);
    for (nmoves=iii=0; iii<nvtxs; iii++) {
      i = perm[iii];

      if (rinfo[i].ed >= rinfo[i].id) { /* Total ED is too high */
        degrees = rinfo[i].degrees;
        from = where[i];
        vwgt = graph->vwgt[i];

        if (lpwgts[from]-vwgt < badminpwgt)
          continue;   /* This cannot be moved! */

        for (k=0; k<rinfo[i].ndegrees; k++) {
          if (lpwgts[degrees[k].edge]+vwgt <= badmaxpwgt)
            break;
        }

        if (k < rinfo[i].ndegrees) { /* You actually found one */
          for (j=k+1; j<rinfo[i].ndegrees; j++) {
            if ((degrees[j].ewgt > degrees[k].ewgt && lpwgts[degrees[j].edge]+vwgt <= badmaxpwgt) ||
                (degrees[j].ewgt == degrees[k].ewgt && lpwgts[degrees[j].edge] < lpwgts[degrees[k].edge]))
              k = j;
          }

          other = degrees[k].edge;

          if (degrees[k].ewgt > rinfo[i].id || (degrees[k].ewgt == rinfo[i].id && (lpwgts[from] - lpwgts[other] >= vwgt))) {
            graph->mincut -= (degrees[k].ewgt-rinfo[i].id);

            /* Update where, weight, and ID/ED information of the vertex you moved */
            where[i] = other;
            INC_DEC(lpwgts[other], lpwgts[from], vwgt);
            SWAP(rinfo[i].id, degrees[k].ewgt, j);
            if (degrees[k].ewgt == 0) 
              degrees[k] = degrees[--(rinfo[i].ndegrees)];
            else
              degrees[k].edge = from;


            /* Put the rinfo of adjacent vertices */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              ii = adjncy[j];
              me = where[ii];
              myrinfo = rinfo+ii;

              myrinfo->ndegrees = myrinfo->id = myrinfo->ed = 0;
              degrees = myrinfo->degrees;

              for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
                other = where[adjncy[jj]];
                if (me != other) {
                  myrinfo->ed += adjwgt[jj];

                  for (k=0; k<myrinfo->ndegrees; k++) {
                    if (degrees[k].edge == other) {
                      degrees[k].ewgt += adjwgt[jj];
                      break;
                    }
                  }
                  if (k == myrinfo->ndegrees) {
                    degrees[k].edge = other;
                    degrees[k].ewgt = adjwgt[jj];
                    myrinfo->ndegrees++;
                  }
                }
                else {
                  myrinfo->id += adjwgt[jj];
                }
              }
            }
            nmoves++;
          }
        }
      }
    }

    /*
    printf("\t[%5d, %5d], [%5d, %5d] Cut: %d\n", minpwgt, maxpwgt, badminpwgt, badmaxpwgt, graph->mincut);
    */

    if (graph->mincut == oldcut)
      break;
  }

  free(perm);

}


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ccgraph.c
 *
 * This file contains the functions that create the coarse graph
 *
 * Started 8/11/97
 * George
 *
 * $Id: ccgraph.c,v 1.1 1997/11/04 23:19:07 karypis Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraph(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
  int i, j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, cnedges, v, u, mask;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cadjncy, *cadjwgt, *cadjwgtsum;
  GraphType *cgraph;

  mask = (1<<11)-1;
  if (cnvtxs < 8*mask || graph->nedges/graph->nvtxs > 15) { 
    CreateCoarseGraphNoMask(ctrl, graph, cnvtxs, match, perm);
    return;
  }

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;


  /* Initialize the coarser graph */
  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;

  cgraph->finer = graph;
  graph->coarser = cgraph;


  /*************************************************************
  * Create the coarser graph
  **************************************************************/
  /* Allocate memory for the coarser graph, and fire up coarsening */
  cgraph->gdata = idxmalloc(4*cnvtxs+1 + 2*graph->nedges, "CreateCoarserGraph: gdata");
  cxadj = cgraph->xadj 			= cgraph->gdata;
  cvwgt = cgraph->vwgt 			= cgraph->gdata + cnvtxs+1;
  cadjwgtsum = cgraph->adjwgtsum 	= cgraph->gdata + 2*cnvtxs+1;
  cgraph->cmap 				= cgraph->gdata + 3*cnvtxs+1;
  cadjncy = cgraph->adjncy 		= cgraph->gdata + 4*cnvtxs+1;
  cadjwgt = cgraph->adjwgt 		= cgraph->gdata + 4*cnvtxs+1 + graph->nedges;

  iend = xadj[nvtxs];
  auxadj = (idxtype *)ctrl->wspace.edegrees;
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  htable = idxset(mask+1, -1, idxwspacemalloc(ctrl, mask+1)); 

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    cvwgt[cnvtxs] = vwgt[v];
    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m] += adjwgt[j];
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj] += adjwgt[j];
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges] = k;
          cadjwgt[nedges++] = adjwgt[j];
        }
      }
    }

    if (v != u) { 
      cvwgt[cnvtxs] += vwgt[u];
      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[kk] = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m] += adjwgt[j];
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj] += adjwgt[j];
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges] = k;
            cadjwgt[nedges++] = adjwgt[j];
          }
        }
      }

      /* Remove the contracted adjacency weight */
      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs) 
            break;
        }
      }
      if (jj >= 0 && cadjncy[jj] == cnvtxs) { /* This 2nd check is needed for non-adjacent matchings */
        cadjwgtsum[cnvtxs] -= cadjwgt[jj];
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      }
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d %d %d %d\n", cnvtxs, cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt), adjwgtsum[u], adjwgtsum[v]));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* Zero out the htable */
    htable[cnvtxs&mask] = -1;

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  /* If significant savings, readjust the amount of memory that you have allocated */
  if (cnedges > 100000 && cnedges < 0.7*graph->nedges) {
    idxcopy(cnedges, cgraph->adjwgt, cgraph->adjncy+cnedges);
    cgraph->gdata = realloc(cgraph->gdata, (4*cnvtxs+1 + 2*cnedges)*sizeof(idxtype));

    /* Do this, in case everything was copied into new space */
    cgraph->xadj 	= cgraph->gdata;
    cgraph->vwgt 	= cgraph->gdata + cnvtxs+1;
    cgraph->adjwgtsum 	= cgraph->gdata + 2*cnvtxs+1;
    cgraph->cmap 	= cgraph->gdata + 3*cnvtxs+1;
    cgraph->adjncy 	= cgraph->gdata + 4*cnvtxs+1;
    cgraph->adjwgt 	= cgraph->gdata + 4*cnvtxs+1 + cnedges;
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, mask+1);

}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraphNoMask(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
  int i, j, k, m, istart, iend, nvtxs, nedges, cnedges, v, u;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cadjncy, *cadjwgt, *cadjwgtsum;
  GraphType *cgraph;


  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;


  /* Initialize the coarser graph */
  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;

  cgraph->finer = graph;
  graph->coarser = cgraph;


  /*************************************************************
  * Create the coarser graph
  **************************************************************/
  /* Allocate memory for the coarser graph, and fire up coarsening */
  cgraph->gdata = idxmalloc(4*cnvtxs+3 + 2*graph->nedges, "CreateCoarserGraph: gdata");
  cxadj = cgraph->xadj 			= cgraph->gdata;
  cvwgt = cgraph->vwgt 			= cgraph->gdata + cnvtxs+1;
  cadjwgtsum = cgraph->adjwgtsum 	= cgraph->gdata + 2*cnvtxs+1;
  cgraph->cmap 				= cgraph->gdata + 3*cnvtxs+1;
  cadjncy = cgraph->adjncy 		= cgraph->gdata + 4*cnvtxs+1;
  cadjwgt = cgraph->adjwgt 		= cgraph->gdata + 4*cnvtxs+1 + graph->nedges + 1;

  htable = idxset(cnvtxs, -1, idxwspacemalloc(ctrl, cnvtxs));

  iend = xadj[nvtxs];
  auxadj = (idxtype *)ctrl->wspace.edegrees;
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    cvwgt[cnvtxs] = vwgt[v];
    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      if ((m = htable[k]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[k] = nedges++;
      }
      else {
        cadjwgt[m] += adjwgt[j];
      }
    }

    if (v != u) { 
      cvwgt[cnvtxs] += vwgt[u];
      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        if ((m = htable[k]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[k] = nedges++;
        }
        else {
          cadjwgt[m] += adjwgt[j];
        }
      }

      /* Remove the contracted adjacency weight */
      if ((j = htable[cnvtxs]) != -1) {
        ASSERT(cadjncy[j] == cnvtxs);
        cadjwgtsum[cnvtxs] -= cadjwgt[j];
        cadjncy[j] = cadjncy[--nedges];
        cadjwgt[j] = cadjwgt[nedges];
        htable[cnvtxs] = -1;
      }
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d\n", cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt)));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]] = -1;  /* Zero out the htable */

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  /* If significant savings, readjust the amount of memory that you have allocated */
  if (cnedges > 100000 && cnedges < 0.7*graph->nedges) {
    idxcopy(cnedges, cgraph->adjwgt, cgraph->adjncy+cnedges);
    cgraph->gdata = realloc(cgraph->gdata, (4*cnvtxs+1 + 2*cnedges)*sizeof(idxtype));

    /* Do this, in case everything was copied into new space */
    cgraph->xadj 	= cgraph->gdata;
    cgraph->vwgt 	= cgraph->gdata + cnvtxs+1;
    cgraph->adjwgtsum 	= cgraph->gdata + 2*cnvtxs+1;
    cgraph->cmap 	= cgraph->gdata + 3*cnvtxs+1;
    cgraph->adjncy 	= cgraph->gdata + 4*cnvtxs+1;
    cgraph->adjwgt 	= cgraph->gdata + 4*cnvtxs+1 + cnedges;
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, cnvtxs);
}



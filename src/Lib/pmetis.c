/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * pmetis.c
 *
 * This file contains the top level routines for the multilevel recursive
 * bisection algorithm PMETIS.
 *
 * Started 7/24/97
 * George
 *
 * $Id: pmetis.c,v 1.5 1998/01/28 17:23:06 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function is the entry point for PMETIS
**************************************************************************/
void METIS_PartGraphRecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                              idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, 
                              int *options, int *edgecut, idxtype *part)
{
  int i;
  float *tpwgts;

  tpwgts = fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++) 
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, 
                            tpwgts, options, edgecut, part);

  free(tpwgts);
}



/*************************************************************************
* This function is the entry point for PWMETIS that accepts exact weights
* for the target partitions
**************************************************************************/
void METIS_WPartGraphRecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                               idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, 
                               float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  int i, j;
  GraphType graph;
  CtrlType ctrl;
  float *mytpwgts;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_PMETIS, *nvtxs, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = PMETIS_CTYPE;
    ctrl.IType = PMETIS_ITYPE;
    ctrl.RType = PMETIS_RTYPE;
    ctrl.dbglvl = PMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_PMETIS;
  ctrl.CoarsenTo = 20;
  ctrl.maxvwgt = 1.5*(idxsum(*nvtxs, graph.vwgt)/ctrl.CoarsenTo);

  mytpwgts = fmalloc(*nparts, "PWMETIS: mytpwgts");
  for (i=0; i<*nparts; i++) 
    mytpwgts[i] = tpwgts[i];

  InitRandom();

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MlevelRecursiveBisection(&ctrl, &graph, *nparts, part, mytpwgts, 1.000, 0);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);
  free(mytpwgts);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}



/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
int MlevelRecursiveBisection(CtrlType *ctrl, GraphType *graph, int nparts, idxtype *part, float *tpwgts, float ubfactor, int fpart)
{
  int i, j, nvtxs, cut, tvwgt, tpwgts2[2];
  idxtype *label, *where;
  GraphType lgraph, rgraph;
  float wsum;

  nvtxs = graph->nvtxs;
  if (nvtxs == 0) {
    printf("\t***Cannot bisect a graph with 0 vertices!\n\t***You are trying to partition a graph into too many parts!\n");
    return 0;
  }

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt*ssum(nparts/2, tpwgts);
  tpwgts2[1] = tvwgt-tpwgts2[0];

  /* printf("%d %d\n", tpwgts2[0], tpwgts2[1]);*/

  MlevelEdgeBisection(ctrl, graph, tpwgts2, ubfactor);
  cut = graph->mincut;

  label = graph->label;
  where = graph->where;
  for (i=0; i<nvtxs; i++)
    part[label[i]] = where[i] + fpart;

  if (nparts > 2) {
    SplitGraphPart(ctrl, graph, &lgraph, &rgraph);
    /* printf("%d %d\n", lgraph.nvtxs, rgraph.nvtxs); */
  }


  /* Free the memory of the top level graph */
  GKfree(&graph->gdata, &graph->rdata, &graph->label, -1);

  /* Scale the fractions in the tpwgts according to the true weight */
  wsum = ssum(nparts/2, tpwgts);
  sscale(nparts/2, 1.0/wsum, tpwgts);
  sscale(nparts-nparts/2, 1.0/(1.0-wsum), tpwgts+nparts/2);
  /*
  for (i=0; i<nparts; i++)
    printf("%5.3f ", tpwgts[i]);
  printf("[%5.3f]\n", wsum);
  */

  /* Do the recursive call */
  if (nparts > 3) {
    cut += MlevelRecursiveBisection(ctrl, &lgraph, nparts/2, part, tpwgts, ubfactor, fpart);
    cut += MlevelRecursiveBisection(ctrl, &rgraph, nparts-nparts/2, part, tpwgts+nparts/2, ubfactor, fpart+nparts/2);
  }
  else if (nparts == 3) {
    cut += MlevelRecursiveBisection(ctrl, &rgraph, nparts-nparts/2, part, tpwgts+nparts/2, ubfactor, fpart+nparts/2);
    GKfree(&lgraph.gdata, &lgraph.label, -1);
  }

  return cut;

}


/*************************************************************************
* This function performs multilevel bisection
**************************************************************************/
void MlevelEdgeBisection(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  GraphType *cgraph;

  cgraph = Coarsen2Way(ctrl, graph);

  Init2WayPartition(ctrl, cgraph, tpwgts, ubfactor);

  Refine2Way(ctrl, graph, cgraph, tpwgts, ubfactor);

/*
  IsConnectedSubdomain(ctrl, graph, 0);
  IsConnectedSubdomain(ctrl, graph, 1);
*/
}




/*************************************************************************
* This function takes a graph and a bisection and splits it into two graphs.
**************************************************************************/
void SplitGraphPart(CtrlType *ctrl, GraphType *graph, GraphType *lgraph, GraphType *rgraph)
{
  int i, j, k, l, istart, iend, mypart, nvtxs, snvtxs[2], snedges[2], sum;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum, *label, *where, *bndptr;
  idxtype *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *sadjwgtsum[2], *slabel[2];
  idxtype *rename;
  idxtype *auxadjncy, *auxadjwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SplitTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  label = graph->label;
  where = graph->where;
  bndptr = graph->bndptr;
  ASSERT(bndptr != NULL);

  rename = idxwspacemalloc(ctrl, nvtxs);
  
  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  for (i=0; i<nvtxs; i++) {
    k = where[i];
    rename[i] = snvtxs[k]++;
    snedges[k] += xadj[i+1]-xadj[i];
  }

  InitGraph(lgraph);
  InitGraph(rgraph);

  /* Allocate memory and set the pointers accordingly */
  lgraph->gdata = idxmalloc(4*snvtxs[0]+1+2*snedges[0], "SplitGraphPart: lgdata");
  sxadj[0] = lgraph->xadj 		= lgraph->gdata;
  svwgt[0] = lgraph->vwgt 		= lgraph->gdata + snvtxs[0]+1;
  sadjwgtsum[0] = lgraph->adjwgtsum 	= lgraph->gdata + 2*snvtxs[0]+1;
  lgraph->cmap 				= lgraph->gdata + 3*snvtxs[0]+1;
  sadjncy[0] = lgraph->adjncy 		= lgraph->gdata + 4*snvtxs[0]+1;
  sadjwgt[0] = lgraph->adjwgt 		= lgraph->gdata + 4*snvtxs[0]+1 + snedges[0];
  slabel[0] = lgraph->label = idxmalloc(snvtxs[0], "SplitGraphPart: lgraph->label");

  rgraph->gdata = idxmalloc(4*snvtxs[1]+1+2*snedges[1], "SplitGraphPart: rgdata");
  sxadj[1] = rgraph->xadj 		= rgraph->gdata;
  svwgt[1] = rgraph->vwgt 		= rgraph->gdata + snvtxs[1]+1;
  sadjwgtsum[1] = rgraph->adjwgtsum 	= rgraph->gdata + 2*snvtxs[1]+1;
  rgraph->cmap 				= rgraph->gdata + 3*snvtxs[1]+1;
  sadjncy[1] = rgraph->adjncy 		= rgraph->gdata + 4*snvtxs[1]+1;
  sadjwgt[1] = rgraph->adjwgt 		= rgraph->gdata + 4*snvtxs[1]+1 + snedges[1];
  slabel[1] = rgraph->label = idxmalloc(snvtxs[1], "Ser_SplitGraphPart: rgraph->label");


  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  sxadj[0][0] = sxadj[1][0] = 0;
  for (i=0; i<nvtxs; i++) {
    mypart = where[i];
    sum = adjwgtsum[i];

    istart = xadj[i];
    iend = xadj[i+1];
    if (bndptr[i] == -1) { /* This is an interior vertex */
      auxadjncy = sadjncy[mypart] + snedges[mypart] - istart;
      auxadjwgt = sadjwgt[mypart] + snedges[mypart] - istart;
      for(j=istart; j<iend; j++) {
        auxadjncy[j] = adjncy[j];
        auxadjwgt[j] = adjwgt[j]; 
      }
      snedges[mypart] += iend-istart;
    }
    else {
      auxadjncy = sadjncy[mypart];
      auxadjwgt = sadjwgt[mypart];
      l = snedges[mypart];
      for (j=istart; j<iend; j++) {
        k = adjncy[j];
        if (where[k] == mypart) {
          auxadjncy[l] = k;
          auxadjwgt[l++] = adjwgt[j]; 
        }
        else {
          sum -= adjwgt[j];
        }
      }
      snedges[mypart] = l;
    }

    svwgt[mypart][snvtxs[mypart]] = vwgt[i];
    sadjwgtsum[mypart][snvtxs[mypart]] = sum;
    slabel[mypart][snvtxs[mypart]] = label[i];
    sxadj[mypart][++snvtxs[mypart]] = snedges[mypart];
  }

  for (mypart=0; mypart<2; mypart++) {
    iend = sxadj[mypart][snvtxs[mypart]];
    auxadjncy = sadjncy[mypart];
    for (i=0; i<iend; i++) 
      auxadjncy[i] = rename[auxadjncy[i]];
  }

  lgraph->nvtxs = snvtxs[0];
  lgraph->nedges = snedges[0];
  rgraph->nvtxs = snvtxs[1];
  rgraph->nedges = snedges[1];

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SplitTmr));

  idxwspacefree(ctrl, nvtxs);
}



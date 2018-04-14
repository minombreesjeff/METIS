/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ometis.c
 *
 * This file contains the top level routines for the multilevel recursive
 * bisection algorithm PMETIS.
 *
 * Started 7/24/97
 * George
 *
 * $Id: ometis.c,v 1.5 1997/12/01 16:53:03 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function is the entry point for OEMETIS
**************************************************************************/
void METIS_EdgeND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, 
                  idxtype *perm, idxtype *iperm) 
{
  int i, j;
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_OEMETIS, *nvtxs, xadj, adjncy, NULL, NULL, 0);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = OEMETIS_CTYPE;
    ctrl.IType = OEMETIS_ITYPE;
    ctrl.RType = OEMETIS_RTYPE;
    ctrl.dbglvl = OEMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_OEMETIS;
  ctrl.CoarsenTo = 20;
  ctrl.maxvwgt = 1.5*(idxsum(*nvtxs, graph.vwgt)/ctrl.CoarsenTo);

  InitRandom();

  AllocateWorkSpace(&ctrl, &graph, 2);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, *nvtxs);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);
}


/*************************************************************************
* This function is the entry point for ONCMETIS
**************************************************************************/
void METIS_NodeND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, 
                  idxtype *perm, idxtype *iperm) 
{
  int i, ii, j, l, wflag, nflag;
  GraphType graph;
  CtrlType ctrl;
  idxtype *cptr, *cind, *piperm;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType   = ONMETIS_CTYPE;
    ctrl.IType   = ONMETIS_ITYPE;
    ctrl.RType   = ONMETIS_RTYPE;
    ctrl.dbglvl  = ONMETIS_DBGLVL;
    ctrl.oflags  = ONMETIS_OFLAGS;
    ctrl.pfactor = ONMETIS_PFACTOR;
    ctrl.nseps   = ONMETIS_NSEPS;
  }
  else {
    ctrl.CType   = options[OPTION_CTYPE];
    ctrl.IType   = options[OPTION_ITYPE];
    ctrl.RType   = options[OPTION_RTYPE];
    ctrl.dbglvl  = options[OPTION_DBGLVL];
    ctrl.oflags  = options[OPTION_OFLAGS];
    ctrl.pfactor = options[OPTION_PFACTOR];
    ctrl.nseps   = options[OPTION_NSEPS];
  }
  if (ctrl.nseps < 1)
    ctrl.nseps = 1;

  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  InitRandom();

  if (ctrl.pfactor > 0) { 
    /*============================================================
    * Prune the dense columns
    ==============================================================*/
    piperm = idxmalloc(*nvtxs, "ONMETIS: piperm");

    PruneGraph(&ctrl, &graph, *nvtxs, xadj, adjncy, piperm, (float)(0.1*ctrl.pfactor));
  }
  else if (ctrl.oflags&OFLAG_COMPRESS) {
    /*============================================================
    * Compress the graph 
    ==============================================================*/
    cptr = idxmalloc(*nvtxs+1, "ONMETIS: cptr");
    cind = idxmalloc(*nvtxs, "ONMETIS: cind");

    CompressGraph(&ctrl, &graph, *nvtxs, xadj, adjncy, cptr, cind);

    if (graph.nvtxs >= COMPRESSION_FRACTION*(*nvtxs)) {
      ctrl.oflags--; /* We actually performed no compression */
      GKfree(&cptr, &cind, -1);
    }
    else if (2*graph.nvtxs < *nvtxs && ctrl.nseps == 1)
      ctrl.nseps = 2;
  }
  else {
    SetUpGraph(&graph, OP_ONMETIS, *nvtxs, xadj, adjncy, NULL, NULL, 0);
  }


  /*=============================================================
  * Do the nested dissection ordering 
  --=============================================================*/
  ctrl.maxvwgt = 1.5*(idxsum(graph.nvtxs, graph.vwgt)/ctrl.CoarsenTo);
  AllocateWorkSpace(&ctrl, &graph, 2);

  if (ctrl.oflags&OFLAG_CCMP) 
    MlevelNestedDissectionCC(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, graph.nvtxs);
  else
    MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, graph.nvtxs);

  FreeWorkSpace(&ctrl, &graph);

  if (ctrl.pfactor > 0) { /* Order any prunned vertices */
    if (graph.nvtxs < *nvtxs) { 
      idxcopy(graph.nvtxs, iperm, perm);  /* Use perm as an auxiliary array */
      for (i=0; i<graph.nvtxs; i++)
        iperm[piperm[i]] = perm[i];
      for (i=graph.nvtxs; i<*nvtxs; i++)
        iperm[piperm[i]] = i;
    }

    GKfree(&piperm, -1);
  }
  else if (ctrl.oflags&OFLAG_COMPRESS) { /* Uncompress the ordering */
    if (graph.nvtxs < COMPRESSION_FRACTION*(*nvtxs)) { 
      /* construct perm from iperm */
      for (i=0; i<graph.nvtxs; i++)
        perm[iperm[i]] = i; 
      for (l=ii=0; ii<graph.nvtxs; ii++) {
        i = perm[ii];
        for (j=cptr[i]; j<cptr[i+1]; j++)
          iperm[cind[j]] = l++;
      }
    }

    GKfree(&cptr, &cind, -1);
  }


  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);

}


/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void METIS_NodeWND(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, 
                   int *options, idxtype *perm, idxtype *iperm) 
{
  int i, j, tvwgt;
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, xadj, adjncy, vwgt, NULL, 2);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = ONMETIS_CTYPE;
    ctrl.IType = ONMETIS_ITYPE;
    ctrl.RType = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags  = OFLAG_COMPRESS;
  ctrl.pfactor = 0;
  ctrl.nseps = 2;
  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;
  ctrl.maxvwgt = 1.5*(idxsum(*nvtxs, graph.vwgt)/ctrl.CoarsenTo);

  InitRandom();

  AllocateWorkSpace(&ctrl, &graph, 2);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, *nvtxs);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);
}




/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void MlevelNestedDissection(CtrlType *ctrl, GraphType *graph, idxtype *order, float ubfactor, int lastvtx)
{
  int i, j, nvtxs, nbnd, tvwgt, tpwgts2[2];
  idxtype *label, *bndind;
  GraphType lgraph, rgraph;

  nvtxs = graph->nvtxs;

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];

  switch (ctrl->optype) {
    case OP_OEMETIS:
      MlevelEdgeBisection(ctrl, graph, tpwgts2, ubfactor);

      IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SepTmr));
      ConstructMinCoverSeparator(ctrl, graph, ubfactor);
      IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SepTmr));

      break;
    case OP_ONMETIS:
      MlevelNodeBisectionMultiple(ctrl, graph, tpwgts2, ubfactor);

      IFSET(ctrl->dbglvl, DBG_SEPINFO, printf("Nvtxs: %6d, [%6d %6d %6d]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

      break;
  }

  /* Order the nodes in the separator */
  nbnd = graph->nbnd;
  bndind = graph->bndind;
  label = graph->label;
  for (i=0; i<nbnd; i++) 
    order[label[bndind[i]]] = --lastvtx;

  SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  GKfree(&graph->gdata, &graph->rdata, &graph->label, -1);

  if (rgraph.nvtxs > MMDSWITCH) 
    MlevelNestedDissection(ctrl, &rgraph, order, ubfactor, lastvtx);
  else {
    MMDOrder(ctrl, &rgraph, order, lastvtx); 
    GKfree(&rgraph.gdata, &rgraph.rdata, &rgraph.label, -1);
  }
  if (lgraph.nvtxs > MMDSWITCH) 
    MlevelNestedDissection(ctrl, &lgraph, order, ubfactor, lastvtx-rgraph.nvtxs);
  else {
    MMDOrder(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs); 
    GKfree(&lgraph.gdata, &lgraph.rdata, &lgraph.label, -1);
  }
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void MlevelNestedDissectionCC(CtrlType *ctrl, GraphType *graph, idxtype *order, float ubfactor, int lastvtx)
{
  int i, j, nvtxs, nbnd, tvwgt, tpwgts2[2], nsgraphs, ncmps, rnvtxs;
  idxtype *label, *bndind;
  idxtype *cptr, *cind;
  GraphType *sgraphs;

  nvtxs = graph->nvtxs;

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];

  MlevelNodeBisectionMultiple(ctrl, graph, tpwgts2, ubfactor);
  IFSET(ctrl->dbglvl, DBG_SEPINFO, printf("Nvtxs: %6d, [%6d %6d %6d]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

  /* Order the nodes in the separator */
  nbnd = graph->nbnd;
  bndind = graph->bndind;
  label = graph->label;
  for (i=0; i<nbnd; i++) 
    order[label[bndind[i]]] = --lastvtx;

  cptr = idxmalloc(nvtxs, "MlevelNestedDissectionCC: cptr");
  cind = idxmalloc(nvtxs, "MlevelNestedDissectionCC: cind");
  ncmps = FindComponents(ctrl, graph, cptr, cind);

/*
  if (ncmps > 2)
    printf("[%5d] has %3d components\n", nvtxs, ncmps);
*/

  sgraphs = (GraphType *)GKmalloc(ncmps*sizeof(GraphType), "MlevelNestedDissectionCC: sgraphs");

  nsgraphs = SplitGraphOrderCC(ctrl, graph, sgraphs, ncmps, cptr, cind);

  GKfree(&cptr, &cind, -1);

  /* Free the memory of the top level graph */
  GKfree(&graph->gdata, &graph->rdata, &graph->label, -1);

  /* Go and process the subgraphs */
  for (rnvtxs=i=0; i<nsgraphs; i++) {
    if (sgraphs[i].adjwgt == NULL) {
      MMDOrder(ctrl, sgraphs+i, order, lastvtx-rnvtxs);
      GKfree(&sgraphs[i].gdata, &sgraphs[i].label, -1);
    }
    else {
      MlevelNestedDissectionCC(ctrl, sgraphs+i, order, ubfactor, lastvtx-rnvtxs);
    }
    rnvtxs += sgraphs[i].nvtxs;
  }

  free(sgraphs);
}



/*************************************************************************
* This function performs multilevel bisection. It performs multiple 
* bisections and selects the best.
**************************************************************************/
void MlevelNodeBisectionMultiple(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  int i, nvtxs, cnvtxs, mincut, tmp;
  GraphType *cgraph; 
  idxtype *bestwhere;

  if (ctrl->nseps == 1 || graph->nvtxs < (ctrl->oflags&OFLAG_COMPRESS ? 1000 : 2000)) {
    MlevelNodeBisection(ctrl, graph, tpwgts, ubfactor);
    return;
  }

  nvtxs = graph->nvtxs;

  if (ctrl->oflags&OFLAG_COMPRESS) { /* Multiple separators at the original graph */
    bestwhere = idxmalloc(nvtxs, "MlevelNodeBisection2: bestwhere");
    mincut = nvtxs;

    for (i=ctrl->nseps; i>0; i--) {
      MlevelNodeBisection(ctrl, graph, tpwgts, ubfactor);

      /* printf("%5d ", cgraph->mincut); */

      if (graph->mincut < mincut) {
        mincut = graph->mincut;
        idxcopy(nvtxs, graph->where, bestwhere);
      }

      GKfree(&graph->rdata, -1);
    
      if (mincut == 0)
        break;
    }
    /* printf("[%5d]\n", mincut); */

    Allocate2WayNodePartitionMemory(ctrl, graph);
    idxcopy(nvtxs, bestwhere, graph->where);
    free(bestwhere);

    Compute2WayNodePartitionParams(ctrl, graph);
  }
  else {  /* Coarsen it a bit */
    ctrl->CoarsenTo = nvtxs-1;

    cgraph = Coarsen2Way(ctrl, graph);

    cnvtxs = cgraph->nvtxs;

    bestwhere = idxmalloc(cnvtxs, "MlevelNodeBisection2: bestwhere");
    mincut = nvtxs;

    for (i=ctrl->nseps; i>0; i--) {
      ctrl->CType += 20; /* This is a hack. Look at coarsen.c */
      MlevelNodeBisection(ctrl, cgraph, tpwgts, ubfactor);

      /* printf("%5d ", cgraph->mincut); */

      if (cgraph->mincut < mincut) {
        mincut = cgraph->mincut;
        idxcopy(cnvtxs, cgraph->where, bestwhere);
      }

      GKfree(&cgraph->rdata, -1);
    
      if (mincut == 0)
        break;
    }
    /* printf("[%5d]\n", mincut); */

    Allocate2WayNodePartitionMemory(ctrl, cgraph);
    idxcopy(cnvtxs, bestwhere, cgraph->where);
    free(bestwhere);

    Compute2WayNodePartitionParams(ctrl, cgraph);

    Refine2WayNode(ctrl, graph, cgraph, ubfactor);
  }

}

/*************************************************************************
* This function performs multilevel bisection
**************************************************************************/
void MlevelNodeBisection(CtrlType *ctrl, GraphType *graph, int *tpwgts, float ubfactor)
{
  GraphType *cgraph;

  ctrl->CoarsenTo = graph->nvtxs/8;
  if (ctrl->CoarsenTo > 100)
    ctrl->CoarsenTo = 100;
  else if (ctrl->CoarsenTo < 40)
    ctrl->CoarsenTo = 40;
  ctrl->maxvwgt = 1.5*((tpwgts[0]+tpwgts[1])/ctrl->CoarsenTo);

  cgraph = Coarsen2Way(ctrl, graph);

  switch (ctrl->IType) {
    case IPART_GGPKL:
      Init2WayPartition(ctrl, cgraph, tpwgts, ubfactor);

      IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SepTmr));

      Compute2WayPartitionParams(ctrl, cgraph);
      ConstructSeparator(ctrl, cgraph, ubfactor);

      IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SepTmr));
      break;
    case IPART_GGPKLNODE:
      InitSeparator(ctrl, cgraph, ubfactor);
      break;
  }

  Refine2WayNode(ctrl, graph, cgraph, ubfactor);

}




/*************************************************************************
* This function takes a graph and a bisection and splits it into two graphs.
* This function relies on the fact that adjwgt is all equal to 1.
**************************************************************************/
void SplitGraphOrder(CtrlType *ctrl, GraphType *graph, GraphType *lgraph, GraphType *rgraph)
{
  int i, ii, j, k, l, istart, iend, mypart, nvtxs, snvtxs[3], snedges[3];
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum, *label, *where, *bndptr, *bndind;
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
  bndind = graph->bndind;
  ASSERT(bndptr != NULL);

  rename = idxwspacemalloc(ctrl, nvtxs);
  
  snvtxs[0] = snvtxs[1] = snvtxs[2] = snedges[0] = snedges[1] = snedges[2] = 0;
  for (i=0; i<nvtxs; i++) {
    k = where[i];
    rename[i] = snvtxs[k]++;
    snedges[k] += xadj[i+1]-xadj[i];
  }

  InitGraph(lgraph);
  InitGraph(rgraph);

  /* Allocate memory and set the pointers accordingly */
  lgraph->gdata = idxmalloc(4*snvtxs[0]+1+2*snedges[0], "SplitGraphOrder: lgdata");
  sxadj[0] = lgraph->xadj               = lgraph->gdata;
  svwgt[0] = lgraph->vwgt               = lgraph->gdata + snvtxs[0]+1;
  sadjwgtsum[0] = lgraph->adjwgtsum     = lgraph->gdata + 2*snvtxs[0]+1;
  lgraph->cmap                          = lgraph->gdata + 3*snvtxs[0]+1;
  sadjncy[0] = lgraph->adjncy           = lgraph->gdata + 4*snvtxs[0]+1;
  sadjwgt[0] = lgraph->adjwgt           = lgraph->gdata + 4*snvtxs[0]+1 + snedges[0];
  slabel[0] = lgraph->label = idxmalloc(snvtxs[0], "SplitGraphOrder: lgraph->label");

  rgraph->gdata = idxmalloc(4*snvtxs[1]+1+2*snedges[1], "SplitGraphOrder: rgdata");
  sxadj[1] = rgraph->xadj               = rgraph->gdata;
  svwgt[1] = rgraph->vwgt               = rgraph->gdata + snvtxs[1]+1;
  sadjwgtsum[1] = rgraph->adjwgtsum     = rgraph->gdata + 2*snvtxs[1]+1;
  rgraph->cmap                          = rgraph->gdata + 3*snvtxs[1]+1;
  sadjncy[1] = rgraph->adjncy           = rgraph->gdata + 4*snvtxs[1]+1;
  sadjwgt[1] = rgraph->adjwgt           = rgraph->gdata + 4*snvtxs[1]+1 + snedges[1];
  slabel[1] = rgraph->label = idxmalloc(snvtxs[1], "SplitGraphOrder: rgraph->label");


  /* Go and use bndptr to also mark the boundary nodes in the two partitions */
  for (ii=0; ii<graph->nbnd; ii++) {
    i = bndind[ii];
    for (j=xadj[i]; j<xadj[i+1]; j++)
      bndptr[adjncy[j]] = 1;
  }

  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  sxadj[0][0] = sxadj[1][0] = 0;
  for (i=0; i<nvtxs; i++) {
    if ((mypart = where[i]) == 2)
      continue;

    istart = xadj[i];
    iend = xadj[i+1];
    if (bndptr[i] == -1) { /* This is an interior vertex */
      auxadjncy = sadjncy[mypart] + snedges[mypart] - istart;
      for(j=istart; j<iend; j++) 
        auxadjncy[j] = adjncy[j];
      snedges[mypart] += iend-istart;
    }
    else {
      auxadjncy = sadjncy[mypart];
      l = snedges[mypart];
      for (j=istart; j<iend; j++) {
        k = adjncy[j];
        if (where[k] == mypart) 
          auxadjncy[l++] = k;
      }
      snedges[mypart] = l;
    }

    svwgt[mypart][snvtxs[mypart]] = vwgt[i];
    sadjwgtsum[mypart][snvtxs[mypart]] = snedges[mypart]-sxadj[mypart][snvtxs[mypart]];
    slabel[mypart][snvtxs[mypart]] = label[i];
    sxadj[mypart][++snvtxs[mypart]] = snedges[mypart];
  }

  for (mypart=0; mypart<2; mypart++) {
    iend = snedges[mypart];
    idxset(iend, 1, sadjwgt[mypart]);

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

/*************************************************************************
* This function uses MMD to order the graph. The vertices are numbered
* from lastvtx downwards
**************************************************************************/
void MMDOrder(CtrlType *ctrl, GraphType *graph, idxtype *order, int lastvtx)
{
  int i, j, k, nvtxs, nofsub, firstvtx;
  idxtype *xadj, *adjncy, *label;
  idxtype *perm, *iperm, *head, *qsize, *list, *marker;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  /* Relabel the vertices so that it starts from 1 */
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]++;
  for (i=0; i<nvtxs+1; i++)
    xadj[i]++;

  perm = idxmalloc(6*(nvtxs+5), "MMDOrder: perm");
  iperm = perm + nvtxs + 5;
  head = iperm + nvtxs + 5;
  qsize = head + nvtxs + 5;
  list = qsize + nvtxs + 5;
  marker = list + nvtxs + 5;

  genmmd(nvtxs, xadj, adjncy, iperm, perm, 1, head, qsize, list, marker, 100000, &nofsub);

  label = graph->label;
  firstvtx = lastvtx-nvtxs;
  for (i=0; i<nvtxs; i++)
    order[label[i]] = firstvtx+iperm[i]-1;

  free(perm);

  /* Relabel the vertices so that it starts from 0 */
  for (i=0; i<nvtxs+1; i++)
    xadj[i]--;
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]--;
}


/*************************************************************************
* This function takes a graph and a bisection and splits it into two graphs.
* It relies on the fact that adjwgt is all set to 1.
**************************************************************************/
int SplitGraphOrderCC(CtrlType *ctrl, GraphType *graph, GraphType *sgraphs, int ncmps, idxtype *cptr, idxtype *cind)
{
  int i, ii, iii, j, k, l, istart, iend, mypart, nvtxs, snvtxs, snedges;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum, *label, *where, *bndptr, *bndind;
  idxtype *sxadj, *svwgt, *sadjncy, *sadjwgt, *sadjwgtsum, *slabel;
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
  bndind = graph->bndind;
  ASSERT(bndptr != NULL);

  /* Go and use bndptr to also mark the boundary nodes in the two partitions */
  for (ii=0; ii<graph->nbnd; ii++) {
    i = bndind[ii];
    for (j=xadj[i]; j<xadj[i+1]; j++)
      bndptr[adjncy[j]] = 1;
  }

  rename = idxwspacemalloc(ctrl, nvtxs);
  
  /* Go and split the graph a component at a time */
  for (iii=0; iii<ncmps; iii++) {
    RandomPermute(cptr[iii+1]-cptr[iii], cind+cptr[iii], 0);
    snvtxs = snedges = 0;
    for (j=cptr[iii]; j<cptr[iii+1]; j++) {
      i = cind[j];
      rename[i] = snvtxs++;
      snedges += xadj[i+1]-xadj[i];
    }

    InitGraph(sgraphs+iii);

    /* Allocate memory and set the pointers accordingly */
    sgraphs[iii].gdata = idxmalloc(4*snvtxs+1+2*snedges, "SplitGraphOrder: lgdata");
    sxadj = sgraphs[iii].xadj               = sgraphs[iii].gdata;
    svwgt = sgraphs[iii].vwgt               = sgraphs[iii].gdata + snvtxs+1;
    sadjwgtsum = sgraphs[iii].adjwgtsum     = sgraphs[iii].gdata + 2*snvtxs+1;
    sgraphs[iii].cmap                       = sgraphs[iii].gdata + 3*snvtxs+1;
    sadjncy = sgraphs[iii].adjncy           = sgraphs[iii].gdata + 4*snvtxs+1;
    sadjwgt = sgraphs[iii].adjwgt           = sgraphs[iii].gdata + 4*snvtxs+1 + snedges;
    slabel = sgraphs[iii].label = idxmalloc(snvtxs, "SplitGraphOrder: lgraph->label");

    snvtxs = snedges = sxadj[0] = 0;
    for (ii=cptr[iii]; ii<cptr[iii+1]; ii++) {
      i = cind[ii];

      istart = xadj[i];
      iend = xadj[i+1];
      if (bndptr[i] == -1) { /* This is an interior vertex */
        auxadjncy = sadjncy + snedges - istart;
        auxadjwgt = sadjwgt + snedges - istart;
        for(j=istart; j<iend; j++) 
          auxadjncy[j] = adjncy[j];
        snedges += iend-istart;
      }
      else {
        l = snedges;
        for (j=istart; j<iend; j++) {
          k = adjncy[j];
          if (where[k] != 2) 
            sadjncy[l++] = k;
        }
        snedges = l;
      }

      svwgt[snvtxs] = vwgt[i];
      sadjwgtsum[snvtxs] = snedges-sxadj[snvtxs];
      slabel[snvtxs] = label[i];
      sxadj[++snvtxs] = snedges;
    }

    idxset(snedges, 1, sadjwgt);
    for (i=0; i<snedges; i++) 
      sadjncy[i] = rename[sadjncy[i]];

    sgraphs[iii].nvtxs = snvtxs;
    sgraphs[iii].nedges = snedges;

    if (snvtxs < MMDSWITCH)
      sgraphs[iii].adjwgt = NULL;  /* A marker to call MMD on the driver */
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SplitTmr));

  idxwspacefree(ctrl, nvtxs);

  return ncmps;

}






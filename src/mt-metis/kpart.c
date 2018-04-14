/**
 * @file kpart.c
 * @brief KWay partitioning functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */


#ifndef KPART_C
#define KPART_C

#include "includes.h"

#define DBGLVL 0 /* 15 */

idx_t ParPartGraphKway(dctrl_t * dctrl, dgraph_t * const graph, 
    idx_t * const where) 
{
  INIT_PARALLEL();

  idx_t i,j,rv,sigrval;
  ctrl_t * ctrl;

  ctrl = dctrl->ctrl;
  const idx_t nparts = ctrl->nparts;

  #pragma omp master
  {
    /* set debugging */
    ctrl->dbglvl = DBGLVL;

    dctrl->tptr_void[0] = (void*)pimalloc(nthreads,"X");
  }
  #pragma omp barrier
  idx_t ** const gwhere = (idx_t**)dctrl->tptr_void[0];

  gwhere[myid] = imalloc(graph->mynvtxs[myid],"X");

  /* set up multipliers for making balance computations easier */
  ParSetupKWayBalMultipliers(dctrl, graph);

  startwctimer(ctrl->TotalTmr);
  rv = ParMlevelKWayPartitioning(dctrl,graph,gwhere); 
  stopwctimer(ctrl->TotalTmr);

  if (rv < 0) {
    printf("MlevelKWayPartitioning() returned %"PRIDX"\n",rv);
    return 0;
  }
  #pragma omp master
  {
    PrintTimers(ctrl);
    PrintExtraTimers(dctrl);

    PrintPartInfo(ctrl,graph,(const idx_t * const*)gwhere);
  }
  #pragma omp barrier

  /* copy partition information back */
  if (where) {
    for (i=0;i<graph->mynvtxs[myid];++i) {
        j = (i&BLOCKMASK)+(((LVTX_2_LBID(i)*nthreads)+myid)<<BLOCKSHIFT); 
      where[j] = gwhere[myid][i];
    }
  }

  gk_free((void**)&gwhere[myid],LTERM);
  #pragma omp barrier

  #pragma omp master
  {
    gk_free((void**)&gwhere,LTERM);
  }

  return rv;
}

idx_t ParMlevelKWayPartitioning(dctrl_t * const dctrl, dgraph_t * graph, 
    idx_t *const * const part)
{
  INIT_PARALLEL();

  idx_t i, j, objval=0, curobj=0, bestobj=0;
  real_t curbal=0.0, bestbal=0.0;
  dgraph_t *cgraph;
  int status;

  ctrl_t * ctrl = dctrl->ctrl;
  startcputimer(dctrl->cpuTmr);

  for (i=0; i<ctrl->ncuts; i++) {
    startwctimer(ctrl->CoarsenTmr);
    cgraph = ParCoarsenGraph(dctrl, graph);
    stopwctimer(ctrl->CoarsenTmr);

    pprintf("Graph coarsened to %"PRIDX"/%"PRIDX" vertices\n",cgraph->nvtxs,
        ctrl->CoarsenTo);

    dprintf("Exposed Edge Weight = %"PRIDX"\n", computeEEW(cgraph));

    startwctimer(ctrl->InitPartTmr);
    ParInitKWayPartitioning(dctrl, cgraph);
    stopwctimer(ctrl->InitPartTmr);

    pprintf("Initial cut of %"PRIDX" made\n",cgraph->mincut);
    ASSERT(ParComputeCut(dctrl,cgraph,(const idx_t **)cgraph->where) 
        == cgraph->mincut); 

    ParAllocateRefinementWorkSpace(dctrl, cgraph);

    startwctimer(ctrl->UncoarsenTmr);
    ParRefineKWay(dctrl, graph, cgraph);
    stopwctimer(ctrl->UncoarsenTmr);

    if (graph->mincut < 0) {
      pprintf("Something went horribly wrong partitioning because I have a "
          "cut of %"PRIDX"!!!\n",graph->mincut);
    }

    curobj = graph->mincut;

    curbal = ParComputeLoadImbalanceDiff(dctrl,graph, ctrl->nparts, 
        ctrl->pijbm, ctrl->ubfactors);

    if (i == 0 
        || (curbal <= 0.0005 && bestobj > curobj)
        || (bestbal > 0.0005 && curbal < bestbal)) {
      icopy(graph->mynvtxs[myid], graph->where[myid], part[myid]);
      bestobj = curobj;
      bestbal = curbal;
    }

    ParFreeRData(graph);

    if (bestobj == 0)
      break;
  }

  stopcputimer(dctrl->cpuTmr);

  return bestobj;
}

void ParSetup2WayBalMultipliers(dctrl_t *dctrl, dgraph_t *graph, 
    const real_t *tpwgts)
{
  idx_t i, j;
  ctrl_t * ctrl = dctrl->ctrl;

  #pragma omp master
  {
    for (i=0; i<2; i++) {
      ctrl->pijbm[i] = graph->invtvwgt[0]/tpwgts[i];
    }
  }
}

void ParSetupKWayBalMultipliers(dctrl_t * const dctrl, dgraph_t *const graph)
{
  idx_t i;
  ctrl_t * ctrl = dctrl->ctrl;

  #pragma omp for 
  for (i=0;i<ctrl->nparts;++i) {
    ctrl->pijbm[i] = graph->invtvwgt[0] / ctrl->tpwgts[i];
  }
}

void PrintPartInfo(ctrl_t * const ctrl, const dgraph_t * const graph, 
    const idx_t * const * const where)
{
  idx_t i, ii, j,ej, k, nvtxs, nparts, tvwgt,myid,mynvtxs;
  const idx_t *adjncy, *vwgt, *adjwgt; 
  idx_t*kpwgts;
  real_t *tpwgts, unbalance;
  idx_t pid, ndom, maxndom, minndom, tndom, *pptr, *pind, *pdom;
  idx_t ncmps, nover, *cptr, *cind, *cpwgts;
  const idx_t *mywhere;

  nvtxs  = graph->nvtxs;

  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  printf("****************************************************************\n");
  printf(" size of idx_t: %zubits, real_t: %zubits, idx_t *: %zubits\n", 
      8*sizeof(idx_t), 8*sizeof(real_t), 8*sizeof(idx_t *));
  printf("\n");
  printf("Graph Information ----------------------------------------------\n");
  printf("#Vertices: %"PRIDX", #Edges: %"PRIDX", #Parts: %"PRIDX"\n", 
    graph->nvtxs, graph->nedges/2, ctrl->nparts);

  printf("\n");
  printf("\n");
  printf("Direct k-way Partitioning --------------------------------------\n");


  /* Compute objective-related infomration */
  printf(" - Edgecut: %"PRIDX", communication volume: %"PRIDX".\n\n", 
    SerComputeCut(graph, where), 
    SerComputeVolume(graph, where));


  /* Compute constraint-related information */
  kpwgts = ismalloc(nparts, 0, "ComputePartitionInfo: kpwgts");

  for (myid=0;myid<graph->ndist;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    vwgt = graph->vwgt[myid];
    mywhere = where[myid];
    for (i=0; i<mynvtxs; ++i) {
      kpwgts[mywhere[i]] += vwgt[i];
    }
  }

  /* Report on balance */
  printf(" - Balance:\n");
  tvwgt = isum(nparts, kpwgts, 1);
  for (k=0, unbalance=1.0*kpwgts[k]/(tpwgts[k]*tvwgt), i=1; i<nparts; i++) {
    if (unbalance < 1.0*kpwgts[i]/(tpwgts[i]*tvwgt)) {
      unbalance = 1.0*kpwgts[i]/(tpwgts[i]*tvwgt);
      k = i;
    }
  }
  idx_t mvwgt = 0;
  for (myid=0;myid<graph->ndist;++myid) {
    if ((mynvtxs = graph->mynvtxs[myid]) > 0) {
      vwgt = graph->vwgt[myid];
      k = iargmax_strd(mynvtxs, (idx_t*)vwgt, 1);
      if (vwgt[k] > mvwgt) {
        mvwgt = vwgt[k];
      }
    }
  }

  printf("     constraint #%"PRIDX":  %5.3"PRREAL" out of %5.3"PRREAL"\n", 
      0, unbalance,
       1.0*nparts*mvwgt/
          (1.0*isum(nparts, kpwgts, 1)));
  printf("\n");

  tvwgt = isum(nparts, kpwgts, 1);
  for (k=0, unbalance=kpwgts[k]/(tpwgts[k]*tvwgt), i=1; i<nparts; i++) {
    if (unbalance < kpwgts[i]/(tpwgts[i]*tvwgt)) {
      unbalance = kpwgts[i]/(tpwgts[i]*tvwgt);
      k = i;
    }
  }

  printf(" - Most overweight partition:\n"
         "     pid: %"PRIDX", actual: %"PRIDX", desired: %"PRIDX", ratio: %.2"PRREAL"\n\n",
      k, kpwgts[k], (idx_t)(tvwgt*tpwgts[k]), unbalance);

  gk_free((void **)&kpwgts, LTERM);


  #ifdef XXX
  /* Compute subdomain adjacency information */
  pptr = imalloc(nparts+1, "ComputePartitionInfo: pptr");
  pind = imalloc(nvtxs, "ComputePartitionInfo: pind");
  pdom = imalloc(nparts, "ComputePartitionInfo: pdom");

  iarray2csr(nvtxs, nparts, (idx_t*)where, pptr, pind);

  maxndom = nparts+1;
  minndom = 0;
  for (tndom=0, pid=0; pid<nparts; pid++) {
    iset(nparts, 0, pdom);
    for (ii=pptr[pid]; ii<pptr[pid+1]; ii++) {
      i = pind[ii];
      for (j=xadj[i].start; j<xadj[i].end; j++)
        pdom[where[adjncy[j]]] += adjwgt[j];
    }
    pdom[pid] = 0;
    for (ndom=0, i=0; i<nparts; i++)
      ndom += (pdom[i] > 0 ? 1 : 0);
    tndom += ndom;
    if (pid == 0 || maxndom < ndom)
      maxndom = ndom;
    if (pid == 0 || minndom > ndom)
      minndom = ndom;
  }

  printf(" - Subdomain connectivity: max: %"PRIDX", min: %"PRIDX", avg: %.2"PRREAL"\n\n",
      maxndom, minndom, 1.0*tndom/nparts);
      
  gk_free((void **)&pptr, &pind, &pdom, LTERM);
  #endif


  printf("****************************************************************\n");
}            

void PrintExtraTimers(dctrl_t * dctrl) {
  printf("----------------------------------------------------------------\n");
  printf("Parallel Match Block: %7.5f sec\n",dctrl->matchParPartTmr);
  printf("  Random Matching: %7.5f sec\n",dctrl->rmMatchTmr);
  printf("    Prep Work: %7.5f sec\n",dctrl->rmHeaderTmr);
  printf("  SHEM Matching: %7.5f sec\n",dctrl->shemMatchTmr);
  printf("    Prep Work: %7.5f sec\n",dctrl->shemHeaderTmr);
  printf("  Match Loop: %7.5f sec\n",dctrl->matchLoopTmr);
  printf("  CMAP Generation: %7.5f sec\n",dctrl->cmapTmr);
  printf("Parallel Contraction Block: %7.5f sec\n",dctrl->contractParPartTmr);
  printf("  Nedges Count: %7.5f sec\n",dctrl->mmapTmr);
  printf("  Contraction Loop: %7.5f sec\n",dctrl->conLoopTmr);
  printf("Parallel Init Part: %7.5f sec\n",dctrl->ipTmr);
  printf("  Convert DGraph: %7.5f sec\n",dctrl->convertTmr);
  printf("  Make Local Copy: %7.5f sec\n",dctrl->lCopyTmr);
  printf("  Parallel Bisection: %7.5f sec\n",dctrl->parBisTmr);
  printf("  Rename Where: %7.5f sec\n",dctrl->renParTmr);
  printf("Parallel Project: %7.5f sec\n",dctrl->projectParPartTmr);
  printf("  Project Where: %7.5f sec\n",dctrl->projectWhereTmr);
  printf("  Project Ckrinfo: %7.5f sec\n",dctrl->projectCKRTmr);
  printf("    Project Internal: %7.5f sec\n",dctrl->projectIdTmr);
  printf("    Project External: %7.5f sec\n",dctrl->projectEdTmr);
  printf("Parallel Refinement Block: %7.3f sec\n",dctrl->refineParPartTmr);
  printf("  Pre-Queue Work: %7.5f sec\n",dctrl->refHdrTmr);
  printf("  Filling Queue: %7.5f sec\n",dctrl->rpqTmr);
  printf("  Emptying Queue: %7.5f sec\n",dctrl->emptyQueueTmr);
  printf("    Checking Move: %7.5f sec\n",dctrl->checkMovTmr);
  printf("    Making Move: %7.5f sec\n",dctrl->makMovTmr);
  printf("      Sync Partition Weights: %7.5f sec\n",dctrl->pwgtsTmr);
  printf("      Undoing Moves: %7.5f sec\n",dctrl->undoMovTmr);
  printf("      Finalizing Moves: %7.5f sec\n",dctrl->finMovTmr);
  printf("    Updating Neighbors: %7.5f sec\n",dctrl->updNeiTmr);
  printf("Total CPU Time: %7.5f sec\n",dctrl->cpuTmr);
  printf("----------------------------------------------------------------\n");
}

#endif


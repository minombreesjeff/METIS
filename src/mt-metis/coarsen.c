/**
 * @file coarsen.c
 * @brief Coarsening functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */


/**
* @brief Coarsen a graph to BSTART vertices
*
* @param dctrl The control structure
* @param graph The fine graph to coarsen
* @param nthreads the number of threads to use
*
* @return 
*/
dgraph_t *ParCoarsenGraph(dctrl_t * dctrl, dgraph_t *graph)
{
  INIT_PARALLEL();
  idx_t i, neqewgts, level=0;
  ctrl_t * ctrl = dctrl->ctrl;

  idx_t * const fcmap = imalloc(graph->mynvtxs[myid],"X");
  idx_t * const match = dctrl->tptr_idx[myid] = 
      ismalloc(graph->mynvtxs[myid],UNMATCHED,"X");
  idx_t *const * const gmatch = dctrl->tptr_idx;

  ASSERT(checkGraph(graph));

  /* determine if the weights on the edges are all the same */
  const idx_t base = graph->adjwgt[0][0];
  neqewgts=0;
  
  if (graph->mynvtxs[myid] > 0) {
    for (i=0;i<graph->xadj[myid][graph->mynvtxs[myid]];++i) {
      if (graph->adjwgt[myid][i] != base) {
        neqewgts = 1;
        break;
      }
    }
  }
  dl_omp_sumreduce(myid,neqewgts,dctrl->buffer1,nthreads);

  /* set the maximum allowed coarsest vertex weight */
  ctrl->maxvwgt[0] = 1.5*graph->tvwgt[0]/ctrl->CoarsenTo;

  ASSERT(gmatch[myid] == match);

  do {
    /* allocate memory for cmap, if it has not already been done due to
       multiple cuts */
    #pragma omp master
    {
      if (graph->cmap == NULL) {
        /* threads need to allocate their own chunk inside the matching 
         * functions */
        graph->cmap = pimalloc(nthreads, "X");
      }
    }
    #pragma omp barrier

    /* coarsening scheme selection used to go here */
    #ifndef NOSHEM
    if (!neqewgts) {
      ParMatch_RM(dctrl,graph,gmatch,fcmap);
    } else {
      ParMatch_SHEM(dctrl,graph,gmatch,fcmap);
    }
    #else
    ParMatch_RM(dctrl,graph,gmatch,fcmap);
    #endif

    graph = graph->coarser;

    neqewgts = 1;
    #ifdef DEBUG
    printf("[%"PRIDX"] NVTXS = %"PRIDX", NEDGES = %"PRIDX", EE = %"PRIDX"\n",level,graph->nvtxs,
        graph->nedges,computeEEW(graph));
    #endif

    #ifdef TRACK_BAD_MATCHES
    #pragma omp master
    {
      printf("END OF LEVEL nvtxs = %"PRIDX"\n",graph->nvtxs);
    }
    #endif


  } while (graph->nvtxs > ctrl->CoarsenTo && 
           graph->nvtxs < PAR_COARSEN_FRACTION*graph->finer->nvtxs && 
           graph->nedges > graph->nvtxs/2); 
  gk_free((void**)&fcmap,&match,LTERM);

  ASSERT(checkGraph(graph));

  return graph;
}

/**
* @brief Randomly match vertices in parallel and create a coarse graph
*
* @param dctrl The control structure
* @param graph The graph who's vertices to matche
* @param nthreads The number of threads to use
*
* @return The number of coarse vertices after matching 
*/
idx_t ParMatch_RM(dctrl_t * const dctrl, dgraph_t * const graph,
    idx_t * const * const gmatch, idx_t * const fcmap) 
{
  INIT_PARALLEL();

  ctrl_t * const ctrl = dctrl->ctrl;

  startwctimer(ctrl->MatchTmr);
  startwctimer(dctrl->rmMatchTmr);

  const idx_t nvtxs  = graph->nvtxs;

  const idx_t * const * const gxadj = (const idx_t **)graph->xadj;
  const idx_t * const * const gvwgt = (const idx_t **)graph->vwgt;
  const idx_t * const * const gadjncy = (const idx_t **)graph->adjncy;
  const idx_t * const * const gadjwgt = (const idx_t **)graph->adjwgt;
  idx_t ** const gcmap   = graph->cmap;

  const idx_t * const maxvwgt  = ctrl->maxvwgt;

  idx_t * const partial = dctrl->buffer1;
  idx_t * const buffer = dctrl->buffer2;

  startwctimer(dctrl->matchParPartTmr);

  idx_t cnvtxs = 0;

  startwctimer(dctrl->rmHeaderTmr);

  idx_t i, pi, ii, j, jj, k, maxidx,
      last_unmatched, seed,nbrid,lvtx,gvtx;

  /* local graph pointers */
  const idx_t mynvtxs = graph->mynvtxs[myid];
  const idx_t dshift = graph->dshift;
  const idx_t dmask = graph->dmask;
  const idx_t * const xadj = gxadj[myid];
  const idx_t * const vwgt = gvwgt[myid];
  const idx_t * const adjncy = gadjncy[myid];
  const idx_t * const adjwgt = gadjwgt[myid];

  idx_t * const perm = imalloc(mynvtxs,"X");
  idx_t * const cmap = gcmap[myid] = perm;  /* re-use perm */

  k = 0;
  seed = ctrl->seed + myid;

  irandArrayPermute_r(mynvtxs,perm,mynvtxs/8,1,&seed);

  last_unmatched=0; /* this is private but ... */

  idx_t * const match = gmatch[myid];

  stopwctimer(dctrl->rmHeaderTmr);
  startwctimer(dctrl->matchLoopTmr);
  for (pi=0; pi<mynvtxs;++pi) {
    /* request my matches */
    i = perm[pi];
    i = pi;

    if (match[i] == UNMATCHED) {  /* Unmatched */
      gvtx = LVTX_2_GVTX(i,myid,dshift);
      maxidx = gvtx;
      
      if (vwgt[i] < maxvwgt[0]) {
        /* Deal with island vertices. Match locally */
        if (xadj[i+1] == xadj[i]) { 
          last_unmatched = gk_max(pi, last_unmatched)+1;
          for (; last_unmatched<mynvtxs; last_unmatched++) {
            j = perm[last_unmatched];
            if (match[j] == UNMATCHED) {
              maxidx = LVTX_2_GVTX(j,myid,dshift);
              break;
            }
          }
        } else {
          /* Find a heavy-edge matching, subject to maxvwgt constraints */
          for (j=xadj[i]; j<xadj[i+1]; ++j) {
            k = adjncy[j];
            nbrid = GVTX_2_THRID(k,dmask);
            lvtx = GVTX_2_LVTX(k,dshift);
            if (vwgt[i]+gvwgt[nbrid][lvtx] <= maxvwgt[0] && 
                (gmatch[nbrid][lvtx] == UNMATCHED)) {
              #if defined(LWGT) && LWGT != 1
              if (nbrid == myid) { /* I own it, lets go with it */
                maxidx = k;
                break;
              } else if (maxidx == gvtx) {
                maxidx = k;
              }
              #else
              maxidx = k;
              break;
              #endif
            }
          }
        }
      }
      nbrid = GVTX_2_THRID(maxidx,dmask);
      lvtx = GVTX_2_LVTX(maxidx,dshift);
      /* handle unmatched vertices later -- and orphaned vertices */
      if (gvtx < maxidx) {
        match[i] = maxidx;
        gmatch[nbrid][lvtx] = gvtx;
      } else if (gvtx > maxidx) {
        gmatch[nbrid][lvtx] = gvtx;
        match[i] = maxidx;
      }
    }
  } /* outer match loop */
  #pragma omp barrier
  stopwctimer(dctrl->matchLoopTmr);
  startwctimer(dctrl->cmapTmr);

  #ifdef TRACK_BAD_MATCHES
  idx_t nbroken = 0;
  #endif
  cnvtxs = 0;
  /* match any unmatched vertices with themselves */
  for (i=0;i<mynvtxs;++i) {
    gvtx = LVTX_2_GVTX(i,myid,dshift);
    ASSERT(GVTX_2_LVTX(gvtx,dshift) == i);
    nbrid = GVTX_2_THRID(match[i],dmask);
    lvtx = GVTX_2_LVTX(match[i],dshift);
    if (match[i] == UNMATCHED) {
      match[i] = gvtx;
    } else if (gmatch[nbrid][lvtx] != gvtx) {
      match[i] = gvtx;
      #ifdef TRACK_BAD_MATCHES
      ++nbroken;
      #endif
    }

    if (MY_CVTX(gvtx,match[i],xadj[i+1]-xadj[i],
          gxadj[nbrid][lvtx+1]-gxadj[nbrid][lvtx])) {
      /* use global cmap[i] id */
      cmap[i] = LVTX_2_GVTX(cnvtxs,myid,dshift);
      fcmap[cnvtxs++] = gvtx;
    }
  }

  #ifdef TRACK_BAD_MATCHES
  printf("[%"PRIDX"] %"PRIDX"/%"PRIDX" vertices in broken matches %0.05f\n",myid,nbroken,
      mynvtxs,(nbroken*100.0f)/mynvtxs);
  #endif

  #pragma omp barrier
  /* tries to avoid false sharing */
  for (i=mynvtxs-1;i>=0;--i) {
    gvtx = LVTX_2_GVTX(i,myid,dshift);
    nbrid = GVTX_2_THRID(match[i],dmask);
    lvtx = GVTX_2_LVTX(match[i],dshift);
    if (!MY_CVTX(gvtx,match[i],xadj[i+1]-xadj[i],
          gxadj[nbrid][lvtx+1]-gxadj[nbrid][lvtx])) { 
      cmap[i] = gcmap[nbrid][lvtx];
    }
  }

  stopwctimer(dctrl->cmapTmr);
  stopwctimer(dctrl->matchParPartTmr);
  stopwctimer(ctrl->MatchTmr);
  stopwctimer(dctrl->rmMatchTmr);

  #pragma omp barrier

  #ifdef NOMASK
  ParCreateCoarseGraphNoMask(dctrl,graph,cnvtxs,(const idx_t **)gmatch,fcmap);
  #else
  ParCreateCoarseGraph(dctrl, graph, cnvtxs, (const idx_t **)gmatch,fcmap);
  #endif

  iset(cnvtxs,UNMATCHED,match);

  return cnvtxs;
}

#ifndef NOSHEM
/**
* @brief Match the vertices in a graph using the SHEM coarsening scheme in
*   parallel and create a coarse graph and store it in graph->coarse
*
* @param dctrl THe control structure to use
* @param graph The fine graph to coarsen
* @param nthreads The number of vertices to use
*
* @return The number of vertices in the coarse graph 
*/
idx_t ParMatch_SHEM(dctrl_t * const dctrl, dgraph_t * const graph,
    idx_t * const * const gmatch, idx_t * const fcmap) 
{
  INIT_PARALLEL();

  ctrl_t * const ctrl = dctrl->ctrl;

  startwctimer(ctrl->MatchTmr);
  startwctimer(dctrl->shemMatchTmr);

  const idx_t * const maxvwgt  = ctrl->maxvwgt;
  const idx_t nvtxs  = graph->nvtxs;
  const idx_t dshift = graph->dshift;
  const idx_t dmask = graph->dmask;

  const idx_t * const * const gxadj = (const idx_t **)graph->xadj;
  const idx_t * const * const gvwgt = (const idx_t **)graph->vwgt;
  const idx_t * const * const gadjncy = (const idx_t **)graph->adjncy;
  const idx_t * const * const gadjwgt = (const idx_t **)graph->adjwgt;

  idx_t ** const gcmap   = graph->cmap;

  idx_t * const partial = dctrl->buffer1;
  idx_t * const buffer = dctrl->buffer2;

  startwctimer(dctrl->matchParPartTmr);

  idx_t cnvtxs = 0;

  idx_t i, pi, ii, j, jj, t, k, maxidx, maxwgt, 
      last_unmatched, avgdegree, seed, nbrid,lvtx,gvtx,lcl;

  startwctimer(dctrl->shemHeaderTmr);

  /* local graph pointers */
  const idx_t mynvtxs = graph->mynvtxs[myid];
  const idx_t * const xadj = gxadj[myid];
  const idx_t * const vwgt = gvwgt[myid];
  const idx_t * const adjncy = gadjncy[myid];
  const idx_t * const adjwgt = gadjwgt[myid];

  idx_t * const perm = imalloc(mynvtxs,"X");
  idx_t * const tperm = imalloc(mynvtxs,"X");
  idx_t * const cmap = gcmap[myid] = perm;  /* re-use perm */
  
  idx_t * const degrees = imalloc(mynvtxs,"X");

  k = 0;

  seed = ctrl->seed + myid;

  avgdegree = 0.7*(graph->nedges/nvtxs);
  for (i=0;i<mynvtxs;++i) {
    j = xadj[i+1] - xadj[i];
    degrees[i] = (j > avgdegree ? avgdegree : j);
  }

  idx_t * const counts = ismalloc(avgdegree+2,0,"X");

  irandArrayPermute_r(mynvtxs,tperm,mynvtxs/8,1,&seed);
  BucketSortKeysInc_t(counts,mynvtxs,avgdegree,degrees,tperm,perm);

  idx_t * const match = gmatch[myid];

  last_unmatched=0; 

  stopwctimer(dctrl->shemHeaderTmr);
  startwctimer(dctrl->matchLoopTmr);

  for (pi=0; pi<mynvtxs;++pi) {
    /* request my matches */
    i = perm[pi];
      
    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxwgt = -1;
      gvtx = LVTX_2_GVTX(i,myid,dshift);
      maxidx = gvtx;
      lcl = 0;

      if (vwgt[i] < maxvwgt[0]) {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i+1] == xadj[i]) { 
          last_unmatched = gk_max(pi, last_unmatched)+1;
          for (; last_unmatched<mynvtxs; last_unmatched++) {
            j = perm[last_unmatched];
            if (match[j] == UNMATCHED) {
              maxidx = LVTX_2_GVTX(j,myid,dshift);
              break;
            }
          }
        } else {
          /* Find a heavy-edge matching, subject to maxvwgt constraints */
          for (j=xadj[i]; j<xadj[i+1]; ++j) {
            k = adjncy[j];
            nbrid = GVTX_2_THRID(k,dmask);
            lvtx = GVTX_2_LVTX(k,dshift);
            #if defined(LWGT) && LWGT != 1
            if (nbrid == myid) {
              lcl = adjwgt[j] * LWGT;
            } else {
              lcl = adjwgt[j];
            }
            #else
            lcl = adjwgt[j];
            #endif
            if (maxwgt < lcl && vwgt[i]+gvwgt[nbrid][lvtx] <= maxvwgt[0] && 
                gmatch[nbrid][lvtx] == UNMATCHED) {
              maxidx = k;
              maxwgt = lcl;
            }
          }
        }
      }
      /* handle unmatched vertices later */
      nbrid = GVTX_2_THRID(maxidx,dmask);
      lvtx = GVTX_2_LVTX(maxidx,dshift);
      if (gvtx < maxidx) {
        match[i] = maxidx;
        gmatch[nbrid][lvtx] = gvtx;
      } else if (gvtx > maxidx) {
        gmatch[nbrid][lvtx] = gvtx;
        match[i] = maxidx;
      }
    }
  }/* outer match loop */
  #pragma omp barrier
  stopwctimer(dctrl->matchLoopTmr);
  startwctimer(dctrl->cmapTmr);


  cnvtxs = 0;
  /* match any unmatched vertices with themselves */
  for (i=0;i<mynvtxs;++i) {
    gvtx = LVTX_2_GVTX(i,myid,dshift);
    ASSERT(GVTX_2_LVTX(gvtx,dshift) == i);
    nbrid = GVTX_2_THRID(match[i],dmask);
    lvtx = GVTX_2_LVTX(match[i],dshift);
    if (match[i] == UNMATCHED || gmatch[nbrid][lvtx] != gvtx) {
      match[i] = gvtx;
    }
    if (MY_CVTX(gvtx,match[i],xadj[i+1]-xadj[i],
          gxadj[nbrid][lvtx+1]-gxadj[nbrid][lvtx])) {
      /* use global cmap[i] id */
      cmap[i] = LVTX_2_GVTX(cnvtxs,myid,dshift);
      fcmap[cnvtxs++] = gvtx;
    }
  }

  #pragma omp barrier
  /* tries to avoid false sharing */
  for (i=mynvtxs-1;i>=0;--i) {
    gvtx = LVTX_2_GVTX(i,myid,dshift);
    nbrid = GVTX_2_THRID(match[i],dmask);
    lvtx = GVTX_2_LVTX(match[i],dshift);
    if (!MY_CVTX(gvtx,match[i],xadj[i+1]-xadj[i],
          gxadj[nbrid][lvtx+1]-gxadj[nbrid][lvtx])) { 
      cmap[i] = gcmap[nbrid][lvtx];
    }
  }

  gk_free((void**)&degrees,&counts,&tperm,LTERM);

  stopwctimer(dctrl->cmapTmr);
  stopwctimer(dctrl->matchParPartTmr);
  stopwctimer(dctrl->shemMatchTmr);
  stopwctimer(ctrl->MatchTmr);

  #pragma omp barrier

  #ifdef NOMASK
  ParCreateCoarseGraphNoMask(dctrl,graph,cnvtxs,(const idx_t **)gmatch,fcmap);
  #else
  ParCreateCoarseGraph(dctrl,graph,cnvtxs,(const idx_t **)gmatch,fcmap);
  #endif

  iset(cnvtxs,UNMATCHED,match);

  return cnvtxs;
}
#endif

#ifndef NOMASK
/**
* @brief Create a coarse graph given a matching using hash table to identify
*   edges to be merged
*
* @param dctrl The control structure
* @param graph The fine graph
* @param cnvtxs The number of coarse vertices in the coarse graph
* @param match The matchings of vertices (match[match[v]] = v)
* @param fcmap The mapping of the coarse vertex number to the lowest number
*   fine vertex in the coarse vertex
* @param nthreads The number of threads to be used
*/
void ParCreateCoarseGraph(dctrl_t * const dctrl, dgraph_t * const graph, 
    const idx_t cnvtxs, const idx_t * const * const gmatch, 
    const idx_t * const fcmap)
{
  INIT_PARALLEL();

  ctrl_t * const ctrl = dctrl->ctrl;

  startwctimer(ctrl->ContractTmr);

  const idx_t nvtxs = graph->nvtxs;
  const idx_t dshift = graph->dshift;
  const idx_t dmask = graph->dmask;

  const idx_t * const * const gcmap = (const idx_t **)graph->cmap;

  dgraph_t * const cgraph = ParSetupCoarseGraph(dctrl->buffer1,&graph->coarser,
      graph, cnvtxs);

  ASSERT(graph->coarser == cgraph);

  /* fine graph pointers */
  const idx_t * const * const gxadj = (const idx_t **)graph->xadj;
  const idx_t * const * const gvwgt = (const idx_t **)graph->vwgt;
  const idx_t * const * const gadjncy = (const idx_t **)graph->adjncy;
  const idx_t * const * const gadjwgt = (const idx_t **)graph->adjwgt;

  /* coarse graph pointers */
  idx_t ** const gcxadj = cgraph->xadj;
  idx_t ** const gcvwgt = cgraph->vwgt;
  idx_t ** const gcadjncy = cgraph->adjncy;
  idx_t ** const gcadjwgt = cgraph->adjwgt;

  #pragma omp master
  {
    cgraph->nedges = 0;

    cgraph->tvwgt[0] = graph->tvwgt[0];
    cgraph->invtvwgt[0] = graph->invtvwgt[0];
  }

  startwctimer(dctrl->contractParPartTmr);

  idx_t i, j, jj, k, kk, l, m, istart, iend, nedges, 
      pv, pi, v, u, adjsize,lv,utid,lu,lvtx,nbrid,deg;

  idx_t maxdeg = -1;

  const idx_t mynvtxs = cgraph->mynvtxs[myid];
  const idx_t * match = gmatch[myid];

  /* fine graph parts */
  const idx_t * xadj = gxadj[myid];
  const idx_t * vwgt = gvwgt[myid];
  const idx_t * adjncy = gadjncy[myid];
  const idx_t * adjwgt = gadjwgt[myid];
  const idx_t * cmap = gcmap[myid];

  /* coarse graph parts */
  idx_t * const cxadj = gcxadj[myid];
  idx_t * const cvwgt = gcvwgt[myid];
  
  idx_t * const partial = dctrl->buffer1;
  idx_t * const buffer = dctrl->buffer2;

  startwctimer(dctrl->mmapTmr);

  adjsize = deg = 0;
  /* initialize the mmap */
  for (pi=0;pi<mynvtxs;++pi) {
    pv = LVTX_2_GVTX(pi,myid,dshift);
    v = fcmap[pi];
    lv = GVTX_2_LVTX(v,dshift);
    ASSERT(myid == GVTX_2_THRID(v,dmask));
    u = gmatch[myid][lv];
    if (v != u) {
      utid = GVTX_2_THRID(u,dmask);
      lu = GVTX_2_LVTX(u,dshift);
      ASSERT(gmatch[utid][lu] == v);
      ASSERT(gcmap[utid][lu] == pv);
      deg =xadj[lv+1] - xadj[lv] + gxadj[utid][lu+1] - gxadj[utid][lu];
    } else {
      deg = xadj[lv+1] - xadj[lv];
    }
    if (deg > maxdeg) {
      maxdeg = deg;
    }
    adjsize += deg;
  }

  stopwctimer(dctrl->mmapTmr);
  startwctimer(dctrl->conLoopTmr);

  idx_t mask = 0x1000;

  /* adaptive */
  while (maxdeg > (mask >> 3) || graph->nedges/graph->nvtxs > mask/20) {
    if (cgraph->nvtxs < mask*2) {
      mask = 0;
      break;
    } else {
      mask <<= 1;
    }
  }

  idx_t * const cadjncy = gcadjncy[myid] = imalloc(adjsize,"X");
  idx_t * const cadjwgt = gcadjwgt[myid] = imalloc(adjsize,"X");

  maxdeg = nedges = 0; 
  if (mask) {
    idx_t * const htable = ismalloc(mask,-1,"X"); 
    --mask; /* make 1000 -> 0FFF */
    for (pi=0; pi<mynvtxs; ++pi) {
      pv = LVTX_2_GVTX(pi,myid,dshift);

      cxadj[pi] = nedges;
      v = fcmap[pi];

      i = GVTX_2_LVTX(v,dshift);

      ASSERT(myid == GVTX_2_THRID(v,dmask));
      ASSERT(LVTX_2_GVTX(i,myid,dshift) == v);
      ASSERT(cmap[i] == pv);
      ASSERT(i<graph->mynvtxs[myid]);
      ASSERT(i>=0);

      istart = xadj[i];
      iend   = xadj[i+1];
      for (j=istart; j<iend; ++j) {
        nbrid = GVTX_2_THRID(adjncy[j],dmask); 
        lvtx = GVTX_2_LVTX(adjncy[j],dshift);
        k = gcmap[nbrid][lvtx];
        if (k == pv) {
          continue;
        }
        ASSERT(k < cgraph->gnvtxs);
        cEdgeMask(cxadj[pi],j,jj,k,kk,m,nedges,mask,htable,adjncy,adjwgt,cadjncy,
            cadjwgt);
      }
      cvwgt[pi] = vwgt[i];

      if ((u = gmatch[myid][i]) != v) {
        lu = GVTX_2_LVTX(u,dshift);
        utid = GVTX_2_THRID(u,dmask);
        ASSERT(gcmap[utid][lu] == gcmap[myid][i]);
        ASSERT(gmatch[utid][lu] == v);
        /* combine u and v stuff */
        cvwgt[pi] += gvwgt[utid][lu];
        istart = gxadj[utid][lu];
        iend   = gxadj[utid][lu+1];
        for (j=istart; j<iend; ++j) {
          nbrid = GVTX_2_THRID(gadjncy[utid][j],dmask); 
          lvtx = GVTX_2_LVTX(gadjncy[utid][j],dshift);
          k = gcmap[nbrid][lvtx];
          if (k == pv) {
            continue;
          }
          ASSERT(k < cgraph->gnvtxs);
          cEdgeMask(cxadj[pi],j,jj,k,kk,m,nedges,mask,htable,gadjncy[utid],
              gadjwgt[utid],cadjncy,cadjwgt);
        }
      }

      /* Zero out the htable */
      for (j=cxadj[pi]; j<nedges; j++) {
        k = cadjncy[j];
        htable[k&mask] = -1;  
      }

      cxadj[pi+1] = nedges;
      ASSERT(nedges <= adjsize);

      if (cxadj[pi+1] - cxadj[pi] > maxdeg) {
         maxdeg = cxadj[pi+1] - cxadj[pi];
      }
    }
    gk_free((void**)&htable,LTERM);
  } else {
    idx_t * const htable = ismalloc(cgraph->gnvtxs,-1,"X"); 
    for (pi=0;pi<mynvtxs;++pi) {
      pv = LVTX_2_GVTX(pi,myid,dshift);

      cxadj[pi] = nedges;
      v = fcmap[pi];

      i = GVTX_2_LVTX(v,dshift);

      ASSERT(myid == GVTX_2_THRID(v,dmask));
      ASSERT(LVTX_2_GVTX(i,myid,dshift) == v);
      ASSERT(cmap[i] == pv);
      ASSERT(i<graph->mynvtxs[myid]);
      ASSERT(i>=0);

      istart = xadj[i];
      iend   = xadj[i+1];
      for (j=istart; j<iend; ++j) {
        nbrid = GVTX_2_THRID(adjncy[j],dmask); 
        lvtx = GVTX_2_LVTX(adjncy[j],dshift);
        k = gcmap[nbrid][lvtx];
        if (k == pv) {
          continue;
        }
        ASSERT(k < cgraph->gnvtxs);
        cEdgeNoMask(j,k,m,nedges,htable,adjncy,adjwgt,cadjncy,cadjwgt);
      }
      cvwgt[pi] = vwgt[i];

      if ((u = gmatch[myid][i]) != v) {
        lu = GVTX_2_LVTX(u,dshift);
        utid = GVTX_2_THRID(u,dmask);
        ASSERT(gcmap[utid][lu] == gcmap[myid][i]);
        ASSERT(gmatch[utid][lu] == v);
        /* combine u and v stuff */
        cvwgt[pi] += gvwgt[utid][lu];
        istart = gxadj[utid][lu];
        iend   = gxadj[utid][lu+1];
        for (j=istart; j<iend; ++j) {
          nbrid = GVTX_2_THRID(gadjncy[utid][j],dmask); 
          lvtx = GVTX_2_LVTX(gadjncy[utid][j],dshift);
          k = gcmap[nbrid][lvtx];
          if (k == pv) {
            continue;
          }
          cEdgeNoMask(j,k,m,nedges,htable,gadjncy[utid],gadjwgt[utid],cadjncy,
              cadjwgt);
        }
      }

      /* Zero out the htable and find the heaveist edge */
      for (j=cxadj[pi]; j<nedges; j++) {
        k = cadjncy[j];
        htable[k] = -1;
      }

      cxadj[pi+1] = nedges;
      ASSERT(nedges <= adjsize);
      
      if (cxadj[pi+1] - cxadj[pi] > maxdeg) {
         maxdeg = cxadj[pi+1] - cxadj[pi];
      }
    }
    gk_free((void**)&htable,LTERM);
  }
  stopwctimer(dctrl->conLoopTmr);

  cgraph->mymaxdeg[myid] = maxdeg;
  nParReAdjustMemory(myid,dctrl,cgraph,adjsize);

  dl_omp_sumreduce(myid,nedges,dctrl->buffer1,nthreads);
  dl_omp_maxreduce(myid,maxdeg,dctrl->buffer1,nthreads);

  #pragma omp master
  {
    cgraph->nedges = nedges;
    cgraph->maxdeg = maxdeg;
    ASSERT(checkGraph(cgraph) == 1);
  }

  #pragma omp barrier
  stopwctimer(dctrl->contractParPartTmr);
  stopwctimer(ctrl->ContractTmr);
}

#else

/**
* @brief Create a coarse graph given a matching using a dense array to identify
*   edges to be merged
*
* @param dctrl The control structure
* @param graph The fine graph
* @param cnvtxs The number of coarse vertices in the coarse graph
* @param match The matchings of vertices (match[match[v]] = v)
* @param fcmap The mapping of the coarse vertex number to the lowest number
*   fine vertex in the coarse vertex
* @param nthreads The number of threads to be used
*/
void ParCreateCoarseGraphNoMask(dctrl_t * const dctrl, dgraph_t * const graph, 
    const idx_t cnvtxs, const idx_t * const * const gmatch, 
    const idx_t * const * fcmap, const idx_t nthreads)
{
  INIT_PARALLEL();

  ctrl_t * const ctrl = dctrl->ctrl;

  startwctimer(ctrl->ContractTmr);

  const idx_t nvtxs = graph->nvtxs;
  const idx_t dshift = graph->dshift;
  const idx_t dmask = graph->dmask;

  const idx_t * const * const gcmap = (const idx_t **)graph->cmap;

  dgraph_t * const cgraph = ParSetupCoarseGraph(&graph->coarser,dctrl->buffer1,
      graph, cnvtxs);

  /* fine graph pointers */
  const idx_t * const * const gxadj = (const idx_t **)graph->xadj;
  const idx_t * const * const gvwgt = (const idx_t **)graph->vwgt;
  const idx_t * const * const gadjncy = (const idx_t **)graph->adjncy;
  const idx_t * const * const gadjwgt = (const idx_t **)graph->adjwgt;

  /* coarse graph pointers */
  idx_t ** const gcxadj = cgraph->xadj;
  idx_t ** const gcvwgt = cgraph->vwgt;
  idx_t ** const gcadjncy = cgraph->adjncy;
  idx_t ** const gcadjwgt = cgraph->adjwgt;

  #pragma omp master
  {
    cgraph->nedges = 0;
    cgraph->tvwgt[0] = graph->tvwgt[0];
    cgraph->invtvwgt[0] = graph->invtvwgt[0];
  }

  startwctimer(dctrl->contractParPartTmr);

  idx_t i,j, k, l, m, istart, iend, nedges, pv, pi, v, u, adjsize,
        lvtx,nbrid,lv,lu,utid;
  idx_t maxdeg = -1;

  const idx_t mynvtxs = cgraph->mynvtxs[myid];
  const idx_t * match = gmatch[myid];

  /* fine graph parts */
  const idx_t * xadj = gxadj[myid];
  const idx_t * vwgt = gvwgt[myid];
  const idx_t * adjncy = gadjncy[myid];
  const idx_t * adjwgt = gadjwgt[myid];
  const idx_t * cmap = gcmap[myid];

  /* coarse graph parts */
  idx_t * const cxadj = gcxadj[myid];
  idx_t * const cvwgt = gcvwgt[myid];
  
  idx_t * const partial = dctrl->buffer1;
  idx_t * const buffer = dctrl->buffer2;

  istart = 0;
  iend = 0;

  startwctimer(dctrl->mmapTmr);

  /* initialize the mmap */
  for (pi=0;pi<mynvtxs;++pi) {
    pv = LVTX_2_GVTX(pi,myid,dshift);
    v = fcmap[pi];
    lv = GVTX_2_LVTX(v,dshift);
    ASSERT(myid == GVTX_2_THRID(v,dmask));
    istart += xadj[lv];
    iend += xadj[lv+1];
    u = gmatch[myid][lv];
    if (v != u) {
      utid = GVTX_2_THRID(u,dmask);
      lu = GVTX_2_LVTX(u,dshift);
      ASSERT(gmatch[utid][lu] == v);
      ASSERT(gcmap[utid][lu] == pv);
      istart += gxadj[utid][lu]; 
      iend += gxadj[utid][lu+1];
    }
  }

  stopwctimer(dctrl->mmapTmr);
  startwctimer(dctrl->conLoopTmr);

  adjsize = iend-istart;

  idx_t * const cadjncy = gcadjncy[myid] = imalloc(adjsize,"X");
  idx_t * const cadjwgt = gcadjwgt[myid] = imalloc(adjsize,"X");

  idx_t * const htable = ismalloc(cgraph->gnvtxs,-1,"X"); 

  nedges = 0; 
  for (pi=0;pi<mynvtxs;++pi) {
    pv = LVTX_2_GVTX(pi,myid,dshift);

    cxadj[pi] = nedges;
    v = fcmap[pi];

    i = GVTX_2_LVTX(v,dshift);

    ASSERT(myid == GVTX_2_THRID(v,dmask));
    ASSERT(LVTX_2_GVTX(i,myid,dshift) == v);
    ASSERT(cmap[i] == pv);
    ASSERT(i<graph->mynvtxs[myid]);
    ASSERT(i>=0);

    istart = xadj[i];
    iend   = xadj[i+1];
    for (j=istart; j<iend; ++j) {
      nbrid = GVTX_2_THRID(adjncy[j],dmask); 
      lvtx = GVTX_2_LVTX(adjncy[j],dshift);
      k = gcmap[nbrid][lvtx];
      if (k == pv) {
        continue;
      }
      ASSERT(k < cgraph->gnvtxs);
      cEdgeNoMask(j,k,m,nedges,htable,adjncy,adjwgt,cadjncy,cadjwgt);
    }
    cvwgt[pi] = vwgt[i];

    if ((u = gmatch[myid][i]) != v) {
      lu = GVTX_2_LVTX(u,dshift);
      utid = GVTX_2_THRID(u,dmask);
      ASSERT(gcmap[utid][lu] == gcmap[myid][i]);
      ASSERT(gmatch[utid][lu] == v);
      /* combine u and v stuff */
      cvwgt[pi] += gvwgt[utid][lu];
      istart = gxadj[utid][lu];
      iend   = gxadj[utid][lu+1];
      for (j=istart; j<iend; ++j) {
        nbrid = GVTX_2_THRID(gadjncy[utid][j],dmask); 
        lvtx = GVTX_2_LVTX(gadjncy[utid][j],dshift);
        k = gcmap[nbrid][lvtx];
        if (k == pv) {
          continue;
        }
        cEdgeNoMask(j,k,m,nedges,htable,gadjncy[utid],gadjwgt[utid],cadjncy,
            cadjwgt);
      }
    }

    /* Zero out the htable and find the heaveist edge */
    for (j=cxadj[pi]; j<nedges; j++) {
      k = cadjncy[j];
      htable[k] = -1;
    }

    cxadj[pi+1] = nedges;
    ASSERT(nedges <= adjsize);
    
    if (cxadj[pi+1] - cxadj[pi] > maxdeg) {
       maxdeg = cxadj[pi+1] - cxadj[pi];
    }
  }
  gk_free((void**)&htable,LTERM);
  stopwctimer(dctrl->conLoopTmr);

  nParReAdjustMemory(myid,dctrl,cgraph,adjsize);

  dl_omp_sumreduce(myid,nedges,dctrl->buffer1,rnthreads);
  dl_omp_maxreduce(myid,maxdeg,dctrl->buffer2,rnthreads);

  #pragma omp master
  {
    cgraph->nedges = nedges;
    cgraph->maxdeg = maxdeg;
    ASSERT(checkGraph(cgraph) == 1);
  }

  stopwctimer(dctrl->contractParPartTmr);
  stopwctimer(ctrl->ContractTmr);
}
#endif



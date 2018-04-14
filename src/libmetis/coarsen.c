/*!
\file  
\brief Functions for computing matchings during graph coarsening

\date Started 7/23/97
\author George  
\author Copyright 1997-2011, Regents of the University of Minnesota 
\version\verbatim $Id: coarsen.c 10105 2011-06-07 20:36:11Z karypis $ \endverbatim
*/


#include "metislib.h"


/*************************************************************************/
/*! This function takes a graph and creates a sequence of coarser graphs.
    It implements the coarsening phase of the multilevel paradigm. */
/*************************************************************************/
graph_t *CoarsenGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, eqewgts;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->CoarsenTmr));

  /* determine if the weights on the edges are all the same */
  for (eqewgts=1, i=1; i<graph->nedges; i++) {
    if (graph->adjwgt[0] != graph->adjwgt[i]) {
      eqewgts = 0;
      break;
    }
  }

  /* set the maximum allowed coarsest vertex weight */
  for (i=0; i<graph->ncon; i++)
    ctrl->maxvwgt[i] = 1.5*graph->tvwgt[i]/ctrl->CoarsenTo;

  do {
    IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));

    /* determine which matching scheme you will use */
    switch (ctrl->ctype) {
      case METIS_CTYPE_RM:
        Match_RM(ctrl, graph);
        break;
      case METIS_CTYPE_SHEM:
        if (eqewgts || graph->nedges == 0)
          Match_RM(ctrl, graph);
        else
          Match_SHEM(ctrl, graph);
        break;
      default:
        errexit("Unknown ctype: %d\n", ctrl->ctype);
    }

    graph = graph->coarser;
    eqewgts = 0;

  } while (graph->nvtxs > ctrl->CoarsenTo && 
           graph->nvtxs < COARSEN_FRACTION2*graph->finer->nvtxs && 
           graph->nedges > graph->nvtxs/2); 

  IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->CoarsenTmr));

  ASSERT(CheckGraph(graph, 0, 1));

  return graph;
}


/*************************************************************************/
/*! This function finds a matching by randomly selecting one of the 
    unmatched adjacent vertices. */
/**************************************************************************/
idx_t Match_RM(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ii, j, k, nvtxs, ncon, cnvtxs, maxidx, last_unmatched;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
  idx_t *match, *cmap, *perm;

  WCOREPUSH;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->MatchTmr));

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  cmap   = graph->cmap;

  maxvwgt  = ctrl->maxvwgt;

  match = iwspacemalloc(ctrl, nvtxs);
  perm  = iwspacemalloc(ctrl, nvtxs);

  irandArrayPermute(nvtxs, perm, nvtxs/5, 1);
  iset(nvtxs, UNMATCHED, match); 

  for (cnvtxs=0, last_unmatched=0, ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      if ((ncon == 1 ? vwgt[i] < maxvwgt[0] : ivecle(ncon, vwgt+i*ncon, maxvwgt))) {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i] == xadj[i+1]) { 
          last_unmatched = gk_max(ii, last_unmatched)+1;
          for (; last_unmatched<nvtxs; last_unmatched++) {
            j = perm[last_unmatched];
            if (match[j] == UNMATCHED) {
              maxidx = j;
              break;
            }
          }
        }
        else {
          /* Find a random matching, subject to maxvwgt constraints */
          if (ncon == 1) {
            /* single constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && vwgt[i]+vwgt[k] <= maxvwgt[0]) {
                maxidx = k;
                break;
              }
            }
          }
          else {
            /* multi-constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && 
                  ivecaxpylez(ncon, 1, vwgt+i*ncon, vwgt+k*ncon, maxvwgt)) {
                maxidx = k;
                break;
              }
            }
          }
        }
      }

      cmap[i]  = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  WCOREPOP;

  return cnvtxs;
}


/**************************************************************************/
/*! This function finds a matching using the HEM heuristic. The vertices 
 are visited based on increasing degree to ensure that all vertices are 
 given a chance to match with something. */
/**************************************************************************/
idx_t Match_SHEM(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ii, j, k, nvtxs, ncon, cnvtxs, maxidx, maxwgt, 
        last_unmatched, avgdegree;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
  idx_t *match, *cmap, *degrees, *perm, *tperm;

  WCOREPUSH;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->MatchTmr));

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  cmap   = graph->cmap;

  maxvwgt  = ctrl->maxvwgt;

  match   = iwspacemalloc(ctrl, nvtxs);
  perm    = iwspacemalloc(ctrl, nvtxs);
  tperm   = iwspacemalloc(ctrl, nvtxs);
  degrees = iwspacemalloc(ctrl, nvtxs);

  iset(nvtxs, UNMATCHED, match); 
  irandArrayPermute(nvtxs, tperm, nvtxs/10, 1);

  avgdegree = 0.7*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++) 
    degrees[i] = (xadj[i+1]-xadj[i] > avgdegree ? avgdegree : xadj[i+1]-xadj[i]);
  BucketSortKeysInc(ctrl, nvtxs, avgdegree, degrees, tperm, perm);

  for (cnvtxs=0, last_unmatched=0, ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = -1;

      if ((ncon == 1 ? vwgt[i] < maxvwgt[0] : ivecle(ncon, vwgt+i*ncon, maxvwgt))) {
        /* Deal with island vertices. Find a non-island and match it with. 
           The matching ignores ctrl->maxvwgt requirements */
        if (xadj[i] == xadj[i+1]) { 
          last_unmatched = gk_max(ii, last_unmatched)+1;
          for (; last_unmatched<nvtxs; last_unmatched++) {
            j = perm[last_unmatched];
            if (match[j] == UNMATCHED) {
              maxidx = j;
              break;
            }
          }
        }
        else {
          /* Find a heavy-edge matching, subject to maxvwgt constraints */
          if (ncon == 1) {
            /* single constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && 
                  maxwgt < adjwgt[j] && vwgt[i]+vwgt[k] <= maxvwgt[0]) {
                maxidx = k;
                maxwgt = adjwgt[j];
              }
            }
          }
          else {
            /* multi-constraint version */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = adjncy[j];
              if (match[k] == UNMATCHED && 
                  ivecaxpylez(ncon, 1, vwgt+i*ncon, vwgt+k*ncon, maxvwgt) &&
                  (maxwgt < adjwgt[j] || 
                   (maxwgt == adjwgt[j] && 
                    BetterVBalance(ncon, graph->invtvwgt, vwgt+i*ncon, 
                        vwgt+maxidx*ncon, vwgt+k*ncon)))) {
                maxidx = k;
                maxwgt = adjwgt[j];
              }
            }
          }
        }
      }

      cmap[i]  = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }


  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  WCOREPOP;

  return cnvtxs;
}


/*************************************************************************/
/*! This function prints various stats for each graph during coarsening */
/*************************************************************************/
void PrintCGraphStats(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i;

  printf("%6"PRIDX" %7"PRIDX" %10"PRIDX" [%"PRIDX"] [", 
      graph->nvtxs, graph->nedges, isum(graph->nvtxs, graph->adjrsum, 1), ctrl->CoarsenTo);

  for (i=0; i<graph->ncon; i++)
    printf(" %6"PRIDX":%6"PRIDX, ctrl->maxvwgt[i], graph->tvwgt[i]);
  printf(" ]\n");
}

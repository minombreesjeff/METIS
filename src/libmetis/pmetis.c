/**
\file
\brief This file contains the top level routines for the multilevel recursive bisection 
       algorithm PMETIS.

\date   Started 7/24/1997
\author George  
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version\verbatim $Id: pmetis.c 10237 2011-06-14 15:22:13Z karypis $ \endverbatim
*/


#include "metislib.h"


/*************************************************************************/
/*! \ingroup api 
    \brief Recursive partitioning routine.

    This function computes a partitioning of a graph based on multilevel
    recursive bisection. It can be used to partition a graph into \e k 
    parts. The objective of the partitioning is to minimize the edgecut
    subject to one or more balancing constraints.

    \param[in] nvtxs is the number of vertices in the graph.

    \param[in] ncon is the number of balancing constraints. For the standard
           partitioning problem in which each vertex is either unweighted
           or has a single weight, ncon should be 1.

    \param[in] xadj is an array of size nvtxs+1 used to specify the starting
           positions of the adjacency structure of the vertices in the
           adjncy array.

    \param[in] adjncy is an array of size to the sum of the degrees of the
           graph that stores for each vertex the set of vertices that
           is adjancent to.

    \param[in] vwgt is an array of size nvtxs*ncon that stores the weights
           of the vertices for each constraint. The ncon weights for the
           ith vertex are stored in the ncon consecutive locations starting
           at vwgt[i*ncon]. When ncon==1, a NULL value can be passed indicating
           that all the vertices in the graph have the same weight.

    \param[in] adjwgt is an array of size equal to adjncy, specifying the weight
           for each edge (i.e., adjwgt[j] corresponds to the weight of the
           edge stored in adjncy[j]). 
           A NULL value can be passed indicating that all the edges in the 
           graph have the same weight.

    \param[in] nparts is the number of desired partitions.

    \param[in] tpwgts is an array of size nparts*ncon that specifies the
           desired weight for each part and constraint. The \e{target partition
           weight} for the ith part and jth constraint is specified
           at tpwgts[i*ncon+j] (the numbering of i and j starts from 0).
           For each constraint, the sum of the tpwgts[] entries must be
           1.0 (i.e., \f$ \sum_i tpwgts[i*ncon+j] = 1.0 \f$). 
           A NULL value can be passed indicating that the graph should
           be equally divided among the parts.

    \param[in] ubvec is an array of size ncon that specifies the allowed 
           load imbalance tolerance for each constraint. 
           For the ith part and jth constraint the allowed weight is the 
           ubvec[j]*tpwgts[i*ncon+j] fraction of the jth's constraint total
           weight. The load imbalances must be greater than 1.0. 
           A NULL value can be passed indicating that the load imbalance
           tolerance for each constraint should be 1.001 (for ncon==1)
           or 1.01 (for ncon>1).

    \params[in] options is the array for passing additional parameters
           in order to customize the behaviour of the partitioning
           algorithm.

    \params[out] edgecut stores the cut of the partitioning.

    \params[out] part is an array of size nvtxs used to store the 
           computed partitioning. The partition number for the ith
           vertex is stored in part[i]. Based on the numflag parameter,
           the numbering of the parts starts from either 0 or 1.


    \returns 
      \retval METIS_OK  indicates that the function returned normally.
      \retval METIS_ERROR_INPUT indicates an input error.
      \retval METIS_ERROR_MEMORY indicates that it could not allocate 
              the required memory.
           
*/
/*************************************************************************/
int METIS_PartGraphRecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, 
          idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, 
          idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
          idx_t *objval, idx_t *part)
{
  graph_t *graph;
  ctrl_t *ctrl;

  /* set up the run parameters */
  ctrl = SetupCtrl(METIS_OP_PMETIS, options, *ncon, *nparts, tpwgts, ubvec);
  if (!ctrl) return METIS_ERROR_INPUT;

  /* set up the graph */
  graph = SetupGraph(ctrl, *nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);

  /* allocate workspace memory */
  AllocateWorkSpace(ctrl, graph);

  /* start the partitioning */
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->TotalTmr));

  *objval = MlevelRecursiveBisection(ctrl, graph, *nparts, part, ctrl->tpwgts, 0);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->TotalTmr));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

  if (ctrl->numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);

  /* clean up */
  FreeCtrl(&ctrl);

  return METIS_OK;
}


/*************************************************************************/
/*! This function is the top-level driver of the recursive bisection 
    routine. */
/*************************************************************************/
idx_t MlevelRecursiveBisection(ctrl_t *ctrl, graph_t *graph, idx_t nparts, 
          idx_t *part, real_t *tpwgts, idx_t fpart)
{
  idx_t i, j, nvtxs, ncon, objval;
  idx_t *label, *where;
  graph_t *lgraph, *rgraph;
  real_t wsum, *tpwgts2;

  if ((nvtxs = graph->nvtxs) == 0) {
    printf("\t***Cannot bisect a graph with 0 vertices!\n"
           "\t***You are trying to partition a graph into too many parts!\n");
    return 0;
  }

  ncon = graph->ncon;

  /* determine the weights of the two partitions as a function of the weight of the
     target partition weights */
  WCOREPUSH;
  tpwgts2 = rwspacemalloc(ctrl, 2*ncon);
  for (i=0; i<ncon; i++) {
    tpwgts2[i]      = rsum((nparts>>1), tpwgts+i, ncon);
    tpwgts2[ncon+i] = 1.0 - tpwgts2[i];
  }

  /* perform the bisection */
  objval = MultilevelBisect(ctrl, graph, tpwgts2);

  WCOREPOP;

  label = graph->label;
  where = graph->where;
  for (i=0; i<nvtxs; i++)
    part[label[i]] = where[i] + fpart;

  if (nparts > 2) 
    SplitGraphPart(ctrl, graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  FreeGraph(&graph);

  /* Scale the fractions in the tpwgts according to the true weight */
  for (i=0; i<ncon; i++) {
    wsum = rsum((nparts>>1), tpwgts+i, ncon);
    rscale((nparts>>1), 1.0/wsum, tpwgts+i, ncon);
    rscale(nparts-(nparts>>1), 1.0/(1.0-wsum), tpwgts+(nparts>>1)*ncon+i, ncon);
  }

  /* Do the recursive call */
  if (nparts > 3) {
    objval += MlevelRecursiveBisection(ctrl, lgraph, (nparts>>1), part, 
               tpwgts, fpart);
    objval += MlevelRecursiveBisection(ctrl, rgraph, nparts-(nparts>>1), part, 
               tpwgts+(nparts>>1)*ncon, fpart+(nparts>>1));
  }
  else if (nparts == 3) {
    FreeGraph(&lgraph);
    objval += MlevelRecursiveBisection(ctrl, rgraph, nparts-(nparts>>1), part, 
               tpwgts+(nparts>>1)*ncon, fpart+(nparts>>1));
  }


  return objval;
}


/*************************************************************************/
/*! This function performs a multilevel bisection */
/*************************************************************************/
idx_t MultilevelBisect(ctrl_t *ctrl, graph_t *graph, real_t *tpwgts)
{
  idx_t i, bestobj=0, curobj, *bestwhere;
  graph_t *cgraph;

  Setup2WayBalMultipliers(ctrl, graph, tpwgts);

  if (ctrl->ncuts > 1)
    bestwhere = imalloc(graph->nvtxs, "MultilevelBisect: bestwhere");

  for (i=0; i<ctrl->ncuts; i++) {
    cgraph = CoarsenGraph(ctrl, graph);

    Init2WayPartition(ctrl, cgraph, tpwgts);

    Refine2Way(ctrl, graph, cgraph, tpwgts);

    curobj = graph->mincut;

    if (ctrl->ncuts > 1) {
      if (i == 0 || bestobj > curobj) {
        icopy(graph->nvtxs, graph->where, bestwhere);
        bestobj = curobj;
      }

      if (i < ctrl->ncuts-1)
        FreeRData(graph);
      else if (bestobj != curobj) {
        icopy(graph->nvtxs, bestwhere, graph->where);
        Compute2WayPartitionParams(ctrl, graph);
      }
    }
    else /* ctrl->ncuts == 1 */
      bestobj = curobj;
  }

  if (ctrl->ncuts > 1)
    gk_free((void **)&bestwhere, LTERM);

  return bestobj;
}


/*************************************************************************/
/*! This function splits a graph into two based on its bisection */
/*************************************************************************/
void SplitGraphPart(ctrl_t *ctrl, graph_t *graph, graph_t **r_lgraph, 
         graph_t **r_rgraph)
{
  idx_t i, j, k, l, istart, iend, mypart, nvtxs, ncon, snvtxs[2], snedges[2], sum;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *adjrsum, *label, *where, *bndptr;
  idx_t *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *sadjrsum[2], 
        *slabel[2];
  idx_t *rename;
  idx_t *auxadjncy, *auxadjwgt;
  graph_t *lgraph, *rgraph;

  WCOREPUSH;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->SplitTmr));

  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  xadj    = graph->xadj;
  vwgt    = graph->vwgt;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  adjrsum = graph->adjrsum;
  label   = graph->label;
  where   = graph->where;
  bndptr  = graph->bndptr;

  ASSERT(bndptr != NULL);

  rename = iwspacemalloc(ctrl, nvtxs);
  
  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  for (i=0; i<nvtxs; i++) {
    k = where[i];
    rename[i] = snvtxs[k]++;
    snedges[k] += xadj[i+1]-xadj[i];
  }

  lgraph      = SetupSplitGraph(graph, snvtxs[0], snedges[0]);
  sxadj[0]    = lgraph->xadj;
  svwgt[0]    = lgraph->vwgt;
  sadjrsum[0] = lgraph->adjrsum;
  sadjncy[0]  = lgraph->adjncy; 	
  sadjwgt[0]  = lgraph->adjwgt; 
  slabel[0]   = lgraph->label;

  rgraph      = SetupSplitGraph(graph, snvtxs[1], snedges[1]);
  sxadj[1]    = rgraph->xadj;
  svwgt[1]    = rgraph->vwgt;
  sadjrsum[1] = rgraph->adjrsum;
  sadjncy[1]  = rgraph->adjncy; 	
  sadjwgt[1]  = rgraph->adjwgt; 
  slabel[1]   = rgraph->label;

  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  sxadj[0][0] = sxadj[1][0] = 0;
  for (i=0; i<nvtxs; i++) {
    mypart = where[i];
    sum = adjrsum[i];

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

    /* copy vertex weights */
    for (k=0; k<ncon; k++)
      svwgt[mypart][snvtxs[mypart]*ncon+k] = vwgt[i*ncon+k];

    sadjrsum[mypart][snvtxs[mypart]] = sum;
    slabel[mypart][snvtxs[mypart]]   = label[i];
    sxadj[mypart][++snvtxs[mypart]]  = snedges[mypart];
  }

  for (mypart=0; mypart<2; mypart++) {
    iend = sxadj[mypart][snvtxs[mypart]];
    auxadjncy = sadjncy[mypart];
    for (i=0; i<iend; i++) 
      auxadjncy[i] = rename[auxadjncy[i]];
  }

  lgraph->nedges = snedges[0];
  rgraph->nedges = snedges[1];

  SetupGraph_tvwgt(lgraph);
  SetupGraph_tvwgt(rgraph);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->SplitTmr));

  *r_lgraph = lgraph;
  *r_rgraph = rgraph;

  WCOREPOP;
}


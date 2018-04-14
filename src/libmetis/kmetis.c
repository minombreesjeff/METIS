/*!
\file  
\brief The top-level routines for  multilevel k-way partitioning that minimizes
       the edge cut.

\date   Started 7/28/1997
\author George  
\author Copyright 1997-2011, Regents of the University of Minnesota 
\version\verbatim $Id: kmetis.c 10495 2011-07-06 16:04:45Z karypis $ \endverbatim
*/

#include "metislib.h"


/*************************************************************************/
/*! This function is the entry point for MCKMETIS */
/*************************************************************************/
int METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, 
          idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts, 
          real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, 
          idx_t *part)
{
  int sigrval=0, renumber=0;
  graph_t *graph;
  ctrl_t *ctrl;

  /* set up malloc cleaning code and signal catchers */
  if (!gk_malloc_init()) 
    return METIS_ERROR_MEMORY;

  gk_sigtrap();

  if ((sigrval = gk_sigcatch()) != 0)
    goto SIGTHROW;


  /* set up the run parameters */
  ctrl = SetupCtrl(METIS_OP_KMETIS, options, *ncon, *nparts, tpwgts, ubvec);
  if (!ctrl) {
    gk_siguntrap();
    return METIS_ERROR_INPUT;
  }

  /* if required, change the numbering to 0 */
  if (ctrl->numflag == 1) {
    Change2CNumbering(*nvtxs, xadj, adjncy);
    renumber = 1;
  }

  /* set up the graph */
  graph = SetupGraph(ctrl, *nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);

  /* set up multipliers for making balance computations easier */
  SetupKWayBalMultipliers(ctrl, graph);

  /* set various run parameters that depend on the graph */
  ctrl->CoarsenTo = gk_max((*nvtxs)/(20*gk_log2(*nparts)), 30*(*nparts));
  ctrl->nIparts   = (ctrl->CoarsenTo == 30*(*nparts) ? 1 : 3);

  /* take care contiguity requests for disconnected graphs */
  if (ctrl->contig && !IsConnected(graph, 0)) {
    printf("***The graph is not contiguous. The contig option is ignored.\n");
    ctrl->contig = 0;
  }
    
  /* allocate workspace memory */  
  AllocateWorkSpace(ctrl, graph);

  /* start the partitioning */
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->TotalTmr));

  *objval = MlevelKWayPartitioning(ctrl, graph, part);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->TotalTmr));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

  /* clean up */
  FreeCtrl(&ctrl);

SIGTHROW:
  /* if required, change the numbering back to 1 */
  if (renumber)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);

  gk_siguntrap();
  gk_malloc_cleanup(0);

  return metis_rcode(sigrval);
}


/*************************************************************************/
/*! This function computes a k-way partitioning of a graph that minimizes
    the specified objective function.

    \param ctrl is the control structure
    \param graph is the graph to be partitioned
    \param part is the vector that on return will store the partitioning

    \returns the objective value of the partitoning. The partitioning 
             itself is stored in the part vector.
*/
/*************************************************************************/
idx_t MlevelKWayPartitioning(ctrl_t *ctrl, graph_t *graph, idx_t *part)
{
  idx_t i, j, options[METIS_NOPTIONS], objval=0, curobj=0, bestobj=0;
  real_t *ubvec, curbal=0.0, bestbal=0.0;
  graph_t *cgraph;
  int status;

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NCUTS]   = ctrl->nIparts;
  options[METIS_OPTION_NITER]   = 10;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

  for (i=0; i<ctrl->ncuts; i++) {
    cgraph = CoarsenGraph(ctrl, graph);

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->InitPartTmr));
    AllocateKWayPartitionMemory(ctrl, cgraph);

    /* Release the work space */
    FreeWorkSpace(ctrl);

    ubvec = rmalloc(graph->ncon, "MlevelKWayPartitioning: ubvec");

    for (j=0; j<graph->ncon; j++) 
      ubvec[j] = (real_t)pow(ctrl->ubfactors[j], 1.0/log(ctrl->nparts));

    status = METIS_PartGraphRecursive(&cgraph->nvtxs, &cgraph->ncon, cgraph->xadj, 
                 cgraph->adjncy, cgraph->vwgt, cgraph->vsize, cgraph->adjwgt, 
                 &ctrl->nparts, ctrl->tpwgts, ubvec, options, &objval, cgraph->where);

    if (status != METIS_OK)
      gk_errexit(SIGERR, "Failed during initial partitioning\n");

    gk_free((void **)&ubvec, LTERM);

    /* Re-allocate the work space */
    AllocateWorkSpace(ctrl, graph);
    AllocateRefinementWorkSpace(ctrl, 2*cgraph->nedges);

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->InitPartTmr));
    IFSET(ctrl->dbglvl, METIS_DBG_IPART, 
        printf("Initial %"PRIDX"-way partitioning cut: %"PRIDX"\n", ctrl->nparts, objval));

    RefineKWay(ctrl, graph, cgraph);

    switch (ctrl->objtype) {
      case METIS_OBJTYPE_CUT:
        curobj = graph->mincut;
        break;

      case METIS_OBJTYPE_VOL:
        curobj = graph->minvol;
        break;

      default:
        gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
    }

    curbal = ComputeLoadImbalanceDiff(graph, ctrl->nparts, ctrl->pijbm, ctrl->ubfactors);

    if (i == 0 
        || (curbal <= 0.0005 && bestobj > curobj)
        || (bestbal > 0.0005 && curbal < bestbal)) {
      icopy(graph->nvtxs, graph->where, part);
      bestobj = curobj;
      bestbal = curbal;
    }

    FreeRData(graph);

    if (bestobj == 0)
      break;
  }

  FreeGraph(&graph);

  return bestobj;
}



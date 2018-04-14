/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mgrsetup.c
 *
 * This file contain various graph setting up routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: grsetup.c 10076 2011-06-03 15:36:39Z karypis $
 *
 */

#include <parmetislib.h>



/*************************************************************************/
/*! This function creates the graph from the user's inputs */
/*************************************************************************/
graph_t *SetUpGraph(ctrl_t *ctrl, idx_t ncon, idx_t *vtxdist, idx_t *xadj, 
               idx_t *vwgt, idx_t *adjncy, idx_t *adjwgt, idx_t *wgtflag)
{
  idx_t i, j;
  graph_t *graph;
  idx_t ltvwgts[MAXNCON];

  graph          = CreateGraph();
  graph->level   = 0;
  graph->gnvtxs  = vtxdist[ctrl->npes];
  graph->nvtxs   = vtxdist[ctrl->mype+1]-vtxdist[ctrl->mype];
  graph->ncon    = ncon;
  graph->nedges  = xadj[graph->nvtxs];
  graph->xadj    = xadj;
  graph->vwgt    = vwgt;
  graph->adjncy  = adjncy;
  graph->adjwgt  = adjwgt;
  graph->vtxdist = vtxdist;


  if (((*wgtflag)&2) == 0) 
    graph->vwgt = ismalloc(graph->nvtxs*ncon, 1, "Par_KMetis: vwgt");

  if (((*wgtflag)&1) == 0) 
    graph->adjwgt = ismalloc(graph->nedges, 1, "Par_KMetis: adjwgt");

  /* compute tvwgts */
  for (j=0; j<ncon; j++)
    ltvwgts[j] = 0;

  for (i=0; i<graph->nvtxs; i++)
    for (j=0; j<ncon; j++)
      ltvwgts[j] += graph->vwgt[i*ncon+j];

  for (j=0; j<ncon; j++) {
    ctrl->tvwgts[j] = GlobalSESum(ctrl, ltvwgts[j]);
    ctrl->invtvwgts[j] = 1.0/ctrl->tvwgts[j];
  }

  /* check for zero wgt constraints */
  for (i=0; i<ncon; i++) {
    /* ADD: take care of the case in which tvwgts is zero */
    if (ctrl->tvwgts[i] == 0) {
      rprintf(ctrl, "ERROR: sum weight for constraint %"PRIDX" is zero\n", i);
      gkMPI_Finalize();
      exit(-1);
    }
  }

  /* compute nvwgts */
  graph->nvwgt = rmalloc(graph->nvtxs*ncon, "graph->nvwgt");
  for (i=0; i<graph->nvtxs; i++) {
    for (j=0; j<ncon; j++)
      graph->nvwgt[i*ncon+j] = ctrl->invtvwgts[j]*graph->vwgt[i*ncon+j];
  }

  srand(ctrl->seed);

  return graph;
}


/*************************************************************************/
/*! This function sets the ctrl_t structure */
/*************************************************************************/
void SetUpCtrl(ctrl_t *ctrl, idx_t nparts, idx_t dbglvl, MPI_Comm comm)
{
  memset(ctrl, 0, sizeof(ctrl_t));

  MPI_Comm_dup(comm, &(ctrl->gcomm));
  gkMPI_Comm_rank(ctrl->gcomm, &ctrl->mype);
  gkMPI_Comm_size(ctrl->gcomm, &ctrl->npes);

  ctrl->dbglvl  = dbglvl;
  ctrl->nparts  = nparts;    /* Set the # of partitions is de-coupled from the # of domains */
  ctrl->comm    = ctrl->gcomm;
  ctrl->xyztype = XYZ_SPFILL;

  srand(ctrl->mype);
}


/*************************************************************************/
/*! Setups the global communicator and related info */
/*************************************************************************/
void SetUpComm(ctrl_t *ctrl, MPI_Comm comm)
{

  MPI_Comm_dup(comm, &(ctrl->gcomm));
  gkMPI_Comm_rank(ctrl->gcomm, &ctrl->mype);
  gkMPI_Comm_size(ctrl->gcomm, &ctrl->npes);

  ctrl->comm    = ctrl->gcomm;
}


/*************************************************************************
* This function changes the numbering from 1 to 0 or 0 to 1
**************************************************************************/
void ChangeNumbering(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *part, idx_t npes, idx_t mype, idx_t from)
{
  idx_t i, nvtxs;

  nvtxs = vtxdist[mype+1]-vtxdist[mype];

  if (from == 1) {  /* Change it from 1 to 0 */
    for (i=0; i<npes+1; i++)
      vtxdist[i]--;

    for (i=0; i<nvtxs+1; i++) 
      xadj[i]--;
    for (i=0; i<xadj[nvtxs]; i++) 
      adjncy[i]--;
  }
  else {  /* Change it from 0 to 1 */
    for (i=0; i<npes+1; i++) 
      vtxdist[i]++;

    for (i=0; i<xadj[nvtxs]; i++) 
      adjncy[i]++; 
    for (i=0; i<nvtxs+1; i++) 
      xadj[i]++; 

    for (i=0; i<nvtxs; i++)
      part[i]++;

  }
}


/*************************************************************************
* This function changes the numbering from 1 to 0 or 0 to 1
**************************************************************************/
void ChangeNumberingMesh(idx_t *elmdist, idx_t *eptr, idx_t *eind, 
                         idx_t *xadj, idx_t *adjncy, idx_t *part, 
			 idx_t npes, idx_t mype, idx_t from)
{
  idx_t i, nelms;

  nelms = elmdist[mype+1]-elmdist[mype];

  if (from == 1) {  /* Change it from 1 to 0 */
    for (i=0; i<npes+1; i++)
      elmdist[i]--;

    for (i=0; i<nelms+1; i++) 
      eptr[i]--;
    for (i=0; i<eptr[nelms]; i++) 
      eind[i]--;
  }
  else {  /* Change it from 0 to 1 */
    for (i=0; i<npes+1; i++) 
      elmdist[i]++;

    for (i=0; i<eptr[nelms]; i++) 
      eind[i]++;
    for (i=0; i<nelms+1; i++) 
      eptr[i]++;

    for (i=0; i<xadj[nelms]; i++) 
      adjncy[i]++; 
    for (i=0; i<nelms+1; i++) 
      xadj[i]++; 

    if (part != NULL)
      for (i=0; i<nelms; i++)
        part[i]++;
  }
}




/*************************************************************************
* This function randomly permutes the locally stored adjacency lists
**************************************************************************/
void GraphRandomPermute(graph_t *graph) 
{
  idx_t i, j, k, tmp;

  for (i=0; i<graph->nvtxs; i++) {
    for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++) {
      k = graph->xadj[i] + RandomInRange(graph->xadj[i+1]-graph->xadj[i]);
     gk_SWAP(graph->adjncy[j], graph->adjncy[k], tmp);
     gk_SWAP(graph->adjwgt[j], graph->adjwgt[k], tmp);
    }
  }
}


/*************************************************************************
* This function computes movement statistics for adaptive refinement
* schemes
**************************************************************************/
void ComputeMoveStatistics(ctrl_t *ctrl, graph_t *graph, idx_t *nmoved, idx_t *maxin, idx_t *maxout)
{
  idx_t i, j, nvtxs;
  idx_t *vwgt, *where;
  idx_t *lpvtxs, *gpvtxs;

  nvtxs = graph->nvtxs;
  vwgt = graph->vwgt;
  where = graph->where;

  lpvtxs = ismalloc(ctrl->nparts, 0, "ComputeMoveStatistics: lpvtxs");
  gpvtxs = ismalloc(ctrl->nparts, 0, "ComputeMoveStatistics: gpvtxs");

  for (j=i=0; i<nvtxs; i++) {
    lpvtxs[where[i]]++;
    if (where[i] != ctrl->mype)
      j++;
  }

  /* PrintVector(ctrl, ctrl->npes, 0, lpvtxs, "Lpvtxs: "); */

  gkMPI_Allreduce((void *)lpvtxs, (void *)gpvtxs, ctrl->nparts, IDX_T, MPI_SUM, ctrl->comm);

  *nmoved = GlobalSESum(ctrl, j);
  *maxout = GlobalSEMax(ctrl, j);
  *maxin = GlobalSEMax(ctrl, gpvtxs[ctrl->mype]-(nvtxs-j));

  gk_free((void **)&lpvtxs, (void **)&gpvtxs, LTERM);
}

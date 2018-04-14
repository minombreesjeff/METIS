/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * weird.c
 *
 * This file contain various graph setting up routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: weird.c 10361 2011-06-21 19:16:22Z karypis $
 *
 */

#include <parmetislib.h>



/*************************************************************************
* This function computes a partitioning of a small graph
**************************************************************************/
void PartitionSmallGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, h, ncon, nparts, npes, mype;
  idx_t moptions[METIS_NOPTIONS];
  idx_t me;
  idx_t *mypart;
  int lpecut[2], gpecut[2];
  graph_t *agraph;
  idx_t *sendcounts, *displs;
  real_t *gnpwgts, *lnpwgts;

  ncon   = graph->ncon;
  nparts = ctrl->nparts;

  gkMPI_Comm_size(ctrl->comm, &npes);
  gkMPI_Comm_rank(ctrl->comm, &mype);

  CommSetup(ctrl, graph);
  graph->where = imalloc(graph->nvtxs+graph->nrecv, "PartitionSmallGraph: where");
  agraph       = AssembleAdaptiveGraph(ctrl, graph);
  mypart       = imalloc(agraph->nvtxs, "mypart");

  METIS_SetDefaultOptions(moptions);
  moptions[METIS_OPTION_SEED] = ctrl->sync + mype;

  METIS_PartGraphKway(&agraph->nvtxs, &ncon, agraph->xadj, agraph->adjncy, 
        agraph->vwgt, NULL, agraph->adjwgt, &nparts, ctrl->tpwgts, NULL,
	moptions, &graph->mincut, mypart);

  lpecut[0] = graph->mincut;
  lpecut[1] = mype;
  gkMPI_Allreduce(lpecut, gpecut, 1, MPI_2INT, MPI_MINLOC, ctrl->comm);
  graph->mincut = gpecut[0];

  if (lpecut[1] == gpecut[1] && gpecut[1] != 0)
    gkMPI_Send((void *)mypart, agraph->nvtxs, IDX_T, 0, 1, ctrl->comm);
  if (lpecut[1] == 0 && gpecut[1] != 0)
    gkMPI_Recv((void *)mypart, agraph->nvtxs, IDX_T, gpecut[1], 1, ctrl->comm, &ctrl->status);

  sendcounts = imalloc(npes, "sendcounts");
  displs     = imalloc(npes, "displs");

  for (i=0; i<npes; i++) {
    sendcounts[i] = graph->vtxdist[i+1]-graph->vtxdist[i];
    displs[i] = graph->vtxdist[i];
  }

  gkMPI_Scatterv((void *)mypart, sendcounts, displs, IDX_T,
               (void *)graph->where, graph->nvtxs, IDX_T, 0, ctrl->comm);

  lnpwgts = graph->lnpwgts = rmalloc(nparts*ncon, "lnpwgts");
  gnpwgts = graph->gnpwgts = rmalloc(nparts*ncon, "gnpwgts");
  rset(nparts*ncon, 0, lnpwgts);
  for (i=0; i<graph->nvtxs; i++) {
    me = graph->where[i];
    for (h=0; h<ncon; h++)
      lnpwgts[me*ncon+h] += graph->nvwgt[i*ncon+h];
  }
  gkMPI_Allreduce((void *)lnpwgts, (void *)gnpwgts, nparts*ncon, REAL_T, MPI_SUM, ctrl->comm);
  gk_free((void**)&mypart, (void**)&sendcounts, (void**)&displs, LTERM);
  FreeGraph(agraph);

  return;
}



/*************************************************************************
* This function checks the inputs for the partitioning routines
**************************************************************************/
void CheckInputs(idx_t partType, idx_t npes, idx_t dbglvl, idx_t *wgtflag, idx_t *iwgtflag,
                 idx_t *numflag, idx_t *inumflag, idx_t *ncon, idx_t *incon, idx_t *nparts, 
		 idx_t *inparts, real_t *tpwgts, real_t **itpwgts, real_t *ubvec, 
		 real_t *iubvec, real_t *ipc2redist, real_t *iipc2redist, idx_t *options, 
		 idx_t *ioptions, idx_t *part, MPI_Comm *comm)
{
  idx_t i, j;
  idx_t doweabort, doiabort = 0;
  real_t tsum, *myitpwgts;
  idx_t mgcnums[5] = {-1, 2, 3, 4, 2};

  /**************************************/
  if (part == NULL) {
    doiabort = 1;
    IFSET(dbglvl, DBG_INFO, printf("ERROR: part array is set to NULL.\n"));
  }
  /**************************************/


  /**************************************/
  if (wgtflag == NULL) {
    *iwgtflag = 0;
    IFSET(dbglvl, DBG_INFO, printf("WARNING: wgtflag is NULL.  Using a value of 0.\n"));
  }
  else {
    *iwgtflag = *wgtflag;
  }
  /**************************************/


  /**************************************/
  if (numflag == NULL) {
    *inumflag = 0;
    IFSET(dbglvl, DBG_INFO, printf("WARNING: numflag is NULL.  Using a value of 0.\n"));
  }
  else {
    if (*numflag != 0 && *numflag != 1) {
      IFSET(dbglvl, DBG_INFO, printf("WARNING: bad value for numflag %"PRIDX".  Using a value of 0.\n", *numflag));
      *inumflag = 0;
    }
    else {
      *inumflag = *numflag;
    }
  }
  /**************************************/


  /**************************************/
  if (ncon == NULL) {
    *incon = 1;
    IFSET(dbglvl, DBG_INFO, printf("WARNING: ncon is NULL.  Using a value of 1.\n"));
  }
  else {
    if (*ncon < 1 || *ncon > MAXNCON) {
      IFSET(dbglvl, DBG_INFO, printf("WARNING: bad value for ncon %"PRIDX".  Using a value of 1.\n", *ncon));
      *incon = 1;
    }
    else {
      *incon = *ncon;
    }
  }
  /**************************************/


  /**************************************/
  if (nparts == NULL) {
    *inparts = npes;
    IFSET(dbglvl, DBG_INFO, printf("WARNING: nparts is NULL.  Using a value of %"PRIDX".\n", npes));
  }
  else {
    if (*nparts < 1 || *nparts > MAX_NPARTS) {
      IFSET(dbglvl, DBG_INFO, printf("WARNING: bad value for nparts %"PRIDX".  Using a value of %"PRIDX".\n", *nparts, npes));
      *inparts = npes;
    }
    else {
      *inparts = *nparts;
    }
  }
  /**************************************/


  /**************************************/
  myitpwgts = *itpwgts = rmalloc((*inparts)*(*incon), "CheckInputs: itpwgts");
  if (tpwgts == NULL) {
    rset((*inparts)*(*incon), 1.0/(real_t)(*inparts), myitpwgts);
    IFSET(dbglvl, DBG_INFO, printf("WARNING: tpwgts is NULL.  Setting all array elements to %.3"PRREAL".\n", 1.0/(real_t)(*inparts)));
  }
  else {
    for (i=0; i<*incon; i++) {
      tsum = 0.0;
      for (j=0; j<*inparts; j++) {
        tsum += tpwgts[j*(*incon)+i];
      } 
      if (fabs(1.0-tsum) < SMALLFLOAT)
        tsum = 1.0;
      for (j=0; j<*inparts; j++)
       myitpwgts[j*(*incon)+i] = tpwgts[j*(*incon)+i] / tsum;
    }
  }
  /**************************************/


  /**************************************/
  if (ubvec == NULL) {
    rset(*incon, 1.05, iubvec);
    IFSET(dbglvl, DBG_INFO, printf("WARNING: ubvec is NULL.  Setting all array elements to 1.05.\n"));
  }
  else {
    for (i=0; i<*incon; i++) {
      if (ubvec[i] < 1.0 || ubvec[i] > (real_t)(*inparts)) {
        iubvec[i] = 1.05;
        IFSET(dbglvl, DBG_INFO, printf("WARNING: bad value for ubvec[%"PRIDX"]: %.3"PRREAL".  Setting value to 1.05.[%"PRIDX"]\n", i, ubvec[i], *inparts));
      }
      else {
        iubvec[i] = ubvec[i];
      }
    }
  }
  /**************************************/


  /**************************************/
  if (partType == ADAPTIVE_PARTITION) {
    if (ipc2redist != NULL) {
      if (*ipc2redist < SMALLFLOAT || *ipc2redist > 1000000.0) {
        IFSET(dbglvl, DBG_INFO, printf("WARNING: bad value for ipc2redist %.3"PRREAL".  Using a value of 1000.\n", *ipc2redist));
        *iipc2redist = 1000.0;
      }
      else {
        *iipc2redist = *ipc2redist;
      }
    }
    else {
      IFSET(dbglvl, DBG_INFO, printf("WARNING: ipc2redist is NULL.  Using a value of 1000.\n"));
      *iipc2redist = 1000.0;
    }
  }
  /**************************************/


  /**************************************/
  if (options == NULL) {
    ioptions[0] = 0;
    IFSET(dbglvl, DBG_INFO, printf("WARNING: options is NULL.  Using defaults\n"));
  }
  else {
    ioptions[0] = options[0];
    ioptions[1] = options[1];
    ioptions[2] = options[2];
    if (partType == ADAPTIVE_PARTITION || partType == REFINE_PARTITION)
      ioptions[3] = options[3];
  }
  /**************************************/


  /**************************************/
  if (comm == NULL) {
    IFSET(dbglvl, DBG_INFO, printf("ERROR: comm is NULL.  Aborting\n"));
    abort();
  }
  else {
    gkMPI_Allreduce((void *)&doiabort, (void *)&doweabort, 1, IDX_T, MPI_MAX, *comm);
    if (doweabort > 0)
      abort();
  }
  /**************************************/

}



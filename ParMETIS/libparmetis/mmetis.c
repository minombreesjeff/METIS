/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mmetis.c
 *
 * This is the entry point of ParMETIS_V3_PartMeshKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: mmetis.c 10052 2011-06-01 22:29:57Z karypis $
 *
 */

#include <parmetislib.h>


/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel mesh partitionioner. 
* This function assumes nothing about the mesh distribution.
* It is the general case.
************************************************************************************/
void ParMETIS_V3_PartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt, 
                 idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommonnodes, idx_t *nparts, 
		 real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
		 MPI_Comm *comm)
{
  idx_t i, nvtxs, nedges, gnedges, npes, mype;
  idx_t *xadj, *adjncy;
  timer TotalTmr, Mesh2DualTmr, ParMETISTmr;
  ctrl_t ctrl;

  /********************************/
  /* Try and take care bad inputs */
  /********************************/
  if (elmdist == NULL || eptr == NULL || eind == NULL || wgtflag == NULL || 
      numflag == NULL || ncon == NULL || ncommonnodes == NULL || nparts == NULL ||
      tpwgts == NULL || ubvec == NULL || options == NULL || edgecut == NULL || 
      part == NULL || comm == NULL) {
    printf("ERROR: One or more required parameters is NULL. Aborting.\n");
    abort();
  }
  if (((*wgtflag)&2) && elmwgt == NULL) {
    printf("ERROR: elmwgt == NULL when vertex weights were specified. Aborting.\n");
    abort();
  }


  SetUpCtrl(&ctrl, *nparts, (options[0] == 1 ? options[PMV3_OPTION_DBGLVL] : 0), *comm);

  npes = ctrl.npes;
  mype = ctrl.mype;

  cleartimer(TotalTmr);
  cleartimer(Mesh2DualTmr);
  cleartimer(ParMETISTmr);

  gkMPI_Barrier(ctrl.comm);
  starttimer(TotalTmr);
  starttimer(Mesh2DualTmr);

  ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, numflag, ncommonnodes, &xadj, &adjncy, 
      &(ctrl.comm));

  if (ctrl.dbglvl&DBG_INFO) {
    nvtxs = elmdist[mype+1]-elmdist[mype];
    nedges = xadj[nvtxs] + (*numflag == 0 ? 0 : -1);
    rprintf(&ctrl, "Completed Dual Graph -- Nvtxs: %"PRIDX", Nedges: %"PRIDX" \n", 
            elmdist[npes], GlobalSESum(&ctrl, nedges));
  }

  gkMPI_Barrier(ctrl.comm);
  stoptimer(Mesh2DualTmr);


  /***********************/
  /* Partition the graph */
  /***********************/
  starttimer(ParMETISTmr);

  ParMETIS_V3_PartKway(elmdist, xadj, adjncy, elmwgt, NULL, wgtflag, numflag, ncon, 
      nparts, tpwgts, ubvec, options, edgecut, part, &(ctrl.comm));

  gkMPI_Barrier(ctrl.comm);
  stoptimer(ParMETISTmr);
  stoptimer(TotalTmr);

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimer(&ctrl, Mesh2DualTmr,	"   Mesh2Dual"));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimer(&ctrl, ParMETISTmr,	"    ParMETIS"));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimer(&ctrl, TotalTmr,	"       Total"));

  gk_free((void **)&xadj, (void **)&adjncy, LTERM);

  FreeCtrl(&ctrl);

  return;
}

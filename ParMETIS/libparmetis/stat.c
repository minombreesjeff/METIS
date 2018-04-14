/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * stat.c
 *
 * This file computes various statistics
 *
 * Started 7/25/97
 * George
 *
 * $Id: stat.c 10367 2011-06-22 13:45:44Z karypis $
 *
 */

#include <parmetislib.h>



/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
void ComputeSerialBalance(ctrl_t *ctrl, graph_t *graph, idx_t *where, real_t *ubvec)
{
  idx_t i, j, nvtxs, ncon, nparts;
  idx_t *pwgts, *tvwgts, *vwgt;
  real_t *tpwgts, maximb;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  vwgt   = graph->vwgt;
  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  pwgts = ismalloc(nparts*ncon, 0, "pwgts");
  tvwgts = ismalloc(ncon, 0, "tvwgts");

  for (i=0; i<graph->nvtxs; i++) {
    for (j=0; j<ncon; j++) {
      pwgts[where[i]*ncon+j] += vwgt[i*ncon+j];
      tvwgts[j] += vwgt[i*ncon+j];
    }
  }

  /* The +1 in the following code is to deal with bad cases of tpwgts[i*ncon+j] == 0 */
  for (j=0; j<ncon; j++) {
    maximb = 0.0;
    for (i=0; i<nparts; i++)
      maximb =gk_max(maximb, (1.0+(real_t)pwgts[i*ncon+j])/(1.0+(tpwgts[i*ncon+j]*(real_t)tvwgts[j])));
    ubvec[j] = maximb;
  }

  gk_free((void **)&pwgts, (void **)&tvwgts, LTERM);
}


/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
void ComputeParallelBalance(ctrl_t *ctrl, graph_t *graph, idx_t *where, real_t *ubvec)
{
  idx_t i, j, nvtxs, ncon, nparts;
  real_t *nvwgt, *lnpwgts, *gnpwgts;
  real_t *tpwgts, maximb;
  real_t lminvwgts[MAXNCON], gminvwgts[MAXNCON];

  ncon   = graph->ncon;
  nvtxs  = graph->nvtxs;
  nvwgt  = graph->nvwgt;
  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  lnpwgts = rmalloc(nparts*ncon, "CPB: lnpwgts");
  gnpwgts = rmalloc(nparts*ncon, "CPB: gnpwgts");
  rset(nparts*ncon, 0.0, lnpwgts);
  rset(ncon, 1.0, lminvwgts);

  for (i=0; i<nvtxs; i++) {
    for (j=0; j<ncon; j++) {
      lnpwgts[where[i]*ncon+j] += nvwgt[i*ncon+j];

      /* The following is to deal with tpwgts[] that are 0.0 for certain partitions/constraints */
      lminvwgts[j] = (nvwgt[i*ncon+j] > 0.0 && lminvwgts[j] > nvwgt[i*ncon+j] ? nvwgt[i*ncon+j] : lminvwgts[j]);
    }
  }

  gkMPI_Allreduce((void *)(lnpwgts), (void *)(gnpwgts), nparts*ncon, REAL_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce((void *)(lminvwgts), (void *)(gminvwgts), ncon, REAL_T, MPI_MIN, ctrl->comm);

  /* The +gminvwgts[j] in the following code is to deal with bad cases of tpwgts[i*ncon+j] == 0 */
  for (j=0; j<ncon; j++) {
    maximb = 0.0;
    for (i=0; i<nparts; i++)
      maximb =gk_max(maximb, (gminvwgts[j]+gnpwgts[i*ncon+j])/(gminvwgts[j]+tpwgts[i*ncon+j]));
    ubvec[j] = maximb;
  }

  gk_free((void **)&lnpwgts, (void **)&gnpwgts, LTERM);

  return;
}


/*************************************************************************
* This function prints a matrix
**************************************************************************/
void Mc_PrintThrottleMatrix(ctrl_t *ctrl, graph_t *graph, real_t *matrix)
{
  idx_t i, j;

  for (i=0; i<ctrl->npes; i++) {
    if (i == ctrl->mype) {
      for (j=0; j<ctrl->npes; j++)
        printf("%.3"PRREAL" ", matrix[j]);
      printf("\n");
      fflush(stdout);
    }
    gkMPI_Barrier(ctrl->comm);
  }

  if (ctrl->mype == 0) {
    printf("****************************\n");
    fflush(stdout);
  }
  gkMPI_Barrier(ctrl->comm);

  return;
}


/*************************************************************************
*  This function computes stats for refinement
**************************************************************************/
void Mc_ComputeRefineStats(ctrl_t *ctrl, graph_t *graph, real_t *ubvec)
{
  idx_t h, i, j, k;
  idx_t nvtxs, ncon;
  idx_t *xadj, *adjncy, *adjwgt, *where;
  real_t *nvwgt, *lnpwgts, *gnpwgts;
  idx_t mype = ctrl->mype, nparts = ctrl->nparts;
  idx_t *gborder, *border, *gfrom, *from, *gto, *to, *connect, *gconnect;
  idx_t gain[20] = {0}, ggain[20];
  idx_t lnborders, gnborders;
  idx_t bestgain, pmoves, gpmoves, other;
  real_t tpwgts[MAXNCON], badmaxpwgt[MAXNCON];
  idx_t HIST_FACTOR = graph->level + 1;
  ckrinfo_t *rinfo;
  cnbr_t *mynbrs;

  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  where   = graph->where;
  lnpwgts = graph->lnpwgts;
  gnpwgts = graph->gnpwgts;
  rinfo   = graph->ckrinfo;

  connect  = ismalloc(nparts*nparts, 0, "CRS: connect");
  gconnect = imalloc(nparts*nparts, "CRS: gconnect");
  border   = ismalloc(nparts, 0, "CRS: border");
  gborder  = imalloc(nparts, "CRS: gborder");
  from     = ismalloc(nparts, 0, "CRS: from");
  gfrom    = imalloc(nparts, "CRS: gfrom");
  to       = ismalloc(nparts, 0, "CRS: to");
  gto      = imalloc(nparts, "CRS: gto");

  for (h=0; h<ncon; h++) {
    tpwgts[h] = rsum(nparts, gnpwgts+h, ncon)/nparts;
    badmaxpwgt[h] = ubvec[h]*tpwgts[h];
  }

  if (mype == 0) printf("******************************\n");
  if (mype == 0) printf("******************************\n");

  /***************************************/
  if (mype == 0) {
    printf("subdomain weights:\n");
    for (h=0; h<ncon; h++) {
      for (i=0; i<nparts; i++)
        printf("%9.3"PRREAL" ", gnpwgts[i*ncon+h]);
      printf("\n");
    }
    printf("\n");
  }

  /***************************************/
  if (mype == 0) {
    printf("subdomain imbalance:\n");
    for (h=0; h<ncon; h++) {
      for (i=0; i<nparts; i++)
        printf("%9.3"PRREAL" ", gnpwgts[i*ncon+h] * (real_t)(nparts));
      printf("\n");
    }
    printf("\n");
  }

  /***************************************/
  for (i=0; i<nparts; i++)
    connect[i*nparts+i] = -1;

  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) {
        connect[where[i]*nparts+where[adjncy[j]]] = 1;
        connect[where[adjncy[j]]*nparts+where[i]] = 1;
      }
    }
  }

  gkMPI_Reduce((void *)connect, (void *)gconnect, nparts*nparts, IDX_T, MPI_MAX, 0, ctrl->comm);
  if (mype == 0) { 
    printf("connectivity\n");
    for (i=0; i<nparts; i++) {
      printf("%"PRIDX": ", i);
      for (j=0; j<nparts; j++)
        printf("%9"PRIDX" ", gconnect[i*nparts+j]);
      printf("\n");
    }
    printf("\n");
  }

  /***************************************/
  lnborders = 0;
  for (i=0; i<nvtxs; i++) {
    if (rinfo[i].nnbrs > 0) {
      lnborders++;
      border[where[i]]++;
    } 
  }

  gkMPI_Reduce((void *)border, (void *)gborder, nparts, IDX_T, MPI_SUM, 0, ctrl->comm);
  gnborders = GlobalSESum(ctrl, lnborders);
  if (mype == 0) {
    printf("number of borders: %"PRIDX"\n", gnborders);
    for (i=0; i<nparts; i++)
      printf("%9"PRIDX" ", gborder[i]);
    printf("\n\n");
  }

  /***************************************/
  pmoves = 0;
  for (i=0; i<nvtxs; i++) {
    nvwgt  = graph->nvwgt+i*ncon;
    mynbrs = ctrl->cnbrpool + rinfo[i].inbr;

    for (j=0; j<rinfo[i].nnbrs; j++) {
      other = mynbrs[j].pid;
      for (h=0; h<ncon; h++) {
        if (gnpwgts[other*ncon+h]+nvwgt[h] > badmaxpwgt[h])
          break;
      }

      if (h == ncon)
        break;
    }

    if (j < rinfo[i].nnbrs) {
      pmoves++;
      from[where[i]]++;
      to[other]++;
      for (k=j+1; k<rinfo[i].nnbrs; k++) {
        other = mynbrs[j].pid;
        for (h=0; h<ncon; h++)
          if (gnpwgts[other*ncon+h]+nvwgt[h] > badmaxpwgt[h])
            break;

        if (h == ncon) {
          pmoves++;
          from[where[i]]++;
          to[other]++;
        }
      }
    }
  }

  gpmoves = GlobalSESum(ctrl, pmoves);
  gkMPI_Reduce((void *)from, (void *)gfrom, nparts, IDX_T, MPI_SUM, 0, ctrl->comm);
  gkMPI_Reduce((void *)to, (void *)gto, nparts, IDX_T, MPI_SUM, 0, ctrl->comm);

  if (mype == 0) {
    printf("possible moves: %"PRIDX"\n", gpmoves);
    printf("from   ");
    for (i=0; i<nparts; i++) {
      printf("%9"PRIDX" ", gfrom[i]);
    }
    printf("\n");
    printf("to     ");
    for (i=0; i<nparts; i++) {
      printf("%9"PRIDX" ", gto[i]);
    }
    printf("\n\n");
  }

  /***************************************/
  for (i=0; i<nvtxs; i++) {
    if (rinfo[i].nnbrs > 0) {
      mynbrs = ctrl->cnbrpool + rinfo[i].inbr;

      bestgain = mynbrs[0].ed-rinfo[i].id;
      for (j=1; j<rinfo[i].nnbrs; j++)
        bestgain = gk_max(bestgain, mynbrs[j].ed-rinfo[i].id);

      if (bestgain / HIST_FACTOR >= 10) {
        gain[19]++;
        continue;
      }

      if (bestgain / HIST_FACTOR < -10) {
        gain[0]++;
        continue;
      }

      gain[(bestgain/HIST_FACTOR)+10]++;
    }
  }

  gkMPI_Reduce((void *)gain, (void *)ggain, 20, IDX_T, MPI_SUM, 0, ctrl->comm);
  if (mype == 0) {
    printf("gain histogram (buckets of %"PRIDX")\n", HIST_FACTOR);
    for (i=0; i<20; i++) {
      if (i == 10 || i == 11)
        printf("    ");
      printf("%"PRIDX" ", ggain[i]);
    }
    printf("\n\n");
  }


  /***************************************/
  if (mype == 0) printf("******************************\n");
  if (mype == 0) printf("******************************\n");

  gk_free((void **)&gconnect, (void **)&connect, (void **)&gborder, (void **)&border, (void **)&gfrom, (void **)&from, (void **)&gto, (void **)&to, LTERM);
  return;
}

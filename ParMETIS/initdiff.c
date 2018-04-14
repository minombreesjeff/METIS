/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initdiff.c 
 *
 * This file contains code for the initial directed diffusion at the coarsest
 * graph
 *
 * Started 5/19/97, Kirk, George
 *
 * $Id: initdiff.c,v 1.2 1997/07/18 00:32:06 karypis Exp $
 *
 */

#include <par_kmetis.h>

#define MAXPDEGREE	50

/*************************************************************************
* This function is the entry point of the directed diffusion algorithm
* for the coarsest graph.  This function assembles the graph to all the 
* processors and preceed serially.
**************************************************************************/
void InitDiffusion(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, j, lpecut[2], gpecut[2];
  idxtype *vtxdist;
  GraphType *agraph;
  int *sendcounts, *displs;

  ASSERT(ctrl, ctrl->npes == ctrl->nparts);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  agraph = AssembleAdaptiveGraph(ctrl, graph, wspace);

  /* myprintf(ctrl, "Assembled graph: %d %d\n", agraph->nvtxs, agraph->nedges); */
  MPI_Barrier(ctrl->comm);

  vtxdist = graph->vtxdist;

  agraph->where = idxmalloc(agraph->nvtxs, "InitDiffusion: agraph->where");
  for (i=0; i<ctrl->npes; i++) {
    for (j=vtxdist[i]; j<vtxdist[i+1]; j++)
      agraph->where[j] = i;
  }

  KWay_InitialDiffuser(ctrl, agraph, ctrl->nparts, .8);

  /* Determine which PE got the minimum cut */
  lpecut[0] = agraph->mincut;
  MPI_Comm_rank(ctrl->comm, lpecut+1);
  MPI_Allreduce(lpecut, gpecut, 1, MPI_2INT, MPI_MINLOC, ctrl->gcomm);

  /* myprintf(ctrl, "Mincut: %d, GMincut: %d\n", agraph->mincut, gpecut[0]);  */

  graph->where = idxmalloc(graph->nvtxs+graph->nrecv, "InitPartition: where");
  sendcounts = imalloc(ctrl->npes, "InitPartitionNew: sendcounts");
  displs = imalloc(ctrl->npes, "InitPartitionNew: displs");

  for (i=0; i<ctrl->npes; i++) {
    sendcounts[i] = vtxdist[i+1]-vtxdist[i];
    displs[i] = vtxdist[i];
  }

  MPI_Scatterv((void *)agraph->where, sendcounts, displs, IDX_DATATYPE, 
               (void *)graph->where, graph->nvtxs, IDX_DATATYPE, gpecut[1], ctrl->comm);

  GKfree(&sendcounts, &displs, -1);

  FreeGraph(agraph);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

}



/*************************************************************************
* This function assembles the graph into a single processor
**************************************************************************/
GraphType *AssembleAdaptiveGraph(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, j, k, l, gnvtxs, nvtxs, gnedges, nedges, firstvtx, gsize;
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *vtxdist, *edgedist, *imap;
  idxtype *axadj, *aadjncy, *aadjwgt, *avwgt, *avsize, *alabel;
  idxtype *mygraph, *ggraph;
  int *recvcounts, *displs, mysize;
  GraphType *agraph;

  gnvtxs = graph->gnvtxs;
  nvtxs = graph->nvtxs;
  nedges = graph->xadj[nvtxs];
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  vtxdist = graph->vtxdist;
  imap = graph->imap;

  /* Determine the # of idxtype to receive from each processor */
  recvcounts = imalloc(ctrl->npes, "AssembleGraph: recvcounts");
  mysize = 3*nvtxs + 2*nedges;
  MPI_Allgather((void *)(&mysize), 1, MPI_INT, (void *)recvcounts, 1, MPI_INT, ctrl->comm);
  
  displs = imalloc(ctrl->npes+1, "AssembleGraph: displs");
  displs[0] = 0;
  for (i=1; i<ctrl->npes+1; i++) 
    displs[i] = displs[i-1] + recvcounts[i-1];

  /* Construct the one-array storage format of the assembled graph */
  mygraph = (mysize <= wspace->maxcore ? wspace->core : idxmalloc(mysize, "AssembleGraph: mygraph"));
  for (k=i=0; i<nvtxs; i++) {
    mygraph[k++] = xadj[i+1]-xadj[i];
    mygraph[k++] = vwgt[i];
    mygraph[k++] = vsize[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      mygraph[k++] = imap[adjncy[j]];
      mygraph[k++] = adjwgt[j];
    }
  }
  ASSERT(ctrl, mysize == k);

  /* Assemble the entire graph */
  gsize = displs[ctrl->npes];
  ggraph = (gsize <= wspace->maxcore-mysize ? wspace->core+mysize : idxmalloc(gsize, "AssembleGraph: ggraph"));
  MPI_Allgatherv((void *)mygraph, mysize, IDX_DATATYPE, (void *)ggraph, recvcounts, displs, IDX_DATATYPE, ctrl->comm);

  /* MPI_Bcast((void *)ggraph, gsize, IDX_DATATYPE, 0, ctrl->comm); */

  GKfree(&recvcounts, &displs, -1);
  if (mysize > wspace->maxcore)
    free(mygraph);

  agraph = CreateGraph();
  agraph->maxvwgt = graph->maxvwgt;
  agraph->nvtxs = gnvtxs;
  agraph->nedges = gnedges = (gsize-3*gnvtxs)/2;

  /* Allocate memory for the assembled graph */
  axadj = agraph->xadj = idxmalloc(gnvtxs+1, "AssembleGraph: axadj");
  avwgt = agraph->vwgt = idxmalloc(gnvtxs, "AssembleGraph: avwgt");
  avsize = agraph->vsize = idxmalloc(gnvtxs, "AssembleGraph: avsize");
  aadjncy = agraph->adjncy = idxmalloc(gnedges, "AssembleGraph: adjncy");
  aadjwgt = agraph->adjwgt = idxmalloc(gnedges, "AssembleGraph: adjwgt");
  alabel = agraph->label = idxmalloc(gnvtxs, "AssembleGraph: alabel");

  for (k=j=i=0; i<gnvtxs; i++) {
    axadj[i] = ggraph[k++];
    avwgt[i] = ggraph[k++];
    avsize[i] = ggraph[k++];
    for (l=0; l<axadj[i]; l++) {
      aadjncy[j] = ggraph[k++];
      aadjwgt[j] = ggraph[k++];
      j++;
    }
  }

  /* Now fix up the received graph */
  MAKECSR(i, gnvtxs, axadj);

  for (i=0; i<gnvtxs; i++)
    alabel[i] = i;

  if (gsize > wspace->maxcore-mysize)
    free(ggraph);

/*
  if (mysize > wspace->nlarge || displs[ctrl->npes] > 2*wspace->nlarge)
    myprintf(ctrl, "Had to allocate memory %d %d %d\n", mysize, displs[ctrl->npes], wspace->nlarge);
*/

  return agraph;
}



/*************************************************************************
* This function performs a k-way directed diffusion
**************************************************************************/
void KWay_InitialDiffuser(CtrlType *ctrl, GraphType *graph, int nparts, float SuppressFactor)
{
  int ii, i, j, k, l, nvtxs, from, other, me, minpwgt, maxpwgt, oldcut, done;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *vsize;
  idxtype *where, *pwgts, *perm;
  idxtype *transfer, *rowptr, *colind;
  int supped;
  float suppression;
  RInfoType *rinfo, *myrinfo;
  EdgeType *alldegrees, *degrees;
  float sum, mean, *solution, *load, *value;
  int npasses = 20;
  int tvsize, tvwgt;
  float balance, oldbalance;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  /* Compute the partition parameters */
  rinfo = graph->rinfo = (RInfoType *)GKmalloc(sizeof(RInfoType)*nvtxs, "KWay_InitialDiffuser: rinfo");
  pwgts = graph->lpwgts = idxsmalloc(nparts, 0, "KWay_InitialDiffuser: pwgts");
  alldegrees = (EdgeType *)GKmalloc(sizeof(EdgeType)*graph->nedges, "KWay_InitialDiffuser: alldegrees");

  graph->mincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    myrinfo = rinfo+i;

    pwgts[me] += vwgt[i];

    myrinfo->degrees = alldegrees + xadj[i];
    myrinfo->ndegrees = myrinfo->id = myrinfo->ed = 0;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]])
        myrinfo->id += adjwgt[j];
      else
        myrinfo->ed += adjwgt[j];
    }

    if (myrinfo->ed > 0) {  /* Time to do some serious work */
      graph->mincut += myrinfo->ed;
      degrees = myrinfo->degrees;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (degrees[k].edge == other) {
              degrees[k].ewgt += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            degrees[k].edge = other;
            degrees[k].ewgt = adjwgt[j];
            myrinfo->ndegrees++;
          }
        }
      }
    }
  }
  graph->mincut = graph->mincut/2;

  /* printf("Cut: %d, %d %d\n", graph->mincut, pwgts[0], pwgts[1]); */

  tvwgt = idxsum(nvtxs, vwgt);
  tvsize = idxsum(nvtxs, vsize);
  mean = tvwgt/nparts;
  suppression = ((float)(tvwgt) / (float)(tvsize)) * SuppressFactor;

  rowptr = idxmalloc(nparts+1, "KWay_InitialDiffuser: rowptr");
  colind = idxmalloc(MAXPDEGREE*nparts, "KWay_InitialDiffuser: colind");
  value = fmalloc(MAXPDEGREE*nparts, "KWay_InitialDiffuser: value");

  transfer = idxmalloc(nparts*nparts, "KWay_InitialDiffuser: transfer");
  load = fmalloc(nparts, "KWay_InitialDiffuser: load");
  solution = fmalloc(nparts, "KWay_InitialDiffuser: solution");

  perm = idxmalloc(nvtxs, "KWay_InitialDiffuser: perm");
  for (i=0; i<nvtxs; i++)
    perm[i] = i;

  minpwgt = pwgts[idxamin(nparts, pwgts)];
  maxpwgt = pwgts[idxamax(nparts, pwgts)];
  balance = 1.0*nparts*maxpwgt/(1.0*tvwgt);

  done = 0;
  for (l=0; l<npasses; l++) {
    oldcut = graph->mincut;
    oldbalance = balance;

    /* Set-up and solve the diffusion equations */
    for (j=0; j<nparts; j++) 
      load[j] = pwgts[j] - mean;

    setupLaplace(graph, nparts, rowptr, colind, value, transfer);
    ConjGrad(nparts, rowptr, colind, value, load, solution, .001);

    /* printf("Difussion solution: %lf %lf %lf %lf\n", solution[0], solution[1], solution[2], solution[3]); */

    idxset(nparts*nparts, 0, transfer);
    for (j=0; j<nparts; j++) {
      for (k=rowptr[j]+1; k<rowptr[j+1]; k++)
        transfer[j*nparts+colind[k]] = (int) (solution[j] - solution[colind[k]]);
    }

    FastRandomPermute(nvtxs, perm, 0);

    supped = 0;
    for (ii=0; ii<nvtxs; ii++) {
      i = perm[ii];
      if ((float)vwgt[i] / (float)(vsize[i]) > suppression) {
        if (rinfo[i].ndegrees > 0) {
          degrees = rinfo[i].degrees;
          from = where[i];
          if (pwgts[from] == vwgt[i])
            continue;  /* We do not want to transfer everything out of a domain */

          for (k=0; k<rinfo[i].ndegrees; k++)
            if (transfer[from*nparts+degrees[k].edge] > (.9 * vwgt[i]))
              break;

          if (k < rinfo[i].ndegrees) {
            for (j=k+1; j<rinfo[i].ndegrees; j++) {
              if (transfer[from*nparts+degrees[j].edge] > (.9 * vwgt[i]))
                if (degrees[j].ewgt > degrees[k].ewgt ||
                   (degrees[j].ewgt == degrees[k].ewgt &&
                      pwgts[degrees[j].edge] < pwgts[degrees[k].edge]))
                  k = j;
            }
            graph->mincut -= (degrees[k].ewgt-rinfo[i].id);
            transfer[from*nparts+degrees[k].edge] -= vwgt[i];
            Ser_KWayUpdateDegrees(graph, i, k);
          }
        }
      }
      else
        supped++;
    }

    minpwgt = pwgts[idxamin(nparts, pwgts)];
    maxpwgt = pwgts[idxamax(nparts, pwgts)];
    balance = 1.0*nparts*maxpwgt/(1.0*tvwgt);

    if (balance < UNBALANCE_FRACTION)
      done = 1;

    if (supped > (int) (.8 * (float)(nvtxs))) {
      printf("Too much suppression\n");
      suppression /= 1.2;
    }
    else {
      if (fabs(balance-oldbalance) < 0.01)
        done = 1;
    }

    if (GlobalSESum(ctrl, done) > 0)
      break;
  }
  if (!done)
    graph->mincut = 10000000;

  GKfree(&transfer, &load, &solution, &perm, &alldegrees, &value, &colind, &rowptr, -1);

  /* 
  printf("Directed Forced edge-cut is %d   Max pwgt is %d   Min pwgt is %d, Balance: %5.4f, Iter: %d\n", 
         graph->mincut, pwgts[idxamax(nparts, pwgts)], pwgts[idxamin(nparts, pwgts)], balance, l);
  {
    int i, j, cut = 0;
    for (i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        if (where[i] != where[adjncy[j]])
          cut += adjwgt[j];
    }
    printf("Actual cut: %d\n", cut/2);
  }
  */

  return;
}


/*************************************************************************
* This function sets up the Laplacian matrix
**************************************************************************/
void setupLaplace(GraphType *graph, int nparts, idxtype *rowptr, idxtype *colind, float *value, idxtype *tempgraph)
{
  int i, j, nvtxs, nedges, me, other;
  idxtype *where, *tptr;
  RInfoType *myrinfo;

  nvtxs = graph->nvtxs;
  where = graph->where;

  idxset(nparts*nparts, 0, tempgraph);

  for (i=0; i<nparts; i++) 
    rowptr[i] = i*MAXPDEGREE+1;

  for (i=0; i<nvtxs; i++) {
    myrinfo = graph->rinfo + i;
    me = where[i];
    tptr = tempgraph + me*nparts;

    for (j=0; j<myrinfo->ndegrees; j++) {
      other = myrinfo->degrees[j].edge;
      if (tptr[other] != 1) {
        tptr[other] = 1;
        colind[rowptr[me]++] = other;
      }
    }
  }

  for (i=0; i<nparts; i++) {
    if (rowptr[i]-i*MAXPDEGREE > MAXPDEGREE)
      printf("Preallocated memory is not enough!\n");
  }

  nedges=0;
  tptr = tempgraph;
  for (i=0; i<nparts; i++) {
    colind[nedges] = i;
    value[nedges++] = rowptr[i]-(i*MAXPDEGREE+1);
    for (j=i*MAXPDEGREE+1; j<rowptr[i]; j++) {
      tptr[colind[j]] = 0;  /* Zero tempgraph out */
      colind[nedges] = colind[j];
      value[nedges++] = -1;
    }
    rowptr[i] = nedges;
    tptr += nparts;
  }
  for (i=nparts; i>0; i--)
    rowptr[i] = rowptr[i-1];
  rowptr[0] = 0;

} 


/*************************************************************************
* This function implements the CG solver used during the directed diffusion
**************************************************************************/
int ConjGrad(int n, idxtype *rowptr, idxtype *colind, float *value, float *b, float *x, float tol)
{
  int i, j, k;
  float *p, *r, *q, *z, *M; 
  float alpha, beta, rho, rho_1, error, bnrm2;

  /* Initial Setup */   
  p = fmalloc(n, "ConjGrad: p");
  r = fmalloc(n, "ConjGrad: r");
  q = fmalloc(n, "ConjGrad: q");
  z = fmalloc(n, "ConjGrad: z");
  M = fmalloc(n, "ConjGrad: M");

  for (i=0; i<n; i++) {
    x[i] = 0.0;
    if (value[rowptr[i]] == 0.0)
      printf("%f\n", value[rowptr[i]]);
    M[i] = 1.0/value[rowptr[i]];  /* This is by construction */
  }

  /* r = b - Ax */
  mvMult(n, rowptr, colind, value, x, r);
  for (i=0; i<n; i++)
    r[i] = b[i]-r[i];

  bnrm2 = snorm2(n, b);
  if (bnrm2 > 0.0) {
    error = snorm2(n, r) / bnrm2;

    if (error > tol) { 
      /* Begin Iterations */
      for (k=0; k<n; k++) {
        for (i=0; i<n; i++)
          z[i] = r[i]*M[i];

        rho = sdot(n, r, z);

        if (k == 0)
          scopy(n, z, p);
        else {
          beta = rho/rho_1;
          for (i=0; i<n; i++)
            p[i] = z[i] + beta*p[i]; 
        }

        mvMult(n, rowptr, colind, value, p, q); /* q = A*p */

        alpha = rho/sdot(n, p, q);
        saxpy(n, alpha, p, x);    /* x = x + alpha*p */
        saxpy(n, -alpha, q, r);   /* r = r - alpha*q */
        error = snorm2(n, r) / bnrm2;
        if (error < tol)
          break;

        rho_1 = rho;
      }
    }
  }

  /* Free memory */
  GKfree(&p, &r, &q, &z, &M, -1);
}


/*************************************************************************
* This function implements a sparse matrix-vector multiplication
**************************************************************************/
void mvMult(int nvtxs, idxtype *rowptr, idxtype *colind, float *value, float *x, float *y)
{
  int i, j;

  for (i=0; i<nvtxs; i++) {
    y[i] = 0.0;
    for (j=rowptr[i]; j<rowptr[i+1]; j++) 
      y[i] += value[j]*x[colind[j]];
  }
}




/*************************************************************************
* This function updates the degrees of a k-way partition
**************************************************************************/
void Ser_KWayUpdateDegrees(GraphType *graph, int v, int k)
{
  int ii, i, u, tmp, from, to, ewgt, nedges;
  RInfoType *myrinfo;
  EdgeType *degrees;

  myrinfo = graph->rinfo + v;
  degrees = myrinfo->degrees;

  myrinfo->ed += (myrinfo->id - degrees[k].ewgt);

  SWAP(myrinfo->id, degrees[k].ewgt, tmp);

  from = graph->where[v];
  graph->where[v] = to = degrees[k].edge;

  degrees[k].edge = from;

  if (degrees[k].ewgt == 0) 
    degrees[k] = degrees[--myrinfo->ndegrees];

  INC_DEC(graph->lpwgts[to], graph->lpwgts[from], graph->vwgt[v]);

  /* Update the degrees of the adjacent vertices */
  for (ii=graph->xadj[v]; ii<graph->xadj[v+1]; ii++) {
    u = graph->adjncy[ii];
    ewgt = graph->adjwgt[ii];
    nedges = graph->xadj[u+1] - graph->xadj[u];
    myrinfo = graph->rinfo + u;
    degrees = myrinfo->degrees;

    /* If adjacent vertex is from original partition... */
    if (graph->where[u] == from) {
      INC_DEC(myrinfo->ed, myrinfo->id, ewgt);

      for (i=0; i<myrinfo->ndegrees; i++) {
        if (degrees[i].edge == to) {
          degrees[i].ewgt += ewgt;
          break;
        }
      }
      if (i == myrinfo->ndegrees) {
        degrees[i].ewgt = ewgt;
        degrees[i].edge = to;
        myrinfo->ndegrees++;
      }
    }
    /* If adjacent vertex is from new partition... */
    else if (graph->where[u] == to) {
      INC_DEC(myrinfo->id, myrinfo->ed, ewgt);

      for (i=0; i<myrinfo->ndegrees; i++) {
        if (degrees[i].edge == from) {
          degrees[i].ewgt -= ewgt;
          break;
        }
      }

      if (degrees[i].ewgt == 0) 
        degrees[i] = degrees[--myrinfo->ndegrees];
    }
    /* Adjacent vertex is from a different partition... */
    else {
      for (i=0; i<myrinfo->ndegrees; i++) {
        if (degrees[i].edge == from) {
          degrees[i].ewgt -= ewgt;
          break;
        }
      }

      if (degrees[i].ewgt == 0) 
        degrees[i] = degrees[--myrinfo->ndegrees];

      for (i=0; i<myrinfo->ndegrees; i++) {
        if (degrees[i].edge == to) {
          degrees[i].ewgt += ewgt;
          break;
        }
      }
      if (i == myrinfo->ndegrees) {
        degrees[i].ewgt = ewgt;
        degrees[i].edge = to;
        myrinfo->ndegrees++;
      }
    }
  }

}


#ifdef XYZ
/*************************************************************************
* This function performs a k-way directed diffusion
**************************************************************************/
void KWay_InitialDiffuser2(GraphType *graph, int nparts, float SuppressFactor)
{
  int ii, i, j, k, l, nvtxs, from, other, me, minpwgt, maxpwgt, oldcut;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *vsize;
  idxtype *where, *pwgts, *perm;
  int **transfer;
  int supped, elgble;
  float suppression;
  RInfoType *rinfo, *myrinfo;
  EdgeType *alldegrees, *degrees;
  int *rowptr, *colind;
  double sum, mean, *solution, *load, *value;
  int npasses = 40;
  int tvsize, tvwgt;
  float balance, oldbalance;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  vsize = vwgt;  /* This is temporary. GEORGE */

  /* Compute the partition parameters */
  rinfo = graph->rinfo = (RInfoType *)GKmalloc(sizeof(RInfoType)*nvtxs, "KWay_InitialDiffuser: rinfo");
  pwgts = graph->lpwgts = idxsmalloc(nparts, 0, "KWay_InitialDiffuser: pwgts");
  alldegrees = (EdgeType *)GKmalloc(sizeof(EdgeType)*graph->nedges, "KWay_InitialDiffuser: alldegrees");

  graph->mincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    myrinfo = rinfo+i;

    pwgts[me] += vwgt[i];

    myrinfo->degrees = alldegrees + xadj[i];
    myrinfo->ndegrees = myrinfo->id = myrinfo->ed = 0;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]])
        myrinfo->id += adjwgt[j];
      else
        myrinfo->ed += adjwgt[j];
    }

    if (myrinfo->ed > 0) {  /* Time to do some serious work */
      graph->mincut += myrinfo->ed;
      degrees = myrinfo->degrees;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (degrees[k].edge == other) {
              degrees[k].ewgt += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            degrees[k].edge = other;
            degrees[k].ewgt = adjwgt[j];
            myrinfo->ndegrees++;
          }
        }
      }
    }
  }
  graph->mincut = graph->mincut/2;

  /* printf("Cut: %d, %d %d\n", graph->mincut, pwgts[0], pwgts[1]); */

  tvwgt = idxsum(nvtxs, vwgt);
  tvsize = idxsum(nvtxs, vsize);

  suppression = 0.5;

  transfer = (int **) (malloc(sizeof(int *)*nparts));
  for(j=0; j<nparts; j++)
    transfer[j] = (int *) (malloc(sizeof(int)*nparts));

  perm = idxmalloc(nvtxs, "KWay_InitialDiffuser: perm");
  load = fmalloc(nparts, "KWay_InitialDiffuser: load");
  solution = fmalloc(nparts, "KWay_InitialDiffuser: solution");

  minpwgt = pwgts[idxamin(nparts, pwgts)];
  maxpwgt = pwgts[idxamax(nparts, pwgts)];
  balance = 1.0*nparts*maxpwgt/(1.0*tvwgt);

  for (l=0; l<npasses; l++) {
    oldcut = graph->mincut;
    oldbalance = balance;

    /* Set-up Diffuse */
    setupLaplace(graph, nparts, &rowptr, &colind, &value);

    sum = 0.0;
    for (j=0; j<nparts; j++) {
      load[j] = (double) (pwgts[j]);
      sum += load[j];
    }
    mean = sum/nparts;
    for (j = 0; j < nparts; j++)
      load[j] -= mean;

    for(j=0; j<nparts; j++) {
      solution[j] = 0.0;
      for (k = 0; k < nparts; k++)
        (transfer[j])[k] = 0;
    }

    /* Diffusion Solution */
    ConjGrad(nparts, rowptr, colind, value, load, solution, 0.001);

    /* printf("Difussion solution: %lf %lf %lf %lf\n", solution[0], solution[1], solution[2], solution[3]); */

    for (j = 0; j < nparts; j++)
      for (k = rowptr[j]; k < rowptr[j + 1]; k++)
        if (colind[k] != j)
          (transfer[j])[colind[k]] = (int) (solution[j] - solution[colind[k]]);

    FastRandomPermute(nvtxs, perm, 1);

    for (ii=0; ii<nvtxs; ii++) {
      i = perm[ii];
      if (rinfo[i].ndegrees > 0) {
        if (rinfo[i].ed-rinfo[i].id > -suppression*vwgt[i]) {
          degrees = rinfo[i].degrees;
          from = where[i];

          for (k=0; k<rinfo[i].ndegrees; k++)
            if ((transfer[from])[degrees[k].edge] > (.9 * vwgt[i]))
              break;

          if (k < rinfo[i].ndegrees) {
            for (j=k+1; j<rinfo[i].ndegrees; j++) {
              if ((transfer[from])[degrees[j].edge] > (.9 * vwgt[i]))
                if (degrees[j].ewgt > degrees[k].ewgt ||
                   (degrees[j].ewgt == degrees[k].ewgt &&
                      pwgts[degrees[j].edge] < pwgts[degrees[k].edge]))
                  k = j;
            }
            if (degrees[k].ewgt-rinfo[i].id > -suppression*vwgt[i]) {
              graph->mincut -= (degrees[k].ewgt-rinfo[i].id);
              (transfer[from])[degrees[k].edge] -= vwgt[i];
              Ser_KWayUpdateDegrees(graph, i, k);
            }
          }
        }
      }
    }

    GKfree(&value, &colind, &rowptr, -1);

    minpwgt = pwgts[idxamin(nparts, pwgts)];
    maxpwgt = pwgts[idxamax(nparts, pwgts)];
    balance = 1.0*nparts*maxpwgt/(1.0*tvwgt);

    if (balance > 0.95*oldbalance) 
      suppression *= 2;

    if (balance < UNBALANCE_FRACTION)
      break;
    if (graph->mincut == oldcut && minpwgt == pwgts[idxamin(nparts, pwgts)] && maxpwgt == pwgts[idxamax(nparts, pwgts)])
      break; 
  }

  GKfree(&load, &solution, &perm, &alldegrees, -1);
  for (i = 0; i < nparts; i++)
    free (transfer[i]);
  free (transfer);

  /* 
  printf("Directed Forced edge-cut is %d   Max pwgt is %d   Min pwgt is %d\n", graph->mincut, pwgts[idxamax(nparts, pwgts)], pwgts[idxamin(nparts, pwgts)]);
  {
    int i, j, cut = 0;
    for (i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        if (where[i] != where[adjncy[j]])
          cut += adjwgt[j];
    }
    printf("Actual cut: %d\n", cut/2);
  }
  */

  return;
}
#endif



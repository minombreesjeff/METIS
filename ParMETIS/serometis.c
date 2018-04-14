/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * serometis.c
 *
 * This file contains the code that finds the vertex separator of the coarsest 
 * graph using serial multilevel node-base bisection.
 *
 * Started 6/5/97
 * George
 *
 * $Id: serometis.c,v 1.3 1997/07/18 00:32:13 karypis Exp $
 *
 */

#include <par_kmetis.h>


#define BetterBalance(ubfactor, pwgts, gpwgts, from, vwgt) \
        (ubfactor*abs(pwgts[0]-gpwgts[0])+abs(pwgts[1]-gpwgts[1]) >= abs(pwgts[from]-vwgt-gpwgts[from]) + abs(pwgts[(from+1)%2]+vwgt-gpwgts[(from+1)%2]))


/*************************************************************************
* Let the game begin
**************************************************************************/
void Ser_NodeOMetis(GraphType *graph, float ubfactor)
{
  int i, j, k, maxvwgt;
  GraphType *cgraph;
  EdgeType *degrees;
  float balance;

  maxvwgt = idxsum(graph->nvtxs, graph->vwgt) / 20;

  cgraph = Ser_Coarsen(graph, 50, maxvwgt);
  
  Ser_NodeBisection(cgraph, ubfactor);

  Ser_NodeRefine(graph, cgraph, ubfactor);

}



/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void Ser_NodeBisection(GraphType *graph, float ubfactor)
{
  int i, j, k, nvtxs, nleft, first, last, tvwgt, pwgts[2], minpwgt[2], maxpwgt[2], from, bestcut, icut, mincut, me, pass, nbfs;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *where, *queue, *touched, *gain, *bestwhere;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where = idxmalloc(nvtxs, "BisectGraph: where");

  bestwhere = idxmalloc(nvtxs, "BisectGraph: bestwhere");
  queue = idxmalloc(nvtxs, "BisectGraph: queue");
  touched = idxmalloc(nvtxs, "BisectGraph: touched");

  tvwgt = idxsum(nvtxs, vwgt);
  minpwgt[0] = minpwgt[1] = (1.0/ubfactor)*tvwgt/2;
  maxpwgt[0] = maxpwgt[1] = ubfactor*tvwgt/2;

  for (nbfs=0; nbfs<NIPARTS; nbfs++) {
    idxset(nvtxs, 0, touched);

    pwgts[1] = tvwgt;
    pwgts[0] = 0;

    idxset(nvtxs, 1, where);

    queue[0] = RandomInRange(nvtxs);
    touched[queue[0]] = 1;
    first = 0; last = 1;
    nleft = nvtxs-1;

    /* Start the BFS from queue to get a partition */
    for (;;) {
      if (first == last) { /* Empty. Disconnected graph! */
        if (nleft == 0)
          break;

        k = RandomInRange(nleft);
        for (i=0; i<nvtxs; i++) {
          if (touched[i] == 0) {
            if (k == 0)
              break;
            else
              k--;
          }
        }

        queue[0] = i;
        touched[i] = 1;
        first = 0; last = 1;;
        nleft--;
      }

      i = queue[first++];
      if (pwgts[1]-vwgt[i] < minpwgt[1])
        continue;

      where[i] = 0;
      INC_DEC(pwgts[0], pwgts[1], vwgt[i]);
      if (pwgts[1] <= maxpwgt[1])
        break;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (touched[k] == 0) {
          queue[last++] = k;
          touched[k] = 1;
          nleft--;
        }
      }
    }

    /*************************************************************
    * Refine the partition using Edge-Based FM
    **************************************************************/
    Ser_EdgeFM(graph, minpwgt, maxpwgt, 2);
    Ser_ConstructSeparator(graph);
    Ser_NodeComputePartitionParams(graph);
    Ser_NodeFM(graph, ubfactor, 2);
    GKfree(&graph->nrinfo, &graph->lpwgts, -1);

    if (nbfs == 0 || bestcut > graph->mincut) {
      bestcut = graph->mincut;
      idxcopy(nvtxs, where, bestwhere);
    }
  }

  graph->mincut = bestcut;
  idxcopy(nvtxs, bestwhere, where);

  /* printf("InitBestCut: %d\n", bestcut); */

  GKfree(&bestwhere, &queue, &touched, -1);
}



/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void Ser_NodeRefine(GraphType *orggraph, GraphType *graph, float ubfactor)
{
  int i;

  for (i=0; ; i++) {
    Ser_NodeComputePartitionParams(graph);
    Ser_NodeFM(graph, ubfactor, 6); 

    if (graph == orggraph)
      break;

    graph = graph->finer;
    Ser_ProjectPartition(graph);
  }
}



/*************************************************************************
* This function computes the initial id/ed 
**************************************************************************/
void Ser_NodeComputePartitionParams(GraphType *graph)
{
  int i, j, k, l, nvtxs;
  idxtype *xadj, *adjncy, *adjwgt, *vwgt, *lpwgts;
  idxtype *where;
  NRInfoType *rinfo, *myrinfo;
  int me, other;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  rinfo = graph->nrinfo = (NRInfoType *)GKmalloc(sizeof(NRInfoType)*nvtxs, "NodeComputePartitionParams: rinfo");
  lpwgts = graph->lpwgts = idxsmalloc(3, 0, "ComputePartitionParams: lpwgts");


  /*------------------------------------------------------------
  / Compute now the separator external degrees
  /------------------------------------------------------------*/
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    lpwgts[me] += vwgt[i];

    if (me == 2) { /* If it is on the separator do some computations */
      myrinfo = rinfo+i;
      myrinfo->edegrees[0] = myrinfo->edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (other != 2)
          myrinfo->edegrees[other] += vwgt[adjncy[j]];
      }
    }
  }

  graph->mincut = lpwgts[2];
}



/*************************************************************************
* This function performs a node-based Kernighan-Lin refinement 
**************************************************************************/
void Ser_NodeFM(GraphType *graph, float ubfactor, int npasses)
{
  int i, j, k, jj, kk, nvtxs, nswaps, nmind;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *mptr, *mind, *moved, *swaps;
  PQueueType parts[2]; 
  NRInfoType *rinfo, *myrinfo;
  int higain, oldgain, mincut, initcut, mincutorder;	
  int pass=0;
  int to, other, limit;
  int badminpwgt, badmaxpwgt;
  int u[2], g[2];

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;
  where = graph->where;
  pwgts = graph->lpwgts;
  rinfo = graph->nrinfo;

  limit = amax(0.001*nvtxs, 15);
  limit = amin(limit, 50);

  PQueueInit(&parts[0], graph->nvtxs);
  PQueueInit(&parts[1], graph->nvtxs);

  moved = idxmalloc(nvtxs, "Ser_NodeFM: moved");
  swaps = idxmalloc(nvtxs, "Ser_NodeFM: swaps");
  mptr = idxmalloc(nvtxs+1, "Ser_NodeFM: mptr");
  mind = idxmalloc(2*nvtxs, "Ser_NodeFM: mind");

  for (pass=0; pass<npasses; pass++) {
    PQueueReset(&parts[0]);
    PQueueReset(&parts[1]);

    mincutorder = -1;
    initcut = mincut = graph->mincut;
    idxset(nvtxs, -1, moved);

    for (i=0; i<nvtxs; i++) {
      if (where[i] == 2) {
        PQueueInsert(&parts[0], i, vwgt[i]-rinfo[i].edegrees[1]);
        PQueueInsert(&parts[1], i, vwgt[i]-rinfo[i].edegrees[0]);
      }
    }

/*
    if (pass == 0)
      printf("Partitions: [%6d, %6d, %6d, %6d] Initial Cut: %8d\n",
              pwgts[0], pwgts[1], pwgts[2], parts[0].nnodes, initcut);
*/

    /******************************************************
    * Get into the FM loop
    *******************************************************/
    mptr[0] = nmind = 0;
    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      badminpwgt = (1.0/ubfactor)*(pwgts[0]+pwgts[1])/2;
      badmaxpwgt = ubfactor*(pwgts[0]+pwgts[1])/2;

      ASSERTS(CheckPartitionParams(graph));

      /* to = (pwgts[0] < pwgts[1] ? 0 : 1); */

      u[0] = PQueueSeeMax(&parts[0]);  
      u[1] = PQueueSeeMax(&parts[1]);
      if (u[0] == -1 || u[1] == -1)
        break;

      g[0] = vwgt[u[0]]-rinfo[u[0]].edegrees[1];
      g[1] = vwgt[u[1]]-rinfo[u[1]].edegrees[0];

      to = (g[0] > g[1] ? 0 : 1);
      other = (to+1)%2;

      if (pwgts[to]+vwgt[u[to]] > badmaxpwgt) {
        to = (to+1)%2;
        other = (to+1)%2;
      }
      if (pwgts[other]-rinfo[u[to]].edegrees[other] < badminpwgt) {
        to = (to+1)%2;
        other = (to+1)%2;
      }
      

      higain = PQueueGetMax(&parts[to]);
      PQueueDelete(&parts[other], higain);

      pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);

      if (pwgts[2] <= mincut) {
        mincut = pwgts[2];
        mincutorder = nswaps;
      }
      else {
        if (nswaps - mincutorder > limit) {
          pwgts[2] += (vwgt[higain]-rinfo[higain].edegrees[other]);
          break; /* No further improvement, break out */
        }
      }

      pwgts[to] += vwgt[higain];
      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;  

/*
printf("Selected %3d to move to %3d, Gain: %3d, \t[%5d %5d %5d]\n", higain, to, vwgt[higain]-rinfo[higain].edegrees[other], pwgts[0], pwgts[1], pwgts[2]);
*/

      /**********************************************************
      * Update the degrees of the affected nodes
      ***********************************************************/
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        /* printf("Working with %d %d\n", k, where[k]); */
        if (where[k] == 2) { /* For the in-separator vertices modify their edegree[to] */
          oldgain = vwgt[k]-rinfo[k].edegrees[to];
          rinfo[k].edegrees[to] += vwgt[higain];
          if (moved[k] == -1)
            PQueueUpdate(&parts[other], k, oldgain, oldgain-vwgt[higain]);
        }
        else if (where[k] == other) { /* This vertex is pulled into the separator */
          mind[nmind++] = k;  /* Keep track for rollback */
          where[k] = 2;
          pwgts[other] -= vwgt[k];

          rinfo[k].edegrees[0] = rinfo[k].edegrees[1] = 0;
          for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
            kk = adjncy[jj];
            /* printf("\tWorking with %d %d\n", kk, where[kk]); */
            if (where[kk] != 2) {
              rinfo[k].edegrees[where[kk]] += vwgt[kk];
            }
            else {
              oldgain = vwgt[kk]-rinfo[kk].edegrees[other];
              rinfo[kk].edegrees[other] -= vwgt[k];
              if (moved[kk] == -1)
                PQueueUpdate(&parts[to], kk, oldgain, oldgain+vwgt[k]);
            }
          }

          /* Insert the new vertex into the priority queue */
          if (moved[k] == -1) {
            PQueueInsert(&parts[0], k, vwgt[k]-rinfo[k].edegrees[1]);
            PQueueInsert(&parts[1], k, vwgt[k]-rinfo[k].edegrees[0]);
          }
        }
      }
      mptr[nswaps+1] = nmind;
    }

/*
    printf("\tMinimum Cut: %8d at %5d [%6d %6d]\n",mincut, mincutorder, pwgts[0], pwgts[1]);
*/

    /****************************************************************
    * Roll back computation 
    *****************************************************************/
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];
      ASSERTS(moved[higain] > mincutorder);

      ASSERTS(CheckPartitionParams(graph));
      /*
      printf("Rolling back vertex %3d, \t[%5d %5d %5d]\n", higain, pwgts[0], pwgts[1], pwgts[2]);
      */

      to = where[higain];
      other = (to+1)%2;
      INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
      where[higain] = 2;

      rinfo[higain].edegrees[0] = rinfo[higain].edegrees[1] = 0;
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) 
          rinfo[k].edegrees[to] -= vwgt[higain];
        else
          rinfo[higain].edegrees[where[k]] += vwgt[k];
      }

      /* Push nodes out of the separator */
      for (j=mptr[nswaps]; j<mptr[nswaps+1]; j++) {
        k = mind[j];
        ASSERTS(where[k] == 2);
        where[k] = other;
        INC_DEC(pwgts[other], pwgts[2], vwgt[k]);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          kk = adjncy[jj];
          if (where[kk] == 2) 
            rinfo[kk].edegrees[other] += vwgt[k];
        }
      }
    }

    ASSERTS(mincut == pwgts[2]);

    graph->mincut = mincut;

    if (mincutorder == -1 || mincut >= initcut)
      break;
  }

  PQueueFree(&parts[0]);
  PQueueFree(&parts[1]);

  GKfree(&moved, &swaps, &mptr, &mind, -1);

}



/*************************************************************************
* This function is used to construct a separator from a partitioning
**************************************************************************/
void Ser_ConstructSeparator(GraphType *graph)
{
  int i, ii, j, nvtxs, nparts, nsep, sepwgt;
  idxtype *xadj, *adjncy, *mark, *where;
  KeyValueType *vwpair;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  mark = idxsmalloc(nvtxs, 0, "ConstructSeparator: mark");
  vwpair = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*nvtxs, "ConstructSeparator: vwpair");

  for (nsep=i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) {
        vwpair[nsep].key = graph->vwgt[i];
        vwpair[nsep++].val = i;
        break;
      }
    }
  }

  ikeysort(nsep, vwpair);
  sepwgt = 0;
  for (ii=0; ii<nsep; ii++) {
    i = vwpair[ii].val;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]] && !mark[adjncy[j]]) {
        sepwgt += graph->vwgt[i];
        mark[i] = 1;
        break;
      }
    }
  }

  for (i=0; i<nvtxs; i++) {
    if (mark[i]) 
      where[i] = 2;
  }

  GKfree(&mark, &vwpair, -1);
}


/*************************************************************************
* This function checks the correctness of the NodeFM data structures
**************************************************************************/
int CheckPartitionParams(GraphType *graph)
{
  int i, j, k, l, nvtxs, me, other;
  idxtype *xadj, *adjncy, *adjwgt, *vwgt, *where;
  idxtype edegrees[2], lpwgts[3];

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;

  /*------------------------------------------------------------
  / Compute now the separator external degrees
  /------------------------------------------------------------*/
  lpwgts[0] = lpwgts[1] = lpwgts[2] = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    lpwgts[me] += vwgt[i];

    if (me == 2) { /* If it is on the separator do some computations */
      edegrees[0] = edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (other != 2)
          edegrees[other] += vwgt[adjncy[j]];
      }
      if (edegrees[0] != graph->nrinfo[i].edegrees[0] || edegrees[1] != graph->nrinfo[i].edegrees[1])
        printf("Something wrong with edegrees: %d %d %d %d %d\n", i, edegrees[0], edegrees[1], graph->nrinfo[i].edegrees[0], graph->nrinfo[i].edegrees[1]);
    }
  }

  if (lpwgts[0] != graph->lpwgts[0] || lpwgts[1] != graph->lpwgts[1] || lpwgts[2] != graph->lpwgts[2])
    printf("Something wrong with part-weights: %d %d %d %d %d %d\n", lpwgts[0], lpwgts[1], lpwgts[2], graph->lpwgts[0], graph->lpwgts[1], graph->lpwgts[2]);

  return 1;
}


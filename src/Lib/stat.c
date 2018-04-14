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
 * $Id: stat.c,v 1.1 1997/11/04 23:19:35 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function computes cuts and balance information
**************************************************************************/
void ComputePartitionInfo(GraphType *graph, int nparts, idxtype *where)
{
  int i, j, k, nvtxs;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt, *kpwgts, *tmpptr;
  idxtype *padjncy, *padjwgt, *padjcut;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  vwgt = graph->vwgt;
  adjwgt = graph->adjwgt;

  printf("%d-way Edge-cut: %5d, ", nparts, ComputeCut(graph, where));

  /* Compute balance information */
  kpwgts = idxsmalloc(nparts, 0, "ComputePartitionInfo: kpwgts");

  for (i=0; i<nvtxs; i++)
    kpwgts[where[i]] += vwgt[i];
  printf("\tBalance: %5.3f out of %5.3f\n", 
          1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)),
          1.0*nparts*vwgt[idxamax(nvtxs, vwgt)]/(1.0*idxsum(nparts, kpwgts)));


  /* Compute p-adjncy information */
  padjncy = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjncy");
  padjwgt = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjwgt");
  padjcut = idxsmalloc(nparts*nparts, 0, "ComputePartitionInfo: padjwgt");

  idxset(nparts, 0, kpwgts);
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) {
        padjncy[where[i]*nparts+where[adjncy[j]]] = 1;
        padjcut[where[i]*nparts+where[adjncy[j]]] += adjwgt[j];
        if (kpwgts[where[adjncy[j]]] == 0) {
          padjwgt[where[i]*nparts+where[adjncy[j]]]++;
          kpwgts[where[adjncy[j]]] = 1;
        }
      }
    }
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      kpwgts[where[adjncy[j]]] = 0;
  }

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjncy+i*nparts);
  printf("Min/Max/Avg/Bal # of adjacent     subdomains: %5d %5d %5d %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)], idxsum(nparts, kpwgts)/nparts, 
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)));

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjcut+i*nparts);
  printf("Min/Max/Avg/Bal # of adjacent subdomain cuts: %5d %5d %5d %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)], idxsum(nparts, kpwgts)/nparts, 
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)));

  for (i=0; i<nparts; i++)
    kpwgts[i] = idxsum(nparts, padjwgt+i*nparts);
  printf("Min/Max/Avg/Bal/Frac # of interface    nodes: %5d %5d %5d %7.3f %7.3f\n",
    kpwgts[idxamin(nparts, kpwgts)], kpwgts[idxamax(nparts, kpwgts)], idxsum(nparts, kpwgts)/nparts, 
    1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts)), 1.0*idxsum(nparts, kpwgts)/(1.0*nvtxs));

  tmpptr = graph->where;
  graph->where = where;
  for (i=0; i<nparts; i++)
    IsConnectedSubdomain(NULL, graph, i, 1);
  graph->where = tmpptr;

  GKfree(&kpwgts, &padjncy, &padjwgt, &padjcut, -1);
}



/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
float ComputePartitionBalance(GraphType *graph, int nparts, idxtype *where)
{
  int i;
  idxtype *kpwgts;
  float balance;


  kpwgts = idxsmalloc(nparts, 0, "ComputePartitionInfo: kpwgts");

  for (i=0; i<graph->nvtxs; i++)
    kpwgts[where[i]] += graph->vwgt[i];

  balance = 1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts));

  free(kpwgts);

  return balance;

}


/*************************************************************************
* This function computes the balance of the element partitioning
**************************************************************************/
float ComputeElementBalance(int ne, int nparts, idxtype *where)
{
  int i;
  idxtype *kpwgts;
  float balance;

  kpwgts = idxsmalloc(nparts, 0, "ComputeElementBalance: kpwgts");

  for (i=0; i<ne; i++)
    kpwgts[where[i]]++;

  balance = 1.0*nparts*kpwgts[idxamax(nparts, kpwgts)]/(1.0*idxsum(nparts, kpwgts));

  free(kpwgts);

  return balance;

}

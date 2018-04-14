/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * xyzpart.c
 *
 * This file contains code that implements a coordinate based partitioning
 *
 * Started 7/11/97
 * George
 *
 * $Id: xyzpart.c 10361 2011-06-21 19:16:22Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************
* This function implements a simple coordinate based partitioning
**************************************************************************/
void Coordinate_Partition(ctrl_t *ctrl, graph_t *graph, idx_t ndims, 
         real_t *xyz, idx_t setup)
{
  idx_t i, j, k, nvtxs, firstvtx, icoord, coords[3];
  idx_t *vtxdist;
  real_t max[3], min[3], gmin[3], gmax[3], shift[3], scale[3];
  ikv_t *cand;

  WCOREPUSH;

  if (setup)
    CommSetup(ctrl, graph);
  else
    graph->nrecv = 0;

  nvtxs    = graph->nvtxs;
  vtxdist  = graph->vtxdist;
  firstvtx = vtxdist[ctrl->mype];

  cand = ikvwspacemalloc(ctrl, nvtxs);

  /* Compute parameters for coordinate transformation */
  for (k=0; k<ndims; k++) {
    min[k] = +10000000;
    max[k] = -10000000;
  }
  for (i=0; i<nvtxs; i++) {
    for (k=0; k<ndims; k++) {
      if (xyz[i*ndims+k] < min[k])
        min[k] = xyz[i*ndims+k];
      if (xyz[i*ndims+k] > max[k])
        max[k] = xyz[i*ndims+k];
    }
  }

  /* Compute global min and max */
  gkMPI_Allreduce((void *)min, (void *)gmin, ndims, REAL_T, MPI_MIN, ctrl->comm);
  gkMPI_Allreduce((void *)max, (void *)gmax, ndims, REAL_T, MPI_MAX, ctrl->comm);

  /* myprintf(ctrl, "Coordinate Range: %e %e, Global %e %e\n", min[0], max[0], gmin[0], gmax[0]); */

  for (k=0; k<ndims; k++) {
    /* rprintf(ctrl, "Dim#%"PRIDX": %e %e, span: %e\n", k, gmin[k], gmax[k], gmax[k]-gmin[k]); */
    shift[k] = -gmin[k];
    if (gmax[k] != gmin[k])
      scale[k] = 1.0/(gmax[k]-gmin[k]);
    else
      scale[k] = 1.0;
  }

  switch (ctrl->xyztype) {
    case XYZ_XCOORD:
      for (i=0; i<nvtxs; i++) {
        cand[i].key = 1000000*((xyz[i*ndims]+shift[0])*scale[0]);
        PASSERT(ctrl, cand[i].key>=0 && cand[i].key<=1000000);
        cand[i].val = firstvtx+i;
      }
      break;
    case XYZ_SPFILL:
      for (i=0; i<nvtxs; i++) {
        /* make the coordinates to be ints in the 0..1023 range (i.e., 10 bits) */
        for (k=0; k<ndims; k++)
          coords[k] = 1023*((xyz[i*ndims+k]+shift[k])*scale[k]); 
        for (icoord=0, j=9; j>=0; j--) {
          for (k=0; k<ndims; k++)
            icoord = (icoord<<1) + (coords[k]&(1<<j) ? 1 : 0);
        }
        cand[i].key = icoord;
        cand[i].val = firstvtx+i;
      }
      break;
    default:
      errexit("Unknown XYZ_Type type!\n");
  }


  /* Partition using sorting */
  PartSort(ctrl, graph, cand);

  WCOREPOP;
}



/**************************************************************************/
/*! This function sorts a distributed list of ikv_t in increasing 
    order, and uses it to compute a partition. It uses samplesort. 

    This function is poorly implemented and makes the assumption that the
    number of vertices in each processor is greater than npes. 
    This constraint is currently enforced by the calling functions. 
    \todo fix it in 4.0.
*/
/**************************************************************************/
void PartSort(ctrl_t *ctrl, graph_t *graph, ikv_t *elmnts)
{
  idx_t i, j, k, nvtxs, nrecv, npes=ctrl->npes, mype=ctrl->mype, 
        firstvtx, lastvtx;
  idx_t *scounts, *rcounts, *vtxdist, *perm;
  ikv_t *relmnts, *mypicks, *allpicks;

  WCOREPUSH;

  nvtxs   = graph->nvtxs;
  vtxdist = graph->vtxdist;

  /* get memory for the counts */
  scounts = iwspacemalloc(ctrl, npes+1);
  rcounts = iwspacemalloc(ctrl, npes+1);

  /* get memory for the splitters */
  mypicks  = ikvwspacemalloc(ctrl, npes+1);
  WCOREPUSH; /* for freeing allpicks */
  allpicks = ikvwspacemalloc(ctrl, npes*npes);

  /* Sort the local elements */
  ikvsorti(nvtxs, elmnts);

  /* Select the local npes-1 equally spaced elements */
  for (i=1; i<npes; i++) { 
    mypicks[i-1].key = elmnts[i*(nvtxs/npes)].key;
    mypicks[i-1].val = elmnts[i*(nvtxs/npes)].val;
  }

  /* PrintPairs(ctrl, npes-1, mypicks, "Mypicks"); */

  /* Gather the picks to all the processors */
  gkMPI_Allgather((void *)mypicks, 2*(npes-1), IDX_T, (void *)allpicks, 
      2*(npes-1), IDX_T, ctrl->comm);

  /* PrintPairs(ctrl, npes*(npes-1), allpicks, "Allpicks"); */

  /* Sort all the picks */
  ikvsortii(npes*(npes-1), allpicks);

  /* PrintPairs(ctrl, npes*(npes-1), allpicks, "Allpicks"); */

  /* Select the final splitters. Set the boundaries to simplify coding */
  for (i=1; i<npes; i++)
    mypicks[i] = allpicks[i*(npes-1)];
  mypicks[0].key    = IDX_MIN;
  mypicks[npes].key = IDX_MAX;

  /* PrintPairs(ctrl, npes+1, mypicks, "Mypicks"); */

  WCOREPOP;  /* free allpicks */

  /* Compute the number of elements that belong to each bucket */
  iset(npes, 0, scounts);
  for (j=i=0; i<nvtxs; i++) {
    if (elmnts[i].key < mypicks[j+1].key || 
        (elmnts[i].key == mypicks[j+1].key && elmnts[i].val < mypicks[j+1].val))
      scounts[j]++;
    else
      scounts[++j]++;
  }
  gkMPI_Alltoall(scounts, 1, IDX_T, rcounts, 1, IDX_T, ctrl->comm);

  MAKECSR(i, npes, scounts);
  MAKECSR(i, npes, rcounts);

/*
  PrintVector(ctrl, npes+1, 0, scounts, "Scounts");
  PrintVector(ctrl, npes+1, 0, rcounts, "Rcounts");
*/

  /* Allocate memory for sorted elements and receive them */
  nrecv   = rcounts[npes];
  relmnts = ikvwspacemalloc(ctrl, nrecv);

  /* Issue the receives first */
  for (i=0; i<npes; i++) 
    gkMPI_Irecv((void *)(relmnts+rcounts[i]), 2*(rcounts[i+1]-rcounts[i]), 
        IDX_T, i, 1, ctrl->comm, ctrl->rreq+i);

  /* Issue the sends next */
  for (i=0; i<npes; i++) 
    gkMPI_Isend((void *)(elmnts+scounts[i]), 2*(scounts[i+1]-scounts[i]), 
        IDX_T, i, 1, ctrl->comm, ctrl->sreq+i);

  gkMPI_Waitall(npes, ctrl->rreq, ctrl->statuses);
  gkMPI_Waitall(npes, ctrl->sreq, ctrl->statuses);


  /* OK, now do the local sort of the relmnts. Use perm to keep track original order */
  perm = iwspacemalloc(ctrl, nrecv);
  for (i=0; i<nrecv; i++) {
    perm[i]        = relmnts[i].val;
    relmnts[i].val = i;
  }
  ikvsorti(nrecv, relmnts);


  /* Compute what needs to be shifted */
  gkMPI_Scan((void *)(&nrecv), (void *)(&lastvtx), 1, IDX_T, MPI_SUM, ctrl->comm);
  firstvtx = lastvtx-nrecv;  

  /*myprintf(ctrl, "first, last: %"PRIDX" %"PRIDX"\n", firstvtx, lastvtx); */

  for (j=0, i=0; i<npes; i++) {
    if (vtxdist[i+1] > firstvtx) {  /* Found the first PE that is passed me */
      if (vtxdist[i+1] >= lastvtx) {
        /* myprintf(ctrl, "Shifting %"PRIDX" elements to processor %"PRIDX"\n", lastvtx-firstvtx, i); */
        for (k=0; k<lastvtx-firstvtx; k++, j++) 
          relmnts[relmnts[j].val].key = i;
      }
      else {
        /* myprintf(ctrl, "Shifting %"PRIDX" elements to processor %"PRIDX"\n", vtxdist[i+1]-firstvtx, i); */
        for (k=0; k<vtxdist[i+1]-firstvtx; k++, j++) 
          relmnts[relmnts[j].val].key = i;

        firstvtx = vtxdist[i+1];
      }
    }
    if (vtxdist[i+1] >= lastvtx)
      break;
  }

  /* Reverse the ordering on the relmnts[].val */
  for (i=0; i<nrecv; i++) {
    PASSERTP(ctrl, relmnts[i].key>=0 && relmnts[i].key<npes, 
            (ctrl, "%"PRIDX" %"PRIDX"\n", i, relmnts[i].key));
    relmnts[i].val = perm[i];
  }

  /* OK, now sent it back */
  /* Issue the receives first */
  for (i=0; i<npes; i++) 
    gkMPI_Irecv((void *)(elmnts+scounts[i]), 2*(scounts[i+1]-scounts[i]), IDX_T, 
        i, 1, ctrl->comm, ctrl->rreq+i);

  /* Issue the sends next */
  for (i=0; i<npes; i++) 
    gkMPI_Isend((void *)(relmnts+rcounts[i]), 2*(rcounts[i+1]-rcounts[i]), IDX_T, 
        i, 1, ctrl->comm, ctrl->sreq+i);

  gkMPI_Waitall(npes, ctrl->rreq, ctrl->statuses);
  gkMPI_Waitall(npes, ctrl->sreq, ctrl->statuses);


  /* Construct a partition for the graph */
  graph->where = imalloc(graph->nvtxs+graph->nrecv, "PartSort: graph->where");
  firstvtx = vtxdist[mype];
  for (i=0; i<nvtxs; i++) {
    PASSERTP(ctrl, elmnts[i].key>=0 && elmnts[i].key<npes, 
        (ctrl, "%"PRIDX" %"PRIDX"\n", i, elmnts[i].key));
    PASSERTP(ctrl, elmnts[i].val>=vtxdist[mype] && elmnts[i].val<vtxdist[mype+1], 
        (ctrl, "%"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", i, vtxdist[mype], vtxdist[mype+1], elmnts[i].val));
    graph->where[elmnts[i].val-firstvtx] = elmnts[i].key;
  }

  WCOREPOP;
}


/**
 * @file pushrelabel.c
 * @brief MPI implementation of the push-relabel algorithm for max-flow
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-07-01
 */





#ifndef BOWSTRING_MPI_PUSHRELABEL_C
#define BOWSTRING_MPI_PUSHRELABEL_C




#include <mpi.h>
#include <stdlib.h>
#include <domlib.h>
#include "base.h"
#include "io/io.h"
#include "graph.h"




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define MPI_VTX_T MPI_UNSIGNED
#define MPI_ADJ_T MPI_UNSIGNED
#define MPI_WGT_T MPI_FLOAT




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef enum __tags {
  TAG_XADJ,
  TAG_ADJNCY,
  TAG_ADJWGT,
  TAG_VREDGE,
  TAG_UREDGE,
  TAG_FLOW,
  TAG_HEIGHT,
  TAG_MAXFLOW
} tags;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLHT_PREFIX adj
#define DLHT_KEY_T uint64_t
#define DLHT_VAL_T adj_t
#define DLHT_STATIC
#include "dlht_headers.h"
#undef DLHT_STATIC
#undef DLHT_VAL_T
#undef DLHT_KEY_T
#undef DLHT_PREFIX




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static int __check_graph(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    adj_t const * const radj,
    adj_t const nredges)
{
  vtx_t i, k;
  adj_t j, nedges;

  nedges = xadj[nvtxs];

  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k >= nvtxs + nredges) {
        printf("Invalid edge on vertex "PF_VTX_T"/"PF_VTX_T" to vertex " \
            PF_VTX_T"/"PF_VTX_T"+"PF_ADJ_T"\n",i,nvtxs,k,nvtxs,nredges);
        return 0;
      }
      if (radj[j] >= nedges + nredges) {
        printf("Invalid reverse adjancy index on edge "PF_ADJ_T" to " \
            PF_ADJ_T"/"PF_ADJ_T"+"PF_ADJ_T"\n",j,radj[j],nedges,nredges);
        return 0;
      }
    }
  }

  for (j=nedges;j<nedges+nredges;++j) {
    k = adjncy[j];
    if (k >= nvtxs) {
      printf("Invalid edge on remote vertex "PF_ADJ_T"/"PF_VTX_T"+"PF_ADJ_T \
          " to vertex "PF_VTX_T"/"PF_VTX_T"\n",j,nvtxs,nredges,k, \
          nvtxs);
      return 0;
    }
    if (radj[j] >= nedges) {
      printf("Invalid reverse adjancy index on edge "PF_ADJ_T" to " \
          PF_ADJ_T"/"PF_ADJ_T"+"PF_ADJ_T"\n",j,radj[j],nedges,nredges);
      return 0;
    }
  }

  return 1;
}


static void __build_chunk(
    vlbl_t const p,
    adj_t const * const gxadj,
    vtx_t const * const gadjncy,
    wgt_t const * const gadjwgt,
    vtx_t const * const dist,
    vtx_t const * const rename,
    vtx_t const * const grename,
    vlbl_t const * const map,
    adj_t ** const r_xadj,
    vtx_t ** const r_adjncy,
    wgt_t ** const r_adjwgt)
{
  vtx_t nvtxs, i, k, l;
  adj_t j, nedges;
  adj_t * xadj;
  vtx_t * adjncy;
  wgt_t * adjwgt;
  vlbl_t o;

  nvtxs = dist[p+1]-dist[p];
  xadj = adj_alloc(nvtxs+1);
  nedges = 0;
  for (i=0;i<nvtxs;++i) {
    xadj[i] = nedges;
    k = grename[i+dist[p]];
    nedges += gxadj[k+1] - gxadj[k];
  }
  xadj[nvtxs] = nedges;

  adjncy = vtx_alloc(nedges);
  adjwgt = wgt_alloc(nedges);

  nedges = 0;
  for (i=0;i<nvtxs;++i) {
    k = grename[i+dist[p]];
    for (j=gxadj[k];j<gxadj[k+1];++j) {
      l = gadjncy[j];
      o = map[l];
      if (o == p) {
        adjncy[nedges] = rename[l];
      } else {
        adjncy[nedges] = rename[l] + dist[o] + nvtxs;
      }
      adjwgt[nedges] = gadjwgt[j];
      ++nedges;
    }
  }

  *r_xadj = xadj;
  *r_adjncy = adjncy;
  *r_adjwgt = adjwgt;
}


static int __distribute_graph(
    char const * const graphfile,
    char const * const partfile,
    adj_t ** const r_xadj,
    vtx_t ** const r_adjncy,
    wgt_t ** const r_adjwgt,
    vtx_t ** const r_dist,
    vtx_t ** const r_rename,
    vlbl_t ** const r_map)
{
  int rank, size, err;
  vlbl_t nparts, myid, p;
  vtx_t i, k, nvtxs, gnvtxs;
  adj_t nedges;
  vlbl_t * map = NULL;
  vtx_t * gadjncy = NULL, * adjncy = NULL, * dist = NULL, * rename = NULL, \
        * grename = NULL;
  adj_t * gxadj = NULL, * xadj = NULL;
  wgt_t * gadjwgt = NULL, * adjwgt = NULL;
  MPI_Status status;

  err = BOWSTRING_SUCCESS;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  nparts = (vlbl_t)size;
  myid = (vlbl_t)rank;
  dist = vtx_alloc(nparts+1);

  if (myid == 0) {
    /* first process reads in the graph and partition and sends it */
    if ((err = bowstring_read_graph(graphfile,BOWSTRING_FORMAT_AUTO,&gnvtxs, \
        &gxadj,&gadjncy,NULL,&gadjwgt)) != BOWSTRING_SUCCESS) {
      goto END;
    }

    if (gadjwgt == NULL) {
      gadjwgt = wgt_init_alloc(1.0,gxadj[gnvtxs]);
    }

    /* read in partition */
    if ((err = read_vertex_labels(partfile,&gnvtxs,&map)) != \
        BOWSTRING_SUCCESS) {
      goto END;
    }

    /* renumber the vertices locally */
    rename = vtx_alloc(gnvtxs);
    vlbl_set_min(map,0,gnvtxs);
    vtx_set(dist,0,nparts+1);
    for (i=0;i<gnvtxs;++i) {
      p = map[i];
      if (p >= nparts) {
        eprintf("Partition has too many parts: "PF_VLBL_T" parts and " \
            PF_VLBL_T" processors.\n",p,nparts);
        err = BOWSTRING_ERROR_INVALIDINPUT;
        goto END;
      }
      rename[i] = dist[p]++;
    }
    vtx_prefixsum_exc(dist,nparts+1);

    /* renumber the vertices globally */
    grename = vtx_alloc(gnvtxs);
    for (i=0;i<gnvtxs;++i) {
      p = map[i];
      k = rename[i] + dist[p];
      grename[k] = i;
    }

    /* send vertex distribution */
    MPI_Bcast(dist,nparts+1,MPI_VTX_T,0,MPI_COMM_WORLD);

    /* build each part of the graph and send it */
    for (p=1;p<nparts;++p) {
      __build_chunk(p,gxadj,gadjncy,gadjwgt,dist,rename,grename,map,&xadj, \
          &adjncy,&adjwgt);

      nvtxs = dist[p+1] - dist[p];
      nedges = xadj[nvtxs];

      MPI_Send(xadj,nvtxs+1,MPI_ADJ_T,p,TAG_XADJ,MPI_COMM_WORLD);
      MPI_Send(adjncy,nedges,MPI_VTX_T,p,TAG_ADJNCY,MPI_COMM_WORLD);
      MPI_Send(adjwgt,nedges,MPI_WGT_T,p,TAG_ADJWGT,MPI_COMM_WORLD);

      dl_free(xadj);
      xadj = NULL;
      dl_free(adjncy);
      adjncy = NULL;
      dl_free(adjwgt);
      adjwgt = NULL;
    }
    /* build my own part of the graph and free the rest */
    __build_chunk(0,gxadj,gadjncy,gadjwgt,dist,rename,grename,map,&xadj, \
        &adjncy,&adjwgt);
  } else {
    /* other processes recieve the graph */
    MPI_Bcast(dist,nparts+1,MPI_VTX_T,0,MPI_COMM_WORLD);

    nvtxs = dist[myid+1] - dist[myid];

    xadj = adj_alloc(nvtxs+1);
    MPI_Recv(xadj,nvtxs+1,MPI_ADJ_T,0,TAG_XADJ,MPI_COMM_WORLD,&status);

    adjncy = vtx_alloc(xadj[nvtxs]);
    adjwgt = wgt_alloc(xadj[nvtxs]);
    MPI_Recv(adjncy,xadj[nvtxs],MPI_VTX_T,0,TAG_ADJNCY,MPI_COMM_WORLD,&status);
    MPI_Recv(adjwgt,xadj[nvtxs],MPI_WGT_T,0,TAG_ADJWGT,MPI_COMM_WORLD,&status);
  }

  /* assign output variables and keep them from being cleared during cleanup */
  *r_xadj = xadj;
  xadj = NULL;
  *r_adjncy = adjncy;
  adjncy = NULL;
  *r_adjwgt = adjwgt;
  adjwgt = NULL;
  *r_dist = dist;
  dist = NULL;
  *r_rename = rename;
  rename = NULL;
  *r_map = map;
  map = NULL;

  END:

  if (gxadj) {
    dl_free(gxadj);
  }
  if (gadjncy) {
    dl_free(gadjncy);
  }
  if (gadjwgt) {
    dl_free(gadjwgt);
  }
  if (rename) {
    dl_free(rename);
  }
  if (grename) {
    dl_free(grename);
  }
  if (xadj) {
    dl_free(xadj);
  }
  if (adjncy) {
    dl_free(adjncy);
  }
  if (adjwgt) {
    dl_free(adjwgt);
  }
  if (dist) {
    dl_free(dist);
  }
  if (map) {
    dl_free(map);
  }

  return err;
}


static int __choose_source_sink(
    vtx_t const gsource,
    vtx_t const gsink,
    vtx_t const * const dist,
    vtx_t const * const rename,
    vlbl_t const * const map,
    vtx_t * const r_source,
    vtx_t * const r_sink)
{
  int rank;
  vtx_t source_sink[2];

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank == 0) {
    /* adjust for partition */
    source_sink[0] = dist[map[gsource]] + rename[gsource];
    source_sink[1] = dist[map[gsink]] + rename[gsink];
  }

  if (MPI_Bcast(source_sink,2,MPI_VTX_T,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
    return BOWSTRING_ERROR_MPICALL;
  }

  *r_source = source_sink[0];
  *r_sink = source_sink[1];

  return BOWSTRING_SUCCESS;
}


static int __build_comm_structures(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t ** const r_adjncy,
    vtx_t const * const dist,
    adj_t ** const r_radj,
    adj_t ** const r_ncedges,
    wgt_t ** const r_flow)

{
  int rank, size, err;
  uint64_t eid;
  vtx_t i, k;
  adj_t j, l, nredges, nedges;
  vlbl_t nparts, p, myid;
  adj_t * radj = NULL, * ncedges = NULL, * trans = NULL, * adjncy;
  vtx_t * vredges = NULL, * uredges = NULL, * vsedges = NULL, * usedges = NULL; 
  adj_ht_t * ht = NULL;
  wgt_t * flow;
  MPI_Status status;
  MPI_Request * reqv = NULL, * requ = NULL;

  err = BOWSTRING_SUCCESS;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  myid = (vlbl_t)rank;
  nparts = (vlbl_t)size;

  reqv = malloc(sizeof(MPI_Request)*nparts);
  requ = malloc(sizeof(MPI_Request)*nparts);

  nedges = xadj[nvtxs];

  adjncy = *r_adjncy;

  /* count remote vertices and edges */
  ncedges = adj_calloc(nparts+1);
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k >= nvtxs) {
        /* a remote edge */
        p = 0;  
        while (dist[p+1] <= k-nvtxs) {
          ++p;
        }
        DL_ASSERT(p != myid,"Detected remote local edge");
        ++ncedges[p];
      }
    }
  }

  vtx_prefixsum_exc(ncedges,nparts+1);

  radj = adj_alloc(nedges+ncedges[nparts]);
  build_adjncy_index_rem(nvtxs,xadj,adjncy,radj);

  vredges = vtx_alloc(ncedges[nparts]);
  uredges = vtx_alloc(ncedges[nparts]);

  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k >= nvtxs) {
        /* remote edge */
        p = 0;  
        while (dist[p+1] <= k-nvtxs) {
          ++p;
        }
        /* will point to extra area in the flows array */
        radj[j] = ncedges[p] + nedges;
        radj[radj[j]] = j;

        /* send gloabal vertex pairs */
        vredges[ncedges[p]] = k - nvtxs - dist[p];
        uredges[ncedges[p]] = i;
        ++ncedges[p];
      }
      DL_ASSERT_EQUALS(radj[radj[j]],j,PF_ADJ_T);
    }
  }
  /* restore the ncedges array */
  vtx_prefixshift(ncedges,nparts+1);

  DL_ASSERT_EQUALS(ncedges[nparts],ncedges[nparts],PF_ADJ_T);

  /* keep extra room for the reverse side of remote edges */
  flow = wgt_calloc(nedges + ncedges[nparts]);

  vsedges = vtx_alloc(ncedges[nparts]);
  usedges = vtx_alloc(ncedges[nparts]);

  /* asynchrounously recieve edge orderings */
  for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
    nredges = ncedges[p+1] - ncedges[p];
    if (nredges > 0) {
      if (MPI_Irecv(vsedges+ncedges[p],nredges,MPI_VTX_T,p, \
            TAG_VREDGE,MPI_COMM_WORLD,reqv+p) != MPI_SUCCESS || \
          MPI_Irecv(usedges+ncedges[p],nredges,MPI_VTX_T,p, \
            TAG_UREDGE,MPI_COMM_WORLD,requ+p) != MPI_SUCCESS) {
        eprintf("%d: Async-Recieving "PF_ADJ_T" remote edge indexes from " \
            PF_VLBL_T" failed\n",rank,nredges,p);
        err = BOWSTRING_ERROR_MPICALL;
        goto END;
      }
    }
  }

  /* synchronously send edge orderings */
  for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
    nredges = ncedges[p+1] - ncedges[p];
    if (nredges > 0) {
      if (MPI_Send(vredges+ncedges[p],nredges,MPI_VTX_T,p,TAG_VREDGE, \
            MPI_COMM_WORLD) != MPI_SUCCESS || \
          MPI_Send(uredges+ncedges[p],nredges,MPI_VTX_T,p,TAG_UREDGE, \
            MPI_COMM_WORLD) != MPI_SUCCESS) {
        eprintf("%d: Sync-Sending "PF_ADJ_T" local edge indexes to " \
            PF_VLBL_T" failed\n",rank,nredges,p);
        err = BOWSTRING_ERROR_MPICALL;
        goto END;
      }
    }
  }

  /* wait for communication to finish */
  for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
    nredges = ncedges[p+1] - ncedges[p];
    if (nredges > 0) {
      if (MPI_Wait(reqv+p,&status) != MPI_SUCCESS || \
          MPI_Wait(requ+p,&status) != MPI_SUCCESS) {
        eprintf("%d: Failed waiting to recieve "PF_ADJ_T" remote edges from " \
            PF_VLBL_T"\n",rank,nredges,p);
        err = BOWSTRING_ERROR_MPICALL;
        goto END;
      }
    }
  }

  /* build organize my remote edgse */
  trans = adj_alloc(ncedges[nparts]);
  for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
    nredges = ncedges[p+1] - ncedges[p];
    if (nredges > 0) {
      /* lower number procs use the ordering of the high number procs */
      if (rank < p) {
        /* use a hash map for creating the translation table */
        ht = adj_ht_create(nredges*4,nredges);

        /* insert my end of the edges into the hash map */
        for (j=ncedges[p];j<ncedges[p+1];++j) {
          i = uredges[j]; /* local vertex */
          k = vredges[j]; /* remote vertex */

          /* the key needs to be unique amongst this group of edges */
          eid = (((uint64_t)i) << dl_bitsize(uint32_t)) | ((uint64_t)k);

          if ((l = adj_ht_put(eid,j,ht)) != adj_ht_null_value) {
            eprintf("Attempting to overwrite "PF_ADJ_T" with "PF_ADJ_T" at " \
                "index %"PRIu64", "PF_VTX_T":"PF_VTX_T"\n",l,j,eid,i,k);
            abort();
          }
        }

        /* save my current radj */
        adj_copy(trans+ncedges[p],radj+ncedges[p]+nedges,nredges);

        /* use p's end of the edges to order my local edges */
        for (j=ncedges[p];j<ncedges[p+1];++j) {
          i = vsedges[j]; /* local vertex */
          k = usedges[j]; /* remote vertex */

          /* lookup the old ordering */
          eid = (((uint64_t)i) << dl_bitsize(uint32_t)) | ((uint64_t)k);
          l = adj_ht_get(eid,ht);

          /* update the radj ordering */
          radj[j+nedges] = trans[l];
          radj[radj[j+nedges]] = j+nedges; 

          /* save the new ordering */
          uredges[j] = i;
        }

        adj_ht_free(ht);
      }

      /* modify the adjncy array to point to the extra spaces */
      for (j=ncedges[p];j<ncedges[p+1];++j) {
        l = radj[j+nedges];
        DL_ASSERT_EQUALS(radj[l],j+nedges,PF_ADJ_T);
        adjncy[l] = nvtxs + j;
      }
    }
  }

  /* extend adjncy to look back at local vertices */
  *r_adjncy = adjncy = adj_realloc(adjncy,nedges+ncedges[nparts]);
  for (j=0;j<ncedges[nparts];++j) {
    i = uredges[j];
    adjncy[j+nedges] = i;
    DL_ASSERT_EQUALS(adjncy[radj[j+nedges]],nvtxs+j,PF_VTX_T);
  }

  *r_radj = radj;
  radj = NULL;
  *r_ncedges = ncedges;
  ncedges = NULL;
  *r_flow = flow;
  flow = NULL;

  END:

  if (reqv) {
    dl_free(reqv);
  }
  if (requ) {
    dl_free(requ);
  }
  if (radj) {
    dl_free(radj);
  }
  if (ncedges) {
    dl_free(ncedges);
  }
  if (flow) {
    dl_free(flow);
  }
  if (trans) {
    dl_free(trans);
  }
  if (vsedges) {
    dl_free(vsedges);
  }
  if (usedges) {
    dl_free(usedges);
  }
  if (vredges) {
    dl_free(vredges);
  }
  if (uredges) {
    dl_free(uredges);
  }

  return err;
}


static int __usage(
    const char * const name, 
    FILE * fout)
{
  fprintf(fout,"USAGE:\n");
  fprintf(fout,"%s <graphfile> <partfile> <source> <sink>\n",name);
  fprintf(fout,"\n");

  return 1;
}



/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


/**
 * @brief The main function. Should be called:
 * pushrelabel graphfile seed partfile
 *
 * @param argc
 * @param argv
 *
 * @return 
 */
int main(
    int argc,
    char ** argv)
{
  int rank, size, err;
  size_t iter;
  adj_t j, l, nedges, rnedges;
  vlbl_t nparts, p, myid;
  vtx_t gnvtxs, nvtxs, i, k, nops, gnops, sink, source, minh; 
  dl_timer_t pre_tmr, com_tmr, cpu_tmr;
  wgt_t f, maxflow;
  vlbl_t * map = NULL;
  adj_t * xadj = NULL, * radj = NULL, * ncedges = NULL;
  vtx_t * adjncy = NULL, * height = NULL, * sheight = NULL, * dist = NULL, \
        * rename = NULL;
  wgt_t * adjwgt = NULL, * flow = NULL, * rflow = NULL, * excess = NULL;
  MPI_Status status;
  MPI_Request * req = NULL;

  /* setup mpi */
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  myid = (vlbl_t)rank;
  nparts = (vlbl_t)size;

  if (argc != 5) {
    if (rank == 0) {
      eprintf("Invalid number of arguments.\n");
      __usage(argv[0],stderr);
    }
    err = BOWSTRING_ERROR_INVALIDINPUT;
    goto END;
  }

  /* init timers */
  dl_init_timer(&pre_tmr);
  dl_init_timer(&cpu_tmr);
  dl_init_timer(&com_tmr);

  dl_start_timer(&pre_tmr);

  req = malloc(sizeof(MPI_Request)* size);

  if ((err = __distribute_graph(argv[1],argv[2],&xadj,&adjncy,&adjwgt,&dist, \
      &rename,&map)) != BOWSTRING_SUCCESS) {
    goto END;
  }

  /* calc number of vertices */
  nvtxs = dist[myid+1] - dist[myid];
  nedges = xadj[nvtxs];
  gnvtxs = dist[nparts];

  /* communicate source and sink */
  if ((err = __choose_source_sink(strtoull(argv[3],NULL,10)-1,
      strtoull(argv[4],NULL,10)-1,dist,rename,map,&source,&sink)) \
      != BOWSTRING_SUCCESS) {
    goto END;
  }

  /* cleanup */
  if (rank == 0) {
    dl_free(map);
    map = NULL;
    dl_free(rename);
    rename = NULL;
  }

  /* setup communication structures */
  if ((err = __build_comm_structures(nvtxs,xadj,&adjncy,dist,&radj,&ncedges, \
      &flow)) != BOWSTRING_SUCCESS) {
    goto END;
  }

  DL_ASSERT(__check_graph(nvtxs,xadj,adjncy,radj,ncedges[nparts]),"Bad " \
      "graph after preprocessing\n");

  /* allocate excess and height structures */
  excess = wgt_calloc(nvtxs+ncedges[nparts]);
  height = vtx_calloc(nvtxs+ncedges[nparts]);

  /* recieving buffers */
  sheight = vtx_alloc(ncedges[nparts]);
  rflow = wgt_alloc(ncedges[nparts]);

  /* set initial height of source */
  if (source >= dist[myid] && source < dist[myid+1]) {
    source -= dist[myid];
    height[source] = gnvtxs;
    excess[source] = INFINITY;
  } else {
    source = NULL_VTX;
  }

  /* set initial excess of sink */
  if (sink >= dist[myid] && sink < dist[myid+1]) {
    sink -= dist[myid];
    excess[sink] = -INFINITY;
  } else {
    sink = NULL_VTX;
  }

  dl_stop_timer(&pre_tmr);

  iter = 0;

  do {
    nops = 0;
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
      err = BOWSTRING_ERROR_MPICALL;
      eprintf("Barrier failed\n");
      goto END;
    }

    /* communicate height updates */
    dl_start_timer(&com_tmr);

    /* copy my heights to the buffer */
    for (j=0;j<ncedges[nparts];++j) {
      k = adjncy[nedges+j];
      sheight[j] = height[k];
    }

    /* recieve heights asynchronously */
    for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
      rnedges = ncedges[p+1] - ncedges[p];
      if (rnedges > 0) {
        if (MPI_Irecv(height+ncedges[p]+nvtxs,rnedges,MPI_VTX_T,p,TAG_HEIGHT, \
            MPI_COMM_WORLD,req+p) != MPI_SUCCESS) {
          err = BOWSTRING_ERROR_MPICALL;
          goto END;
        }
      }
    }
    /* send heights */
    for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
      rnedges = ncedges[p+1] - ncedges[p];
      if (rnedges > 0) {
        if (MPI_Send(sheight+ncedges[p],rnedges,MPI_VTX_T,p,TAG_HEIGHT, \
            MPI_COMM_WORLD) != MPI_SUCCESS) {
          err = BOWSTRING_ERROR_MPICALL;
          goto END;
        }
      }
    }
    /* wait for asynchronous messages to finish */
    for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
      rnedges = ncedges[p+1] - ncedges[p];
      if (rnedges > 0) {
        if (MPI_Wait(req+p,&status) != MPI_SUCCESS) {
          err = BOWSTRING_ERROR_MPICALL;
          goto END;
        }
      }
    }

    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
      err = BOWSTRING_ERROR_MPICALL;
      eprintf("Barrier failed\n");
      goto END;
    }
    dl_stop_timer(&com_tmr);

    /* compute flows */
    dl_start_timer(&cpu_tmr);
    wgt_set(flow+nedges,0,ncedges[nparts]);
    for (i=0;i<nvtxs;++i) {
      if (excess[i] > 0) {
        /* push my excess */
        for (j=xadj[i];j<xadj[i+1];++j) {
          k = adjncy[j];
          if (height[i] > height[k]) {
            f = dl_min(excess[i],adjwgt[j] - flow[j]);
            if (f > 0) {
              excess[k] += f;
              excess[i] -= f;
              flow[j] += f;
              if (k < nvtxs) {
                flow[radj[j]] -= f;
              } else {
                flow[radj[j]] = f;
              }
              ++nops;
            }
          }
        }
      }
    }

    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
      err = BOWSTRING_ERROR_MPICALL;
      eprintf("Barrier failed\n");
      goto END;
    }
    dl_stop_timer(&cpu_tmr);

    /* communicate flow updates */
    dl_start_timer(&com_tmr);

    /* recieve flows asynchronously */
    for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
      rnedges = ncedges[p+1] - ncedges[p];
      if (rnedges > 0) {
        if (MPI_Irecv(rflow+ncedges[p],rnedges,MPI_WGT_T,p,TAG_FLOW, \
            MPI_COMM_WORLD,req+p) != MPI_SUCCESS) {
          err = BOWSTRING_ERROR_MPICALL;
          goto END;
        }
      }
    }
    /* send flows */
    for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
      rnedges = ncedges[p+1] - ncedges[p];
      if (rnedges > 0) {
        if (MPI_Send(flow+nedges+ncedges[p],rnedges,MPI_WGT_T,p,TAG_FLOW, \
            MPI_COMM_WORLD) != MPI_SUCCESS) {
          err = BOWSTRING_ERROR_MPICALL;
          goto END;
        }
      }
    }
    /* wait for asynchronous messages to finish */
    for (p=(myid+1)%nparts;p!=myid;p=(p+1)%nparts) {
      rnedges = ncedges[p+1] - ncedges[p];
      if (rnedges > 0) {
        if (MPI_Wait(req+p,&status) != MPI_SUCCESS) {
          err = BOWSTRING_ERROR_MPICALL;
          goto END;
        }
      }
    }

    /* update flows */
    for (j=0;j<ncedges[nparts];++j) {
      if (rflow[j] > 0) {
        f = rflow[j];
        l = radj[nedges+j];
        i = adjncy[nedges+j]; /* local vertex */
        k = adjncy[l]; /* remote vertex */

        DL_ASSERT(height[k] > height[i],"Invalid flow from "PF_VTX_T"/" \
            PF_VTX_T" :" PF_VTX_T" to "PF_VTX_T"/"PF_VTX_T" :"PF_VTX_T" " \
            "of "PF_WGT_T" ("PF_WGT_T":"PF_WGT_T")\n",k,nvtxs,height[k],i, \
            nvtxs,height[i],f,flow[nedges+j],rflow[j]);

        flow[nedges+j] += f;
        flow[l] -= f;
        excess[i] += f;
      }
    }

    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
      err = BOWSTRING_ERROR_MPICALL;
      eprintf("Barrier failed\n");
      goto END;
    }
    dl_stop_timer(&com_tmr);

    dl_start_timer(&cpu_tmr);
    for (i=0;i<nvtxs;++i) {
      if (excess[i] > 0) {
        minh = gnvtxs;
        for (j=xadj[i];j<xadj[i+1];++j) {
          k = adjncy[j];
          DL_ASSERT(k<nvtxs+ncedges[nparts],"The value of k is outside of "
              "range: "PF_VTX_T"/"PF_VTX_T"\n",k,nvtxs+ncedges[nparts]);
          f = dl_min(excess[i],adjwgt[j] - flow[j]);
          if (f > 0) {
            if (height[i] > height[k]) {
              if (f > 0) {
                goto FLOW_AVAILABLE;
              }
            }
            minh = dl_min(minh,height[k]);
          }
        }
        /* relabel if necessary */
        if (i != source && excess[i] > 0 && height[i] < minh+1) {
          height[i] = minh + 1;
          ++nops;
        }
        FLOW_AVAILABLE:;
      }
    }

    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
      err = BOWSTRING_ERROR_MPICALL;
      eprintf("Barrier failed\n");
      goto END;
    }

    dl_stop_timer(&cpu_tmr);

    MPI_Allreduce(&nops,&gnops,1,MPI_VTX_T,MPI_SUM,MPI_COMM_WORLD);

    ++iter;
  } while (gnops > 0);

  MPI_Barrier(MPI_COMM_WORLD);

  /* calcuate max flow, and send it to the root process if needed */
  if (sink != NULL_VTX) {
    maxflow = 0;
    for (j=xadj[sink];j<xadj[sink+1];++j) {
      maxflow += -flow[j];
    }
    if (rank != 0) {
      MPI_Send(&maxflow,1,MPI_WGT_T,0,TAG_MAXFLOW,MPI_COMM_WORLD);
    }
  } else if (rank == 0) {
    MPI_Recv(&maxflow,1,MPI_WGT_T,MPI_ANY_SOURCE,TAG_MAXFLOW,MPI_COMM_WORLD, \
        &status);
  }

  if (rank == 0) {
    printf("Maximum flow of "PF_WGT_T"\n",maxflow); 
    printf("Completed %zu iterations.\n",iter);
    printf("Preprocessing: %0.05lf\n",dl_poll_timer(&pre_tmr));
    printf("Computation: %0.05lf\n",dl_poll_timer(&cpu_tmr));
    printf("Communication: %0.05lf\n",dl_poll_timer(&com_tmr));
  }

  END:

  if (req) {
    dl_free(req);
  }
  if (map) {
    dl_free(map);
  } 
  if (rename) {
    dl_free(rename);
  }
  if (sheight) {
    dl_free(sheight);
  }
  if (rflow) {
    dl_free(rflow);
  }
  if (height) {
    dl_free(height);
  }
  if (dist) {
    dl_free(dist);
  }
  if (ncedges) {
    dl_free(ncedges);
  }
  if (flow) {
    dl_free(flow);
  }
  if (excess) {
    dl_free(excess);
  }
  if (radj) {
    dl_free(radj);
  }
  if (xadj) {
    dl_free(xadj);
  }
  if (adjncy) {
    dl_free(adjncy);
  }
  if (adjwgt) {
    dl_free(adjwgt);
  }

  MPI_Finalize();

  if (err != BOWSTRING_SUCCESS) {
    return 1;
  } else { 
    return 0;
  }
}




#endif

/**
 * @file kwinfo.c
 * @brief Functions for allocating an manipulating kwinfos. 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-19
 */





#ifndef MTMETIS_KWINFO_C
#define MTMETIS_KWINFO_C




#include "kwinfo.h"




/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


#define DLMEM_PREFIX kwnbrinfo
#define DLMEM_TYPE_T kwnbrinfo_t
#include <dlmem_funcs.h>
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX adjinfo
#define DLMEM_TYPE_T adjinfo_t
#include <dlmem_funcs.h>
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const MIN_NNBRPOOL = 1024;
static size_t const NPOOLS = 64;




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void par_kwinfo_create(
    ctrl_t * const ctrl,
    graph_t * const graph)
{
  size_t pool;
  kwinfo_t * kwinfo;

  tid_t const myid = dlthread_get_id(ctrl->comm);
  tid_t const nthreads = dlthread_get_nthreads(ctrl->comm);

  kwinfo = dlthread_get_shmem(sizeof(kwinfo_t)*nthreads,ctrl->comm);

  kwinfo += myid;

  kwinfo->nparts = ctrl->nparts;
  kwinfo->bnd = vtx_iset_create(0,graph->mynvtxs[myid]);
  kwinfo->nbrinfo = kwnbrinfo_alloc(graph->mynvtxs[myid]);

  /* setup initial pool */
  kwinfo->nnbrpool = 0;
  kwinfo->basebits = size_uplog2(dl_max(graph->mynedges[myid]*4,MIN_NNBRPOOL));
  kwinfo->basennbrs = 1 << kwinfo->basebits;
  kwinfo->npools = NPOOLS;
  kwinfo->nbrpools = malloc(sizeof(adjinfo_t*)*kwinfo->npools);

  kwinfo->nbrpools[0] = malloc(sizeof(adjinfo_t)*kwinfo->basennbrs);
  for (pool=1;pool<kwinfo->npools;++pool) {
    kwinfo->nbrpools[pool] = NULL;
  }

  dlthread_init_lock(&kwinfo->lock);

  if (myid == 0) {
    graph->kwinfo = kwinfo;
  }

  dlthread_barrier(ctrl->comm);
}


void par_kwinfo_free(
    graph_t * const graph)
{
  size_t n;
  tid_t pid;
  kwinfo_t * kwinfo;

  tid_t const myid = dlthread_get_id(graph->comm);
  tid_t const nthreads = dlthread_get_nthreads(graph->comm);

  for (pid=myid;pid<graph->dist.nthreads;pid+=nthreads) {
    kwinfo = graph->kwinfo+pid;

    if (kwinfo->bnd) {
      vtx_iset_free(kwinfo->bnd);
    }
    if (kwinfo->nbrinfo) {
      dl_free(kwinfo->nbrinfo);
    }
    if (kwinfo->nbrpools) {
      for (n=0;n<kwinfo->npools;++n) {
        if (kwinfo->nbrpools[n] != NULL) {
          dl_free(kwinfo->nbrpools[n]);
        } else {
          /* nothing else is allocated, stop looking */
          break;
        }
      }
      dl_free(kwinfo->nbrpools);
    }
    dlthread_free_lock(&kwinfo->lock);

    dlthread_free_shmem(graph->kwinfo,graph->comm);
    if (pid == 0) {
      graph->kwinfo = NULL;
    }
  }
}




#endif

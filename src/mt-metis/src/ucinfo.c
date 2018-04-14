/**
 * @file ucinfo.c
 * @brief Functions for allocating an manipulating ucinfos. 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-19
 */





#ifndef MTMETIS_UCINFO_C
#define MTMETIS_UCINFO_C




#include "ucinfo.h"




/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


#define DLMEM_PREFIX nbrinfo
#define DLMEM_TYPE_T nbrinfo_t
#include <dlmem_funcs.h>
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX adjinfo
#define DLMEM_TYPE_T adjinfo_t
#include <dlmem_funcs.h>
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


ucinfo_t * ucinfo_create(
    ctrl_t const * const ctrl,
    graph_t const * const graph)
{
  ucinfo_t * ucinfo;

  tid_t const myid = omp_get_thread_num();

  ucinfo = (ucinfo_t*)malloc(sizeof(ucinfo_t));

  ucinfo->nparts = ctrl->nparts;
  ucinfo->bnd = vtx_iset_create(0,graph->mynvtxs[myid]);
  ucinfo->nbrinfo = nbrinfo_alloc(graph->mynvtxs[myid]);
  ucinfo->nnbrpool = 0;
  ucinfo->maxnnbrpool = graph->mynedges[myid]*4;
  ucinfo->nbrpool = adjinfo_alloc(ucinfo->maxnnbrpool);

  return ucinfo;
}


void ucinfo_free(
    ucinfo_t * ucinfo)
{
  if (ucinfo->bnd) {
    vtx_iset_free(ucinfo->bnd);
  }
  if (ucinfo->nbrinfo) {
    dl_free(ucinfo->nbrinfo);
  }
  if (ucinfo->nbrpool) {
    dl_free(ucinfo->nbrpool);
  }
  dl_free(ucinfo);
}

#endif

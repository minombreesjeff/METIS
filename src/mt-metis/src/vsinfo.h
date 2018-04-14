/**
 * @file vsinfo.h
 * @brief Tupes and function prototypes for vertex separator info. 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-11-11
 */




#ifndef MTMETIS_VSINFO_H
#define MTMETIS_VSINFO_H




#include "base.h"
#include "graph.h"
#include "ctrl.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct vsnbrinfo_t {
  /* no need to store information about connectivity to separator */
  wgt_t con[2];
} vsnbrinfo_t;


typedef struct vsinfo_t {
  vtx_iset_t * bnd;
  vsnbrinfo_t * nbrinfo;
} vsinfo_t;




/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


#define DLMEM_PREFIX vsnbrinfo
#define DLMEM_TYPE_T vsnbrinfo_t
#include <dlmem_headers.h>
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define vsinfo_free __mtmetis_vsinfo_free
/**
 * @brief Free a vsinfo and its associate memory.
 *
 * @param graph The graph to free the vsinfo of.
 */
void vsinfo_free(
    graph_t * graph);


#define par_vsinfo_create __mtmetis_par_vsinfo_create
/**
 * @brief Allocate the memory arrays for refinement of a vertex separator. 
 *
 * @param ctrl The control structure.
 * @param graph The graph to create vsinfo for.
 */
void par_vsinfo_create(
    ctrl_t * ctrl,
    graph_t * graph);


#define par_vsinfo_free __mtmetis_par_vsinfo_free
/**
 * @brief Free a vsinfo and its associate memory.
 *
 * @param graph The graph to free the vsinfo of.
 */
void par_vsinfo_free(
    graph_t * graph);




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


static inline void __calc_conn(
    vtx_t const v,
    tid_t const myid,
    vtx_t const mynvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const * const gvwgt,
    pid_t const * const * const gwhere,
    graphdist_t const dist,
    wgt_t * const con)
{
  vtx_t k, lvtx;
  adj_t j;
  tid_t nbrid;
  pid_t nbr;
  wgt_t a, b, w;

  a = 0;
  b = 0;

  for (j=xadj[v];j<xadj[v+1];++j) {
    k = adjncy[j];
    if (k < mynvtxs) {
      lvtx = k;
      nbrid = myid;
    } else {
      lvtx = gvtx_to_lvtx(k,dist);
      nbrid = gvtx_to_tid(k,dist);
    }
    nbr = gwhere[nbrid][lvtx];
    w = gvwgt[nbrid][lvtx];
    switch (nbr) {
      case 0:
        a += w; 
        break;
      case 1:
        b += w;
        break;
    }
  }

  con[0] = a;
  con[1] = b;
}




#endif

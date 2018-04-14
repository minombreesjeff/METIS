/**
 * @file ucinfo.h
 * @brief Types and function for uncoarsening information.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-19
 */





#ifndef MTMETIS_UCINFO_H
#define MTMETIS_UCINFO_H




#include "base.h"
#include "graph.h"
#include "ctrl.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct adjinfo_t {
  wgt_t ed;
  pid_t pid;
} adjinfo_t;


typedef struct nbrinfo_t {
  wgt_t id;
  wgt_t ed;
  pid_t nnbrs;
  adj_t nbrstart;
} nbrinfo_t;


typedef struct ucinfo_t {
  pid_t nparts;
  vtx_iset_t * bnd;
  nbrinfo_t * nbrinfo;
  vtx_t  nnbrpool;
  vtx_t  maxnnbrpool;
  adjinfo_t * nbrpool;
} ucinfo_t;



/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


#define DLMEM_PREFIX nbrinfo
#define DLMEM_TYPE_T nbrinfo_t
#include <dlmem_headers.h>
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX adjinfo
#define DLMEM_TYPE_T adjinfo_t
#include <dlmem_headers.h>
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define ucinfo_create __mtmetis_ucinfo_create
/**
 * @brief Allocate the memory arrays for refinement. 
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 *
 * @return The thread's newly allocated ucinfo
 */
ucinfo_t * ucinfo_create(
    ctrl_t const * ctrl,
    graph_t const * graph);


#define ucinfo_free __mtmetis_ucinfo_free
/**
 * @brief Free a ucinfo and its associate memory.
 *
 * @param ucinfo The ucinfo to free.
 */
void ucinfo_free(
    ucinfo_t * ucinfo);



#endif

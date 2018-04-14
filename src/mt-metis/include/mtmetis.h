/**
 * @file mtmetis.h
 * @brief Library entry points
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-01
 */





#ifndef MTMETIS_H
#define MTMETIS_H




/******************************************************************************
* INCLUDES ********************************************************************
******************************************************************************/


#include <unistd.h>
#include <stdint.h>
#include <float.h>




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define MTMETIS_VER_MAJOR 0
#define MTMETIS_VER_MINOR 4
#define MTMETIS_VER_SUBMINOR 3




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


#ifndef MTMETIS_GRAPH_TYPES_DEFINED
/* define these types for internal use */
#ifdef MTMETIS_64BIT_VERTICES
typedef uint64_t mtmetis_vtx_t;
#else
typedef uint32_t mtmetis_vtx_t;
#endif

#ifdef MTMETIS_64BIT_EDGES
typedef uint64_t mtmetis_adj_t;
#else
typedef uint32_t mtmetis_adj_t;
#endif

/* this must be signed for refinement to work */
#ifdef MTMETIS_64BIT_WEIGHTS
typedef int64_t mtmetis_wgt_t;
#else
typedef int32_t mtmetis_wgt_t;
#endif
#endif /* MTMETIS_GRAPH_TYPES_DEFINED */


#ifdef MTMETIS_64BIT_PARTITIONS
typedef uint64_t mtmetis_pid_t;
#else
typedef uint32_t mtmetis_pid_t;
#endif




/* enums *********************************************************************/

typedef enum mtmetis_error_t {
  MTMETIS_SUCCESS = 1,
  MTMETIS_ERROR_INVALIDINPUT,
  MTMETIS_ERROR_NOTENOUGHMEMORY,
  MTMETIS_ERROR_THREADING
} mtmetis_error_t;


typedef enum mtmetis_option_t {
  MTMETIS_OPTION_HELP,
  MTMETIS_OPTION_TIME,
  MTMETIS_OPTION_NPARTS,
  MTMETIS_OPTION_NTHREADS,
  MTMETIS_OPTION_SEED,
  MTMETIS_OPTION_NCUTS,
  MTMETIS_OPTION_NRUNS,
  MTMETIS_OPTION_NINITSOLUTIONS,
  MTMETIS_OPTION_NITER,
  MTMETIS_OPTION_UBFACTOR,
  MTMETIS_OPTION_CTYPE,
  MTMETIS_OPTION_CONTYPE,
  MTMETIS_OPTION_LEAFMATCH,
  MTMETIS_OPTION_RTYPE,
  MTMETIS_OPTION_PTYPE,
  MTMETIS_OPTION_VERBOSITY,
  MTMETIS_OPTION_DISTRIBUTION,
  MTMETIS_OPTION_RUNSTATS,
  MTMETIS_OPTION_METIS,
  MTMETIS_OPTION_REMOVEISLANDS,
  MTMETIS_OPTION_VWGTDEGREE,
  MTMETIS_OPTION_IGNORE,
  __MTMETIS_OPTION_TERM
} mtmetis_option_t;


typedef enum mtmetis_ctype_t {
  MTMETIS_CTYPE_RM,
  MTMETIS_CTYPE_SHEM,
  MTMETIS_CTYPE_SLEM,
  MTMETIS_CTYPE_MDM,
  MTMETIS_CTYPE_FC
} mtmetis_ctype_t;


typedef enum mtmetis_contype_t {
  MTMETIS_CONTYPE_CLS,
  MTMETIS_CONTYPE_DENSE,
  MTMETIS_CONTYPE_SORT
} mtmetis_contype_t;


typedef enum mtmetis_rtype_t {
  MTMETIS_RTYPE_GREEDY,
  MTMETIS_RTYPE_FM,
  MTMETIS_RTYPE_SFM,
  MTMETIS_RTYPE_SFG,
} mtmetis_rtype_t;


typedef enum mtmetis_ptype_t {
  MTMETIS_PTYPE_KWAY,
  MTMETIS_PTYPE_ESEP,
  MTMETIS_PTYPE_RB,
  MTMETIS_PTYPE_VSEP,
  MTMETIS_PTYPE_ND
} mtmetis_ptype_t;


typedef enum mtmetis_verbosity_t {
  MTMETIS_VERBOSITY_NONE,
  MTMETIS_VERBOSITY_LOW,
  MTMETIS_VERBOSITY_MEDIUM,
  MTMETIS_VERBOSITY_HIGH,
  MTMETIS_VERBOSITY_MAXIMUM
} mtmetis_verbosity_t;


typedef enum mtmetis_dtype_t {
  MTMETIS_DISTRIBUTION_BLOCK,
  MTMETIS_DISTRIBUTION_CYCLIC,
  MTMETIS_DISTRIBUTION_BLOCKCYCLIC
} mtmetis_dtype_t;


typedef enum mtmetis_part_t {
  MTMETIS_VSEP_NULL = -1,
  MTMETIS_VSEP_PARTA = 0,
  MTMETIS_VSEP_PARTB = 1,
  MTMETIS_VSEP_SEP = 2,
  MTMETIS_VSEP_NPARTS = 3,
  MTMETIS_ESEP_PARTA = 0,
  MTMETIS_ESEP_PARTB = 1,
  MTMETIS_ESEP_NPARTS = 2
} mtmetis_part_t;


typedef enum mtmetis_ignore_t {
  MTMETIS_IGNORE_NONE = 0x00,
  MTMETIS_IGNORE_VERTEXWEIGHTS = 0x01,
  MTMETIS_IGNORE_EDGEWEIGHTS = 0x02
} mtmetis_ignore_t;




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const MTMETIS_NOPTIONS = __MTMETIS_OPTION_TERM;
static double const MTMETIS_VAL_OFF = -DBL_MAX;




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Allocate and initialize a set of options for use with the
 * mtmetis_partkway_explicit() function.
 *
 * @return The allocated and initialized options. 
 */
double * mtmetis_init_options(void);


/**
 * @brief Partition a graph using default options. 
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 * @param nparts The number of partitions.
 * @param where The partition ID of each vertex (can be NULL, or of length
 * nvtxs)
 * @param r_edgecut A reference to the weight of cut edges (can be NULL).
 *
 * @return MTMETIS_SUCCESS unless an error is encountered. 
 */
int mtmetis_partkway(
    mtmetis_vtx_t nvtxs, 
    mtmetis_adj_t const * xadj, 
    mtmetis_vtx_t const * adjncy, 
    mtmetis_wgt_t const * vwgt, 
    mtmetis_wgt_t const * adjwgt,
    mtmetis_pid_t nparts, 
    mtmetis_pid_t * where, 
    mtmetis_wgt_t * r_edgecut);


/**
 * @brief Partition a graph using an explicit set of options.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 * @param options The set of options.
 * @param where The partition ID of each vertex (can be NULL, or of length
 * nvtxs)
 * @param r_edgecut A reference to the weight of cut edges (can be NULL).
 *
 * @return MTMETIS_SUCCESS unless an error was encountered.
 */
int mtmetis_partition_explicit(
    mtmetis_vtx_t nvtxs,
    mtmetis_adj_t const * xadj,
    mtmetis_vtx_t const * adjncy,
    mtmetis_wgt_t const * vwgt,
    mtmetis_wgt_t const * adjwgt,
    double const * options,
    mtmetis_pid_t * where,
    mtmetis_wgt_t * r_edgecut);



/**
 * @brief Perform a nested disection on a graph, and output a permutation. 
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 * @param options The set of options.
 * @param perm The permutation of the vertices (rows+columns).
 *
 * @return MTMETIS_SUCCESS unless an error is encountered. 
 */
int mtmetis_nd(
    mtmetis_vtx_t const nvtxs,
    mtmetis_adj_t const * xadj,
    mtmetis_vtx_t const * adjncy,
    mtmetis_wgt_t const * vwgt,
    mtmetis_wgt_t const * adjwgt,
    mtmetis_pid_t * perm);


#ifdef __cplusplus
}
#endif



#endif

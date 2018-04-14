/**
 * @file ctrl.h
 * @brief Type and function prototypes for the ctrl structure.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-17
 */




#ifndef MTMETIS_CTRL_H
#define MTMETIS_CTRL_H




#include "base.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct timers_t {
  dl_timer_t total;
  dl_timer_t io;
  dl_timer_t partitioning;
  dl_timer_t coarsening;
  dl_timer_t matching;
  dl_timer_t contraction;
  dl_timer_t initpart;
  dl_timer_t uncoarsening;
  dl_timer_t projection;
  dl_timer_t refinement;
} timers_t;


typedef struct ctrl_t {
  /* runtime parameters */
  unsigned int seed;
  tid_t nthreads;
  int verbosity;
  int time;
  timers_t timers;
  /* partitioning parameters */
  pid_t nparts;
  size_t ncuts;
  real_t * tpwgts;
  real_t * pijbm;
  real_t ubfactor;
  /* coarsening parameters */
  int ctype;
  vtx_t coarsen_to;
  wgt_t maxvwgt;
  /* initial partitiong parameters */
  size_t ninitsolutions;
  /* refinement parameters */
  size_t nrefpass;
} ctrl_t;




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define ctrl_create __mtmetis_ctrl_create
/**
 * @brief Allocate and initialize a control structure.
 *
 * @return The new control structure.
 */
ctrl_t * ctrl_create(void);


#define ctrl_setup __mtmetis_ctrl_setup
/**
 * @brief Setup a control structure to partition a graph with a specified
 * number of vertices. The structure should already be configured with nthreads
 * and nparts.
 *
 * @param ctrl  The control structure to configure.
 * @param nvtxs The number of vertices in the graph to partition.
 */
void ctrl_setup(
    ctrl_t * ctrl,
    vtx_t nvtxs);


#define ctrl_parse __mtmetis_ctrl_parse
/**
 * @brief Create a control structure using the specified set of options.
 *
 * @param options The options specifying how to setup the ctrl.
 * @param ctrl A reference to the control structure to allocate and configure.
 *
 * @return MTMETIS_SUCCESS if the options array is valid.
 */
int ctrl_parse(
    double const * options,
    ctrl_t ** ctrl);


#define ctrl_free __mtmetis_ctrl_free
/**
 * @brief Free a control structure and its associated memory.
 *
 * @param ctrl The control structure to free.
 */
void ctrl_free(
    ctrl_t * ctrl);




#endif

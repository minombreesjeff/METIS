/**
 * @file ctrl.c
 * @brief Functions for allocating, freeing, and manipulating control
 * structures.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-17
 */




#ifndef MTMETIS_CTRL_C
#define MTMETIS_CTRL_C




#include <omp.h>
#include "ctrl.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const DEFAULT_NCUTS = 1; 
static size_t const DEFAULT_NREFPASS = 8;
static real_t const DEFAULT_UBFACTOR = 1.03;
static size_t const DEFAULT_NINITSOLUTIONS = 8;
static int const DEFAULT_CTYPE = MTMETIS_CTYPE_SHEM;
static int const DEFAULT_VERBOSITY = MTMETIS_VERBOSITY_NONE;




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


ctrl_t * ctrl_create(void)
{
  ctrl_t * ctrl;

  /* allocate my memory */
  ctrl = (ctrl_t*)calloc(1,sizeof(ctrl_t));

  ctrl->nthreads = omp_get_max_threads();
  ctrl->seed = (unsigned int)time(NULL);
  ctrl->ncuts = DEFAULT_NCUTS;
  ctrl->nrefpass = DEFAULT_NREFPASS;
  ctrl->ubfactor = DEFAULT_UBFACTOR;
  ctrl->ninitsolutions = DEFAULT_NINITSOLUTIONS;
  ctrl->ctype = DEFAULT_CTYPE;
  ctrl->verbosity = DEFAULT_VERBOSITY;
  ctrl->time = 0;

  return ctrl;
}


void ctrl_setup(
    ctrl_t * const ctrl,
    vtx_t const nvtxs)
{
  vtx_t i;
  pid_t nparts;
  real_t * tpwgts;

  nparts = ctrl->nparts;

  DL_ASSERT(nparts > 0,"Setting up ctrl with 0 parts");

  /* initialize tpwgts */
  ctrl->tpwgts = tpwgts = real_alloc(nparts);
  for (i=0;i<nparts;++i) {
    tpwgts[i] = 1.0 / nparts;
  }

  /* set various run parameters that depend on the graph */
  ctrl->coarsen_to = dl_min(PAR_COARSEN_FACTOR*nparts,
      dl_max(nvtxs/(20*pid_downlog2(nparts)),20*nparts));
}


int ctrl_parse(
    double const * const options,
    ctrl_t ** const r_ctrl)
{
  int rv;
  pid_t nparts;
  ctrl_t * ctrl = NULL;

  rv = MTMETIS_SUCCESS;

  ctrl = ctrl_create();

  /* decide how many threads to use */
  if (options[MTMETIS_OPTION_NTHREADS] != MTMETIS_VAL_OFF) {
    if (options[MTMETIS_OPTION_NTHREADS] < 1) {
      eprintf("Invalid number of threads: %"PF_TID_T"\n", \
          (tid_t)options[MTMETIS_OPTION_NTHREADS]);
      rv = MTMETIS_ERROR_INVALIDINPUT;
      goto CLEANUP;
    }
    ctrl->nthreads = (tid_t)options[MTMETIS_OPTION_NTHREADS];
  }

  /* check the number of partitions */
  if (options[MTMETIS_OPTION_NPARTS] == MTMETIS_VAL_OFF) {
    eprintf("The number of partitions must be specified.\n");
    rv = MTMETIS_ERROR_INVALIDINPUT;
    goto CLEANUP;
  } else if (options[MTMETIS_OPTION_NPARTS] < 1) {
    eprintf("The number of partitions must be at least 1.\n");
    rv = MTMETIS_ERROR_INVALIDINPUT;
    goto CLEANUP;
  } 
  nparts = (pid_t)options[MTMETIS_OPTION_NPARTS];

  ctrl->nparts = nparts;

  if (options[MTMETIS_OPTION_SEED] != MTMETIS_VAL_OFF) {
    ctrl->seed = (unsigned int)options[MTMETIS_OPTION_SEED];
  }

  if (options[MTMETIS_OPTION_NCUTS] != MTMETIS_VAL_OFF) {
    ctrl->ncuts = (size_t)options[MTMETIS_OPTION_NCUTS];
  }

  if (options[MTMETIS_OPTION_NITER] != MTMETIS_VAL_OFF) {
    ctrl->nrefpass = (size_t)options[MTMETIS_OPTION_NITER];
  }

  if (options[MTMETIS_OPTION_UBFACTOR] != MTMETIS_VAL_OFF) {
    ctrl->ubfactor = (real_t)options[MTMETIS_OPTION_UBFACTOR];
  }

  if (options[MTMETIS_OPTION_NINITSOLUTIONS] != MTMETIS_VAL_OFF) {
    ctrl->ninitsolutions = (size_t)options[MTMETIS_OPTION_NINITSOLUTIONS];
  }

  if (options[MTMETIS_OPTION_CTYPE] != MTMETIS_VAL_OFF) {
    ctrl->ctype = (int)options[MTMETIS_OPTION_CTYPE];
  }

  if (options[MTMETIS_OPTION_VERBOSITY] != MTMETIS_VAL_OFF) {
    ctrl->verbosity = (int)options[MTMETIS_OPTION_VERBOSITY];
  }

  if (options[MTMETIS_OPTION_TIME] != MTMETIS_VAL_OFF) {
    ctrl->time = 1;
    dl_init_timer(&ctrl->timers.total); 
    dl_init_timer(&ctrl->timers.io); 
    dl_init_timer(&ctrl->timers.partitioning); 
    dl_init_timer(&ctrl->timers.coarsening); 
    dl_init_timer(&ctrl->timers.matching); 
    dl_init_timer(&ctrl->timers.contraction); 
    dl_init_timer(&ctrl->timers.initpart); 
    dl_init_timer(&ctrl->timers.uncoarsening); 
    dl_init_timer(&ctrl->timers.projection); 
    dl_init_timer(&ctrl->timers.refinement); 
  }

  *r_ctrl = ctrl;
  ctrl = NULL;

  CLEANUP:

  if (ctrl) {
    ctrl_free(ctrl);
  }

  return rv;
}


void ctrl_free(
    ctrl_t * ctrl)
{
  #pragma omp barrier
  #pragma omp master
  {
    if (ctrl->tpwgts) {
      dl_free(ctrl->tpwgts);
    }
    if (ctrl->pijbm) {
      dl_free(ctrl->pijbm);
    }
    dl_free(ctrl);
  }
}



#endif

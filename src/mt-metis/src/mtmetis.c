/**
 * @file mtmetis.c
 * @brief Library entry points
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-07-01
 */

#include "includes.h"

int mtmetis_partkway(const idx_t * const r_nvtxs, 
    const idx_t * const xadj, const idx_t * const adjncy, 
    const idx_t * const vwgt, const idx_t * const adjwgt,
    const idx_t * const r_nparts, idx_t * const where, idx_t * const r_edgecut)
{
  dgraph_t * graph;
  idx_t nthreads = omp_get_max_threads();

  idx_t buffer[nthreads+1];
  dctrl_t * dctrl = AllocateDCtrl(*r_nvtxs,*r_nparts,nthreads); 

  #pragma omp parallel shared(graph,dctrl,buffer) num_threads(nthreads) \
    default(none)
  {
    idx_t cut;
    ParSetupGraph(&graph,buffer,*r_nvtxs,xadj,adjncy,adjwgt,vwgt);

    cut = ParPartGraphKway(dctrl,graph,where);

    ParFreeGraph(&graph);

    if (r_edgecut) {
      #pragma omp master
      {
        *r_edgecut = cut;
      }
    }
  }

  FreeDCtrl(&dctrl,nthreads);

  return 1;
}



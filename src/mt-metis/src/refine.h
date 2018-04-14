/**
 * @file refine.h
 * @brief Refinement function prototypes.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-19
 */




#ifndef MTMETIS_REFINE_H
#define MTMETIS_REFINE_H




#include "base.h"
#include "graph.h"
#include "ctrl.h"
#include "ucinfo.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/



/**
* @brief Parallel kway-refinement
*
* @param ctrl control strucutre
* @param graph the graph who's partition to refine
* @param niter number of iterations of refinement to perform
* @param ffactor Fudge factor (allow positive gains at teh cost of balance)
* @param ucinfo The uncoarsening informatino struct. 
*
* @return Total of moved vertices.
*/
vtx_t refine_kway(
    ctrl_t * const ctrl, 
    graph_t * const graph,
    size_t niter, 
    real_t ffactor,
    ucinfo_t * ucinfo);





#endif

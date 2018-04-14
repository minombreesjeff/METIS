/**
 * @file kwayrefine.h
 * @brief KWay refinement function prototypes.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-19
 */




#ifndef MTMETIS_KWAYREFINE_H
#define MTMETIS_KWAYREFINE_H




#include "base.h"
#include "graph.h"
#include "ctrl.h"
#include "kwinfo.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define par_kwawyrefine __mtmetis_par_kwayrefine
/**
* @brief Parallel kway-refinement.
*
* @param ctrl The control strucutre.
* @param graph The graph who's partition to refine.
* @param kwinfo The uncoarsening information structure. 
*
* @return Total of moved vertices.
*/
vtx_t par_kwayrefine(
    ctrl_t * const ctrl, 
    graph_t * const graph,
    kwinfo_t * kwinfo);




#endif

/**
 * @file vseprefine.h
 * @brief Vertex separator refinement prototypes.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-12-09
 */




#ifndef MTMETIS_VSEPREFINE_H
#define MTMETIS_VSEPREFINE_H




#include "base.h"
#include "graph.h"
#include "ctrl.h"
#include "vsinfo.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define par_vseprefine __mtmetis_par_vseprefine
/**
* @brief Parallel vertex separator refinement.
*
* @param ctrl The control strucutre
* @param graph The graph who's partition to refine
* @param kwinfo The uncoarsening information struct. 
*
* @return Total of moved vertices.
*/
vtx_t par_vseprefine(
    ctrl_t * const ctrl, 
    graph_t * const graph,
    vsinfo_t * vsinfo);




#endif

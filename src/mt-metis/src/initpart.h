/**
 * @file initpart.h
 * @brief Parallel initial partitioning function prototypes.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-17
 */





#ifndef MTMETIS_INITPART_H
#define MTMETIS_INITPART_H




#include "base.h"
#include "graph.h"
#include "ctrl.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define initpart_kway __mtmetis_initpart_kway
/**
 * @brief Create a kway partitioning of a coarsened graph.
 *
 * @param ctrl The control structure with runtime parameters.
 * @param graph The coarse graph to partition.
 *
 * @return The edgecut of the new partitioning.  
 */
vtx_t initpart_kway(
    ctrl_t * ctrl,
    graph_t * graph);




#endif

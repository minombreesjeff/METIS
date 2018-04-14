/**
 * @file uncoarsen.h
 * @brief Refinement functions and structs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-18
 */




#ifndef MTMETIS_UNCOARSEN_H
#define MTMETIS_UNCOARSEN_H




#include "base.h"
#include "graph.h"
#include "ctrl.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


/**
 * @brief Uncoarsen a graph with a kway partition.
 *
 * @param ctrl The control structure containing runtime parameters.
 * @param ograph The original graph to uncoarsen the partition to.
 * @param cgraph The partitioned coarse graph.
 */
void uncoarsen_kway(
    ctrl_t * ctrl,
    graph_t * ograph,
    graph_t * cgraph);




#endif

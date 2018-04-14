/**
 * @file internal.h
 * @brief Top level internal functions for mtmetis
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2015-01-17
 */




#ifndef MTMETIS_INTERNAL_H
#define MTMETIS_INTERNAL_H




#include "base.h"
#include "ctrl.h"
#include "graph.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


/**
 * @brief Innermost function for spawning parallelism for kway partitions.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * @param where The output partition IDs for each vertex.
 */
void mtmetis_partition_kway_int(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * where);


/**
 * @brief Innermost function for spawning parallelism for recursive bisection.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * @param where The output partition IDs for each vertex.
 */
void mtmetis_partition_rb_int(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * where);


/**
 * @brief Innermost function for spawning parallelism for edge separators.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * @param where The output partition IDs for each edge.
 */
void mtmetis_partition_esep_int(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * where);


/**
 * @brief Innermost function for spawning parallelism for vertex separators.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * @param where The output partition IDs for each vertex.
 */
void mtmetis_partition_vsep_int(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * where);




#endif

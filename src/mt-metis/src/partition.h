/**
 * @file partition.h
 * @brief Partitioning function prototypes
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-17
 */




#ifndef MTMETIS_PARTITION_H
#define MTMETIS_PARTITION_H




#include "base.h"
#include "ctrl.h"
#include "graph.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define partition_kway __mtmetis_partition_kway
/**
 * @brief Entry level function for multithreaded kway partitioning. Should be
 * called by all threads in a  parallel region.
 *
 * @param ctrl The control structure.
 * @param graph The graph to partition.
 * @param where The allocated where vector.
 */
void partition_kway(
    ctrl_t * ctrl, 
    graph_t * graph, 
    pid_t ** where);


#define partition_print_info __mtmetis_print_info
/**
 * @brief Print information about a partition.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param where The where vector.
 */
void partition_print_info(
    ctrl_t const * ctrl,
    graph_t const * graph,
    pid_t const * const * where);




#endif

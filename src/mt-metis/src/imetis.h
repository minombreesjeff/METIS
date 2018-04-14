/**
 * @file imetis.h
 * @brief Metis wrappers. 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 * @date 2015-06-08
 */




#ifndef MTMETIS_IMETIS_H
#define MTMETIS_IMETIS_H



#include "ctrl.h"
#include "graph.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


/**
 * @brief Generate an initial partitioning.
 *
 * @param ctrl The control structure with parameters.
 * @param nparts The number of parittions in the partitioning.
 * @param tpwgts The target partition weights.
 * @param ncuts The number of partitionings to make.
 * @param rb Use recursive bisection to generate k-way partitionings.
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The adjacecny weights.
 * @param where The partition ids of each vertex (output).
 *
 * @return The weight of the separator.
 */
wgt_t metis_initcut(
    ctrl_t * const ctrl,
    pid_t const nparts,
    real_t * tpwgts,
    size_t const ncuts,
    int const rb,
    vtx_t nvtxs,
    adj_t * const xadj,
    vtx_t * const adjncy,
    wgt_t * const vwgt,
    wgt_t * const adjwgt,
    pid_t * const where);


/**
 * @brief Generate an initial vertex separator using metis.
 *
 * @param ctrl The control structure with parameters.
 * @param nseps The number of separators to generate.
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The adjacecny weights.
 * @param where The partition ids of each vertex (output).
 *
 * @return The weight of the separator.
 */
wgt_t metis_initsep(
    ctrl_t * ctrl,
    size_t nseps,
    vtx_t nvtxs,
    adj_t * xadj,
    vtx_t * adjncy,
    wgt_t * vwgt,
    wgt_t * adjwgt,
    pid_t * where);


/**
 * @brief Serially generate a k-way partition using metis (direct k-way
 * parittioning).
 *
 * @param ctrl The control structure containing partitioning parameters.
 * @param graph The graph to partition.
 * @param where The partition id for each vertex (output).
 * @param rb Whether or not to use recursive bisection. 
 *
 * @return The weight of the edgecut. 
 */
wgt_t metis_kway(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * const * where,
    int rb);


/**
 * @brief Serially generate a 2-way edge separator using metis.
 *
 * @param ctrl The control structure containing partitioning parameters.
 * @param graph The graph to partition.
 * @param where The partition id for each vertex (output).
 *
 * @return The weight of the edge separator.
 */
wgt_t metis_esep(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * const * where);


/**
 * @brief Serially generate a 2-way vertex separator using metis.
 *
 * @param ctrl The control structure containing partitioning parameters.
 * @param graph The graph to partition.
 * @param where The partition id for each vertex (output).
 *
 * @return The weight of the vertex separator.
 */
wgt_t metis_vsep(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * const * where);




#endif

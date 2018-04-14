/**
 * @file flow.h
 * @brief Function for computing maximum-flows / minimum cuts
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Dominique LaSalle
 * @version 1
 * @date 2014-10-14
 */




#ifndef BOWSTRING_FLOW_H
#define BOWSTRING_FLOW_H



#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define maxflow_edge __bowstring_maxflow_edge
/**
 * @brief Perform a maximum flow calculation (assuming vertices have infinite
 * capacity and edges have capacities equal to their weights).
 *
 * @param src The vertex to serve as the source.
 * @param dst The vertex to serve as teh destination (sink).
 * @param nvtxs THe number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights (capacity).
 * @param flow The flow on each edge (output).
 *
 * @return The value of the maximum flow.
 */
wgt_t maxflow_edge(
    vtx_t src,
    vtx_t dst,
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt,
    wgt_t * flow);


#define maxflow_vertex __bowstring_maxflow_vertex
/**
 * @brief Perform a maximum flow calculation (assuming edges have infinite
 * capacity and vertices have capacities equal to their weights).
 *
 * @param src The vertex to serve as the source.
 * @param dst The vertex to serve as teh destination (sink).
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy THe adjacency list.
 * @param vwgt The vertex weights (capacity).
 * @param flow THe flow on each vertex (output). 
 *
 * @return The value of the maximum flow.
 */
wgt_t maxflow_vertex(
    vtx_t src,
    vtx_t dst,
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * vwgt,
    wgt_t * flow);




#endif

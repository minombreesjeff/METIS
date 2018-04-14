/**
 * @file order.h
 * @brief Graph re-ording utilities
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-07-18
 */




#ifndef PERMUTE_H
#define PERMUTE_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define order_graph __bowstring_order_graph
/**
 * @brief Re-order a graph given a permutateion.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer array.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 * @param perm The permuation.
 */
void order_graph(
    vtx_t nvtxs, 
    adj_t * xadj, 
    vtx_t * adjncy, 
    wgt_t * vwgt, 
    wgt_t * adjwgt, 
    const vtx_t * perm);


#define order_permutation __bowstring_order_permutation
/**
 * @brief Generate a vertex permutation vector of a graph using the given
 * method. 
 *
 * @param ordering The type of ordering to create.
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 * @param perm The vertex permutation.
 */
void order_permutation(
    int ordering,
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * vwgt, 
    wgt_t const * adjwgt,
    vtx_t * perm);






#endif

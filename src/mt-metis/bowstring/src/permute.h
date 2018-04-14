/**
 * @file permute.h
 * @brief Graph re-ording utilities
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-07-18
 */




#ifndef PERMUTE_H
#define PERMUTE_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define permute_graph __bowstring_permute_graph
/**
 * @brief Permute a graph using a random permutation. 
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 */
void permute_graph(
    vtx_t nvtxs, 
    adj_t * xadj, 
    vtx_t * adjncy, 
    wgt_t * vwgt, 
    wgt_t * adjwgt);


void reorder_graph(
    vtx_t nvtxs, 
    adj_t * xadj, 
    vtx_t * adjncy, 
    wgt_t * vwgt, 
    wgt_t * adjwgt, 
    const vtx_t * perm);




#endif

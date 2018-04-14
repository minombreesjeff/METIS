/**
 * @file generate.h
 * @brief Function prototypes for generating graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-10-19
 */




#ifndef BOWSTRING_GENERATE_H
#define BOWSTRING_GENERATE_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define generate_complete_graph __bowstring_generate_complete_graph
/**
 * @brief Generate a complete graph of specified size.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param r_xadj A reference to the adjacency pointer.
 * @param r_adjncy A reference to the adjacecny list.
 * @param r_vwgt A reference to the vertex weights.
 * @param r_adjwgt A reference to the edge weights.
 *
 */
void generate_complete_graph(
    vtx_t nvtxs, 
    adj_t ** r_xadj, 
    vtx_t ** r_adjncy,
    wgt_t ** r_vwgt, 
    wgt_t ** r_adjwgt);


#define generate_grid_graph __bowstring_generate_grid_graph
/**
 * @brief Generate a grid graph in up to three dimensions.
 *
 * @param nvx The number of vertices in the 1st dimension.
 * @param nvy The number of vertices in the 2nd dimension.
 * @param nvz The number of vertices in the 3rd dimension.
 * @param r_xadj A reference to the adjacency list pointer.
 * @param r_adjncy A reference to the adjacency list.
 * @param r_vwgt A reference to the vertex weights.
 * @param r_adjwgt A reference to the edge weights.
 *
 */
void generate_grid_graph(
    vtx_t nvx, 
    vtx_t nvy, 
    vtx_t nvz, 
    adj_t ** r_xadj, 
    vtx_t ** r_adjncy, 
    wgt_t ** r_vwgt, 
    wgt_t ** r_adjwgt);


#ifdef XXX
int generate_powerlaw_graph(
    vtx_t nvtxs, 
    adj_t nedges,);
#endif




#endif


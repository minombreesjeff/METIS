/**
 * @file cut.h
 * @brief Function prototypes for cutting graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2014-01-26
 */




#ifndef BOWSTRING_CUT_H
#define BOWSTRING_CUT_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define coordinate_bisection __bowstring_coordinate_bisection
/**
 * @brief Performa a coordinate-based bisection of vertices. 
 *
 * @param nvtxs The number of vertices in the graph.
 * @param vwgt The weight of the vertices (may be NULL).
 * @param x The coordinate of the vertices to bisect along.
 * @param lwgt The weight of the left partition.
 * @param where The vector where partition ID's will be placed.
 *
 * @return 
 */
void coordinate_bisection(
    vtx_t nvtxs, 
    const wgt_t * vwgt, 
    const coord_t * x, 
    wgt_t lwgt, 
    vlbl_t * where);


#define recursive_coordinate_bisection \
    __bowstring_recursive_coordinate_bisection
/**
 * @brief Perform RCB.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param vwgt The weight of the vertices.
 * @param x The coordinate information.
 * @param ndim The number of dimensions.
 * @param twgts The target partition weights.
 * @param nparts The number of partitions.
 * @param where The vector where partition ID's will be placed.
 *
 * @return  
 */
void recursive_coordinate_bisection(
    vtx_t nvtxs, 
    const wgt_t * vwgt, 
    const coord_t ** x, 
    size_t ndim, 
    const wgt_t * twgts, 
    vlbl_t nparts,
    vlbl_t * where);




#endif

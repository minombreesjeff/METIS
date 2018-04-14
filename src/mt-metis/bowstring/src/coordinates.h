/**
 * @file coordinates.h
 * @brief Function prototypes for assigning coordinates
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2014-01-26
 */




#ifndef BOWSTRING_COORDINATES_H
#define BOWSTRING_COORDINATES_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int bfs_coordinates(vtx_t nvtxs, const adj_t * xadj, const vtx_t * adjncy,
    const wgt_t * adjwgt, size_t ndim, coord_t ** coords);




#endif

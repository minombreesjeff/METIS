/**
 * @file dimacs.h
 * @brief Function prototypes for reading and writing graphs in DIMACS format
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-02-09
 */




#ifndef BOWSTRING_IO_DIMACS_H
#define BOWSTRING_IO_DIMACS_H




#include "io.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int read_dimacs_graph(const char * filename, vtx_t * nvtxs, adj_t ** xadj, 
    vtx_t ** adjncy, wgt_t ** vwgt, wgt_t ** adjwgt);


int write_dimacs_graph(const char * filename, vtx_t nvtxs, const adj_t * xadj,
    const vtx_t * adjncy, const wgt_t * vwgt, const wgt_t * adjwgt);




#endif

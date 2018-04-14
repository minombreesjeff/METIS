/**
 * @file metis.h
 * @brief Function prototypes for reading and writing graphs in METIS format
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-02-09
 */




#ifndef BOWSTRING_IO_METIS_H
#define BOWSTRING_IO_METIS_H




#include "io.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int read_metis_graph(const char * filename, vtx_t * nvtxs, adj_t ** xadj, 
    vtx_t ** adjncy, wgt_t ** vwgt, wgt_t ** adjwgt);


int write_metis_graph(const char * filename, vtx_t nvtxs, const adj_t * xadj,
    const vtx_t * adjncy, const wgt_t * vwgt, const wgt_t * adjwgt);


int read_metis_partition(const char * filename, vtx_t * nvtxs, 
    vlbl_t ** labels);


int write_metis_partition(const char * filename, vtx_t nvtxs, 
    const vlbl_t * labels);




#endif

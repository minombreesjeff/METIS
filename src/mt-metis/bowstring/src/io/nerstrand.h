/**
 * @file nerstrand.h
 * @brief Function prototypes for reading and writing graphs in NERSTRAND 
 * format.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-02-09
 */




#ifndef BOWSTRING_IO_NERSTRAND_H
#define BOWSTRING_IO_NERSTRAND_H




#include "io.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int read_nerstrand_graph(
    const char * filename, 
    vtx_t * r_nvtxs, 
    adj_t ** xadj,
    vtx_t ** adjncy, 
    wgt_t ** vwgt, 
    wgt_t ** adjwgt, 
    int binary);


int write_nerstrand_graph(
    const char * filename, 
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    const wgt_t * vwgt, 
    const wgt_t * adjwgt, 
    int binary);


int read_nerstrand_clustering(
    const char * filename, 
    vtx_t * nvtxs, 
    vlbl_t ** labels,
    int binary);


int write_nerstrand_clustering(
    const char * filename, 
    vtx_t nvtxs, 
    const vlbl_t * labels,
    int binary);




#endif

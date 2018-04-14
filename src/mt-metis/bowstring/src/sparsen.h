/**
 * @file sparsen.h
 * @brief Edge removal function prototypes
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-09-25
 */




#ifndef SPARSEN_H
#define SPARSEN_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


adj_t prune_ranked_edges(
    vtx_t nvtxs, 
    const adj_t * gxadj, 
    const vtx_t * gadjncy, 
    const wgt_t * gadjwgt, 
    const elbl_t * rank, 
    elbl_t maxrank, 
    double frac, 
    adj_t ** r_xadj, 
    vtx_t ** r_adjncy, 
    wgt_t ** r_adjwgt, 
    int reweight);


elbl_t build_nirank(
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    const wgt_t * adjwgt, 
    elbl_t * rank);


elbl_t build_mstrank(
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    const wgt_t * adjwgt, 
    elbl_t * rank);


elbl_t build_astrank(
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    const wgt_t * adjwgt, 
    elbl_t * rank);


elbl_t build_lstrank(
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    const wgt_t * adjwgt, 
    elbl_t * rank);




#endif

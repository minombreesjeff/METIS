/**
 * @file contract.h
 * @brief Functions for contracting a graph.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2015-01-21
 */




#ifndef MTMETIS_CONTRACT_H
#define MTMETIS_CONTRACT_H



#include "base.h"
#include "graph.h"
#include "ctrl.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


void par_contract_graph(
    ctrl_t * ctrl, 
    graph_t * graph, 
    vtx_t mycnvtxs, 
    vtx_t const * const * gmatch, 
    vtx_t const * fcmap);




#endif

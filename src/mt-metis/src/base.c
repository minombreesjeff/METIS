/**
 * @file base.c
 * @brief Base type specific functions.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-17
 */




#ifndef MTMETIS_BASE_C
#define MTMETIS_BASE_C




#include "base.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


DL_MK_OMP_REDUCTION_FUNCS(vtx,vtx_t)
DL_MK_OMP_REDUCTION_FUNCS(adj,adj_t)
DL_MK_OMP_REDUCTION_FUNCS(wgt,wgt_t)
DL_MK_OMP_REDUCTION_FUNCS(double,double)


#endif

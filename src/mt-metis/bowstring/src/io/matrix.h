/**
 * @file matrix.h
 * @brief Functions for reading and writing dense matrices
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2014-02-13
 */




#ifndef BOWSTRING_IO_MATRIX_H
#define BOWSTRING_IO_MATRIX_H




#include "io.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


int write_dense_matrix(const char * filename, const coord_t * const * coords, 
    vtx_t n, size_t m);




#endif

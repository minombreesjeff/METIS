/**
 * @file matrix.c
 * @brief Functions for reading and writing dense matrices
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014
 * @version 1
 * @date 2014-02-13
 */




#ifndef BOWSTRING_IO_MATRIX_C
#define BOWSTRING_IO_MATRIX_C



#include "matrix.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int write_dense_matrix(const char * const filename, 
    const coord_t * const * const coords, const vtx_t n, const size_t m)
{
  int rv;
  file_t * file;
  vtx_t i;
  size_t d;


  if ((rv = __open_file(filename,"w",&file)) != BOWSTRING_SUCCESS) {
    return rv;
  }

  for (i=0;i<n;++i) {
    for (d=0;d<m;++d) {
      dl_fprintf(file,PF_COORD_T" ",coords[d][i]);
    }
    dl_fprintf(file,"\n");
  }

  dl_close_file(file);

  return BOWSTRING_SUCCESS; 
}




#endif


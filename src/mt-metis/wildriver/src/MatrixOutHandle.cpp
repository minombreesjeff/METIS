/**
 * @file MatrixOutHandle.cpp
 * @brief Class for writing all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */





#ifndef SRC_MATRIXOUTHANDLE_CPP
#define SRC_MATRIXOUTHANDLE_CPP




#include "MatrixOutHandle.hpp"
#include "MatrixFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


MatrixOutHandle::MatrixOutHandle(
    std::string const & fname)
{
  writer = MatrixFactory::OpenFile(fname);
}


MatrixOutHandle::~MatrixOutHandle()
{
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixOutHandle::setInfo(
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz)
{
  writer->setInfo(nrows,ncols,nnz);
}


void MatrixOutHandle::writeSparse(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  writer->write(rowptr,rowind,rowval);
}



}



#endif

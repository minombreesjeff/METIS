/**
 * @file MatrixHandle.cpp
 * @brief Class for reading all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_MATRIXINHANDLE_CPP
#define SRC_MATRIXINHANDLE_CPP




#include "MatrixInHandle.hpp"
#include "MatrixFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


MatrixInHandle::MatrixInHandle(
    std::string const & fname)
{

  reader = MatrixFactory::OpenFile(fname);
}


MatrixInHandle::~MatrixInHandle()
{
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixInHandle::getInfo(
    dim_t & nrows,
    dim_t & ncols,
    ind_t & nnz)
{
  reader->getInfo(nrows,ncols,nnz);
}


void MatrixInHandle::readSparse(
    ind_t * rowptr,
    dim_t * rowind,
    val_t * rowval)
{
  reader->read(rowptr,rowind,rowval);
}


dim_t MatrixInHandle::getNextRow(
        std::vector<MatrixEntry> & next)
{
  return reader->getNextRow(next);
}




}




#endif


/**
 * @file MatrixFile.cpp
 * @brief Implementation of the base abstract class for matrix files. 
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-05
 */




#include "MatrixFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MatrixFile::MatrixFile() :
  nnz(NULL_IND),
  info_set(false)
{
  // do nothing
}


MatrixFile::~MatrixFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixFile::getInfo(
    dim_t & nrows,
    dim_t & ncols,
    ind_t & nnz)
{
  // see if need to read the header
  if (!info_set) {
    readHeader();
  }

  // set values
  nrows = getNumRows();
  ncols = getNumCols();
  nnz = getNNZ();

  info_set = true;
}


void MatrixFile::setInfo(
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz)
{
  setNumRows(nrows);
  setNumCols(ncols);
  setNNZ(nnz);

  info_set = true;

  writeHeader();
}




}

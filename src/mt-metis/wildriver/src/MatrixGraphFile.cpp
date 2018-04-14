/**
 * @file MatrixGraphFile.cpp
 * @brief An adapter class for treating matrices as graphs.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 *
 */




#ifndef WILDRIVER_MATRIXGRAPHFILE_CPP
#define WILDRIVER_MATRIXGRAPHFILE_CPP




#include "MatrixGraphFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MatrixGraphFile::MatrixGraphFile(
    IMatrixFile * const mfile) :
  file(mfile)
{
}


MatrixGraphFile::~MatrixGraphFile()
{
  delete file;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/

void MatrixGraphFile::getInfo(
    dim_t & nvtxs,
    ind_t & nedges,
    int & nvwgt,
    bool & ewgts)
{
  dim_t nrows, ncols;
  ind_t nnz;

  file->getInfo(nrows,ncols,nnz);

  nvtxs = nrows;
  nedges = nnz;
  
  // matrices never have vertex weights and always have edge weights
  nvwgt = 0;
  ewgts = true;
}


void MatrixGraphFile::setInfo(
    dim_t const nvtxs,
    ind_t const nedges,
    int const nvwgt,
    bool const ewgts)
{
  if (nvwgt != 0) {
    throw BadParameterException("MatrixGraphFile cannot handle vertex " \
        "weights");
  }

  // TODO: add an internal state for edge weights, such that they can be
  // surpessed via this function.

  file->setInfo(nvtxs,nvtxs,nedges);
}


void MatrixGraphFile::read(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt)
{
  dim_t nrows, ncols;
  ind_t nnz;

  file->getInfo(nrows,ncols,nnz);

  // fill vertex weights with ones if requested
  if (vwgt) {
    for (dim_t i=0; i<nrows; ++i) {
      vwgt[i] = 1;
    }
  }

  file->read(xadj,adjncy,adjwgt);
}


void MatrixGraphFile::write(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  file->write(xadj,adjncy,adjwgt);
}


void MatrixGraphFile::firstVertex()
{
  file->firstRow();
}


bool MatrixGraphFile::getNextVertex(
    std::vector<val_t> & vwgt,
    std::vector<MatrixEntry> & list)
{
  return file->getNextRow(list);
}


void MatrixGraphFile::setNextVertex(
    std::vector<val_t> const & vwgts,
    std::vector<MatrixEntry> const & list)
{
  file->setNextRow(list);
}


std::string const & MatrixGraphFile::getName() const noexcept
{
  return file->getName();
}


}




#endif

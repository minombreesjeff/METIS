/**
 * @file MatrixTextFile.cpp
 * @brief Abstract providing common code for matrix reading/writing classes.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_MATRIXFILE_CPP
#define SRC_MATRIXFILE_CPP




#include "MatrixTextFile.hpp"
#include "MetisFile.hpp"
#include "CSRFile.hpp"




namespace WildRiver {


/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


dim_t MatrixTextFile::incRow()
{
  if (current_row >= getNumRows()) {
    throw BadFileStateException(std::string("Attempt to increase current " \
          "row beyond the number of rows: ") + std::to_string(getNumRows()));
  }

  ++current_row;

  return current_row;
}




/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MatrixTextFile::MatrixTextFile(
    std::string const & name) :
  TextFile(name),
  MatrixFile(),
  current_row(0),
  read_nnz(0)
{
  // do nothing
}


MatrixTextFile::~MatrixTextFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixTextFile::read(
    ind_t * const rowptr,
    dim_t * const rowind,
    val_t * const rowval)
{
  std::vector<MatrixEntry> row;

  dim_t nrows = getNumRows();

  // read in the rows the matrix into our csr
  dim_t j = 0;
  rowptr[0] = j;
  for (dim_t i=0; i<nrows; ++i) {
    if (!getNextRow(row)) {
      // fewer rows than expected
      throw EOFException(std::string("Failed to read row ") + \
          std::to_string(i) + std::string("/") + std::to_string(nrows)); 
    }

    // populate the row
    for (MatrixEntry entry : row) {
      if (j >= getNNZ()) {
        throw BadFileException(std::string("Found more than ") + \
            std::to_string(getNNZ()) + std::string(" non-zeroes in file."));
      }

      rowind[j] = entry.ind;
      if (rowval) {
        rowval[j] = entry.val;
      }
      ++j;
    }
    rowptr[i+1] = j;
  }

  if (j != getNNZ()) {
    // we read in the wrong number of non-zeroes
    throw EOFException(std::string("Only found ") + std::to_string(j) + \
        std::string("/") + std::to_string(getNNZ()) + \
        std::string(" non-zeroes in file"));
  }
}


void MatrixTextFile::write(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  std::vector<MatrixEntry> row;

  dim_t const nrows = getNumRows();

  for (dim_t i=0;i<nrows;++i) {
    // build and insert a new row
    row.clear();
    for (ind_t j=rowptr[i];j<rowptr[i+1];++j) {
      MatrixEntry entry;

      entry.ind = rowind[j];
      if (rowval) {
        entry.val = rowval[j];
      } else {
        entry.val = 1;
      }
      row.push_back(entry);
    }
    setNextRow(row);
  }
}




}




#endif

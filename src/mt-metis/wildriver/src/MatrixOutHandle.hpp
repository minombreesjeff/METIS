/**
 * @file MatrixOutHandle.hpp
 * @brief Class for writing all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_MATRIXOUTHANDLE_HPP
#define SRC_MATRIXOUTHANDLE_HPP




#include <vector>
#include <memory>

#include "IMatrixWriter.hpp"




namespace WildRiver 
{


class MatrixOutHandle
{
  private:
    std::shared_ptr<IMatrixWriter> writer;

    // disable copying
    MatrixOutHandle(
        MatrixOutHandle const & handle);
    MatrixOutHandle & operator=(
        MatrixOutHandle const & handle);


  public:
    MatrixOutHandle(
        std::string const & fname);


    ~MatrixOutHandle();


    void setInfo(
        dim_t nrows,
        dim_t ncols,
        ind_t nnz);


    void writeSparse(
        ind_t const * rowptr,
        dim_t const * rowind,
        val_t const * rowval);


    dim_t setNextRow(
        std::vector<MatrixEntry> const & next);


};




}




#endif

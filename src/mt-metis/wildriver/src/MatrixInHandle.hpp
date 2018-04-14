/**
 * @file MatrixInHandle.hpp
 * @brief Class for reading all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_MATRIXINHANDLE_HPP
#define SRC_MATRIXINHANDLE_HPP




#include <vector>
#include <memory>

#include "IMatrixReader.hpp"




namespace WildRiver 
{


class MatrixInHandle
{
  private:
    std::shared_ptr<IMatrixReader> reader;

    // disable copying
    MatrixInHandle(
        MatrixInHandle const & handle);
    MatrixInHandle & operator=(
        MatrixInHandle const & handle);


  public:
    /**
     * @brief Create a new file handle for reading matrices.
     *
     * @param fname The filename/path of the file to read.
     */
    MatrixInHandle(
        std::string const & fname);


    /**
     * @brief Get the number of rows, columns, and non-zeros in the matrix.
     *
     * @param nrows The number of rows.
     * @param ncols The number of columns.
     * @param nnz THe number of non-zeros.
     */
    void getInfo(
        dim_t & nrows,
        dim_t & ncols,
        ind_t & nnz);


    /**
     * @brief Get the sparse matrix in CSR form. The pointers must be
     * pre-allocated to the sizes required by the info of the matrix 
     *
     * |rowptr| = nrows + 1
     * |rowind| = nnz
     * |rowval| = nnz
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     */
    void readSparse(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval);


    /**
     * @brief Set the next row in the matrix file.
     *
     * @param next The row to set.
     *
     * @return The number of rows in the matrix. 
     */
    dim_t getNextRow(
        std::vector<MatrixEntry> & next);


    /**
     * @brief Close the handle.
     */
    ~MatrixInHandle();
};




}




#endif

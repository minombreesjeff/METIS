/**
 * @file IMatrixReader.hpp
 * @brief Interface for reading sparse matrices.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef WILDRIVER_IMATRIXREADER_HPP
#define WILDRIVER_IMATRIXREADER_HPP




#include <string>
#include <vector>

#include "base.h"
#include "MatrixEntry.hpp"
#include "Exception.hpp"




namespace WildRiver 
{


class IMatrixReader
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IMatrixReader() 
    {
    }


    /**
     * @brief Get the matrix in CSR form. The pointers must be
     * pre-allocated to the sizes required by the info of the matrix.
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
    virtual void read(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval) = 0;


    /**
     * @brief Get the number of rows, columns, and non-zeros in the matrix.
     *
     * @param nrows The number of rows.
     * @param ncols The number of columns.
     * @param nnz THe number of non-zeros.
     */
    virtual void getInfo(
        dim_t & nrows,
        dim_t & ncols,
        ind_t & nnz) = 0;


    /**
     * @brief Reset the current position in the matrix file to the first row.
     */
    virtual void firstRow() = 0;


    /**
     * @brief Get the next row in the matrix.
     *
     * @param next The row to set.
     *
     * @return True if the row was successfully read, and false if EOF is 
     * reached. 
     */
    virtual bool getNextRow(
        std::vector<MatrixEntry> & next) = 0;


    /**
     * @brief Get the name of this matrix file type.
     *
     * @return The matrix file type name.
     */
    virtual std::string const & getName() const noexcept = 0;


};




}




#endif

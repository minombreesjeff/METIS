/**
 * @file Matrix.hpp
 * @brief Base abstract class for matrices.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-06
 */




#ifndef WILDRIVER_MATRIX_HPP
#define WILDRIVER_MATRIX_HPP




#include "base.h"




namespace WildRiver
{


class Matrix
{
  private:
    /**
     * @brief The number of rows in the matrix.
     */
    dim_t nrows;


    /**
     * @brief The number of columns in the matrix.
     */
    dim_t ncols;


  public:
    /**
     * @brief Default contructor. Initializes the number of rows and number of
     * columns to invalid values (NULL_DIM).
     */
    Matrix() :
      nrows(NULL_DIM),
      ncols(NULL_DIM)
    {
    }


    /**
     * @brief Virtual destructor.
     */
    virtual ~Matrix()
    {
    }


    /**
     * @brief Get the number of rows in this matrix.
     *
     * @return The number of rows in this matrix.
     */
    inline dim_t getNumRows() const noexcept
    {
      return nrows;
    }


    /**
     * @brief Get the number of columns in this matrix.
     *
     * @return The number of columns in the matrix.
     */
    inline dim_t getNumCols() const noexcept
    {
      return ncols;
    }


    /**
     * @brief Set the number of rows in this matrix.
     *
     * @param n The new number of row sin tha matrix.
     */
    inline void setNumRows(
        dim_t const n) noexcept
    {
      nrows = n;
    }


    /**
     * @brief Set the number of columns in this matrix.
     *
     * @param n The new number of columns.
     */
    inline void setNumCols(
        dim_t const n) noexcept
    {
      ncols = n;
    }


};




}




#endif

/**
 * @file MatrixEntry.hpp
 * @brief Structure a matrix entry.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_MATRIXENTRY_HPP
#define SRC_MATRIXENTRY_HPP




#include "base.h"




namespace WildRiver
{


/**
 * @brief Structure for holding the column index and value of an entry in a
 * sparse matrix.
 */
struct MatrixEntry 
{
  /**
   * @brief Column index of the entry.
   */
  dim_t ind;


  /**
   * @brief Numerical value of the entry.
   */
  val_t val;


};


}




#endif

/**
 * @file IVectorWriter.hpp
 * @brief Interface for writing out vectors.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_IVECTORWRITER_HPP
#define WILDRIVER_IVECTORWRITER_HPP




#include <vector>

#include "base.h"




namespace WildRiver
{


class IVectorWriter
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IVectorWriter()
    {
    }


    /**
     * @brief Set the size of the vector.  
     *
     * @param size The new size of the vector.
     */
    virtual void setSize(
        size_t size) = 0;


    /**
     * @brief Write the vector to the underlying medium. 
     *
     * @param vals The dense array of values in the vector.
     */
    virtual void write(
        val_t const * vals) = 0;


    /**
     * @brief Write the vector to the underlying medium. This function ignores
     * previous calls to setSize().
     *
     * @param vals The values to write as the vector.
     */
    virtual void write(
        std::vector<val_t> const & vals) = 0;


};




}




#endif

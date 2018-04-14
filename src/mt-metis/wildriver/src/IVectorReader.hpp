/**
 * @file IVectorReader.hpp
 * @brief Interface for reading in vectors.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_IVECTORREADER_HPP
#define WILDRIVER_IVECTORREADER_HPP




#include <vector>

#include "base.h"




namespace WildRiver
{


class IVectorReader
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IVectorReader()
    {
    }


    /**
     * @brief Get the size of the vector. May alter the internal state of the
     * reader.
     *
     * @return The size of the vector.
     */
    virtual size_t getSize() = 0;


    /**
     * @brief Read the values of the vector. 
     *
     * @param vals The values in the vector (output).
     */
    virtual void read(
        val_t * vals) = 0;


    /**
     * @brief Read the values of the vector.
     *
     * @param vals The structure to fill with values.
     */
    virtual void read(
        std::vector<val_t> & vals) = 0;


};




}




#endif

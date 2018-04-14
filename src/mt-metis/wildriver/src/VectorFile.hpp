/**
 * @file VectorFile.hpp
 * @brief Base abstract class for vector files.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_VECTORFILE_HPP
#define WILDRIVER_VECTORFILE_HPP




#include <string>
#include <memory>

#include "IVectorFile.hpp"
#include "Vector.hpp"




namespace WildRiver
{


class VectorFile :
  public IVectorFile

{
  private:
    /**
     * @brief Underlying vector structure.
     */
    Vector vector;

  public:
    /**
     * @brief Default constructor which intializes the vector properties. 
     */
    VectorFile();


    /**
     * @brief Virtual destructor.
     */
    virtual ~VectorFile();


    /**
     * @brief Get the size of the vector. May alter the internal state of the
     * reader.
     *
     * @return The size of the vector.
     */
    virtual size_t getSize() override
    {
      return vector.getSize();
    }


    /**
     * @brief Set the size of the vector.  
     *
     * @param size The new size of the vector.
     */
    virtual void setSize(
        const size_t size) override
    {
      vector.setSize(size);
    }


    /**
     * @brief Check to see if the size of the vector has been set.
     *
     * @return True if the size has been set, false otherwise. 
     */
    bool isSizeSet() const noexcept
    {
      return vector.isSizeSet();
    }



};




}




#endif

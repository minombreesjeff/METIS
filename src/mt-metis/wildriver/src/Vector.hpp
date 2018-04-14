/**
 * @file Vector.hpp
 * @brief Class for vectors.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_VECTOR_HPP
#define WILDRIVER_VECTOR_HPP



namespace WildRiver
{


class Vector
{
  private:
    /**
     * @brief Constant value for an unset size.
     */
    static const size_t UNSET_SIZE = (size_t)-1;


    /**
     * @brief The size of the vector.
     */
    size_t size;


  public:
    /**
     * @brief Create a new vector with an unset size.
     */
    Vector() noexcept :
      size(UNSET_SIZE)
    {
      // do nothing
    }


    /**
     * @brief Virtual destructor.
     */
    ~Vector()
    {
      // do nothing
    }


    /**
     * @brief Get the size of the vector. May alter the internal state of the
     * reader.
     *
     * @return The size of the vector.
     */
    size_t getSize() const noexcept
    {
      return size;
    }


    /**
     * @brief Set the size of the vector.  
     *
     * @param size The new size of the vector.
     */
    void setSize(
        const size_t size) noexcept
    {
      this->size = size;
    }


    /**
     * @brief Check to see if the size of the vector has been set.
     *
     * @return True if the size has been set, false otherwise. 
     */
    bool isSizeSet() const noexcept
    {
      return size != UNSET_SIZE; 
    }


};




}




#endif

/**
 * @file VectorTextFile.hpp
 * @brief Abstract class for text based vector files.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_VECTORTEXTFILE_HPP
#define WILDRIVER_VECTORTEXTFILE_HPP




#include "TextFile.hpp"
#include "VectorFile.hpp"




namespace WildRiver
{


class VectorTextFile :
  public TextFile,
  public VectorFile
{
  protected:
    /**
     * @brief get the next value from the file.
     *
     * @return the next value.
     */
    virtual val_t getNextValue() = 0;


    /**
     * @brief Set the next value in the file.
     *
     * @param val The next value.
     */
    virtual void setNextValue(
        val_t val) = 0;


  public:
    /**
     * @brief Open a vector file.
     *
     * @param name
     */
    VectorTextFile(
        std::string const & name);


    /**
     * @brief Virtual destructor.
     */
    virtual ~VectorTextFile();


    /**
     * @brief Read the values of the vector. 
     *
     * @param vals The values in the vector (output).
     */
    virtual void read(
        val_t * vals) override;


    /**
     * @brief Read the values of the vector.
     *
     * @param vals The structure to fill with values.
     */
    virtual void read(
        std::vector<val_t> & vals) override;


    /**
     * @brief Write the vector to the underlying medium. 
     *
     * @param vals The dense array of values in the vector.
     */
    virtual void write(
        val_t const * vals) override;


    /**
     * @brief Write the vector to the underlying medium. This function ignores
     * previous calls to setSize().
     *
     * @param vals The values to write as the vector.
     */
    virtual void write(
        std::vector<val_t> const & vals) override;


    /**
     * @brief Get the filename/path of this file.
     *
     * @return The filename/path.
     */
    virtual std::string const & getFilename() const noexcept
    {
      return TextFile::getFilename();
    }




};




}




#endif

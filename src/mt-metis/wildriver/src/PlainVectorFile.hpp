/**
 * @file PlainVectorFile.hpp
 * @brief Concrete class text vector files.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-08
 */





#ifndef WILDIRVER_PLAINVECTORFILE_HPP
#define WILDIRVER_PLAINVECTORFILE_HPP



#include "VectorTextFile.hpp"




namespace WildRiver
{


class PlainVectorFile :
  public VectorTextFile
{
  private:
    /**
     * @brief The string for buffer read lines into.
     */
    std::string buffer;


  protected:
    /**
     * @brief get the next value from the file.
     *
     * @return the next value.
     */
    val_t getNextValue() override;


    /**
     * @brief Set the next value in the file.
     *
     * @param val The next value.
     */
    void setNextValue(
        val_t val) override;

    /**
     * @brief Determine the given line is a comment.
     *
     * @param line The line.
     *
     * @return True if the line is a comment.
     */
    virtual bool isComment(
        std::string const & line) const noexcept override;


  public:
    /**
     * @brief Open a plain vector file for reading or writing.
     *
     * @param name The name of the file.
     */
    PlainVectorFile(
        std::string const & name);


    /**
     * @brief Close file and associated memory.
     */
    ~PlainVectorFile();


    /**
     * @brief Get the size of the vector. May alter the internal state of the
     * reader.
     *
     * @return The size of the vector.
     */
    virtual size_t getSize() override;




};




}




#endif

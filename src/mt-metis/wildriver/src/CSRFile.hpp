/**
 * @file CSRFile.hpp
 * @brief Class for reading/writing metis files.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_CSRFILE_HPP
#define SRC_CSRFILE_HPP




#include "MatrixTextFile.hpp"




namespace WildRiver {


class CSRFile : 
  public MatrixTextFile
{
  private:
    /**
     * @brief Name of this filetype.
     */
    static std::string const name;


  protected:
    /**
     * @brief Read the header of this matrix file. Populates internal fields
     * with the header information.
     */
    virtual void readHeader() override;


    /**
     * @brief Write the header of this matrix file. The header consists of
     * internal fields set by "setInfo()".
     */
    virtual void writeHeader() override; 


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
     * @brief Check if the given filename matches an extension for this 
     * filetype.
     *
     * @param f The filename.
     *
     * @return True if the extension matches this filetype.
     */
    static bool hasExtension(
        std::string const & f);


    /**
     * @brief Create a new CSRFile for reading and writing.
     *
     * @param fname The filename/path.
     */
    CSRFile(
        std::string const & fname);


    /**
     * @brief Close file and free any memory.
     */
    virtual ~CSRFile();


    /**
     * @brief Reset the current position in the matrix file to the first row.
     */
    virtual void firstRow() override;


    /**
     * @brief Set the next row in the matrix file.
     *
     * @param next The row to set.
     *
     * @return The true if the row was read, false if EOF was encountered. 
     */
    virtual bool getNextRow(
        std::vector<MatrixEntry> & next) override;


    /**
     * @brief Set the next row in the matrix file.
     *
     * @param next The row to set.
     */
    virtual void setNextRow(
        std::vector<MatrixEntry> const & next) override;


    /**
     * @brief Get the name of this matrix file type.
     *
     * @return The matrix file type name.
     */
    virtual std::string const & getName() const noexcept override
    {
      return name;
    } 


};




}




#endif

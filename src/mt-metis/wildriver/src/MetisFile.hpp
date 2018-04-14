/**
 * @file MetisFile.hpp
 * @brief Class for reading/writing metis files.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_METISFILE_HPP
#define SRC_METISFILE_HPP




#include "GraphTextFile.hpp"




namespace WildRiver {


class MetisFile : 
  public GraphTextFile
{
  private:
    /**
     * @brief Name of this filetype.
     */
    static std::string const name;


    /**
     * @brief Get the flags representing the weights associated with this
     * graph.
     *
     * @return 
     */
    int getWeightFlags();


  protected:
    /**
     * @brief Determine if teh given line is a comment.
     *
     * @param line The line.
     *
     * @return True if the line is a comment.
     */
    virtual bool isComment(
        std::string const & line) const noexcept override;


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


  public:
    /**
     * @brief Close file and free any memory.
     */
    virtual ~MetisFile();


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
     * @brief Create a new MetisFile for reading and writing.
     *
     * @param fname The filename/path.
     */
    MetisFile(
        std::string const & fname);


    /**
     * @brief Reset the current position in the matrix file to the first row.
     */
    virtual void firstRow() override;


    /**
     * @brief Reset the current position in the graph to the first vertex.
     */
    virtual void firstVertex() override;


    /**
     * @brief Get the name of this matrix file type.
     *
     * @return The matrix file type name.
     */
    virtual std::string const & getName() const noexcept override
    {
      return name;
    } 


    /**
     * @brief Set the next row in the matrix (adjacecny list in the graph).
     *
     * @param row The next row in the matrix.
     *
     * @return True if another row was found in the file.
     */
    virtual bool getNextRow(
        std::vector<MatrixEntry> & row) override;


    /**
     * @brief Get the information of the next vertex.
     *
     * @param vwgt The vertex weight(s).
     * @param next The adjacency list of the vertex.
     *
     * @return True if another vertex was found in the file.
     */
    virtual bool getNextVertex(
        std::vector<val_t> & vwgts,
        std::vector<MatrixEntry> & list) override;


    /**
     * @brief Set the next row in the matrix file.
     *
     * @param next The row to set.
     */
    virtual void setNextRow(
        std::vector<MatrixEntry> const & next) override;


    /**
     * @brief Set the adjacency list and vertex weight of the next vertex.
     *
     * @param vwgts The vertex weights for this vertex.
     * @param list The adjacecny list.
     */
    virtual void setNextVertex(
        std::vector<val_t> const & vwgts,
        std::vector<MatrixEntry> const & list) override;


};




}




#endif

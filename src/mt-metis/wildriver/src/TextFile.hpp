/**
 * @file TextFile.hpp
 * @brief Abstract class for text files (both graphs and matrices).
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 *
 */





#ifndef WILDRIVER_TEXTFILE_HPP
#define WILDRIVER_TEXTFILE_HPP




#include <string>
#include <vector>
#include <fstream>

#include "Exception.hpp"




namespace WildRiver
{


class TextFile
{
  private:
    /**
     * @brief The current state of the file.
     */
    int state;


    /**
     * @brief The current line number of the filestream.
     */
    size_t current_line;


    /**
     * @brief The filename/path of this file.
     */
    std::string fname;


    /**
     * @brief The I/O stream for this file.
     */
    std::fstream stream;


  protected:
    /**
     * @brief Match the given filename with the given extensions.
     *
     * @param f The filename.
     * @param extensions The extensions.
     *
     * @return True if the filename matches one of the extensions.
     */
    static bool matchExtension(
        std::string const & f,
        std::vector<std::string> const & extensions);


    /**
     * @brief Open the underlying file for writing.
     */
    void openWrite();


    /**
     * @brief Open the underlying file for writing.
     */
    void openRead();


    /**
     * @brief Return to the start of the i/o stream.
     */
    void resetStream();


    /**
     * @brief Retrieve the next line in the file.
     *
     * @param line The string to populate with the contents of the next line.
     *
     * @return True if the line was successfully read, false otherwise.
     */
    bool nextLine(std::string & line);


    /**
     * @brief Retrieve the next non-comment line from the file.
     *
     * @param line The string to populate with the contents of the next
     * non-comment line.
     *
     * @return True if the line was successfully read, false otherwise.
     */
    bool nextNoncommentLine(std::string & line);


    /**
     * @brief Get the current line in the file.
     *
     * @return The current line in the file.
     */
    inline size_t getCurrentLine() const
    {
      return current_line;
    }


    /**
     * @brief Retrive the stream associated with this object.
     *
     * @return The iostream.
     */
    inline std::iostream & getStream()
    {
      return stream;
    }


    /**
     * @brief Check if this file is already open for writing.
     *
     * @return True if the file is open for writing.
     */
    bool isOpenWrite() const;


    /**
     * @brief Check if this file is already open for reading.
     *
     * @return True if the file is open for reading.
     */
    bool isOpenRead() const;


    /**
     * @brief Determine if the given line is a comment or not.
     *
     * @param line The line to check.
     *
     * @return True if the line is a comment.
     */
    virtual bool isComment(
        std::string const & line) const noexcept = 0;


  public:
    /**
     * @brief Create a new text file from the given filename/path.
     *
     * @param fname The filename/path.
     */
    TextFile(
        std::string const & fname);


    /**
     * @brief Close the underlying filestream if opened.
     */
    ~TextFile();


    /**
     * @brief Get the filename/path of this file.
     *
     * @return The filename/path.
     */
    virtual std::string const & getFilename() const noexcept
    {
      return fname;
    }


};




}




#endif

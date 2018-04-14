/**
 * @file CSRFile.cpp
 * @brief Class functions for reading and writing metis graph files.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_METISFILE_CPP
#define SRC_METISFILE_CPP




#include "CSRFile.hpp"




namespace WildRiver 
{


/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


std::string const CSRFile::name = "CSR";




/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


bool CSRFile::hasExtension(
    std::string const & f)
{
  std::vector<std::string> extensions;

  extensions.push_back(".csr");

  return TextFile::matchExtension(f,extensions);
}



bool CSRFile::isComment(
    std::string const & line) const noexcept
{
  // awful solution since I can't declare this statically in c++ -- at
  // somepoint generate all 256 entries
  bool comment_chars[256] = {false};
  comment_chars['#'] = true;
  comment_chars['%'] = true;
  comment_chars['"'] = true;
  comment_chars['/'] = true;

  return line.size() > 0 && comment_chars[line[0]];
}



void CSRFile::readHeader()
{
  int flags;
  size_t offset;
  std::string line;

  // open our file for reading
  openRead();

  // go to the start of the file
  firstRow();

  // read the whole file
  std::vector<MatrixEntry> row;

  dim_t nrows = 0;
  dim_t ncols = 0;
  ind_t nnz = 0;
  while (getNextRow(row)) {
    for (MatrixEntry const & e : row) {
      if (e.ind > ncols) {
        ncols = e.ind;
      }
      ++nnz;
    }
    ++nrows;
  }
  ++ncols;

  setNumRows(nrows);
  setNumCols(ncols);
  setNNZ(nnz);

  // go back to the start
  firstRow();
}


void CSRFile::writeHeader()
{
  // open for writing
  openWrite();
}




/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


CSRFile::CSRFile(
    std::string const & fname) : 
  MatrixTextFile(fname)
{
}


CSRFile::~CSRFile()
{
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void CSRFile::firstRow()
{
  // go to the start of the file
  resetStream();
  setCurrentRow(0);
}


bool CSRFile::getNextRow(
    std::vector<MatrixEntry> & next)
{
  MatrixEntry e;
  std::string line;
  size_t offset;

  if (!nextNoncommentLine(line)) {
    return false;
  }

  next.clear();

  try {
    dim_t const ncols = getNumCols();

    // Loop through row until we streamed to the end
    while (!line.empty()) {
      e.ind = std::stoull(line,&offset,10) - 1;
      if (ncols != NULL_DIM && e.ind >= ncols) {
        throw BadFileException(std::string("Entry with column of ") + \
            std::to_string(e.ind) + std::string("/") + \
            std::to_string(ncols));
      }

      line = line.substr(offset);
      e.val = std::stod(line,&offset);
      line = line.substr(offset);

      next.push_back(e);

      if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
        // only whitespace left in this string
        break;
      }
    }

    incRow();

  } catch (std::invalid_argument & e) {
    // maybe make a custom exception to format this better
    throw BadFileException( \
        std::string("Failed to read column index at line:" + \
          std::to_string(getCurrentLine()) + std::string(":'") + \
          line + std::string("' : ") + e.what()));
  }

  return true;
}


void CSRFile::setNextRow(
    std::vector<MatrixEntry> const & next)
{
  size_t const size = next.size();
  if (size > 0) {
    for (size_t i=0; i<size-1; ++i) {
      MatrixEntry const & e = next[i];
      getStream() << (e.ind+1) << " " << e.val << " "; 
    }
    MatrixEntry const & e = next.back();
    getStream() << (e.ind+1) << " " << e.val;
  }
  getStream() << std::endl;

  incRow();
}




}




#endif

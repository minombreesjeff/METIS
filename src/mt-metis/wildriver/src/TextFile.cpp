/**
 * @file TextFile.cpp
 * @brief Abstract class for text files (both graphs and matrices).
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */





#ifndef WILDRIVER_TEXTFILE_CPP
#define WILDRIVER_TEXTFILE_CPP




#include "TextFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


enum {
  FILE_STATE_UNOPENED,
  FILE_STATE_READ,
  FILE_STATE_WRITE
};




/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


bool TextFile::matchExtension(
    std::string const & f,
    std::vector<std::string> const & extensions)
{
  for (std::string const & ext : extensions) {
    size_t len = ext.size();

    if (f.compare(f.size()-len,len,ext) == 0) {
      return true;
    }
  }

  return false;
}


void TextFile::resetStream()
{
  stream.clear();
  stream.seekg(0,std::ifstream::beg);
  current_line = 0;
}


void TextFile::openWrite()
{
  if (state != FILE_STATE_UNOPENED) {
    throw BadFileStateException(std::string("Attempting to re-open file '") + \
        fname + std::string("' for writing."));
  }

  stream.open(getFilename(),std::fstream::out | std::fstream::trunc);

  if (!stream.good()) {
    throw BadFileException(std::string("Failed to open file '") + \
        fname + std::string("'"));
  }

  state = FILE_STATE_WRITE;
}


void TextFile::openRead()
{
  if (state != FILE_STATE_UNOPENED) {
    throw BadFileStateException(std::string("Attempting to re-open file '") + \
        fname + std::string("' for reading."));
  }

  stream.open(getFilename(),std::fstream::in);

  if (!stream.good()) {
    throw BadFileException(std::string("Failed to open file '") + \
        fname + std::string("'"));
  }

  state = FILE_STATE_READ;
}


bool TextFile::nextLine(
    std::string & line)
{
  std::getline(stream,line);

  if (stream.good()) {
    ++current_line;
    return true;
  } else {
    return false;
  }
}


bool TextFile::nextNoncommentLine(
    std::string & line)
{
  do {
    if (!nextLine(line)) { 
      return false;
    }
  } while (isComment(line));

  return true;
}


bool TextFile::isOpenWrite() const
{
  return state == FILE_STATE_WRITE;
}


bool TextFile::isOpenRead() const
{
  return state == FILE_STATE_READ;
}




/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


TextFile::TextFile(
    std::string const & name) :
  fname(name)
{
  state = FILE_STATE_UNOPENED;
  current_line = 0;
}


TextFile::~TextFile()
{
  stream.close();
}




}




#endif

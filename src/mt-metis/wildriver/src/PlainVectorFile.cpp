/**
 * @file PlainVectorFile.cpp
 * @brief Implemenation of PlainVectorFile class.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-08
 */





#include "PlainVectorFile.hpp"




namespace WildRiver
{


/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


val_t PlainVectorFile::getNextValue()
{
  if (!nextNoncommentLine(buffer)) {
    throw EOFException("Hit end of file before getting next value");
  }

  if (typeid(val_t) == typeid(double) || typeid(val_t) == typeid(float)) {
    return static_cast<val_t>(std::stod(buffer));
  } else {
    return static_cast<val_t>(std::stoll(buffer));
  }
}


void PlainVectorFile::setNextValue(
    val_t val)
{
  getStream() << val << std::endl; 
}

bool PlainVectorFile::isComment(
    std::string const & line) const noexcept
{
  if (line.size() > 0) {
    switch (line[0]) {
      case '#':
      case '%':
      case '/':
        return true;
      default:
        return false;
    }
  } else {
    return false;
  }
}


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


PlainVectorFile::PlainVectorFile(
    std::string const & name) :
  VectorTextFile(name)
{
  // do nothing
}


PlainVectorFile::~PlainVectorFile()
{
  // do nothing
}


/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


size_t PlainVectorFile::getSize()
{
  if (!isSizeSet()) {
    size_t nlines = 0;
    std::string buffer;

    if (!isOpenRead()) {
      openRead();
    }

    // cout non-comment lines
    while (nextNoncommentLine(buffer)) {
      ++nlines;
    }

    VectorFile::setSize(nlines);
  }

  return VectorFile::getSize();
}




}

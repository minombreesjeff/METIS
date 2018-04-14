/**
 * @file VectorTextFile.cpp
 * @brief Implementation of abstract class for text based vector files.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-07
 */




#include "VectorTextFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


VectorTextFile::VectorTextFile(
    std::string const & name) :
  TextFile(name)
{
  // do nothing
}


VectorTextFile::~VectorTextFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void VectorTextFile::read(
    val_t * vals)
{
  if (!isOpenRead()) {
    openRead();
  }

  const size_t n = getSize();

  resetStream();
  
  for (size_t i=0;i<n;++i) {
    vals[i] = getNextValue();
  }
}


void VectorTextFile::read(
    std::vector<val_t> & vals)
{
  size_t i;

  const size_t n = getSize();
  vals.resize(n);

  read(vals.data());
}


void VectorTextFile::write(
    val_t const * vals)
{
  if (!isSizeSet()) {
    throw UnsetInfoException("Size of vector is not set before call to " \
        "write()");
  }

  if (!isOpenWrite()) {
    openWrite();
  }

  const size_t n = getSize();
  for (size_t i=0;i<n;++i) {
    setNextValue(vals[i]);
  }
}


void VectorTextFile::write(
    std::vector<val_t> const & vals)
{
  setSize(vals.size());

  write(vals.data());
}




}

/**
 * @file MatrixFactory.cpp
 * @brief Implementation of the MatrixFactory class.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-05
 */




#include "MatrixFactory.hpp"
#include "CSRFile.hpp"
#include "MetisFile.hpp"




namespace WildRiver
{


/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::shared_ptr<IMatrixFile> MatrixFactory::OpenFile(
    std::string const & fname)
{
  std::shared_ptr<IMatrixFile> file;

  // determine what type of reader to instantiate based on extension
  if (MetisFile::hasExtension(fname)) {
    file = std::shared_ptr<IMatrixFile>(new MetisFile(fname));
  } else if (CSRFile::hasExtension(fname)) {
    file = std::shared_ptr<IMatrixFile>(new CSRFile(fname));
  } else {
    throw UnknownExtensionException(std::string("Unknown filetype: ") + fname);
  }
 
  return file;
}




}

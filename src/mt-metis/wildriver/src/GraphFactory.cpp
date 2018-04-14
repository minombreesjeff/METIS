/**
 * @file GraphFactory.cpp
 * @brief Implementation of GraphFactory for instantiating graph files.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-05
 */



#include "GraphFactory.hpp"
#include "MetisFile.hpp"
#include "CSRFile.hpp"
#include "MatrixGraphFile.hpp"




namespace WildRiver
{


/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::shared_ptr<IGraphFile> GraphFactory::OpenFile(
    std::string const & fname)
{
  std::shared_ptr<IGraphFile> file;
  // determine what type of reader to instantiate based on extension
  if (MetisFile::hasExtension(fname)) {
    file = std::shared_ptr<IGraphFile>(new MetisFile(fname));
  } else if (CSRFile::hasExtension(fname)) {
    // need to wrap it with an adapter
    file = \
        std::shared_ptr<IGraphFile>(new MatrixGraphFile(new CSRFile(fname)));
  } else {
    throw UnknownExtensionException(std::string("Unknown filetype: ") + fname);
  }
 
  return file;
}




}

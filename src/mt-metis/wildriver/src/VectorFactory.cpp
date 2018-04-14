/**
 * @file VectorFactory.cpp
 * @brief Implemnetation of class for 
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-07
 */




#include "Exception.hpp"
#include "VectorFactory.hpp"




namespace WildRiver
{



/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::shared_ptr<IVectorFile> OpenFile(
    std::string const & fname)
{
  std::shared_ptr<IVectorFile> file;

  // determine what type of reader to instantiate based on extension
  throw UnknownExtensionException(std::string("Unknown filetype: ") + fname);
 
  return file;

}




}

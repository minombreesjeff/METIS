/**
 * @file IGraphFile.hpp
 * @brief Interface for graph files.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 *
 */




#ifndef WILDRIVER_IGRAPHFILE_HPP
#define WILDRIVER_IGRAPHFILE_HPP




#include <string>

#include "IGraphReader.hpp"
#include "IGraphWriter.hpp"




namespace WildRiver
{


class IGraphFile :
  public IGraphReader,
  public IGraphWriter
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IGraphFile() {}

    /**
     * @brief Set the cursor to the first vertex in the file.
     */
    virtual void firstVertex() = 0;


    /**
     * @brief Get the filename/path of the graph file.
     *
     * @return The filename/path. 
     */
    virtual std::string const & getFilename() const = 0;


};


}




#endif

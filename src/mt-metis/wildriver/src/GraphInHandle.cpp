/**
 * @file GraphHandle.cpp
 * @brief Class for reading all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_MATRIXINHANDLE_CPP
#define SRC_MATRIXINHANDLE_CPP




#include "GraphInHandle.hpp"
#include "GraphFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


GraphInHandle::GraphInHandle(
    std::string const & fname)
{
  reader = GraphFactory::OpenFile(fname);
}


GraphInHandle::~GraphInHandle()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphInHandle::getInfo(
    dim_t & nvtxs,
    ind_t & nedges,
    int & nvwgt,
    bool & ewgts)
{
  reader->getInfo(nvtxs,nedges,nvwgt,ewgts);
}


void GraphInHandle::readGraph(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt)
{
  reader->read(xadj,adjncy,vwgt,adjwgt);
}


dim_t GraphInHandle::getNextVertex(
        std::vector<val_t> & vwgt,
        std::vector<MatrixEntry> & list)
{
  return reader->getNextVertex(vwgt,list);
}


}




#endif


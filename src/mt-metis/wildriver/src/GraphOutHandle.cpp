/**
 * @file GraphOutHandle.cpp
 * @brief Classfor writing all graph types.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef WILDRIVER_GRAPHOUTHANEL_CPP
#define WILDRIVER_GRAPHOUTHANEL_CPP




#include "GraphOutHandle.hpp"
#include "GraphFactory.hpp"




namespace WildRiver
{



/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


GraphOutHandle::GraphOutHandle(
    std::string const & fname)
{
  writer = GraphFactory::OpenFile(fname);
}


GraphOutHandle::~GraphOutHandle()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphOutHandle::writeGraph(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  writer->write(xadj,adjncy,vwgt,adjwgt);
}


void GraphOutHandle::setInfo(
    dim_t const nvtxs,
    ind_t const nedges,
    int const nvwgt,
    bool const ewgts)
{
  writer->setInfo(nvtxs,nedges,nvwgt,ewgts);
}


void GraphOutHandle::setNextVertex(
    std::vector<val_t> const & vwgts,
    std::vector<MatrixEntry> const & list)
{
  writer->setNextVertex(vwgts,list);
}




}




#endif

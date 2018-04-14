/**
 * @file GraphTextFile.cpp
 * @brief Abstract class for reading and writing graphs.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 *
 */




#ifndef WILDRIVER_GRAPHFILE_CPP
#define WILDRIVER_GRAPHFILE_CPP




#include "GraphTextFile.hpp"
#include "MetisFile.hpp"
#include "CSRFile.hpp"
#include "MatrixGraphFile.hpp"



namespace WildRiver
{


/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


dim_t GraphTextFile::incVertex()
{
  if (current_vertex >= getNumVertices()) {
    throw BadFileStateException(std::string("Attempt to increase current " \
        "vertex beyond the number of vertices: ") + \
        std::to_string(getNumRows()));
  }

  ++current_vertex;

  return current_vertex;
}




/******************************************************************************
* CONSTRUCTOR / DESTRUCTOR ****************************************************
******************************************************************************/


GraphTextFile::GraphTextFile(
    std::string const & fname) :
  TextFile(fname),
  current_vertex(0)
{
  // do nothing
}


GraphTextFile::~GraphTextFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphTextFile::read(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt)
{
  dim_t const nvtxs = getNumVertices();
  int const ncon = getNumVertexWeights();
  bool const ewgts = hasEdgeWeights();

  xadj[0] = 0;
  for (dim_t i=0; i<nvtxs; ++i) {
    // retrieve row
    std::vector<val_t> vwgts;
    std::vector<MatrixEntry> list;
    if (!getNextVertex(vwgts,list)) {
      throw BadFileException(std::string("Premature end of file: ") + \
          std::to_string(i) + std::string("/") + std::to_string(nvtxs) + \
          std::string(" vertices found."));
    }

    // handle vertex weights
    if (ncon > 0) {
      for (dim_t k=0; k<ncon; ++k) {
        vwgt[(i*ncon)+k] = vwgts[k];
      }
    } else if (vwgt) {
      // set unit vertex weights
      vwgt[i] = 1;
    }

    // handle edges 
    ind_t const adjstart = xadj[i];
    dim_t j;
    for (j=0; j<list.size(); ++j) {
      adjncy[adjstart+j] = list[j].ind; 
      if (ewgts) {
        adjwgt[adjstart+j] = list[j].val;
      } else if (adjwgt) {
        // set unit edge weights
        adjwgt[adjstart+j] = 1;
      }
    }
    xadj[i+1] = adjstart+j;
  }
}


void GraphTextFile::write(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  dim_t const nvtxs = getNumVertices();
  int const ncon = getNumVertexWeights();
  bool const ewgts = hasEdgeWeights();

  for (dim_t i=0; i<nvtxs; ++i) {
    std::vector<val_t> vwgts;
    std::vector<MatrixEntry> list;

    // build vertex weight vector
    for (dim_t k=0; k<ncon; ++k) {
      vwgts.push_back(vwgt[(i*ncon)+k]);
    }

    // build edges
    for (ind_t j=xadj[i]; j<xadj[i+1]; ++j) {
      MatrixEntry e;
      e.ind = adjncy[j];
      if (ewgts) {
        if (adjwgt) {
          e.val = adjwgt[j];
        } else {
          e.val = 1;
        }
      }
      list.push_back(e);
    }

    // set the vertex
    setNextVertex(vwgts,list);
  }
}




}




#endif

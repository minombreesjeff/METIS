/**
 * @file MetisFile.cpp
 * @brief Class functions for reading and writing metis graph files.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_METISFILE_CPP
#define SRC_METISFILE_CPP




#include "MetisFile.hpp"




namespace WildRiver 
{


/******************************************************************************
* INTERNAL TYPES **************************************************************
******************************************************************************/


enum {
  HAS_NOWEIGHTS = 0,
  HAS_EDGEWEIGHTS = 1,
  HAS_VERTEXWEIGHTS = 10
};




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


std::string const MetisFile::name = "Metis";




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


int MetisFile::getWeightFlags()
{
  int flags = HAS_NOWEIGHTS;

  if (hasEdgeWeights()) {
    flags |= HAS_EDGEWEIGHTS;
  }

  if (getNumVertexWeights() > 0) {
    flags |= HAS_VERTEXWEIGHTS;
  }

  return flags;
}




/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


bool MetisFile::hasExtension(
    std::string const & f)
{
  std::vector<std::string> extensions;

  extensions.push_back(".graph");
  extensions.push_back(".metis");
  extensions.push_back(".chaco");

  return TextFile::matchExtension(f,extensions);
}


bool MetisFile::isComment(
    std::string const & line) const noexcept
{
  // awful solution since I can't declare this statically in c++ -- at
  // somepoint generate all 256 entries using template programming
  bool comment_chars[256] = {false};
  comment_chars['#'] = true;
  comment_chars['%'] = true;
  comment_chars['"'] = true;
  comment_chars['/'] = true;

  return line.size() > 0 && comment_chars[line[0]];
}


void MetisFile::readHeader()
{
  if (!isOpenRead()) {
    // open our file for reading if not already
    openRead();
  }

  // get the first line
  std::string line;
  nextNoncommentLine(line);

  // parse out my header
  size_t offset;
  dim_t nvtxs = std::stoull(line,&offset,10);
  line = line.substr(offset);

  setNumVertices(nvtxs);

  ind_t nedges = std::stoull(line,&offset,10)*2;
  line = line.substr(offset);

  setNumEdges(nedges);

  // handle weights
  if (line.find_first_not_of(" \t") != std::string::npos) {
    int flags = std::stoul(line,&offset,10);  
    line = line.substr(offset);

    if (flags & HAS_EDGEWEIGHTS) {
      hasEdgeWeights(true);
    }

    if (flags & HAS_VERTEXWEIGHTS) {
      setNumVertexWeights(std::stoul(line,&offset,10));
    }
  }
}


void MetisFile::writeHeader()
{
  dim_t nvtxs = getNumVertices();
  ind_t nedges = getNumEdges();

  // make sure we have a valid number of edges/nnz
  if (nedges % 2 != 0) {
    throw BadParameterException("Metis files are required to be symmetric: " \
        "odd number of non-zeroes specified.");
  }

  // set header to write mode if not already
  if (!isOpenWrite()) {
    openWrite();
  }

  // write the header -- edges in metis files are undirected (symmetric).
  getStream() << nvtxs << " " << (nedges/2);

  int weightflag = getWeightFlags();  

  if (weightflag != HAS_NOWEIGHTS) {
    // write weight flags
    getStream() << " " << weightflag;
    if (weightflag & HAS_VERTEXWEIGHTS) {
      // write number of vertex weights
      getStream() << " " << getNumVertexWeights();
    }
  } 

  getStream() << std::endl;
}




/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MetisFile::MetisFile(
    std::string const & fname) : 
  GraphTextFile(fname)
{
}


MetisFile::~MetisFile()
{
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MetisFile::firstRow()
{
  resetStream();

  // go past all comments -- read and discard header
  std::string line;
  nextNoncommentLine(line);
}


void MetisFile::firstVertex()
{
  firstRow();
}


bool MetisFile::getNextRow(
    std::vector<MatrixEntry> & row)
{
  std::string line;

  // get my line
  if (!nextNoncommentLine(line)) {
    return false;
  }

  row.clear();

  dim_t const nvtxs = getNumVertices();
  bool ewgts = hasEdgeWeights();

  // read in edges
  while (!line.empty() && \
      line.find_first_not_of(" \t\r\n") != std::string::npos) {
    size_t offset;
    MatrixEntry e;
    e.ind = std::stoull(line,&offset,10) - 1;
    line = line.substr(offset);

    // make sure this is a valid edge
    if (e.ind >= nvtxs) {
      throw BadFileException(std::string("Edge with destination of ") + \
          std::to_string(e.ind) + std::string("/") + \
          std::to_string(nvtxs));
    }

    if (ewgts) {
      e.val = std::stod(line,&offset);
      line = line.substr(offset);
    }
    
    row.push_back(e);
  }
}


bool MetisFile::getNextVertex(
    std::vector<val_t> & vwgts,
    std::vector<MatrixEntry> & list)
{
  int const ncon = getNumVertexWeights();

  std::string line;

  // get my line
  if (!nextNoncommentLine(line)) {
    return false;
  }

  vwgts.clear();
  list.clear();

  // read in vertex weights
  for (dim_t k=0; k<ncon; ++k) {
    size_t offset;
    vwgts.push_back(std::stod(line,&offset));
    line = line.substr(offset);
  }

  dim_t const nvtxs = getNumVertices();
  bool ewgts = hasEdgeWeights();

  // read in edges
  while (!line.empty() && \
      line.find_first_not_of(" \t\r\n") != std::string::npos) {
    size_t offset;
    MatrixEntry e;
    e.ind = std::stoull(line,&offset,10) - 1;
    line = line.substr(offset);

    // make sure this is a valid edge
    if (e.ind >= nvtxs) {
      throw BadFileException(std::string("Edge with destination of ") + \
          std::to_string(e.ind) + std::string("/") + \
          std::to_string(nvtxs));
    }

    if (ewgts) {
      e.val = std::stod(line,&offset);
      line = line.substr(offset);
    }
    
    list.push_back(e);
  }

  // indicate that we successfully found a vertex
  return true;
}


void MetisFile::setNextRow(
    std::vector<MatrixEntry> const & row)
{
  dim_t const nadj = row.size();
  bool const ewgts = hasEdgeWeights();

  for (dim_t j=0; j<nadj; ++j) {
    MatrixEntry e = row[j];
    getStream() << (e.ind+1);
    if (ewgts) {
      getStream() << " " << e.val;
    }
    if (j < nadj-1) {
      // add space
      getStream() << " ";
    }
  }
  getStream() << std::endl;

  incVertex();
}


void MetisFile::setNextVertex(
    std::vector<val_t> const & vwgts,
    std::vector<MatrixEntry> const & list)
{
  int const ncon = getNumVertexWeights();
  dim_t const nadj = list.size();

  // set vertex weights
  for (dim_t k=0; k<ncon; ++k) {
    getStream() << vwgts[k];
    if (k < ncon-1 || nadj > 0) {
      getStream() << " ";
    }
  }

  // set adjacency list and edge weights
  setNextRow(list);
}




}




#endif

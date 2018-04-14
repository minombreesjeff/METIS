/**
 * @file WildRiver.cpp
 * @brief Main function.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef SRC_WILDRIVER_CPP
#define SRC_WILDRIVER_CPP




#include <cstdint>
#include <iostream>
#include <memory>

#include "MatrixInHandle.hpp"
#include "MatrixOutHandle.hpp"
#include "GraphInHandle.hpp"
#include "GraphOutHandle.hpp"
#include "Exception.hpp"




using namespace WildRiver;



/******************************************************************************
* HELPER FUNCTIONS ************************************************************
******************************************************************************/

namespace
{

/**
 * @brief Deleter function for using malloc with unique_ptr.
 */
struct c_delete
{
  void operator()(void * ptr)
  {
    free(ptr);
  }
};

}


/******************************************************************************
* API FUNCTIONS ***************************************************************
******************************************************************************/


extern "C" int wildriver_read_matrix(
    char const * const fname,
    dim_t * const r_nrows,
    dim_t * const r_ncols,
    ind_t * const r_nnz,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    val_t ** const r_rowval)
{
  try {
    MatrixInHandle handle(fname);

    dim_t nrows, ncols;
    ind_t nnz;

    // read the header
    handle.getInfo(nrows,ncols,nnz);

    // allocate matrix
    size_t nbytes = sizeof(ind_t)*(nrows+1);
    std::unique_ptr<ind_t,c_delete> rowptr((ind_t*)malloc(nbytes));
    if (!rowptr.get()) {
      throw OutOfMemoryException(nbytes);
    }

    nbytes = sizeof(dim_t)*nnz;
    std::unique_ptr<dim_t,c_delete> rowind((dim_t*)malloc(nbytes));
    if (!rowind.get()) {
      throw OutOfMemoryException(nbytes);
    }

    std::unique_ptr<val_t,c_delete> rowval;
    if (r_rowval) {
      // we need to use rowval
      nbytes = sizeof(val_t)*nnz;
      rowval.reset((val_t*)malloc(nbytes));
      if (!rowind.get()) {
        throw OutOfMemoryException(nbytes);
      }
    } else {
      // don't allocate rowval
      nbytes = 0;
    }

    handle.readSparse(rowptr.get(),rowind.get(),rowval.get());

    // we've completely succeed -- assign pointers
    *r_nrows = nrows;
    *r_ncols = ncols;
    *r_nnz = nnz;

    *r_rowptr = rowptr.release();
    *r_rowind = rowind.release();
    if (r_rowval) {
      *r_rowval = rowval.release();
    }
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to read matrix due to: " << e.what() \
        << std::endl;
    return 0;
  }

  return 1;
}


extern "C" int wildriver_write_matrix(
    char const * const fname,
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  try{
    MatrixOutHandle handle(fname);

    handle.setInfo(nrows,ncols,nnz);

    handle.writeSparse(rowptr,rowind,rowval);
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to write matrix due to: " << e.what() \
        << std::endl;
    return 0;
  }

  return 1;
}


extern "C" int wildriver_read_graph(
    char const * const fname,
    dim_t * const r_nvtxs,
    ind_t * const r_nedges,
    int * const r_nvwgts,
    int * const r_ewgts,
    ind_t ** const r_xadj,
    dim_t ** const r_adjncy,
    val_t ** const r_vwgt,
    val_t ** const r_adjwgt)
{
  try {
    GraphInHandle handle(fname);

    dim_t nvtxs;
    ind_t nedges;
    int nvwgts;
    bool ewgts;

    // read the header
    handle.getInfo(nvtxs,nedges,nvwgts,ewgts);

    // allocate matrix
    size_t nbytes = sizeof(ind_t)*(nvtxs+1);
    std::unique_ptr<ind_t,c_delete> xadj((ind_t*)malloc(nbytes));
    if (!xadj.get()) {
      throw OutOfMemoryException(nbytes);
    }

    nbytes = sizeof(dim_t)*nedges;
    std::unique_ptr<dim_t,c_delete> adjncy((dim_t*)malloc(nbytes));
    if (!adjncy.get()) {
      throw OutOfMemoryException(nbytes);
    }

    std::unique_ptr<val_t> vwgt;
    if (r_vwgt && nvwgts > 0) {
      // we need to use rowval
      nbytes = sizeof(val_t)*nvtxs*nvwgts;
      vwgt.reset((val_t*)malloc(nbytes));
      if (!adjncy.get()) {
        throw OutOfMemoryException(nbytes);
      }
    }

    std::unique_ptr<val_t,c_delete> adjwgt;
    if (r_adjwgt) {
      // we need to use rowval
      nbytes = sizeof(val_t)*nedges;
      adjwgt.reset((val_t*)malloc(nbytes));
      if (!adjncy.get()) {
        throw OutOfMemoryException(nbytes);
      }
    }

    handle.readGraph(xadj.get(),adjncy.get(),vwgt.get(),adjwgt.get());

    // we've completed exception possible tasks -- assign pointers
    *r_xadj = xadj.release();
    *r_adjncy = adjncy.release();
    if (r_vwgt) {
      *r_vwgt = vwgt.release();
    }
    if (r_adjwgt) {
      *r_adjwgt = adjwgt.release();
    }

    *r_nvtxs = nvtxs;

    if (r_nedges) {
      *r_nedges = nedges;
    }

    if (r_ewgts) {
      // convert to c style int
      *r_ewgts = (int)ewgts;
    }
    if (r_nvwgts) {
      *r_nvwgts = nvwgts;
    }
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to read graph due to: " << e.what() \
        << std::endl;
    return 0;
  }

  return 1;
}


extern "C" int wildriver_write_graph(
    char const * const fname,
    dim_t const nvtxs,
    ind_t const nedges,
    int nvwgts,
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  try{
    GraphOutHandle handle(fname);

    bool const ewgts = adjwgt != NULL;

    handle.setInfo(nvtxs,nedges,nvwgts,ewgts);

    handle.writeGraph(xadj,adjncy,vwgt,adjwgt);
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to write graph due to: " << e.what() \
        << std::endl;
    return 0;
  }

  return 1;
}




#endif

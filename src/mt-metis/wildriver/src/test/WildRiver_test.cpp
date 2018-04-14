/**
 * @file WildRiver_test.cpp
 * @brief Test for reading matrices
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#include <iostream>
#include <fstream>

#include "wildriver.h"
#include "DomTest.hpp"




namespace DomTest
{


static void writeMatrix(
    std::string const & testFile)
{
  wildriver_dim_t nrows = 6, ncols = 6;
  wildriver_ind_t nnz = 14;
  wildriver_ind_t rowptr[] = {0,2,4,7,10,12,14};
  wildriver_dim_t rowind[] = {1,2,0,2,0,1,3,2,4,5,3,5,3,4};
  wildriver_val_t rowval[] = {1,2,3,4,5,6,7,8,9,1,2,3,4,5};

  int rv = wildriver_write_matrix(testFile.data(),nrows,ncols,nnz,rowptr, \
      rowind,rowval);

  testEquals(rv,1);
}


static void writeGraph(
    std::string const & testFile)
{
  wildriver_dim_t nvtxs = 6;
  wildriver_ind_t nedges = 14;
  wildriver_ind_t xadj[] = {0,2,4,7,10,12,14};
  wildriver_dim_t adjncy[] = {1,2,0,2,0,1,3,2,4,5,3,5,3,4};
  wildriver_val_t adjwgt[] = {1,2,3,4,5,6,7,8,9,1,2,3,4,5};

  int rv = wildriver_write_graph(testFile.data(),nvtxs,nedges,0,xadj,adjncy, \
      NULL,adjwgt);

  testEquals(rv,1);
}


static void readMatrix(
    std::string const & testFile)
{
  wildriver_dim_t nrows, ncols;
  wildriver_ind_t nnz;
  wildriver_ind_t * rowptr;
  wildriver_dim_t * rowind;
  wildriver_val_t * rowval;

  int rv = wildriver_read_matrix(testFile.data(),&nrows,&ncols,&nnz,&rowptr, \
      &rowind,&rowval);

  testEquals(rv,1);

  testEquals(nrows,6);
  testEquals(ncols,6);
  testEquals(nnz,14);

  // test rowptr
  testEquals(rowptr[0],0);
  testEquals(rowptr[1],2);
  testEquals(rowptr[2],4);
  testEquals(rowptr[3],7);
  testEquals(rowptr[4],10);
  testEquals(rowptr[5],12);
  testEquals(rowptr[6],14);

  // test rowind
  testEquals(rowind[0],1);
  testEquals(rowind[1],2);
  testEquals(rowval[0],1);
  testEquals(rowval[1],2);

  testEquals(rowind[2],0);
  testEquals(rowind[3],2);
  testEquals(rowval[2],3);
  testEquals(rowval[3],4);

  testEquals(rowind[4],0);
  testEquals(rowind[5],1);
  testEquals(rowind[6],3);
  testEquals(rowval[4],5);
  testEquals(rowval[5],6);
  testEquals(rowval[6],7);

  testEquals(rowind[7],2);
  testEquals(rowind[8],4);
  testEquals(rowind[9],5);
  testEquals(rowval[7],8);
  testEquals(rowval[8],9);
  testEquals(rowval[9],1);

  testEquals(rowind[10],3);
  testEquals(rowind[11],5);
  testEquals(rowval[10],2);
  testEquals(rowval[11],3);

  testEquals(rowind[12],3);
  testEquals(rowind[13],4);
  testEquals(rowval[12],4);
  testEquals(rowval[13],5);
}


static void readGraph(
    std::string const & testFile)
{
  wildriver_dim_t nvtxs;
  wildriver_ind_t nedges;
  int nvwgts;
  int ewgts;
  wildriver_ind_t * xadj;
  wildriver_dim_t * adjncy;
  wildriver_val_t * adjwgt;

  int rv = wildriver_read_graph(testFile.data(),&nvtxs,&nedges,&nvwgts, \
      &ewgts,&xadj,&adjncy,NULL,&adjwgt);

  testEquals(rv,1);

  testEquals(nvtxs,6);
  testEquals(nedges,14);

  // test rowptr
  testEquals(xadj[0],0);
  testEquals(xadj[1],2);
  testEquals(xadj[2],4);
  testEquals(xadj[3],7);
  testEquals(xadj[4],10);
  testEquals(xadj[5],12);
  testEquals(xadj[6],14);

  // test adjncy
  testEquals(adjncy[0],1);
  testEquals(adjncy[1],2);
  testEquals(adjwgt[0],1);
  testEquals(adjwgt[1],2);

  testEquals(adjncy[2],0);
  testEquals(adjncy[3],2);
  testEquals(adjwgt[2],3);
  testEquals(adjwgt[3],4);

  testEquals(adjncy[4],0);
  testEquals(adjncy[5],1);
  testEquals(adjncy[6],3);
  testEquals(adjwgt[4],5);
  testEquals(adjwgt[5],6);
  testEquals(adjwgt[6],7);

  testEquals(adjncy[7],2);
  testEquals(adjncy[8],4);
  testEquals(adjncy[9],5);
  testEquals(adjwgt[7],8);
  testEquals(adjwgt[8],9);
  testEquals(adjwgt[9],1);

  testEquals(adjncy[10],3);
  testEquals(adjncy[11],5);
  testEquals(adjwgt[10],2);
  testEquals(adjwgt[11],3);

  testEquals(adjncy[12],3);
  testEquals(adjncy[13],4);
  testEquals(adjwgt[12],4);
  testEquals(adjwgt[13],5);
}


void Test::run()
{
  // test metis
  writeMatrix("/tmp/wildriver_test.graph");
  readMatrix("/tmp/wildriver_test.graph");

  writeGraph("/tmp/wildriver_test.graph");
  readGraph("/tmp/wildriver_test.graph");

  // test csr
  writeMatrix("/tmp/wildriver_test.csr");
  readMatrix("/tmp/wildriver_test.csr");

  writeGraph("/tmp/wildriver_test.csr");
  readGraph("/tmp/wildriver_test.csr");
}




}

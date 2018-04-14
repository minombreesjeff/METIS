/**
 * @file wildriver.h
 * @brief Top level C include for WildRiver
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015
 * @version 1
 *
 */




#ifndef WILDRIVER_H
#define WILDRIVER_H




#include <stdlib.h>
#include <stdint.h>



/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/

#define WILDRIVER_VER_MAJOR 0
#define WILDRIVER_VER_MINOR 1
#define WILDRIVER_VER_SUBMINOR 0




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


#ifndef WILDRIVER_DIMENSION_TYPE
#define WILDRIVER_DIMENSION_TYPE uint32_t
#endif
typedef WILDRIVER_DIMENSION_TYPE wildriver_dim_t;


#ifndef WILDRIVER_INDEX_TYPE
#define WILDRIVER_INDEX_TYPE uint32_t
#endif
typedef WILDRIVER_INDEX_TYPE wildriver_ind_t;


#ifndef WILDRIVER_VALUE_TYPE
#define WILDRIVER_VALUE_TYPE float
#endif
typedef WILDRIVER_VALUE_TYPE wildriver_val_t;


typedef enum wildriver_format_t {
  WILDIRVER_FORMAT_AUTO,
  WILDRIVER_FORMAT_METIS,
  WILDRIVER_FORMAT_CHACO,
  WILDRIVER_FORMAT_GRAPH,
  WILDRIVER_FORMAT_MATRIXMARKET,
  WILDRIVER_FORMAT_CSR,
  WILDRIVER_FORMAT_BCSR,
  WILDRIVER_FORMAT_CLUTO,
  WILDRIVER_FORMAT_HMETIS,
  WILDRIVER_FORMAT_PATOH
} wildriver_format_t;




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Read a matrix from the given path into a CSR data-structure.
 *
 * @param fname The filename/path of the matrix file.
 * @param r_nrows The number of rows in the matrix (output).
 * @param r_ncols The number of columns in the matrix (output).
 * @param r_nnz The number of non-zeros in the matrix (output).
 * @param r_rowptr The the starting index for each row * (output).
 * @param r_rowind The column indices for entries in each row (output).
 * @param r_rowval The value of the entries in each row (output).
 *
 * @return 1 on success, 0 otherwise.
 */
int wildriver_read_matrix(
    char const * fname,
    wildriver_dim_t * r_nrows,
    wildriver_dim_t * r_ncols,
    wildriver_ind_t * r_nnz,
    wildriver_ind_t ** r_rowptr,
    wildriver_dim_t ** r_rowind,
    wildriver_val_t ** r_rowval);


/**
 * @brief Write a matrix file to the given path from a CSR data-structure.
 *
 * @param fname The filename/path of the matrix file.
 * @param nrows The number of rows in the matrix.
 * @param ncols The number of columns in the matrix.
 * @param nnz The number of non-zeros in the matrix.
 * @param rowptr The starting index for each row.
 * @param rowind The column indices for entries in each row. 
 * @param rowval The value of the entries in each row.
 *
 * @return 1 on success, 0 otherwise. 
 */
int wildriver_write_matrix(
    char const * fname,
    wildriver_dim_t nrows,
    wildriver_dim_t ncols,
    wildriver_ind_t nnz,
    wildriver_ind_t const * rowptr,
    wildriver_dim_t const * rowind,
    wildriver_val_t const * rowval);


/**
 * @brief Read a graph from the given path into a CSR data-structure.
 *
 * @param fname The filename/path of the graph file.
 * @param r_nvtxs The number of vertices in the graph (output).
 * @param r_nedges The number of edges in the graph (output, optional).
 * @param r_nvwgts The number of vertex weights in the graph (output,
 * optional).
 * @param r_ewgts Whether or not edge weights are present in the graph file
 * (output, optional).
 * @param r_xadj The adjacency list pointer (output).
 * @param r_adjncy The adjacency list (output).
 * @param r_vwgt The vertex weights (output, optional).
 * @param r_adjwgt The edge weights (output, optional).
 *
 * @return 1 on success, 0 otherwise.
 */
int wildriver_read_graph(
    char const * fname,
    wildriver_dim_t * r_nvtxs,
    wildriver_ind_t * r_nedges,
    int * r_nvwgts,
    int * r_ewgts,
    wildriver_ind_t ** r_xadj,
    wildriver_dim_t ** r_adjncy,
    wildriver_val_t ** r_vwgt,
    wildriver_val_t ** r_adjwgt);


/**
 * @brief Write a graph to the given path from a CSR data-structure.
 *
 * @param fname The filename/path of the graph file.
 * @param nvtxs The number of vertices in the graph.
 * @param nedges The number of edges in the graph.
 * @param nvwgts The number of vertex weights.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 *
 * @return 1 on success, 0 otherwise.
 */
int wildriver_write_graph(
    char const * fname,
    wildriver_dim_t nvtxs,
    wildriver_ind_t nedges,
    int nvwgts,
    wildriver_ind_t const * xadj,
    wildriver_dim_t const * adjncy,
    wildriver_val_t const * vwgt,
    wildriver_val_t const * adjwgt);




#ifdef __cplusplus
}
#endif




#endif

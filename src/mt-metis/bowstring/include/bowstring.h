/**
 * @file bowstring.h
 * @brief Library external header file
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2015, Dominique LaSalle
 * @version 1
 * @date 2013-04-09
 */




#ifndef BOWSTRING_H
#define BOWSTRING_H




#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define BOWSTRING_VER_MAJOR 0
#define BOWSTRING_VER_MINOR 2
#define BOWSTRING_VER_SUBMINOR 2




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


/**
 * @brief Return values for status
 */
typedef enum bowstring_error_t {
  /* success messages are 0x0* */
  BOWSTRING_SUCCESS=0x01,
  /* file errors are 0x1* */
  BOWSTRING_ERROR_FILENOTFOUND=0x10,
  BOWSTRING_ERROR_FILEREAD=0x11,
  BOWSTRING_ERROR_FILEWRITE=0x12,
  BOWSTRING_ERROR_FILECLOSE=0x13,
  BOWSTRING_ERROR_FILEPERMISSIONDENIED=0x14,
  /* input errors are 0x2* */
  BOWSTRING_ERROR_INVALIDINPUT=0x20,
  BOWSTRING_ERROR_INVALIDVALUE=0x21,
  /* resource errors */
  BOWSTRING_ERROR_NOTENOUGHMEMORY=0x30,
  /* network errors are 0xA* */
  BOWSTRING_ERROR_MPICALL=0xA0,
  /* special errors are 0xF* */
  BOWSTRING_ERROR_UNKNOWN=0xF0,
  BOWSTRING_ERROR_UNIMPLEMENTED=0xFF
} bowstring_error_t;


typedef enum bowstring_format_t {
  BOWSTRING_FORMAT_AUTO,
  BOWSTRING_FORMAT_METIS,
  BOWSTRING_FORMAT_CSR,
  BOWSTRING_FORMAT_CLOUD9,
  BOWSTRING_FORMAT_SNAP,
  BOWSTRING_FORMAT_DIMACS,
  BOWSTRING_FORMAT_PEGASUS,
  BOWSTRING_FORMAT_MATRIXMARKET,
  BOWSTRING_FORMAT_EDGELIST,
  BOWSTRING_FORMAT_STP,
  BOWSTRING_FORMAT_NERSTRAND,
  BOWSTRING_FORMAT_NBG
} bowstring_graph_type_t;


typedef enum bowstring_tree_t {
  BOWSTRING_TREE_MST,
  BOWSTRING_TREE_RST,
  BOWSTRING_TREE_DFS,
  BOWSTRING_TREE_BFS
} bowstring_tree_type_t;


typedef enum bowstring_order_t {
  BOWSTRING_ORDER_BFS,
  BOWSTRING_ORDER_DFS,
  BOWSTRING_ORDER_RANDOM,
  BOWSTRING_ORDER_INCDEGREE,
  BOWSTRING_ORDER_DECDEGREE,
  BOWSTRING_ORDER_POST,
  BOWSTRING_ORDER_RCM
} bowstring_order_type_t;


typedef enum bowstring_edgerank_t {
  BOWSTRING_EDGERANK_NI,
  BOWSTRING_EDGERANK_MST,
  BOWSTRING_EDGERANK_AST,
  BOWSTRING_EDGERANK_LST
} bowstring_edgerank_t;


typedef enum bowstring_reweight_t {
  BOWSTRING_REWEIGHT_NONE,
  BOWSTRING_REWEIGHT_EXACT,
  BOWSTRING_REWEIGHT_APPROX
} bowstring_reweight_t;


typedef enum bowstring_contract_t {
  BOWSTRING_CONTRACT_SUM,
  BOWSTRING_CONTRACT_MIN,
  BOWSTRING_CONTRACT_MAX
} bowstring_contract_t;


typedef enum bowtring_compression_t {
  BOWSTRING_COMPRESSION_NONE=0x0000,
  BOWSTRING_COMPRESSION_GZIP1=0x0011,
  BOWSTRING_COMPRESSION_GZIP2=0x0012,
  BOWSTRING_COMPRESSION_GZIP3=0x0013,
  BOWSTRING_COMPRESSION_GZIP4=0x0014,
  BOWSTRING_COMPRESSION_GZIP5=0x0015,
  BOWSTRING_COMPRESSION_GZIP6=0x0016,
  BOWSTRING_COMPRESSION_GZIP7=0x0017,
  BOWSTRING_COMPRESSION_GZIP8=0x0018,
  BOWSTRING_COMPRESSION_GZIP9=0x0019
} bowstring_compression_t;


#ifdef BOWSTRING_64BIT_VERTICES
typedef uint64_t bowstring_vtx_t;
#else
typedef uint32_t bowstring_vtx_t;
#endif /* BOWSTRING_64BIT_VERTICES */
#ifdef BOWSTRING_64BIT_EDGES
typedef uint64_t bowstring_adj_t;
#else
typedef uint32_t bowstring_adj_t;
#endif /* BOWSTRING_64BIT_EDGES */
#ifdef BOWSTRING_64BIT_WEIGHTS
#ifdef BOWSTRING_INT_WEIGHTS
typedef int64_t bowstring_wgt_t;
#else
typedef double bowstring_wgt_t;
#endif /* BOWSTRING_INT_WEIGHTS */
#else
#ifdef BOWSTRING_INT_WEIGHTS
typedef int32_t bowstring_wgt_t;
#else
typedef float bowstring_wgt_t;
#endif /* BOWSTRING_INT WEIGHTS */
#endif /* BOWSTRING_DOUBLE_WEIGHTS */


#ifdef BOWSTRING_64BIT_VLABELS
typedef uint64_t bowstring_vlbl_t;
#else
typedef uint32_t bowstring_vlbl_t;
#endif /* BOWSTRING_64BIT_VLABELS */
#ifdef BOWSTRING_64BIT_ELABELS
typedef uint64_t bowstring_elbl_t;
#else
typedef uint32_t bowstring_elbl_t;
#endif /* BOWSTRING_64BIT_ELABELS */


#ifdef BOWSTRING_COORD_DOUBLE
typedef double bowstring_coord_t;
#else
typedef float bowstring_coord_t;
#endif /* BOWSTRING_COORD_DOUBLE */




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Read a graph structure from a file.
 *
 * @param filename The file name to read the graph from.
 * @param graphtype The file format.
 * @param r_nvtxs A reference to the number of vertices in the graph (output).
 * @param r_xadj A reference to the adjacecny list pointer (output).
 * @param r_adjncy A reference to the adjacency list (output).
 * @param r_vwgt A reference to the vertex weights (output). May be NULL.
 * @param r_adjwgt A reference to the edge weights (output). May be NULL.
 *
 * @return BOWSTRING_SUCCESS unless an error occurs. 
 */
int bowstring_read_graph(
    char const * filename, 
    int graphtype, 
    bowstring_vtx_t * r_nvtxs, 
    bowstring_adj_t ** r_xadj, 
    bowstring_vtx_t ** r_adjncy, 
    bowstring_wgt_t ** r_vwgt, 
    bowstring_wgt_t ** r_adjwgt);


/**
 * @brief Write a graph structure to the specified file using the specified
 * format.
 *
 * @param filename The file to write to.
 * @param graphtype The type of file to write.
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 *
 * @return BOWSTRING_SUCCESS unless an error occured.
 */
int bowstring_write_graph(
    char const * filename, 
    int graphtype, 
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t const * xadj, 
    bowstring_vtx_t const * adjncy, 
    bowstring_wgt_t const * vwgt, 
    bowstring_wgt_t const * adjwgt);


/**
 * @brief Write vertex labels to a file.
 *
 * @param filename The name of the file to write to.
 * @param nvtxs The number of vertices.
 * @param labels The labels of the vertices.
 *
 * @return BOWSTRING_SUCCESS unless an error is encountered.
 */
int bowstring_write_vertex_labels(
    char const * filename, 
    bowstring_vtx_t nvtxs, 
    bowstring_vlbl_t const * labels);


/**
 * @brief Construct a tree (forest) structure given a graph and type of tree. 
 * The tree (forest) will contain the same vertex set as the input, thus 
 * vertex numbers and weights remane unchanged.
 *
 * @param treetype The type of tree to build.
 * @param nvtxs The number of vertices in the graph.
 * @param gxadj The adjacency list pointer of the graph.
 * @param gadjncy The adjacecny list of the graph.
 * @param gadjwgt The edge weights of the graph.
 * @param xadj A reference to the adjacency list pointer of the tree (output,
 * may be NULL).
 * @param adjncy A reference to the adjacency list of the tree (output, may be 
 * NULL).
 * @param adjwgt A reference to the edge weights of the tree (output, may be 
 * NULL).
 * @param adjmask A reference to the edge mask; a 1 indicates that an edge is
 * part of the tree (output, may be NULL). 
 */
void bowstring_build_tree(
    int treetype, 
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t const * gxadj, 
    bowstring_vtx_t const * gadjncy, 
    bowstring_wgt_t const * const gadjwgt, 
    bowstring_adj_t ** xadj, 
    bowstring_vtx_t ** adjncy, 
    bowstring_wgt_t ** adjwgt,
    int ** r_adjmask);


/**
 * @brief Build a reverse adjacency index, such that:
 *
 * radj[radj[j]] = j
 *
 * and:
 *
 * for (j=xadj[i];j<xadj[i+1];++j) {
 *   adjncy[radj[j]] == i
 * }
 *
 * @param nvtxs Number of vertices in the graph
 * @param xadj The array pointing to the start of each vertices adjacency list
 * @param adjncy The adjacency list of each vertex
 * @param radj The reverse adjacency list -- output array of length 
 *   xadj[nvtxs].
 */
void bowstring_build_adjncy_index(
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy, 
    bowstring_adj_t * radj);


/**
 * @brief Remove a fraction of the edges from a graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weight (may be NULL).
 * @param type The method for selecting the edges to remove.
 * @param frac The fraction of edges to remove.
 * @param r_xadj A reference to the modified adjacency list pointer.
 * @param r_adjncy A reference to the modified adjacecny list.
 * @param r_adjwgt A reference to the modified edge weights.
 */
void bowstring_remove_edges(
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t const * xadj, 
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt, 
    int type, 
    double frac, 
    bowstring_adj_t ** r_xadj, 
    bowstring_vtx_t ** r_adjncy, 
    bowstring_wgt_t ** r_adjwgt);


/**
 * @brief Induce a subgraph given a set of vertices (based on the corresponding
 * values of the present array).
 *
 * @param nvtxs The total number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param vwgt The vertex weights.
 * @param present The array indicating which vertices are present in the
 *   subgraph.
 * @param r_nvtxs A reference to the number of vertices in the subgraph
 *   (output).
 * @param r_xadj A reference to the subgraph's adjacency list pointer (output).
 * @param r_adjncy A reference to the subgraph's adjacncy list (output).
 * @param r_adjwgt A reference to the subgraph's edge weight (output).
 * @param r_vwgt A reference to the subgraph's vertex weight (output).
 * @param r_alias A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 * @param r_alias A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 */
void bowstring_induce_subgraph(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_wgt_t const * vwgt,
    int const * present,
    bowstring_vtx_t * r_nvtxs,
    bowstring_adj_t ** r_xadj,
    bowstring_vtx_t ** r_adjncy,
    bowstring_wgt_t ** r_adjwgt,
    bowstring_wgt_t ** r_vwgt,
    bowstring_vtx_t ** r_alias,
    bowstring_vtx_t ** r_rename);


/**
 * @brief Induce a set of subgraphs given a vertex coloring/partitioning of 
 *    the graph.
 *
 * @param nvtxs The total number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param vwgt The vertex weights.
 * @param part The vertex colors/partitions.
 * @param npart The number of colors/partitions.
 * @param present The array containing the color/partition of each vertex.
 * @param xnvtxs The number of vertices in the subgraphs (should be of
 * length npart).
 * @param xxadj The array containing the subgraphs' adjacency list pointers
 * (should be of size npart).
 * @param xadjncy The array containing the subgraphs' adjacncy lists (should be
 * of length npart).
 * @param xadjwgt The array containing the subgraphs' edge weights (should be
 * of length npart). May be NULL.
 * @param xvwgt The array containing the subgraphs' vertex weights (should be
 * of length npart). May be NULL.
 * @param xalias The array containing the original vertex numbers (should be of
 * size npart). May be NULL.
 * @param r_alias A reference to the original graph's vertex numbers in the
 * subgraphs (output). May be NULL.
 */
void bowstring_induce_subgraphs(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_wgt_t const * vwgt,
    bowstring_vlbl_t const * const part,
    bowstring_vlbl_t const npart,
    bowstring_vtx_t * xnvtxs,
    bowstring_adj_t ** xxadj,
    bowstring_vtx_t ** xadjncy,
    bowstring_wgt_t ** xadjwgt,
    bowstring_wgt_t ** xvwgt,
    bowstring_vtx_t ** xalias,
    bowstring_vtx_t ** r_rename);




/**
 * @brief Extract a subgraph given a set of vertices and edges (based on the 
 * corresponding values of the vpresent and epresent arrays).
 *
 * @param nvtxs The total number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights (can be NULL).
 * @param vwgt The vertex weights (can be NULL).
 * @param present The array indicating which vertices are present in the
 *   subgraph.
 * @param r_nvtxs A reference to the number of vertices in the subgraph
 *   (output).
 * @param r_xadj A reference to the subgraph's adjacency list pointer (output).
 * @param r_adjncy A reference to the subgraph's adjacncy list (output).
 * @param r_adjwgt A reference to the subgraph's edge weight (output).
 * @param r_vwgt A reference to the subgraph's vertex weights (output).
 * @param r_alias A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 * @param r_alias A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 */
void bowstring_extract_subgraph(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_wgt_t const * vwgt,
    int const * vpresent,
    int const * epresent,
    bowstring_vtx_t * r_nvtxs,
    bowstring_adj_t ** r_xadj,
    bowstring_vtx_t ** r_adjncy,
    bowstring_wgt_t ** r_adjwgt,
    bowstring_wgt_t ** r_vwgt,
    bowstring_vtx_t ** r_alias,
    bowstring_vtx_t ** r_rename);


/**
 * @brief Label the connected components of a graph.
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param adjncy The adjacency list.
 * @param r_lbl A reference to the component labels (output). May be NULL.
 * @param r_nlbl A reference to the number of components (output). May be NULL.
 */
void bowstring_label_components(
    bowstring_vtx_t const nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_vlbl_t ** r_lbl,
    bowstring_vlbl_t * r_nlbl);


/**
 * @brief Label the connected partition components of a graph.
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param adjncy The adjacency list.
 * @param where The partition ID's of each vertex.
 * @param r_lbl A reference to the component labels (output). May be NULL.
 * @param r_nlbl A reference to the number of components (output). May be NULL.
 */
void bowstring_label_partition_components(
    bowstring_vtx_t const nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_vlbl_t const * where,
    bowstring_vlbl_t ** r_lbl,
    bowstring_vlbl_t * r_nlbl);


/**
 * @brief Check that a graph structure is valid.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights (can be NULL)
 *
 * @return 1 if it is a valid graph.
 */
int bowstring_check_graph(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt);


/**
 * @brief Create voronoi regions on a graph for a given set of vertices.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param sources The set of source vertices.
 * @param nsources The number of source vertices.
 * @param where The voronoi region ID of each vertex (output).
 * @param parent The parent of each vertex in the voronoi region (output).
 * @param dist The distance of each vertex to its source vertex (output).
 */
void bowstring_voronoi_regions(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_vtx_t const * sources,
    bowstring_vtx_t nsources,
    bowstring_vtx_t * where,
    bowstring_vtx_t * parent,
    bowstring_wgt_t * dist);


/**
 * @brief Create a voronoi diagram of a graph given a subset of vertices. 
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights (may be NULL).
 * @param nsources The number of source vertices.
 * @param where The voronoi region labels.
 * @param dist The distance of each vertex to the center of its voronoi region.
 * @param r_vxadj A reference to the adjacency list pointer for the voronoi
 * graph (may be NULL, output).
 * @param r_vadjncy A reference to the adjacency list for the voronoi graph
 * (may be NULL, output).
 * @param r_vadjwgt A reference to the edge weight list for the voronoi graph
 * (may be NULL, output).
 * @param r_adjsrc The source vertex in the original graph for the edge in the
 * voronoi diagram (may be NULL, output).
 * @param r_adjori The source edge index in the original graph for the edge in 
 * the voronoi diagram (may be NULL, output).
 */
void bowstring_voronoi_diagram(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_vtx_t nsources,
    bowstring_vtx_t const * where,
    bowstring_wgt_t const * dist,
    bowstring_adj_t ** r_vxadj,
    bowstring_vtx_t ** r_vadjncy,
    bowstring_wgt_t ** r_vadjwgt,
    bowstring_vtx_t ** r_adjsrc,
    bowstring_adj_t ** r_adjori);


/**
 * @brief Re-order a graph given a permutateion.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer array.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 * @param perm The permuation.
 */
void bowstring_order_graph(
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t * xadj, 
    bowstring_vtx_t * adjncy, 
    bowstring_wgt_t * vwgt, 
    bowstring_wgt_t * adjwgt, 
    bowstring_vtx_t const * perm);


/**
 * @brief Generate a permutation of the vertices of the graph.
 *
 * @param nvtxs The number of vertices in the tree.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights (may be NULL).
 * @param adjwgt The edge weights (may be NULL).
 * @param perm The permutation of the vertices.
 */
void bowstring_permutation(
    int ordering,
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * vwgt,
    bowstring_wgt_t const * adjwgt,
    bowstring_vtx_t * perm);




/**
 * @brief Read a graph structure from a file.
 *
 * @param filename The file name to read the graph from.
 * @param graphtype The file format.
 * @param r_nvtxs A reference to the number of vertices in the graph (output).
 * @param r_xadj A reference to the adjacecny list pointer (output).
 * @param r_adjncy A reference to the adjacency list (output).
 * @param r_vwgt A reference to the vertex weights (output). May be NULL.
 * @param r_adjwgt A reference to the edge weights (output). May be NULL.
 *
 * @return BOWSTRING_SUCCESS unless an error occurs. 
 */
int bowstring_read_graph(
    char const * filename, 
    int graphtype, 
    bowstring_vtx_t * r_nvtxs, 
    bowstring_adj_t ** r_xadj, 
    bowstring_vtx_t ** r_adjncy, 
    bowstring_wgt_t ** r_vwgt, 
    bowstring_wgt_t ** r_adjwgt);


/**
 * @brief Write a graph structure to the specified file using the specified
 * format.
 *
 * @param filename The file to write to.
 * @param graphtype The type of file to write.
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 *
 * @return BOWSTRING_SUCCESS unless an error occured.
 */
int bowstring_write_graph(
    char const * filename, 
    int graphtype, 
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t const * xadj, 
    bowstring_vtx_t const * adjncy, 
    bowstring_wgt_t const * vwgt, 
    bowstring_wgt_t const * adjwgt);


/**
 * @brief Write vertex labels to a file.
 *
 * @param filename The name of the file to write to.
 * @param nvtxs The number of vertices.
 * @param labels The labels of the vertices.
 *
 * @return BOWSTRING_SUCCESS unless an error is encountered.
 */
int bowstring_write_vertex_labels(
    char const * filename, 
    bowstring_vtx_t nvtxs, 
    bowstring_vlbl_t const * labels);


/**
 * @brief Construct a tree (forest) structure given a graph and type of tree. 
 * The tree (forest) will contain the same vertex set as the input, thus 
 * vertex numbers and weights remane unchanged.
 *
 * @param treetype The type of tree to build.
 * @param nvtxs The number of vertices in the graph.
 * @param gxadj The adjacency list pointer of the graph.
 * @param gadjncy The adjacecny list of the graph.
 * @param gadjwgt The edge weights of the graph.
 * @param xadj A reference to the adjacency list pointer of the tree (output,
 * may be NULL).
 * @param adjncy A reference to the adjacency list of the tree (output, may be 
 * NULL).
 * @param adjwgt A reference to the edge weights of the tree (output, may be 
 * NULL).
 * @param adjmask A reference to the edge mask; a 1 indicates that an edge is
 * part of the tree (output, may be NULL). 
 */
void bowstring_build_tree(
    int treetype, 
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t const * gxadj, 
    bowstring_vtx_t const * gadjncy, 
    bowstring_wgt_t const * const gadjwgt, 
    bowstring_adj_t ** xadj, 
    bowstring_vtx_t ** adjncy, 
    bowstring_wgt_t ** adjwgt,
    int ** r_adjmask);


/**
 * @brief Build a reverse adjacency index, such that:
 *
 * radj[radj[j]] = j
 *
 * and:
 *
 * for (j=xadj[i];j<xadj[i+1];++j) {
 *   adjncy[radj[j]] == i
 * }
 *
 * @param nvtxs Number of vertices in the graph
 * @param xadj The array pointing to the start of each vertices adjacency list
 * @param adjncy The adjacency list of each vertex
 * @param radj The reverse adjacency list -- output array of length 
 *   xadj[nvtxs].
 */
void bowstring_build_adjncy_index(
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy, 
    bowstring_adj_t * radj);


/**
 * @brief Remove a fraction of the edges from a graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weight (may be NULL).
 * @param type The method for selecting the edges to remove.
 * @param frac The fraction of edges to remove.
 * @param r_xadj A reference to the modified adjacency list pointer.
 * @param r_adjncy A reference to the modified adjacecny list.
 * @param r_adjwgt A reference to the modified edge weights.
 */
void bowstring_remove_edges(
    bowstring_vtx_t nvtxs, 
    bowstring_adj_t const * xadj, 
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt, 
    int type, 
    double frac, 
    bowstring_adj_t ** r_xadj, 
    bowstring_vtx_t ** r_adjncy, 
    bowstring_wgt_t ** r_adjwgt);


/**
 * @brief Induce a subgraph given a set of vertices (based on the corresponding
 * values of the present array).
 *
 * @param nvtxs The total number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param vwgt The vertex weights.
 * @param present The array indicating which vertices are present in the
 *   subgraph.
 * @param r_nvtxs A reference to the number of vertices in the subgraph
 *   (output).
 * @param r_xadj A reference to the subgraph's adjacency list pointer (output).
 * @param r_adjncy A reference to the subgraph's adjacncy list (output).
 * @param r_adjwgt A reference to the subgraph's edge weight (output).
 * @param r_vwgt A reference to the subgraph's vertex weight (output).
 * @param r_alias A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 * @param r_alias A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 */
void bowstring_induce_subgraph(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_wgt_t const * vwgt,
    int const * present,
    bowstring_vtx_t * r_nvtxs,
    bowstring_adj_t ** r_xadj,
    bowstring_vtx_t ** r_adjncy,
    bowstring_wgt_t ** r_adjwgt,
    bowstring_wgt_t ** r_vwgt,
    bowstring_vtx_t ** r_alias,
    bowstring_vtx_t ** r_rename);


/**
 * @brief Induce a set of subgraphs given a vertex coloring/partitioning of 
 *    the graph.
 *
 * @param nvtxs The total number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param vwgt The vertex weights.
 * @param part The vertex colors/partitions.
 * @param npart The number of colors/partitions.
 * @param present The array containing the color/partition of each vertex.
 * @param xnvtxs The number of vertices in the subgraphs (should be of
 * length npart).
 * @param xxadj The array containing the subgraphs' adjacency list pointers
 * (should be of size npart).
 * @param xadjncy The array containing the subgraphs' adjacncy lists (should be
 * of length npart).
 * @param xadjwgt The array containing the subgraphs' edge weights (should be
 * of length npart). May be NULL.
 * @param xvwgt The array containing the subgraphs' vertex weights (should be
 * of length npart). May be NULL.
 * @param xalias The array containing the original vertex numbers (should be of
 * size npart). May be NULL.
 * @param r_alias A reference to the original graph's vertex numbers in the
 * subgraphs (output). May be NULL.
 */
void bowstring_induce_subgraphs(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_wgt_t const * vwgt,
    bowstring_vlbl_t const * const part,
    bowstring_vlbl_t const npart,
    bowstring_vtx_t * xnvtxs,
    bowstring_adj_t ** xxadj,
    bowstring_vtx_t ** xadjncy,
    bowstring_wgt_t ** xadjwgt,
    bowstring_wgt_t ** xvwgt,
    bowstring_vtx_t ** xalias,
    bowstring_vtx_t ** r_rename);


/**
 * @brief Extract a subgraph given a set of vertices and edges (based on the 
 * corresponding values of the vpresent and epresent arrays).
 *
 * @param nvtxs The total number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights (can be NULL).
 * @param vwgt The vertex weights (can be NULL).
 * @param present The array indicating which vertices are present in the
 *   subgraph.
 * @param r_nvtxs A reference to the number of vertices in the subgraph
 *   (output).
 * @param r_xadj A reference to the subgraph's adjacency list pointer (output).
 * @param r_adjncy A reference to the subgraph's adjacncy list (output).
 * @param r_adjwgt A reference to the subgraph's edge weight (output).
 * @param r_vwgt A reference to the subgraph's vertex weights (output).
 * @param r_alias A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 * @param r_alias A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 */
void bowstring_extract_subgraph(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_wgt_t const * vwgt,
    int const * vpresent,
    int const * epresent,
    bowstring_vtx_t * r_nvtxs,
    bowstring_adj_t ** r_xadj,
    bowstring_vtx_t ** r_adjncy,
    bowstring_wgt_t ** r_adjwgt,
    bowstring_wgt_t ** r_vwgt,
    bowstring_vtx_t ** r_alias,
    bowstring_vtx_t ** r_rename);


/**
 * @brief Label the connected components of a graph.
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param adjncy The adjacency list.
 * @param r_lbl A reference to the component labels (output). May be NULL.
 * @param r_nlbl A reference to the number of components (output). May be NULL.
 */
void bowstring_label_components(
    bowstring_vtx_t const nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_vlbl_t ** r_lbl,
    bowstring_vlbl_t * r_nlbl);


/**
 * @brief Label the connected partition components of a graph.
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param adjncy The adjacency list.
 * @param where The partition ID's of each vertex.
 * @param r_lbl A reference to the component labels (output). May be NULL.
 * @param r_nlbl A reference to the number of components (output). May be NULL.
 */
void bowstring_label_partition_components(
    bowstring_vtx_t const nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_vlbl_t const * where,
    bowstring_vlbl_t ** r_lbl,
    bowstring_vlbl_t * r_nlbl);


/**
 * @brief Check that a graph structure is valid.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights (can be NULL)
 *
 * @return 1 if it is a valid graph.
 */
int bowstring_check_graph(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt);




#ifdef __cplusplus
}
#endif




#endif

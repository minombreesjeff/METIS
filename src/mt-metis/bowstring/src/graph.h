/**
 * @file graph.h
 * @brief Misc function prototypes for graph operations 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-08-07
 */




#ifndef GRAPH_H
#define GRAPH_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define check_graph __bowstring_check_graph
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
int check_graph(
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt);


#define apply_edge_mask __bowstring_apply_edge_mask
/**
 * @brief Sparsen a graph such that only edges marked by the adjmask array are
 * kept.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param gxadj The original adjacency list pointer.
 * @param gadjncy The original adjacency list.
 * @param gadjwgt The original edge weights.
 * @param adjmask The edge mask -- edges marked with 0's will be removed.
 * @param r_xadj The output adjacency list pointer.
 * @param r_adjncy The output adjacency list.
 * @param r_adjwgt The output edge weights.
 *
 */
void apply_edge_mask(
    vtx_t nvtxs, 
    adj_t const * gxadj, 
    vtx_t const * gadjncy, 
    wgt_t const * gadjwgt, 
    int const * adjmask, 
    adj_t ** r_xadj, 
    vtx_t ** r_adjncy, 
    wgt_t ** r_adjwgt); 


#define mark_boundary_vertices __bowstring_mark_boundary_vertices
/**
 * @brief Mark all vertices located on the boundary of vertex labels. The bnd
 * parameter should be allocated to a size of nvtxs before being passed into
 * this function. Upon function return, the indecies of vertices on the
 * boundary will be non-zero.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param where The vertex labels.
 * @param bnd The boundary marker array (should be of length nvtxs). 
 *
 * @return The number of vertices on the boundary.
 */
vtx_t mark_boundary_vertices(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    vlbl_t const * where, 
    int * bnd);


#define find_boundary_vertices __bowstring_find_boundary_vertices
/**
 * @brief Build a list of all vertices located on label boundaries.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param where The vertex labels.
 * @param nbnd A reference to the number of vertices on the boundary (output).
 * @param bnd A reference to the vertices on the boundary (output).
 *
 */
void find_boundary_vertices(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    vlbl_t const * where, 
    vtx_t * nbnd,
    vtx_t ** bnd);


#define build_adjncy_index __bowstring_build_adjncy_index
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
 *
 * @return BOWSTRING_SUCCESS on success. 
 */
int build_adjncy_index(
    vtx_t nvtxs, 
    adj_t const * xadj,
    vtx_t const * adjncy, 
    adj_t * radj);


#define build_adjncy_index_rem __bowstring_build_adjncy_index_rem
/**
 * @brief Build a reverse adjacency index supporting remote edges, such that:
 *
 * radj[radj[j]] = j
 *
 * and:
 *
 * for (j=xadj[i];j<xadj[i+1];++j) {
 *   adjncy[radj[j]] == i
 * }
 *
 * except for edges where adjncy[j] > nvtxs, in which case radj[j] = NULL_ADJ.
 *
 * @param nvtxs Number of vertices in the graph
 * @param xadj The array pointing to the start of each vertices adjacency list
 * @param adjncy The adjacency list of each vertex
 * @param radj The reverse adjacency list -- output array of length 
 *   xadj[nvtxs].
 *
 * @return BOWSTRING_SUCCESS on success. 
 */
int build_adjncy_index_rem(
    vtx_t nvtxs, 
    adj_t const * xadj,
    vtx_t const * adjncy, 
    adj_t * radj);


#define neighborhoodify __bowstring_neighborhoodify
/**
 * @brief Assign each vertex a neighborhood ID, by peforming many BFS's until
 * each vertex is assigned a neighborhood, and no neighborhood has more than
 * nbrsize members. The nighborhood assignements are written to nbrhd, which
 * should be of size nvtxs.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param nbrsize The maximum size of a neighborhood.
 * @param nbrhd The neighborhood assignements for each vertex.
 *
 * @return The number of neighborhoods created.
 */
vtx_t neighborhoodify(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy,
    vtx_t nbrsize,
    vtx_t * nbrhd);


#define induce_subgraph __bowstring_induce_subgraph
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
 * @param r_rename A reference to the subgraph's vertex numbers in the original
 *   graph (output). May be NULL.
 */
void induce_subgraph(
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt,
    wgt_t const * vwgt,
    int const * present,
    vtx_t * r_nvtxs,
    adj_t ** r_xadj,
    vtx_t ** r_adjncy,
    wgt_t ** r_adjwgt,
    wgt_t ** r_vwgt,
    vtx_t ** r_alias,
    vtx_t ** r_rename);


#define induce_subgraphs __bowstring_induce_subgraphs
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
void induce_subgraphs(
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt,
    wgt_t const * vwgt,
    vlbl_t const * const part,
    vlbl_t const npart,
    vtx_t * xnvtxs,
    adj_t ** xxadj,
    vtx_t ** xadjncy,
    wgt_t ** xadjwgt,
    wgt_t ** xvwgt,
    vtx_t ** xalias,
    vtx_t ** r_rename);


#define extract_subgraph __bowstring_extract_subgraph
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
void extract_subgraph(
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt,
    wgt_t const * vwgt,
    int const * vpresent,
    int const * epresent,
    vtx_t * r_nvtxs,
    adj_t ** r_xadj,
    vtx_t ** r_adjncy,
    wgt_t ** r_adjwgt,
    wgt_t ** r_vwgt,
    vtx_t ** r_alias,
    vtx_t ** r_rename);


#define label_components __bowstring_label_components
/**
 * @brief Label the connected components of a graph.
 *
 * @param nvtxs The number of vertices.
 * @param xadj The adjacecny list pointer.
 * @param adjncy The adjacency list.
 * @param r_lbl A reference to the component labels (output). May be NULL.
 * @param r_nlbl A reference to the number of components (output). May be NULL.
 */
void label_components(
    vtx_t const nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    vlbl_t ** r_lbl,
    vlbl_t * r_nlbl);


#define label_partition_components __bowstring_label_partition_components
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
void label_partition_components(
    vtx_t const nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    vlbl_t const * where,
    vlbl_t ** r_lbl,
    vlbl_t * r_nlbl);


#define voronoi_regions __bowstring_voronoi_regions
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
void voronoi_regions(
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt,
    vtx_t const * sources,
    vtx_t nsources,
    vtx_t * where,
    vtx_t * parent,
    wgt_t * dist);


#define voronoi_add_sources __bowstring_voronoi_add_sources
/**
 * @brief Add a set of sources to existing voronoi regions.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param offset The number of existing voronoi regions.
 * @param sources The set of source vertices.
 * @param nsources The number of source vertices.
 * @param where The voronoi region ID of each vertex (output).
 * @param parent The parent of each vertex in the voronoi region (output).
 * @param dist The distance of each vertex to its source vertex (output).
 */
void voronoi_add_sources(
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt,
    vtx_t offset,
    vtx_t const * sources,
    vtx_t nsources,
    vtx_t * where,
    vtx_t * parent,
    wgt_t * dist);


#define voronoi_diagram __bowstring_voronoi_diagram
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
void voronoi_diagram(
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * adjwgt,
    vtx_t nsources,
    vtx_t const * where,
    wgt_t const * dist,
    adj_t ** r_vxadj,
    vtx_t ** r_vadjncy,
    wgt_t ** r_vadjwgt,
    vtx_t ** r_adjsrc,
    adj_t ** r_adjori);





#endif

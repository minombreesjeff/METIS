/**
 * @file bowstring_dynamic.h
 * @brief Library external header file for dynamic graphs. 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Dominique LaSalle
 * @version 1
 * @date 2014-10-22
 */




#ifndef BOWSTRING_DYNAMIC_H
#define BOWSTRING_DYNAMIC_H




#include "bowstring.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct bowstring_dynnode_t {
  bowstring_vtx_t * adj;
  bowstring_wgt_t * wgt;
  bowstring_vtx_t nadj;
} bowstring_dynnode_t;


typedef struct bowstring_dyngraph_t {
  bowstring_vtx_t nvtxs;
  bowstring_vtx_t max_nvtxs;
  bowstring_adj_t nedges;
  struct bowstring_dynnode_t * nodes;
  bowstring_vtx_t * htable;
} bowstring_dyngraph_t;





/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


/**
 * @brief Create a dynamic graph structure with a maximum of max_nvtxs nodes.
 *
 * @param max_nvtxs The maximum number of nodes in the graph.
 *
 * @return The newly allocated and initialized graph. 
 */
bowstring_dyngraph_t * bowstring_dg_create(
    bowstring_vtx_t max_nvtxs);


/**
 * @brief Add a vertex to the graph. The adjacency list is not added to
 * neighbor vertices. This does not increment nedges.
 *
 * @param v The new vertex number (must be less than the maximum number of
 * vertices).
 * @param adjncy The adjacency list of the vertex. This is copied and thus may
 * be free'd/modified after this function returns without affecting the dynamic
 * graph.
 * @param adjwgt The weights associated with the vertex's edges. This is copied
 * and thus may be free'd/modified after this function returns without 
 * affecting the dynamic graph.
 * @param nadj The size of the adjacency list.
 * @param graph The graph to add the vertex to.
 */
void bowstring_dg_addvtx(
    bowstring_vtx_t v,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt,
    bowstring_vtx_t nadj,
    bowstring_dyngraph_t * graph);


/**
 * @brief Add an edge to the graph to both vertices.
 *
 * @param v The first endpoint.
 * @param u The second endpoint.
 * @param w The weight of the edge.
 * @param graph The graph.
 */
void bowstring_dg_addedge(
    bowstring_vtx_t v,
    bowstring_vtx_t u,
    bowstring_wgt_t w,
    bowstring_dyngraph_t * graph);


/**
 * @brief Delete a vertex and its associated edges from a graph.
 *
 * @param v The vertex to delete.
 * @param graph The graph.
 */
void bowstring_dg_delvtx(
    bowstring_vtx_t v,
    bowstring_dyngraph_t * graph);


/**
 * @brief Delete an edge from the graph.
 *
 * @param v The first endpoint of the edge.
 * @param u The second endpoint of the edge.
 * @param graph The graph.
 */
void bowstring_dg_deledge(
    bowstring_vtx_t v,
    bowstring_vtx_t u,
    bowstring_dyngraph_t * graph);


/**
 * @brief Update teh weight of an edge in the graph.
 *
 * @param v The first endpoint of the edge.
 * @param u The second endpoint of the edge.
 * @param w The new weight of the edge.
 * @param graph The graph.
 */
void bowstring_dg_updateedge(
    bowstring_vtx_t v,
    bowstring_vtx_t u,
    bowstring_wgt_t w,
    bowstring_dyngraph_t * graph);


/**
 * @brief Check if two vertices are connected by an edge. 
 *
 * @param v The first vertex.
 * @param u The second vertex.
 * @param graph The dynamic graph.
 *
 * @return The weight connected the two edges.
 */
bowstring_wgt_t bowstring_dg_connected(
    bowstring_vtx_t v,
    bowstring_vtx_t u,
    bowstring_dyngraph_t const * graph);


/**
 * @brief Contract a group of vertice together in the graph.
 *
 * @param method The contraction method to use (bowstring_contract_t).
 * @param vtxs The vertices to contract together.
 * @param nvtxs The number of vertices to contract together.
 * @param graph The graph.
 *
 * @return The new vertex number.
 */
bowstring_vtx_t bowstring_dg_contract(
    int method,
    bowstring_vtx_t const * vtxs,
    bowstring_vtx_t nvtxs,
    bowstring_dyngraph_t * graph);


/**
 * @brief Convert a dynamic graph to a static graph.
 *
 * @param graph The dynamic graph to convert.
 * @param r_xadj A reference to the adjacency list poitner.
 * @param r_adjncy A reference to teh adjacency list.
 * @param r_adjwgt A reference to the edge weights.
 * @param r_alias A reference to teh vertex id alias array.
 */
void bowstring_dg_tostatic(
    bowstring_dyngraph_t const * graph,
    bowstring_adj_t ** r_xadj,
    bowstring_vtx_t ** r_adjncy,
    bowstring_wgt_t ** r_adjwgt,
    bowstring_vtx_t ** r_alias);


/**
 * @brief Convert a static graph into a dynamic graph. 
 *
 * @param nvtxs The number of vertices in the static graph.
 * @param xadj The adjacecny list pointer of the static graph.
 * @param adjncy The adjacnecy list of the static graph.
 * @param adjwgt The edge weights of the static graph.
 *
 * @return The dynamic graph.
 */
bowstring_dyngraph_t * bowstring_dg_todynamic(
    bowstring_vtx_t nvtxs,
    bowstring_adj_t const * xadj,
    bowstring_vtx_t const * adjncy,
    bowstring_wgt_t const * adjwgt);


/**
 * @brief Check a dynamic graph for validity.
 *
 * @param graph The graph to check.
 *
 * @return 1 if the graph is valid. 
 */
int bowstring_dg_check(
    bowstring_dyngraph_t const * graph);


bowstring_vtx_t bowstring_dg_ncomp(
    bowstring_dyngraph_t const * graph);


/**
 * @brief Free a dynamic graph and its associated memory.
 *
 * @param graph The graph to free.
 */
void bowstring_dg_free(
    bowstring_dyngraph_t * graph);


#endif

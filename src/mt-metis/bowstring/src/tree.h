/**
 * @file tree.h
 * @brief Functions for generating trees in a graph
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-08-06
 */




#ifndef TREE_H
#define TREE_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define build_dfs_tree __bowstring_build_dfs_tree
/**
 * @brief Build a DFS tree.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param start The vertex to start the search from.
 * @param label The label of each vertex in the DFS. 
 * @param rank The rank of each vertex in the DFS.
 * @param parent The parent of each vertex in the DFS.
 * @param adjmap The binary indicating which edges are part of the DFS.
 *
 */
void build_dfs_tree(
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    vtx_t start, 
    vtx_t * label, 
    vtx_t * rank, adj_t * parent, 
    int * adjmap);


#define build_bfs_tree __bowstring_build_bfs_tree
/**
 * @brief Build a BFS tree.
 *
 * @param nvtxs The number of the vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param start The vertex to start the BFS from.
 * @param label The label of each vertex in the BFS.
 * @param rank The rank of each vertex in the BFS.
 * @param parent The parent of each vertex.
 * @param adjmap The binary arrary indicating which edges are part of the BFS.
 */
void build_bfs_tree(
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    vtx_t start, 
    vtx_t * label, 
    vtx_t * rank, 
    adj_t * parent, 
    int * adjmap);


#define build_mst_tree __bowstring_build_mst_tree
/**
 * @brief Build a minimum spanning tree.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param adjmap The binary array indicating which edges are part of the tree.
 *
 */
void build_mst_tree(
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    const wgt_t * adjwgt, 
    int * adjmap);


#define build_rst_tree __bowstring_build_rst_tree
/**
 * @brief Build a random spanning tree.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjmap The binary array indicating which edges are part of the RST.
 */
void build_rst_tree(
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    int * adjmap);




#endif

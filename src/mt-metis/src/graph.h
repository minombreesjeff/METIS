/**
 * @file graph.h
 * @brief Types and functions for distributed graph objects.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-16
 */




#ifndef MTMETIS_GRAPH_H
#define MTMETIS_GRAPH_H




#include "base.h"
#include "ctrl.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct graphdist_t {
  tid_t nthreads;
  vtx_t mask;
  int shift;
  vtx_t offset;
} graphdist_t;


typedef struct graph_t {
  /* global counts */
  vtx_t nvtxs; 
  vtx_t gnvtxs;
  adj_t nedges;
  /* distribution information */
  graphdist_t dist;
  /* distributed graph structure */
  vtx_t * mynvtxs;
  adj_t * mynedges;
  vtx_t * mymaxdeg;
  adj_t ** xadj;
  wgt_t ** vwgt;
  vtx_t ** adjncy;
  wgt_t ** adjwgt;
  /* coarsening info */
  size_t level;
  vtx_t ** cmap;
  /* partition information */
  wgt_t * pwgts;
  pid_t ** where;
  /* total weight */
  wgt_t tvwgt, tadjwgt;
  real_t invtvwgt;
  /* aliasing */
  vtx_t ** rename;
  vtx_t ** label;
  /* metrics */
  wgt_t mincut;
  vtx_t minvol;
  /* "To free, or not free" */
  int free_xadj, free_vwgt, free_vsize, free_adjncy, free_adjwgt;
  /* graphs in the heirarchy */
  struct graph_t *coarser, *finer;
} graph_t;




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define graph_create __mtmetis_graph_create
/**
 * @brief Allocate and initialize a graph structure.
 *
 * @param nthreads The number of threads the graph will be used by.
 *
 * @return The allocated and initialized graph.
 */
graph_t * graph_create(
    tid_t nthreads);


#define graph_setup __mtmetis_graph_setup
/**
 * @brief Setup a graph structure given this threads parts of the graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param vwgt The vertex weights.
 *
 * @return The setup graph structure.
 */
graph_t * graph_setup(
    vtx_t nvtxs, 
    adj_t * xadj, 
    vtx_t * adjncy, 
    wgt_t * adjwgt, 
    wgt_t * vwgt);


#define graph_distribute __mtmetis_graph_distribute
/**
 * @brief Distribute a csr based graph among threads. 
 *
 * @param dist The type of distribution to use.
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list poitner.
 * @param adjncy The adjacecny list.
 * @param adjwgt The edge weights.
 * @param vwgt The vertex weights.
 * @param nthreads The number of threads to distribute the graph between.
 *
 * @return The distributed graph. 
 */
graph_t * graph_distribute(
    int dist,
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    wgt_t const * vwgt,
    tid_t nthreads);


#define graph_gather __mtmetis_graph_gather
/**
 * @brief Gather a copy of the graph to each thread. Alias is private to each
 * thread where the other three arrays are not.
 *
 * @param graph The distributed graph to gather.
 * @param r_xadj A refernce to the adjacency list pointer.
 * @param r_adjncy A reference to the adjacency list.
 * @param r_vwgt A reference to the vertex weight.
 * @param r_adjwgt A reference to the edge weight.
 * @param r_voff A reference to the offset of the vertices for this thread.
 */
void graph_gather(
  graph_t const * graph,
  adj_t ** r_xadj,
  vtx_t ** r_adjncy,
  wgt_t ** r_vwgt,
  wgt_t ** r_adjwgt,
  vtx_t * r_voff);


#define graph_setup_coarse __mtmetis_graph_setup_coarse
/**
 * @brief Setup a coarse graph given the fine graph and the number of coarse
 * vertices.  
 *
 * @param graph The fine graph.
 * @param cnvtxs The number of coarse vertices.
 *
 * @return The allocated and setup coarse graph. 
 */
graph_t * graph_setup_coarse(
    graph_t * const graph, 
    vtx_t cnvtxs);


#define graph_setup_twgts __mtmetis_graph_setup_twgts
/**
 * @brief Calculate and save the tvwgts of a new graph.
 *
 * @param graph The graph.
 */
void graph_setup_twgts(
    graph_t * graph);


#define graph_alloc_partmemory __mtmetis_graph_alloc_partmemory
/**
 * @brief Allocate memory for partition informatin.
 *
 * @param ctrl The control structure containing nparts.
 * @param graph The graph.
 */
void graph_alloc_partmemory(
    ctrl_t * ctrl,
    graph_t * graph);

#define graph_free __mtmetis_graph_free
/**
 * @brief Free a graph structure and its associated memory.
 *
 * @param graph The graph to free.
 */
void graph_free(
    graph_t * graph);


#define graph_free_rdata __mtmetis_graph_free_rdata
/**
 * @brief Free partition/refinement data associated with a graph.
 *
 * @param graph The graph to free the associated partition/refinement data of.
 */
void graph_free_rdata(
    graph_t * graph);


#define graph_imbalance __mtmetis_graph_imbalance
/**
 * @brief Compute the load imbalance of a partitioning.
 *
 * @param graph The graph.
 * @param nparts The number of partitions.
 * @param pijbm The inverted average partition weight. 
 *
 * @return The imbalance of the partitioning.
 */
double graph_imbalance(
    graph_t const * graph,
    pid_t nparts,
    real_t const * pijbm);


#define graph_imbalance_diff __mtmetis_graph_imbalance_diff
/**
 * @brief Compute the amount the load imbalance of the graph violates the
 * constraint.
 *
 * @param graph The graph.
 * @param nparts The number of partitions.
 * @param pijbm The inverted average partition weight. 
 * @param ubfactor The allowed imbalance constraint.
 *
 * @return The amount of imbalance in excess of the constraint. 
 */
double graph_imbalance_diff(
    graph_t const * const graph,
    pid_t const nparts,
    real_t const * const pijbm,
    real_t const ubfactor);


#define graph_cut __mtmetis_graph_cut
/**
 * @brief Compute the edgecut of a partitioning.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param where The partition labels.
 *
 * @return The total weight of cut edges.
 */
wgt_t graph_cut(
    ctrl_t const * ctrl,
    graph_t const * graph,
    pid_t const * const * where);


#define graph_isbalanced __mtmetis_graph_isbalanced
/**
 * @brief Check if a partitioning of a graph is balanced within the given
 * constraint.
 *
 * @param ctrl The control structure.
 * @param graph The partitioned graph.
 * @param ffactor The balance constraint.
 *
 * @return 1 if the partitioning is balanced.
 */
int graph_isbalanced(
    ctrl_t const * ctrl, 
    graph_t const * graph, 
    real_t ffactor);


#define graph_readjust_memory __mtmetis_graph_readjust_memory
/**
 * @brief Re-adjust the memory used the by edge arrays of a graph. 
 *
 * @param graph The graph structure.
 * @param adjsize The current size of the adjacency arrays.
 */
void graph_readjust_memory(
    graph_t * const graph,
    adj_t adjsize);




#endif


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
  dlthread_comm_t comm;
  graphdist_t dist;
  /* pre partitioning info */
  pid_t ** group;
  pid_t ngroup;
  /* distributed graph structure */
  vtx_t * mynvtxs;
  adj_t * mynedges;
  adj_t ** xadj;
  wgt_t ** vwgt;
  vtx_t ** adjncy;
  wgt_t ** adjwgt;
  /* graph info */
  int uniformvwgt;
  int uniformadjwgt;
  vtx_t * nislands;
  /* coarsening info */
  size_t level;
  vtx_t ** cmap;
  /* partition information */
  wgt_t * pwgts;
  pid_t ** where;
  struct kwinfo_t * kwinfo;
  struct esinfo_t * esinfo;
  struct vsinfo_t * vsinfo;
  /* total weight */
  twgt_t tvwgt, tadjwgt;
  real_t invtvwgt;
  /* aliasing */
  vtx_t ** rename;
  vtx_t ** label;
  /* metrics */
  wgt_t mincut, minsep;
  vtx_t minvol;
  /* "To free, or not free" */
  int free_xadj, free_vwgt, free_vsize, free_adjncy, free_adjwgt;
  /* graphs in the heirarchy */
  struct graph_t *coarser, *finer;
} graph_t;




/******************************************************************************
* SERIAL FUNCTION PROTOTYPES **************************************************
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
 * @param nthreads The number of threads the graph is distributed across. 
 *
 * @return The setup graph structure.
 */
graph_t * graph_setup(
    vtx_t * nvtxs, 
    adj_t ** xadj, 
    vtx_t ** adjncy, 
    wgt_t ** adjwgt, 
    wgt_t ** vwgt,
    tid_t nthreads);


#define graph_distribute __mtmetis_graph_distribute
/**
 * @brief Distribute a csr based graph among threads. 
 *
 * @param dist The type of distribution to use.
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list poitner.
 * @param adjncy The adjacecny list.
 * @param vwgt The vertex weights.
 * @param adjwgt The edge weights.
 * @param nthreads The number of threads to distribute the graph between.
 *
 * @return The distributed graph. 
 */
graph_t * graph_distribute(
    int dist,
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * vwgt,
    wgt_t const * adjwgt, 
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
  vtx_t ** r_voff);


#define graph_setup_coarse __mtmetis_graph_setup_coarse
/**
 * @brief Setup a coarse graph given the fine graph and the number of coarse
 * vertices.  
 *
 * @param graph The fine graph.
 * @param cnvtxs The number of coarse vertices for each thread.
 *
 * @return The allocated and setup coarse graph. 
 */
graph_t * graph_setup_coarse(
    graph_t * graph, 
    vtx_t * cnvtxs);


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
 * @param graph The graph structure.
 * @param where The partition labels.
 *
 * @return The total weight of cut edges.
 */
wgt_t graph_cut(
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


#define graph_extract_halves __mtmetis_graph_extract_halves
/**
 * @brief Pull out partition 0 and partition 1 from a vertex separated graph or
 * edge bisected graph. Vertices with high partition labels than 1 are dropped.
 * This sets the label[][] array of the new graphs to be global vertex numbers
 * in the original graph, and sets the rename[][] array in the original graph
 * to be global vertex numbers in the new graph.
 *
 * @param graph The graph to extract subgraphs from.
 * @param where The partition ID's of each vertex.
 * @param halves The two extracted subgraphs (output).
 */
void graph_extract_halves(
    graph_t * graph,
    pid_t const * const * where,
    graph_t ** halves);


#define graph_size __mtmetis_graph_size
/**
 * @brief Determine the amount of memory required to store the graph.
 *
 * @param graph The graph to calculate the size of.
 *
 * @return The number of bytes required to store the graph.
 */
size_t graph_size(
    graph_t const * graph);


#define ser_graph_extract_halves __mtmetis_ser_graph_extract_halves
/**
 * @brief Extract two subgraph from partitions 0 and 1.
 *
 * @param graph The graph to extract subgraphs from.
 * @param gwhere The partition IDs of each vertex in the graph.
 * @param halves The references to each subgraph (output).
 */
void ser_graph_extract_halves(
    graph_t * graph,
    pid_t const * const * gwhere,
    graph_t ** halves);


#define graph_calc_dist __mtmetis_graph_calc_dist
/**
 * @brief Configure the distribution structure.
 *
 * @param maxnvtxs The maximum number of vertices owned by a thread.
 * @param nthreads The number of threads.
 * @param dist The distribution structure to configure.
 */
void graph_calc_dist(
    vtx_t maxnvtxs, 
    tid_t nthreads,
    graphdist_t * dist);




/******************************************************************************
* PARALLEL FUNCTION PROTOTYPES ************************************************
******************************************************************************/


#define par_graph_create __mtmetis_par_graph_create
/**
 * @brief Allocate and initialize a graph structure.
 *
 * @param nthreads The thread communicator.
 *
 * @return The allocated and initialized graph.
 */
graph_t * par_graph_create(
    dlthread_comm_t comm);


#define par_graph_setup __mtmetis_par_graph_setup
/**
 * @brief Setup a graph structure given this threads parts of the graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param vwgt The vertex weights.
 * @param comm The thread communicator.
 *
 * @return The setup graph structure.
 */
graph_t * par_graph_setup(
    vtx_t nvtxs, 
    adj_t * xadj, 
    vtx_t * adjncy, 
    wgt_t * adjwgt, 
    wgt_t * vwgt,
    dlthread_comm_t comm);


#define par_graph_gather __mtmetis_par_graph_gather
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
void par_graph_gather(
  graph_t const * graph,
  adj_t ** r_xadj,
  vtx_t ** r_adjncy,
  wgt_t ** r_vwgt,
  wgt_t ** r_adjwgt,
  vtx_t * r_voff);


#define par_graph_shuffle __mtmetis_par_graph_shuffle
/**
 * @brief Shuffle the vertices in a graph such that the partition matches the
 * owning thread. 
 *
 * @param ctrl The control structure.
 * @param graph The graph to shuffle.
 * @param where The partition (destination thread) ids of each vertex.
 * @param wgts Preserve the current weight information.
 */
void par_graph_shuffle(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t const * const * where,
    int wgts);


#define par_graph_setup_coarse __mtmetis_par_graph_setup_coarse
/**
 * @brief Setup a coarse graph given the fine graph and the number of coarse
 * vertices.  
 *
 * @param graph The fine graph.
 * @param cnvtxs The number of coarse vertices.
 *
 * @return The allocated and setup coarse graph. 
 */
graph_t * par_graph_setup_coarse(
    graph_t * const graph, 
    vtx_t cnvtxs);


#define par_graph_setup_twgts __mtmetis_par_graph_setup_twgts
/**
 * @brief Calculate and save the twgts of a new graph.
 *
 * @param graph The graph.
 */
void par_graph_setup_twgts(
    graph_t * graph);


#define par_graph_alloc_partmemory __mtmetis_par_graph_alloc_partmemory
/**
 * @brief Allocate memory for partition informatin.
 *
 * @param ctrl The control structure containing nparts.
 * @param graph The graph.
 */
void par_graph_alloc_partmemory(
    ctrl_t * ctrl,
    graph_t * graph);


#define par_graph_free __mtmetis_par_graph_free
/**
 * @brief Free a graph structure and its associated memory.
 *
 * @param graph The graph to free.
 */
void par_graph_free(
    graph_t * graph);


#define par_graph_free_rdata __mtmetis_par_graph_free_rdata
/**
 * @brief Free partition/refinement data associated with a graph.
 *
 * @param graph The graph to free the associated partition/refinement data of.
 */
void par_graph_free_rdata(
    graph_t * graph);


#define par_graph_readjust_memory __mtmetis_par_graph_readjust_memory
/**
 * @brief Re-adjust the memory used the by edge arrays of a graph. 
 *
 * @param graph The graph structure.
 * @param adjsize The current size of the adjacency arrays.
 */
void par_graph_readjust_memory(
    graph_t * const graph,
    adj_t adjsize);


#define par_graph_extract_halves __mtmetis_par_graph_extract_halves
/**
 * @brief Pull out partition 0 and partition 1 from a vertex separated graph or
 * edge bisected graph. Vertices with high partition labels than 1 are dropped.
 * This sets the label[][] array of the new graphs to be global vertex numbers
 * in the original graph, and sets the rename[][] array in the original graph
 * to be global vertex numbers in the new graph.
 *
 * @param graph The graph to extract subgraphs from.
 * @param where The partition ID's of each vertex.
 * @param halves The two extracted subgraphs (output).
 *
 * @return Which half of the graph this thread is assigned.
 */
tid_t par_graph_extract_halves(
    graph_t * graph,
    pid_t const * const * where,
    graph_t ** halves);


#define par_graph_extract_boundary __mtmetis_par_graph_extract_boundary
/**
 * @brief Extract a subgraph exposing the boundary of a bisection, where the
 * size of the boundary is determined by the maximum imbalanced allowed in the
 * partitioning parameters.
 *
 * Each thread will have two vertices represent the core of each partition, at
 * indices 0 and 1.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * @param bnd The set of vertices on the immediate boundary.
 *
 * @return The extract boundary subgraph.
 */
graph_t * par_graph_extract_boundary(
    ctrl_t const * ctrl,
    graph_t const * graph,
    vtx_iset_t const * bnd);


#define par_graph_extract_separator __mtmetis_par_graph_extract_separator
/**
 * @brief Extract a subgraph consisting only the current separator and two
 * super vertices representing each half.
 *
 * Each thread will have two vertices represent the core of each partition, at
 * indices 0 and 1.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * @param bnd The set of vertices on the immediate boundary.
 *
 * @return The extract boundary subgraph.
 */
graph_t * par_graph_extract_separator(
    ctrl_t const * ctrl,
    graph_t const * graph,
    vtx_iset_t const * bnd);


#define par_graph_extract_aseparator __mtmetis_par_graph_extract_aseparator
/**
 * @brief Extract a subgraph consisting only the current separator and two
 * super vertices representing each half.
 *
 * Each thread will have two vertices represent the core of each partition, at
 * indices 0 and 1.
 *
 * @param ctrl The control structure.
 * @param graph The graph.
 * @param bnd The set of vertices on the immediate boundary.
 *
 * @return The extract boundary subgraph.
 */
graph_t * par_graph_extract_aseparator(
    ctrl_t const * ctrl,
    graph_t const * graph,
    vtx_iset_t const * bnd);


#define par_graph_build_radj __mtmetis_par_graph_build_radj
/**
 * @brief Build a reverse adjacency index and store it in graph->radj.
 *
 * @param graph The graph to build the reverse adjacency list of.
 *
 * @return The reverse adjacecny index.
 */
adj_t * par_graph_build_radj(
    graph_t const * graph);


#define par_graph_contract __mtmetis_par_graph_contract
/**
* @brief Create a coarse graph given a matching using hash table to identify
*   edges to be merged
*
* @param ctrl The control structure.
* @param graph The fine graph.
* @param cnvtxs The number of coarse vertices in the coarse graph.
* @param match The matchings of vertices (match[match[v]] = v).
* @param fcmap The mapping of the coarse vertex number to the lowest number
*   fine vertex in the coarse vertex.
*/
void par_graph_contract(
    ctrl_t * ctrl, 
    graph_t * graph, 
    vtx_t const cnvtxs, 
    vtx_t const * const * gmatch, 
    vtx_t const * fcmap);


#define par_graph_intext_vtx __mtmetis_par_graph_intext_vtx
/**
 * @brief Count the number of internal and external (interface) vertices owned
 * by this thread.
 *
 * @param graph The graph.
 * @param r_nint A reference to the number of internal vertices.
 * @param r_next A reference to the number of external vertices.
 */
void par_graph_intext_vtx(
    graph_t const * const graph,
    vtx_t * const r_nint,
    vtx_t * const r_next);


#define par_graph_cut __mtmetis_par_graph_cut
/**
 * @brief Compute the edgecut of a partitioning in parallel.
 *
 * @param graph The graph structure.
 * @param where The partition labels.
 *
 * @return The total weight of cut edges.
 */
wgt_t par_graph_cut(
    graph_t const * graph,
    pid_t const * const * where);


#define par_graph_removeislands __mtmetis_par_graph_removeislands
/**
 * @brief Remove island vertices from a graph and adjust parameters.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 */
void par_graph_removeislands(
    ctrl_t * ctrl,
    graph_t * graph);


#define par_graph_restoreislands __mtmetis_par_graph_restoreislands
/**
 * @brief Restore island vertices and parameters to a graph.
 *
 * @param ctrl The control structure.
 * @param graph The graph structure.
 * @param gwhere The where vector (should be large enough to hold island
 * vertices).
 */
void par_graph_restoreislands(
    ctrl_t * ctrl,
    graph_t * graph,
    pid_t * const * gwhere);


#define par_graph_extract_parts __mtmetis_par_graph_extract_parts
/**
 * @brief Extract the partitions from a graph. If called with more threads than
 * partitions, the number of threads should be a multiple of number of
 * partitions. The supplied 'nparts' can be lower than the actual number of
 * partitions, an donly partitions with IDs lower than 'nparts' will be
 * extracted.
 *
 * @param graph The graph to extract partitions from.
 * @param gwhere The partition id of each vertex.
 * @param nparts The number of partitions to extract.
 * @param parts The extracted partitions (output), should be of lenght nparts.
 */
void par_graph_extract_parts(
    graph_t * const graph,
    pid_t const * const * const gwhere,
    pid_t const nparts,
    graph_t ** const parts);




#endif

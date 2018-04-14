/**
 * @file analytics.h
 * @brief Functions for measuring and estimating characteristics of graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-10-09
 */




#ifndef BOWSTRING_ANALYZE_H
#define BOWSTRING_ANALYZE_H




#include "base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define degree_distribution __bowstring_degree_distribution
/**
 * @brief Calculate the degree distribution of a graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param r_degrees A reference to the array containing the degree counts
 *     (output)
 * @param r_maxdegree A reference to the maximum degree of the graph (output)
 *
 * @return 
 */
void degree_distribution(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t ** r_degrees, 
    vtx_t * r_maxdegree);


#define nhop_degree_distribution __bowstring_nhop_degree_distribution
/**
 * @brief Calculate the distribution of the number of vertices reachable in
 * nhops.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param nhops The number of hops to search.
 * @param r_degrees A reference to the array containing the degree counts
 *     (output)
 * @param maxdegree A reference to the maximum degree of the graph (output)
 */
void nhop_degree_distribution(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    vtx_t nhops, 
    vtx_t ** r_degrees, 
    vtx_t * r_maxdegree);


#define count_triangles __bowstring_count_triangles
/**
 * @brief Count the number of triangles in a graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param radj The reverse adjacency index.
 *
 * @return The number of triangles in the graph.
 */
size_t count_triangles(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy,
    adj_t * radj);


#define atomic_cycle_distribution __bowstring_atomic_cycle_distribution
/**
 * @brief Calculate the distribution of atomic cycles.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param radj The reverse adjacency index.
 * @param r_cycles A reference to the cycle distribution array.
 * @param r_maxcycle A reference to maximum size of the cycle.
 */
void atomic_cycle_distribution(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy,
    adj_t * radj, 
    size_t ** r_cycles, 
    vtx_t * r_maxcycle);


#define star_distribution __bowstring_star_distribution
/**
 * @brief Calculate the distribution of stars in the graph. The size of a star
 * is determined by the number of degree one vertices that it is connected to.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacnecy list.
 * @param r_stars A refernce to the array containing the star distribution.
 * @param r_maxstar A referce to the maximum size of a star.
 */
void star_distribution(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy,
    vtx_t ** r_stars, 
    vtx_t * r_maxstar);


#define calc_domaindegree __bowstring_calc_domaindegree
/**
 * @brief Calcuate the total degree of each domain/partition/cluster.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param nparts The number of domains.
 * @param where The array containing the domain label for each vertex.
 * @param dd The degree of each domain (output, should be of length nparts).
 */
void calc_domaindegree(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    vlbl_t nparts, 
    vlbl_t const * where, 
    wgt_t * dd);


#define calc_domainconn __bowstring_calc_domainconn
/**
 * @brief Calcuate the connectivity of each domain/partition/cluster.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param nparts The number of domains.
 * @param where The array containing the domain label for each vertex.
 * @param dd The connectivity of each domain (output, should be of length 
 *   nparts).
 */
void calc_domainconn(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    vlbl_t nparts, 
    vlbl_t const * where, 
    int * dc);


#define calc_domaincomvol __bowstring_calc_domaincomvol
/**
 * @brief Calcuate the total communication volume of each 
 *    domain/partition/cluster.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param nparts The number of domains.
 * @param where The array containing the domain label for each vertex.
 * @param dcvo The domain communication volume out (output array of length
 *     nparts)
 * @param dcvi The domain communication volumn in (output array of length
 *     nparts)
 */
void calc_domaincomvol(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * vwgt, 
    vlbl_t nparts, 
    vlbl_t const * where, 
    wgt_t * dcvo,
    wgt_t * dcvi);


#define calc_edgecut __bowstring_calc_edgecut
/**
 * @brief Calculate the edge cut of a partitioning.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weight (can be NULL).
 * @param where The array containing the partition id for each vertex.
 *
 * @return The total weight of cut edges.
 */
wgt_t calc_edgecut(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    vlbl_t const * where);


#define calc_communicationvolume __bowstring_calc_communicationvolume
/**
 * @brief Calculate the total communication volume of a partitioning.
 *
 * @param nvtxs The number of vertices in a graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The weight/size of each vertex.
 * @param nparts The number of partitions.
 * @param where The array containing the partition id of each vertex.
 *
 * @return Total communication volume.
 */
wgt_t calc_communicationvolume(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * vwgt, 
    vlbl_t nparts,
    vlbl_t const * where);


#define calc_max_domaincomvol __bowstring_calc_max_domaincomvol
/**
 * @brief Calculate the maximum in or out comminucation volume for a
 * partionting in a partitioning.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param nparts The number of partitions.
 * @param where The array containing the partition id of each vertex.
 *
 * @return The maximum communication volume of a partition.
 */
wgt_t calc_max_domaincomvol(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * vwgt, 
    vlbl_t nparts,
    vlbl_t const * where);


#define calc_max_domaindegree __bowstring_calc_max_domaindegree
/**
 * @brief Calculate the maximum degree of a partition in a partitioning.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param nparts The number of partitions.
 * @param where The array containing the partition id of each vertex.
 *
 * @return The maximum degree of a partition.
 */
wgt_t calc_max_domaindegree(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    vlbl_t nparts,
    vlbl_t const * where);


#define calc_vertex_balance __bowstring_calc_vertex_balance
/**
 * @brief Calculate the balance (imbalance) of vertices among the partitions.
 * That is, the maximum weight of a partition divided by the average weight of
 * a partition. If 1.0 is returned, it means the vertices are perfectly 
 * balanced, and if nparts is returned, it means that all of the vertices are
 * assigned to the same partition.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency list pointer.
 * @param adjncy The adjacency list.
 * @param vwgt The vertex weights.
 * @param nparts The number of partitions in the graph.
 * @param where The partition id of each vertex.
 *
 * @return The vertex balance (imbalance).
 */
double calc_vertex_balance(
    vtx_t nvtxs,
    adj_t const * xadj,
    vtx_t const * adjncy,
    wgt_t const * vwgt,
    vlbl_t nparts,
    vlbl_t const * where);


#define calc_edge_balance __bowstring_calc_edge_balance
/**
 * @brief Calculate the edge (im)balance given a partitioning/clustering of a
 * graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacnecy list poitner.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param nparts The number partitions/clusters.
 * @param where The partition/cluster id of each vertex.
 *
 * @return The imbalance (1.0 - nparts) 
 */
double calc_edge_balance(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    vlbl_t nparts,
    vlbl_t const * where);


#define calc_degree_balance __bowstring_calc_degree_balance
/**
 * @brief Calcuate the (im)balance of external facing edges given a
 * partitioning/clustering of a graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacnecy list poitner.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param nparts The number partitions/clusters.
 * @param where The partition/cluster id of each vertex.
 *
 * @return The imbalance (1.0 - nparts) 
 */
double calc_degree_balance(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    vlbl_t nparts,
    vlbl_t const * where);


#define calc_comvol_balance __bowstring_calc_comvol_balance
/**
 * @brief Calcuate the (im)balance of in and out communication volumes given a
 * partitioning/clustering of a graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacnecy list poitner.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param nparts The number partitions/clusters.
 * @param where The partition/cluster id of each vertex.
 *
 * @return The imbalance (1.0 - nparts) 
 */
double calc_comvol_balance(
    vtx_t nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    vlbl_t nparts,
    vlbl_t const * where);


#define calc_modularity __bowstring_calc_modularity
/**
 * @brief Calcuate the modularity of a clustering/partitioning of a graph.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacnecy list poitner.
 * @param adjncy The adjacency list.
 * @param adjwgt The edge weights.
 * @param nparts The number partitions/clusters.
 * @param where The partition/cluster id of each vertex.
 *
 * @return The modularity (-0.5 - 1.0)
 */
double calc_modularity(
    vtx_t const nvtxs, 
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    wgt_t const * adjwgt, 
    vlbl_t nparts,
    vlbl_t const * where);


#define calc_partition_components __bowstring_calc_partition_components
/**
 * @brief Calculate the average number of connected components per partition.
 *
 * @param nvtxs The number of vertices in the graph.
 * @param xadj The adjacency pointer.
 * @param adjncy The adjacnecy list.
 * @param nparts The number of partitions.
 * @param where The vertex partition ids.
 *
 * @return The average number of connected components per partition. 
 */
double calc_partition_components(
    vtx_t nvtxs,
    adj_t const * xadj, 
    vtx_t const * adjncy, 
    vlbl_t nparts,
    vlbl_t const * where);



#endif

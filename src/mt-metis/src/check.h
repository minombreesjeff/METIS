/**
 * @file check.h
 * @brief Function prototypes for sanity checks.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2014-09-16
 */




#ifndef MTMETIS_CHECK_H
#define MTMETIS_CHECK_H




#include "base.h"
#include "graph.h"
#include "ucinfo.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define check_info __mtmetis_check_info
/**
 * @brief Perform a sanity check on refinement information.
 *
 * @param ucinfo The ucinfo to check.
 * @param graph The graph structure.
 * @param where The partition id of each vertex.
 *
 * @return 1 if the information is sane.
 */
int check_info(
    ucinfo_t const * ucinfo,
    graph_t const * graph,
    pid_t const * const * where);


#define check_graph __mtmetis_check_graph
/**
 * @brief Check the sanity of a graph structure.
 *
 * @param graph The graph structure.
 *
 * @return 1 if the graph is sane.
 */
int check_graph(
    graph_t const * graph);


#define check_bnd __mtmetis_check_bnd
/**
 * @brief Check the sanity of a boundary.
 *
 * @param graph The graph to check the boundary of.
 *
 * @return 1 if the boundary is sane.
 */
int check_bnd(
    vtx_iset_t const * bnd,
    graph_t const * graph);




#endif

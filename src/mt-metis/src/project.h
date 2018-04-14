/**
 * @file project.c
 * @brief Projection functions. 
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2015, Regents of the University of Minnesota
 * @version 1
 * @date 2015-05-20
 */




#ifndef PROJECT_H
#define PROJECT_H




#include "kwinfo.h"
#include "vsinfo.h"
#include "esinfo.h"
#include "check.h"





/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


#define par_project_graph __mtmetis_par_project_graph
/**
 * @brief Project a partitioning from a coarse graph to a fine graph.
 * 
 * @param ctrl The control strucutre with partitioning parameters.
 * @param graph The graph to project the partitioning to, from graph->coarser.
 */
void par_project_graph(
    ctrl_t * ctrl,
    graph_t * graph);




#endif



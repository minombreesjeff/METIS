/**
 * @file bowstring.c
 * @brief Top level driver functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-07-24
 */




#ifndef BOWSTRING_C
#define BOWSTRING_C




#include "base.h"
#include "tree.h"
#include "order.h"
#include "graph.h"
#include "sparsen.h"
#include "io/io.h"




/******************************************************************************
* PRIVATE CONSTANTS ***********************************************************
******************************************************************************/


static const char * FORMAT_SUFFIXES[] = {
  [BOWSTRING_FORMAT_METIS] = ".graph",
  [BOWSTRING_FORMAT_SNAP] = ".snap",
  [BOWSTRING_FORMAT_DIMACS] = ".dimacs",
  [BOWSTRING_FORMAT_CSR] = ".csr",
  [BOWSTRING_FORMAT_PEGASUS] = ".pegasus",
  [BOWSTRING_FORMAT_CLOUD9] = ".cloud9",
  [BOWSTRING_FORMAT_STP] = ".stp",
  [BOWSTRING_FORMAT_NERSTRAND] = ".graph",
  [BOWSTRING_FORMAT_NBG] = ".nbg"
};
static const size_t NSUFFIXES = sizeof(FORMAT_SUFFIXES)/sizeof(char*);




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int bowstring_read_graph(
    const char * const filename, 
    int graphtype, 
    vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  size_t i;

  if (graphtype == BOWSTRING_FORMAT_AUTO) {
    for (i=0;i<NSUFFIXES;++i) {
      if (dl_string_endswith(filename,FORMAT_SUFFIXES[i])) {
        graphtype = (bowstring_graph_type_t)i; 
        break;
      }
    }
    if (i==NSUFFIXES) {
      eprintf("Unable to determine filetype: '%s'\n",filename);
      return BOWSTRING_ERROR_INVALIDINPUT;
    }
  }

  return io_read_graph(filename,graphtype,r_nvtxs,r_xadj,r_adjncy,r_vwgt,
      r_adjwgt);
}


int bowstring_write_graph(
    const char * const filename, 
    int graphtype, 
    vtx_t const nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const vwgt, 
    const wgt_t * const adjwgt)
{
  size_t i;

  if (graphtype == BOWSTRING_FORMAT_AUTO) {
    for (i=0;i<NSUFFIXES;++i) {
      if (dl_string_endswith(filename,FORMAT_SUFFIXES[i])) {
        graphtype = (bowstring_graph_type_t)i; 
        break;
      }
    }
    if (i==NSUFFIXES) {
      eprintf("Unable to determine filetype: '%s'\n",filename);
      return BOWSTRING_ERROR_INVALIDINPUT;
    }
  }

  return io_write_graph(filename,graphtype,nvtxs,xadj,adjncy,vwgt,adjwgt);
}


int bowstring_write_vertex_labels(
    char const * filename, 
    vtx_t nvtxs, 
    vlbl_t const * labels)
{
  return write_vertex_labels(filename,nvtxs,labels);
}


void bowstring_build_tree(
    const int treetype, 
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt,
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_adjwgt,
    int ** const r_adjmask)
{
  int * adjmask;

  adjmask = int_alloc(xadj[nvtxs]);

  switch (treetype) {
    case BOWSTRING_TREE_MST:
      build_mst_tree(nvtxs,xadj,adjncy,adjwgt,adjmask);
      break;
    case BOWSTRING_TREE_RST:
      build_rst_tree(nvtxs,xadj,adjncy,adjmask);
      break;
    case BOWSTRING_TREE_DFS:
      build_dfs_tree(nvtxs,xadj,adjncy,vtx_rand(0,nvtxs),NULL,NULL,NULL,
          adjmask);
    case BOWSTRING_TREE_BFS:
      build_bfs_tree(nvtxs,xadj,adjncy,vtx_rand(0,nvtxs),NULL,NULL,NULL,
          adjmask);
    default:
      dl_error("Unknown tree type '%d'\n",treetype);
  }

  apply_edge_mask(nvtxs,xadj,adjncy,adjwgt,adjmask,r_xadj,r_adjncy,
      r_adjwgt);

  if (r_adjmask) {
    *r_adjmask = adjmask;
  } else {
    dl_free(adjmask);
  }
}


void bowstring_build_adjncy_index(
    vtx_t const nvtxs, 
    adj_t const * const xadj,
    vtx_t const * const adjncy, 
    adj_t * const radj)
{
  build_adjncy_index(nvtxs,xadj,adjncy,radj);
}


void bowstring_remove_edges(
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt, 
    const int type,
    const double frac, 
    adj_t ** const r_xadj,
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_adjwgt)
{
  elbl_t maxrank;
  elbl_t * ranks;

  ranks = elbl_calloc(xadj[nvtxs]);

  maxrank = 0;
  switch (type) {
    case BOWSTRING_EDGERANK_NI:
      maxrank = build_nirank(nvtxs,xadj,adjncy,adjwgt,ranks);
      break;
    case BOWSTRING_EDGERANK_MST:
      maxrank = build_mstrank(nvtxs,xadj,adjncy,adjwgt,ranks);
      break;
    case BOWSTRING_EDGERANK_AST:
      maxrank = build_astrank(nvtxs,xadj,adjncy,adjwgt,ranks);
      break;
    case BOWSTRING_EDGERANK_LST:
      maxrank = build_lstrank(nvtxs,xadj,adjncy,adjwgt,ranks);
      break;
    default:
      dl_error("Unknown edge ranking %d\n",type);
  }

  prune_ranked_edges(nvtxs,xadj,adjncy,adjwgt,ranks,maxrank,frac,r_xadj,
      r_adjncy,r_adjwgt,BOWSTRING_REWEIGHT_NONE);

  dl_free(ranks);
}


void bowstring_induce_subgraph(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt,
    wgt_t const * const vwgt,
    int const * const present,
    vtx_t * const r_nvtxs,
    adj_t ** const r_xadj,
    vtx_t ** const r_adjncy,
    wgt_t ** const r_adjwgt,
    wgt_t ** const r_vwgt,
    vtx_t ** const r_alias,
    vtx_t ** const r_rename)
{
  induce_subgraph(nvtxs,xadj,adjncy,adjwgt,vwgt,present,r_nvtxs,r_xadj, \
      r_adjncy,r_adjwgt,r_vwgt,r_alias,r_rename);
}


void bowstring_induce_subgraphs(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt,
    wgt_t const * const vwgt,
    vlbl_t const * const part,
    vlbl_t const npart,
    vtx_t * const xnvtxs,
    adj_t ** const xxadj,
    vtx_t ** const xadjncy,
    wgt_t ** const xadjwgt,
    wgt_t ** const xvwgt,
    vtx_t ** const xalias,
    vtx_t ** const r_rename)
{
  induce_subgraphs(nvtxs,xadj,adjncy,adjwgt,vwgt,part,npart,xnvtxs,xxadj, \
      xadjncy,xadjwgt,xvwgt,xalias,r_rename);
}



void bowstring_extract_subgraph(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt,
    wgt_t const * const vwgt,
    int const * const vpresent,
    int const * const epresent,
    vtx_t * const r_nvtxs,
    adj_t ** const r_xadj,
    vtx_t ** const r_adjncy,
    wgt_t ** const r_adjwgt,
    wgt_t ** const r_vwgt,
    vtx_t ** const r_alias,
    vtx_t ** const r_rename)
{
  extract_subgraph(nvtxs,xadj,adjncy,adjwgt,vwgt,vpresent,epresent,r_nvtxs, \
      r_xadj,r_adjncy,r_adjwgt,r_vwgt,r_alias,r_rename);
}


void bowstring_label_components(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    vlbl_t ** r_lbl,
    vlbl_t * r_nlbl)
{
  label_components(nvtxs,xadj,adjncy,r_lbl,r_nlbl);
}


void bowstring_label_partition_components(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    vlbl_t const * const where,
    vlbl_t ** r_lbl,
    vlbl_t * r_nlbl)
{
  label_partition_components(nvtxs,xadj,adjncy,where,r_lbl,r_nlbl);
}


int bowstring_check_graph(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt)
{
  return check_graph(nvtxs,xadj,adjncy,adjwgt);
}

void bowstring_voronoi_regions(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt,
    vtx_t const * const sources,
    vtx_t const nsources,
    vtx_t * const where,
    vtx_t * const parent,
    wgt_t * const dist)
{
  voronoi_regions(nvtxs,xadj,adjncy,adjwgt,sources,nsources,where,parent,dist);
}


void bowstring_voronoi_diagram(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt,
    vtx_t const nsources,
    vtx_t const * const where,
    wgt_t const * const dist,
    adj_t ** const r_vxadj,
    vtx_t ** const r_vadjncy,
    wgt_t ** const r_vadjwgt,
    vtx_t ** const r_adjsrc,
    adj_t ** const r_adjori)
{
  voronoi_diagram(nvtxs,xadj,adjncy,adjwgt,nsources,where,dist,r_vxadj, \
      r_vadjncy,r_vadjwgt,r_adjsrc,r_adjori);
}


void bowstring_order_graph(
    vtx_t const nvtxs, 
    adj_t * const xadj, 
    vtx_t * const adjncy,
    wgt_t * const vwgt, 
    wgt_t * const adjwgt, 
    vtx_t const * const perm) 
{
  order_graph(nvtxs,xadj,adjncy,vwgt,adjwgt,perm);
}


void bowstring_permutation(
    int const ordering,
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy,
    wgt_t const * const vwgt,
    wgt_t const * const adjwgt,
    vtx_t * const perm)
{
  order_permutation(ordering,nvtxs,xadj,adjncy,vwgt,adjwgt,perm);
}


#endif

/**
 * @file graph.c
 * @brief Functions for performing misc graph operations
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-08-07
 */


#ifndef GRAPH_C
#define GRAPH_C


#include "graph.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const vtx_t UNASSIGNED = (vtx_t)-1;
static const vtx_t QUEUED = (vtx_t)-2;




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLDJSET_PREFIX vtx
#define DLDJSET_TYPE_T vtx_t
#define DLDJSET_STATIC
#include "dldjset_headers.h"
#undef DLDJSET_STATIC
#undef DLDJSET_TYPE_T
#undef DLDJSET_PREFIX


#define DLPQ_PREFIX vw
#define DLPQ_KEY_T wgt_t
#define DLPQ_VAL_T vtx_t
#define DLPQ_MIN
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_MIN
#undef DLPQ_VAL_T
#undef DLPQ_KEY_T
#undef DLPQ_PREFIX


#define DLISET_PREFIX vtx
#define DLISET_TYPE_T vtx_t
#define DLISET_STATIC
#include "dliset_headers.h"
#undef DLISET_STATIC
#undef DLISET_PREFIX
#undef DLISET_TYPE_T


#define DLSORTKV_PREFIX vtx
#define DLSORTKV_KEY_T vtx_t
#define DLSORTKV_VAL_T vtx_t
#define DLSORTKV_STATIC
#include "dlsortkv_headers.h"
#undef DLSORTKV_STATIC
#undef DLSORTKV_VAL_T
#undef DLSORTKV_KEY_T
#undef DLOSRTKV_PREFIX




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Check the validity of a reverse adjacency index.
 *
 * @param nvtxs The number of vertices in teh graph.
 * @param adjncy The adjacency list.
 * @param nedges The number of edges in the graph.
 * @param radj The reverse adjacecny index to check.
 *
 * @return 1 if it is a valid reverse adjacency index.
 */
static int __check_radj(
    vtx_t const nvtxs,
    vtx_t const * const adjncy,
    adj_t const nedges, 
    adj_t const * const radj)
{
  adj_t j;

  for (j=0;j<nedges;++j) {
    if (adjncy[j] < nvtxs) {
      if (radj[radj[j]] != j) {
        eprintf("Bad reverse adjacency list j = "PF_ADJ_T", radj[j] = "
            PF_ADJ_T", radj[radj[j]] = "PF_ADJ_T"\n",j,radj[j],radj[radj[j]]);
        return 0;
      }
    } else {
      if (radj[j] != NULL_ADJ) {
        eprintf("Bad reverse adjacency list j = "PF_ADJ_T", adjncy[j] = "
            PF_VTX_T", radj[j] = "PF_ADJ_T"\n",j,adjncy[j],radj[j]);
        return 0;
      }
    }
  }

  return 1;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void apply_edge_mask(
    vtx_t const nvtxs, 
    adj_t const * const gxadj, 
    vtx_t const * const gadjncy, 
    wgt_t const * const gadjwgt, 
    int const * const adjmask, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_adjwgt) 
{
  vtx_t i;
  adj_t j, tnedges;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL;

  adj_t const nedges = int_sum(adjmask,gxadj[nvtxs]);

  xadj = adj_alloc(nvtxs+1);
  if (gadjncy) {
    adjncy = vtx_alloc(nedges);
  }
  if (gadjwgt) {
    adjwgt = wgt_alloc(nedges); 
  }
  tnedges = xadj[0] = 0;
  if (adjwgt) {
    if (adjncy) {
      for (i=0;i<nvtxs;++i) {
        for (j=gxadj[i];j<gxadj[i+1];++j) {
          if (adjmask[j]) {
            adjncy[tnedges] = gadjncy[j];
            adjwgt[tnedges] = gadjwgt[j];
            ++tnedges;
          }
        }
        xadj[i+1] = tnedges;
      }
    } else {
      for (i=0;i<nvtxs;++i) {
        for (j=gxadj[i];j<gxadj[i+1];++j) {
          if (adjmask[j]) {
            adjwgt[tnedges++] = gadjwgt[j];
          }
        }
        xadj[i+1] = tnedges;
      }
    }
  } else {
    if (adjncy) {
      for (i=0;i<nvtxs;++i) {
        for (j=gxadj[i];j<gxadj[i+1];++j) {
          if (adjmask[j]) {
            adjncy[tnedges++] = gadjncy[j];
          }
        }
        xadj[i+1] = tnedges;
      }
    } else {
      for (i=0;i<nvtxs;++i) {
        for (j=gxadj[i];j<gxadj[i+1];++j) {
          if (adjmask[j]) {
            ++tnedges;
          }
        }
        xadj[i+1] = tnedges;
      }
    }
  }

  if (r_xadj) {
    *r_xadj = xadj;
  } else if (xadj) {
    dl_free(xadj);
  }
  if (r_adjncy) {
    *r_adjncy = adjncy;
  } else if (adjncy) {
    dl_free(adjncy);
  }
  if (r_adjwgt) {
    *r_adjwgt = adjwgt;
  } else if (adjwgt) {
    dl_free(adjwgt);
  }
}


vtx_t mark_boundary_vertices(
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const vlbl_t * const where, 
    int * const bnd)
{
  vtx_t i, k, nbnd;
  adj_t j;
  vlbl_t w;

  nbnd = 0;
  for (i=0;i<nvtxs;++i) {
    w = where[i];
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (w != where[k]) {
        if (bnd) {
          bnd[i] = 1;
        }         
        ++nbnd;
      } else {
        if (bnd) {
          bnd[i] = 0;
        }
      }
    }
  }

  return nbnd;
}


void find_boundary_vertices(
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const vlbl_t * const where, 
    vtx_t * const r_nbnd,
    vtx_t ** const r_bnd)
{
  vtx_t i, k, nbnd;
  adj_t j;
  vlbl_t w;

  vtx_t * bnd = NULL;

  if (r_bnd) {
    bnd = vtx_alloc(nvtxs);
  }

  nbnd = 0;
  for (i=0;i<nvtxs;++i) {
    w = where[i];
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (w != where[k]) {
        if (r_bnd) {
          bnd[nbnd++] = i;
        } else {
          ++nbnd;
        }
      }
    }
  }

  if (r_nbnd) {
    *r_nbnd = nbnd;
  }
  if (r_bnd) {
    *r_bnd = vtx_realloc(bnd,nbnd);
  }
}


int build_adjncy_index(
    const vtx_t nvtxs, 
    const adj_t * const xadj,
    const vtx_t * const adjncy, 
    adj_t * const radj)
{
  vtx_t i, k;
  adj_t j;
  vtx_t * tadj;
  adj_t * txadj, * padj, * tradj, * trans; 

  const adj_t nedges = xadj[nvtxs];

  /* In this loop the lower number vertex on each edge inserts the edge and its
   * reverse index into the radj array. The txadj array is used to track where
   * into that array it should be inserted. The radj array stores a permuted
   * adjncy list for each vertex, and the tradj array stores the reverse index
   * for the tadj array. */ 
  txadj = adj_duplicate(xadj,nvtxs+1);
  tadj = vtx_alloc(nedges); 
  tradj = adj_alloc(nedges);
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (i < k) {
        tadj[txadj[i]] = k;
        tradj[txadj[i]] = txadj[k];
        tadj[txadj[k]] = i;
        tradj[txadj[k]] = txadj[i];
        ++txadj[k];
        ++txadj[i];
      }
    }
  }
  dl_free(txadj);

  /* In this loop the padj array is populated such that it serves as a
   * permutation array from tadj to adjncy. */
  trans = adj_alloc(nvtxs);
  padj = adj_alloc(nedges);
  for (i=0;i<nvtxs;++i) {
    /* write the location of edges in the original graph to trans */
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      trans[k] = j;
    }
    /* using the new ordering of the edges, read in the location of the edges
     * from the original graph. */
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = tadj[j];
      padj[j] = trans[k];
    }
  }
  dl_free(trans);
  dl_free(tadj);

  /* correct radj -- save it to tadj */
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      radj[padj[j]] = padj[tradj[j]];
    }
  }

  dl_free(padj);
  dl_free(tradj);

  DL_ASSERT(__check_radj(nvtxs,adjncy,nedges,radj) == 1,
      "Bad radj generated\n");

  return BOWSTRING_SUCCESS;
}


int build_adjncy_index_rem(
    const vtx_t nvtxs, 
    const adj_t * const xadj,
    const vtx_t * const adjncy, 
    adj_t * const radj)
{
  vtx_t i, k, nr;
  adj_t j, l;
  vtx_t * tadj;
  adj_t * txadj, * padj, * tradj, * trans; 

  const adj_t nedges = xadj[nvtxs];

  /* In this loop the lower number vertex on each edge inserts the edge and its
   * reverse index into the radj array. The txadj array is used to track where
   * into that array it should be inserted. The radj array stores a permuted
   * adjncy list for each vertex, and the tradj array stores the reverse index
   * for the tadj array. */ 
  txadj = adj_duplicate(xadj,nvtxs+1);
  tadj = vtx_alloc(nedges); 
  tradj = adj_alloc(nedges);
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (i < k) {
        if (k < nvtxs) {
          /* handle local edges normally */
          tadj[txadj[i]] = k;
          tradj[txadj[i]] = txadj[k];
          tadj[txadj[k]] = i;
          tradj[txadj[k]] = txadj[i];
          ++txadj[k];
          ++txadj[i];
        } else {
          /* for incomplete graphs */
          tadj[txadj[i]] = k;
          tradj[txadj[i]] = NULL_ADJ;
          ++txadj[i];
        }
      }
    }
  }
  dl_free(txadj);

  /* In this loop the padj array is populated such that it serves as a
   * permutation array from tadj to adjncy. Remote edges have their ordering
   * preserved. */
  trans = adj_alloc(nvtxs+nedges);
  padj = adj_alloc(nedges);
  for (i=0;i<nvtxs;++i) {
    /* write the location of edges in the original graph to trans */
    nr = nvtxs;
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k < nvtxs) {
        trans[k] = j;
      } else {
        trans[nr++] = j;
      }
    }
    /* using the new ordering of the edges, read in the location of the edges
     * from the original graph. */
    nr = nvtxs;
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = tadj[j];
      if (k < nvtxs) {
        padj[j] = trans[k];
      } else {
        padj[j] = trans[nr++];
      }
    }
  }
  dl_free(tadj);
  dl_free(trans);

  /* correct radj -- save it to tadj */
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      l = tradj[j];
      if (l != NULL_ADJ) {
        radj[padj[j]] = padj[l];
      } else {
        radj[padj[j]] = NULL_ADJ;
      }
    }
  }

  dl_free(padj);
  dl_free(tradj);

  DL_ASSERT(__check_radj(nvtxs,adjncy,nedges,radj) == 1,
      "Bad radj generated\n");

  return BOWSTRING_SUCCESS;
}


vtx_t neighborhoodify(
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const vtx_t nbrsize,
    vtx_t * const nbrhd) 
{
  vtx_t i, p, k, nnbrhd, nnbr, nbr, nq; 
  adj_t j;
  vtx_t * perm, * q;
  
  vtx_set(nbrhd,UNASSIGNED,nvtxs);
  q = vtx_alloc(nvtxs);
  perm = vtx_alloc(nvtxs);
  vtx_incset(perm,0,1,nvtxs);
  vtx_pseudo_shuffle(perm,nvtxs/16,nvtxs);

  nnbrhd = 0;
  for (p=0;p<nvtxs;++p) {
    i = perm[p];
    if (nbrhd[i] == UNASSIGNED) {
      nnbr = nq = 0;
      nbr = nnbrhd++;
      q[nq++] = i;
      while (nnbr < nq && nnbr < nbrsize) {
        i = q[nnbr++];
        nbrhd[i] = nbr;
        for (j=xadj[i];j<xadj[i+1];++j) {
          k = adjncy[j];
          if (nbrhd[k] == UNASSIGNED) {
            q[nq++] = k;
            nbrhd[k] = QUEUED;
          }
          if (nq >= nbrsize) {
            break;
          }
        }
      }
      /* clear nodes left in queue */
      while (nnbr < nq) {
        nbrhd[q[nnbr++]] = UNASSIGNED;
      }
    }
  }

  dl_free(perm);

  return nnbrhd;
}


void induce_subgraph(
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
  vtx_t i, k, snvtxs, v;
  adj_t j, snedges;
  adj_t * sxadj;
  vtx_t * sadjncy, * salias, * rename;
  wgt_t * sadjwgt = NULL, * svwgt = NULL;

  int const do_wgt = r_adjwgt && adjwgt;

  /* count number of vertices in subgraph */
  snvtxs = 0;
  for (i=0;i<nvtxs;++i) {
    if (present[i]) {
      ++snvtxs;
    }
  }

  /* allocate subgraph vertex based arrays */
  sxadj = adj_alloc(snvtxs+1);
  salias = vtx_alloc(snvtxs);
  rename = vtx_alloc(nvtxs);

  /* populate subgraphs vertex based arrays */
  snvtxs = 0;
  snedges = sxadj[0] = 0;
  for (i=0;i<nvtxs;++i) {
    if (present[i]) {
      salias[snvtxs] = i;
      rename[i] = snvtxs;
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        if (present[k]) {
          ++snedges;
        }
      }
      sxadj[snvtxs+1] = snedges;
      ++snvtxs;
    }
  }

  /* allocate subgraph edge based arrays */
  sadjncy = vtx_alloc(snedges);
  if (do_wgt) {
    sadjwgt = wgt_alloc(snedges);
  }

  /* populate edge based arrays */
  snedges = 0;
  for (i=0;i<snvtxs;++i) {
    v = salias[i];
    for (j=xadj[v];j<xadj[v+1];++j) {
      k = adjncy[j];
      if (present[k]) {
        sadjncy[snedges] = rename[k];
        if (do_wgt) {
          sadjwgt[snedges] = adjwgt[j];
        }
        ++snedges;
      }
    }
    sxadj[i+1] = snedges;
  }

  /* fill vwgt */
  if (vwgt && *r_vwgt) {
    svwgt = wgt_alloc(snvtxs);
    for (i=0;i<snvtxs;++i) {
      v = salias[i];
      svwgt[i] = vwgt[v];
    }
  }

  *r_nvtxs = snvtxs;
  *r_xadj = sxadj;
  *r_adjncy = sadjncy;
  if (r_adjwgt) {
    *r_adjwgt = sadjwgt;
  }
  if (r_vwgt) {
    *r_vwgt = svwgt;
  }
  if (r_alias) {
    *r_alias = salias;
  } else {
    dl_free(salias);
  }
  if (r_rename) {
    *r_rename = rename;
  } else {
    dl_free(rename);
  }
}


void induce_subgraphs(
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
  vlbl_t p;
  vtx_t i, k, snvtxs, v;
  adj_t j, snedges; 
  vtx_t * rename;
  vtx_t ** lalias;

  int const do_wgt = xadjwgt && adjwgt;

  /* handle null xalias */
  if (xalias) {
    lalias = xalias;
  } else {
    lalias = r_vtx_alloc(npart);
  }

  /* count number of vertices in subgraph */
  vtx_set(xnvtxs,0,npart);
  for (i=0;i<nvtxs;++i) {
    DL_ASSERT(part[i] < npart,"Bad partition id of "PF_VTX_T"/"PF_VTX_T" for "
        "vertex "PF_VTX_T"\n",part[i],npart,i);
    ++xnvtxs[part[i]];
  }

  /* allocate rename vector */
  rename = vtx_alloc(nvtxs);

  /* allocate subgraph vertex based arrays */
  for (p=0;p<npart;++p) {
    snvtxs = xnvtxs[p];
    lalias[p] = vtx_alloc(snvtxs);
    xxadj[p] = adj_alloc(snvtxs+1);
    xxadj[p][0] = 0;
  }

  /* populate subgraphs vertex based arrays */
  vtx_set(xnvtxs,0,npart);
  for (i=0;i<nvtxs;++i) {
    p = part[i];
    snvtxs = xnvtxs[p];
    lalias[p][snvtxs] = i;
    rename[i] = snvtxs;
    snedges = xxadj[p][snvtxs];
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (part[k] == p) {
        ++snedges;
      }
    }
    xxadj[p][snvtxs+1] = snedges;
    ++xnvtxs[p];
  }

  /* allocate subgraph edge based arrays */
  for (p=0;p<npart;++p) {
    snvtxs = xnvtxs[p];
    xadjncy[p] = vtx_alloc(xxadj[p][snvtxs]);
    if (xadjwgt) {
      if (adjwgt) {
        xadjwgt[p] = wgt_alloc(xxadj[p][snvtxs]);
      } else {
        xadjwgt[p] = NULL;
      }
    }
  }

  /* populate edge based arrays */
  for (p=0;p<npart;++p) {
    snvtxs = xnvtxs[p];
    for (i=0;i<snvtxs;++i) {
      v = lalias[p][i];
      snedges = xxadj[p][i];
      for (j=xadj[v];j<xadj[v+1];++j) {
        k = adjncy[j];
        if (p == part[k]) {
          xadjncy[p][snedges] = rename[k];
          if (do_wgt) {
            xadjwgt[p][snedges] = adjwgt[j];
          }
          ++snedges;
        }
      }
      xxadj[p][i+1] = snedges;
    }
  }

  /* fill vwgt */
  if (xvwgt) {
    if (vwgt) {
      for (p=0;p<npart;++p) {
        snvtxs = xnvtxs[p];
        xvwgt[p] = wgt_alloc(snvtxs);
        for (i=0;i<snvtxs;++i) {
          v = lalias[p][i];
          xvwgt[p][i] = vwgt[v];
        }
      }
    } else {
      for (p=0;p<npart;++p) {
        xvwgt[p] = NULL;
      }
    }
  }

  if (!xalias) {
    r_vtx_free(lalias,npart);
  }
  if (r_rename) {
    *r_rename = rename;
  } else {
    dl_free(rename);
  }
}


void extract_subgraph(
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
  vtx_t i, k, snvtxs, v;
  adj_t j, snedges;
  adj_t * sxadj;
  vtx_t * sadjncy, * salias, * rename;
  wgt_t * sadjwgt = NULL, * svwgt = NULL;

  int const do_ewgt = r_adjwgt && adjwgt;
  int const do_vwgt = r_vwgt && vwgt;

  /* count number of vertices in subgraph */
  snvtxs = 0;
  for (i=0;i<nvtxs;++i) {
    if (vpresent[i]) {
      ++snvtxs;
    }
  }
  
  /* allocate subgraph vertex based arrays */
  sxadj = adj_alloc(snvtxs+1);
  salias = vtx_alloc(snvtxs);
  rename = vtx_alloc(nvtxs);
  if (do_vwgt) {
    svwgt = wgt_alloc(snvtxs);
  }

  /* populate subgraphs vertex based arrays */
  snvtxs = 0;
  snedges = sxadj[0] = 0;
  for (i=0;i<nvtxs;++i) {
    if (vpresent[i]) {
      if (do_vwgt) {
        svwgt[snvtxs] = vwgt[i];
      }
      salias[snvtxs] = i;
      rename[i] = snvtxs;
      for (j=xadj[i];j<xadj[i+1];++j) {
        if (epresent[j]) {
          k = adjncy[j];
          DL_ASSERT_EQUALS(vpresent[k],1,"%d");
          ++snedges;
        }
      }
      sxadj[snvtxs+1] = snedges;
      ++snvtxs;
    }
  }

  /* allocate subgraph edge based arrays */
  sadjncy = vtx_alloc(snedges);
  if (do_ewgt) {
    sadjwgt = wgt_alloc(snedges);
  }

  /* populate edge based arrays */
  snedges = 0;
  for (i=0;i<snvtxs;++i) {
    v = salias[i];
    for (j=xadj[v];j<xadj[v+1];++j) {
      if (epresent[j]) {
        k = adjncy[j];
        sadjncy[snedges] = rename[k];
        if (do_ewgt) {
          sadjwgt[snedges] = adjwgt[j];
        }
        ++snedges;
      }
    }
    sxadj[i+1] = snedges;
  }

  *r_nvtxs = snvtxs;
  *r_xadj = sxadj;
  *r_adjncy = sadjncy;
  if (r_adjwgt) {
    *r_adjwgt = sadjwgt;
  }
  if (r_vwgt) {
    *r_vwgt = svwgt;
  }
  if (r_alias) {
    *r_alias = salias;
  } else {
    dl_free(salias);
  }
  if (r_rename) {
    *r_rename = rename;
  } else {
    dl_free(rename);
  }
}


void label_components(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    vlbl_t ** const r_lbl,
    vlbl_t * const r_nlbl)
{
  vtx_t i, k, sq, nq, start, l;
  adj_t j;
  vtx_t * q;
  vlbl_t * label;

  q = vtx_alloc(nvtxs);
  label = vlbl_init_alloc(NULL_VTX,nvtxs);

  l = 0;
  start = 0;
  sq = 0;
  nq = 0;
  while (start < nvtxs) {
    if (label[start] == NULL_VTX) {
      q[nq++] = start; 
      label[start] = l;
      while (sq < nq) {
        i = q[sq++];
        for (j=xadj[i];j<xadj[i+1];++j) {
          k = adjncy[j];
          if (label[k] == NULL_VTX) {
            label[k] = l;
            q[nq++] = k;
          }
        }
      }
      ++l;
    }
    ++start;
  }

  dl_free(q);

  if (r_nlbl) {
    *r_nlbl = l;
  }
  if (r_lbl) {
    *r_lbl = label;
  } else {
    dl_free(label);
  }
}


void label_partition_components(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    vlbl_t const * const where,
    vlbl_t ** const r_lbl,
    vlbl_t * const r_nlbl)
{
  vtx_t i, k, sq, nq, start, l;
  adj_t j;
  vlbl_t me;
  vtx_t * q;
  vlbl_t * label;

  q = vtx_alloc(nvtxs);
  label = vlbl_init_alloc(NULL_VTX,nvtxs);

  l = 0;
  start = 0;
  sq = 0;
  nq = 0;
  while (start < nvtxs) {
    if (label[start] == NULL_VTX) {
      q[nq++] = start; 
      label[start] = l;
      while (sq < nq) {
        i = q[sq++];
        me = where[i];
        for (j=xadj[i];j<xadj[i+1];++j) {
          k = adjncy[j];
          if (label[k] == NULL_VTX && where[k] == me) {
            label[k] = l;
            q[nq++] = k;
          }
        }
      }
      ++l;
    }
    ++start;
  }

  dl_free(q);

  if (r_nlbl) {
    *r_nlbl = l;
  }
  if (r_lbl) {
    *r_lbl = label;
  } else {
    dl_free(label);
  }
}


int check_graph(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt)
{
  vtx_t i, k, kk;
  adj_t j, jj;

  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k >= nvtxs) {
        eprintf("Invalid connection to vertex "PF_VTX_T"/"PF_VTX_T" from " \
            "vertex "PF_VTX_T"\n",k,nvtxs,i);
        return 0;
      }
      for (jj=xadj[k];jj<xadj[k+1];++jj) {
        kk = adjncy[jj];
        if (kk == i) {
          if (adjwgt && ! dl_near_equal(adjwgt[j],adjwgt[jj])) {
            eprintf("Unbalanced edge weight "PF_WGT_T":"PF_WGT_T" for the " \
                "edge from vertex "PF_VTX_T" to "PF_VTX_T"\n",adjwgt[j], \
                adjwgt[jj],i,k);
            return 0;
          }
          goto NEXT_EDGE;
        }
      }
      eprintf("Could not find reverse of edge going from "PF_VTX_T" to " \
          PF_VTX_T"\n",i,k);
      return 0;
      NEXT_EDGE:;
    }
  }

  return 1;
}


void voronoi_regions(
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
  vtx_t t, i, k;
  adj_t j;
  wgt_t d;
  vw_pq_t * q;

  DL_ASSERT(check_graph(nvtxs,xadj,adjncy,adjwgt),"Bad graph passed " \
      "into voronoi_regions().");

  q = vw_pq_create(0,nvtxs);

  vtx_set(where,nsources,nvtxs);

  /* create voronoi partitions */
  for (t=0;t<nsources;++t) {
    i = sources[t];
    where[i] = t;
    dist[i] = 0;
    parent[i] = NULL_VTX;
    vw_pq_push(0,i,q);
  }
  if (adjwgt) {
    while (q->size > 0) {
      i = vw_pq_pop(q);
      d = dist[i];
      t = where[i];
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        if (where[k] == nsources) {
          where[k] = t;
          parent[k] = i;
          dist[k] = d + adjwgt[j];
          vw_pq_push(dist[k],k,q);
        } else {
          DL_ASSERT(where[k] < nsources,"Bad location of vertex "PF_VTX_T \
              ", "PF_VTX_T"/"PF_VTX_T"\n",k,where[k],nsources);
          if (d + adjwgt[j] < dist[k]) {
            parent[k] = i;
            where[k] = t;
            dist[k] = d + adjwgt[j];
            vw_pq_update(dist[k],k,q);
          }
        }
      }
    }
  } else {
    while (q->size > 0) {
      i = vw_pq_pop(q);
      d = dist[i];
      t = where[i];
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        if (where[k] == nsources) {
          where[k] = t;
          parent[k] = i;
          dist[k] = d + 1;
          vw_pq_push(dist[k],k,q);
        } else {
          DL_ASSERT(where[k] < nsources,"Bad location of vertex "PF_VTX_T \
              ", "PF_VTX_T"/"PF_VTX_T"\n",k,where[k],nsources);
          if (d + 1 < dist[k]) {
            parent[k] = i;
            where[k] = t;
            dist[k] = d + 1;
            vw_pq_update(dist[k],k,q);
          }
        }
      }
    }
  }

  vw_pq_free(q);
}


void voronoi_add_sources(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const adjwgt,
    vtx_t const offset,
    vtx_t const * const sources,
    vtx_t const nsources,
    vtx_t * const where,
    vtx_t * const parent,
    wgt_t * const dist)
{
  vtx_t t, i, k;
  adj_t j;
  wgt_t d, nd;

  vw_pq_t * q;

  q = vw_pq_create(0,nvtxs);

  /* initialize new sources */
  for (t=0;t<nsources;++t) {
    i = sources[t];
    where[i] = t + offset;
    dist[i] = 0;
    parent[i] = NULL_VTX;
    vw_pq_push(0,i,q);
  }
  /* expand sources */
  if (adjwgt) {
    while (q->size > 0) {
      i = vw_pq_pop(q);
      d = dist[i];
      t = where[i];
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        nd = d + adjwgt[j];
        if (nd < dist[k]) {
          parent[k] = i;
          where[k] = t;
          dist[k] = nd;
          vw_pq_update(nd,k,q);
        }
      }
    }
  } else {
    while (q->size > 0) {
      i = vw_pq_pop(q);
      d = dist[i];
      t = where[i];
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        nd = d + 1;
        if (nd < dist[k]) {
          parent[k] = i;
          where[k] = t;
          dist[k] = nd;
          vw_pq_update(nd,k,q);
        }
      }
    }
  }

  vw_pq_free(q);
}


void voronoi_diagram(
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
  vtx_t i, k, me, p, t, v;
  adj_t j, l, m, vnedges;
  vtx_t nborder;
  wgt_t d;

  adj_t * htable;
  vtx_t * perm, * keys, * border;

  /* voronoi graph parts */
  adj_t * vxadj, * adjori;
  vtx_t * vadjncy, * adjsrc;
  wgt_t * vadjwgt;


  DL_ASSERT(check_graph(nvtxs,xadj,adjncy,adjwgt),"Bad graph passed " \
      "into voronoi_diagram().");

  DL_ASSERT(nsources <= nvtxs,"More sources than vertices: "PF_VTX_T" vs " \
      PF_VTX_T"\n",nsources,nvtxs);

  /* maximum sized border has everything */
  border = vtx_alloc(nvtxs+1);

  /* find my border vertices and count maximum possible of edges */
  vnedges = 0;
  nborder = 0;
  for (i=0;i<nvtxs;++i) {
    me = where[i];
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (me != where[k]) {
        border[nborder++] = i;
        vnedges += xadj[i+1] - xadj[i];
        break;
      }
    }
  }

  /* shuffle the border vertices */
  vtx_pseudo_shuffle(border,nborder/8,nborder);

  /* sort boundary vertices based on voronoi region */
  keys = vtx_alloc(nvtxs);
  for (i=0;i<nborder;++i) {
    keys[i] = where[border[i]];
  }
  perm = vtx_alloc(nborder);
  vtx_countingsort_kv(keys,border,0,nsources,nborder,perm,NULL);

  /* create a prefix sum array in 'border' */
  vtx_set(border,0,nsources+1);
  for (i=0;i<nborder;++i) {
    ++border[keys[i]];
  }
  vtx_prefixsum_exc(border,nsources+1);
  dl_free(keys);

  /* allocate mst structure */
  vxadj = adj_calloc(nsources+1);
  vadjncy = vtx_alloc(vnedges);
  vadjwgt = wgt_alloc(vnedges);
  adjori = adj_alloc(vnedges);
  adjsrc = vtx_alloc(vnedges);
  htable = adj_init_alloc(NULL_ADJ,nsources);

  /* create sources graph */ 
  t = 0;
  l = vxadj[0] = 0;
  for (t=0;t<nsources;++t) {
    for (v=border[t];v<border[t+1];++v) {
      i = perm[v];
      DL_ASSERT_EQUALS(where[i],t,PF_VTX_T);
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        p = where[k];
        if (p != t) {
          if (adjwgt) {
            d = adjwgt[j] + (dist[i] + dist[k]);
          } else {
            d = 1 + (dist[i] + dist[k]);
          }
          if ((m = htable[p]) == NULL_ADJ) {
            vadjncy[l] = p;
            vadjwgt[l] = d;
            adjsrc[l] = i;
            adjori[l] = j;
            htable[p] = l++;
          } else if (d < vadjwgt[m] || (d == vadjwgt[m] && \
                (i+k < adjsrc[m] + adjncy[adjori[m]] || \
                  (i+k == adjsrc[m] + adjncy[adjori[m]] && \
                   i*k < adjsrc[m]*adjncy[adjori[m]]) ) ) ) {
            /* update edge -- deterministically regardless of order */
            vadjwgt[m] = d; 
            adjsrc[m] = i;
            adjori[m] = j;
          }
        }
      }
    }
    vxadj[t+1] = l;

    /* clear htable */
    for (j=vxadj[t];j<vxadj[t+1];++j) {
      p = vadjncy[j];
      htable[p] = NULL_ADJ;
    }
  }
  dl_free(border);
  dl_free(perm);
  dl_free(htable);

  DL_ASSERT(check_graph(nsources,vxadj,vadjncy,vadjwgt),"Bad " \
      "voronoi diagram generated.");

  if (r_vxadj) {
    *r_vxadj = vxadj;
  } else {
    dl_free(vxadj);
  }
  if (r_vadjncy) {
    *r_vadjncy = vadjncy;
  } else {
    dl_free(vadjncy);
  }
  if (r_vadjwgt) {
    *r_vadjwgt = vadjwgt;
  } else {
    dl_free(vadjwgt);
  }
  if (r_adjsrc) {
    *r_adjsrc = adjsrc;
  } else {
    dl_free(adjsrc);
  }
  if (r_adjori) {
    *r_adjori = adjori;
  } else {
    dl_free(adjori);
  }
}





#endif

/**
 * @file csr.c
 * @brief Functions for reading and writing graphs stored in CSR format
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2014-02-13
 */




#ifndef BOWSTRING_IO_CSR_C
#define BOWSTRING_IO_CSR_C




#include "csr.h"




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct vertex_t {
  vtx_t deg;
  vtx_t * adjncy;
  wgt_t * adjwgt;
  struct vertex_t * next;
} vertex_t;




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const BUFFERSIZE = 0x1000;
static char const COMMENT_CHARS[256] = {
  ['#']=1,
  ['%']=1,
  ['\'']=1,
  ['"']=1,
  ['/']=1
};




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static inline vertex_t * __create_vertex(
    vtx_t const deg)
{
  vertex_t * vtx;

  vtx = malloc(sizeof(vertex_t)+((sizeof(vtx_t)+sizeof(wgt_t))*deg));

  vtx->deg = deg;
  vtx->next = NULL;
  vtx->adjncy = (vtx_t*)(vtx+1);
  vtx->adjwgt = (wgt_t*)(vtx->adjncy+deg);

  return vtx;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int read_csr_graph(
    char const * const filename, 
    vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  int rv;
  ssize_t ll;
  size_t bufsize;
  adj_t j;
  vtx_t i, n, minvtx, maxvtx, deg;
  wgt_t f;
  file_t * fin = NULL;
  char * line = NULL, * eptr, * sptr;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL;
  vertex_t * oldvtx, * vtx;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(filename,"r",&fin)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  n = 0;
  j = 0;
  minvtx = 1;
  maxvtx = 0;
  vtx = NULL;
  oldvtx = NULL;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* a blank line we'll assume means an empty row */
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }

    /* count the row */
    deg = (vtx_t)(dl_get_ne_str(line) / 2);

    if (vtx == NULL) {
      oldvtx = vtx = __create_vertex(deg);
    } else {
      vtx->next = __create_vertex(deg);
      vtx = vtx->next;
    }

    /* save the row */
    deg = 0;
    sptr = line;
    i = (vtx_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    f = (wgt_t)strtod(sptr,&eptr);
    while (eptr != sptr) {
      if (i > maxvtx) {
        maxvtx = i;
      }
      if (i < minvtx) {
        minvtx = i;
      }

      vtx->adjncy[deg] = i;
      vtx->adjwgt[deg] = f;
      ++deg;

      sptr = eptr;
      i = (vtx_t)strtoull(sptr,&eptr,10);
      sptr = eptr;
      f = (wgt_t)strtod(sptr,&eptr);
    }
    j += deg;

    ++n;
  }

  dl_reset_file(fin);

  *r_xadj = xadj = adj_alloc(n+1); 
  *r_adjncy = adjncy = vtx_alloc(j); 
  if (r_adjwgt) {
    *r_adjwgt = adjwgt = wgt_alloc(j); 
  }

  n = 0;
  j = 0;
  xadj[0] = 0;
  vtx = oldvtx;
  while(vtx != NULL) {
    for (i=0;i<vtx->deg;++i) {
      /* mincol is at most 1 */
      adjncy[j] = vtx->adjncy[i]-minvtx;
      if (r_adjwgt) {
        adjwgt[j] = vtx->adjwgt[i];
      }
      ++j;
    }
    oldvtx = vtx;
    vtx = oldvtx->next;
    dl_free(oldvtx);
    xadj[++n] = j;
  }

  dl_free(line);
  dl_close_file(fin);

  *r_nvtxs = n;

  return BOWSTRING_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (xadj) {
    dl_free(xadj);
  }
  if (adjncy) {
    dl_free(adjncy);
  }
  if (adjwgt) {
    dl_free(adjwgt);
  }

  return rv;
}


int write_csr_graph(
    char const * const filename, 
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    wgt_t const * const adjwgt)
{
  int rv;
  file_t * file;
  vtx_t i;
  adj_t j;

  if ((rv = __open_file(filename,"w",&file)) != BOWSTRING_SUCCESS) {
    return rv;
  }

  /* print out the graph */
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      dl_fprintf(file,"%zu ",(size_t)(adjncy[j]+1));
      if (adjwgt) {
        dl_fprintf(file,"%zu ",(size_t)adjwgt[j]);
      } else {
        dl_fprintf(file,"1.0 ");
      }
    }
    dl_fprintf(file,"\n");
  }
  dl_close_file(file);

  return BOWSTRING_SUCCESS;
}




#endif

/**
 * @file snap.c
 * @brief Functions for reading and writing SNAP graph files
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2014-02-13
 */




#ifndef BOWSTRING_IO_SNAP_C
#define BOWSTRING_IO_SNAP_C




#include "snap.h"




/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


#define DLSORT_PREFIX vtx
#define DLSORT_TYPE_T vtx_t
#define DLSORT_DLTYPE DLTYPE_INTEGRAL
#define DLSORT_STATIC
#include "dlsort_headers.h"
#undef DLSORT_STATIC
#undef DLSORT_DLTYPE
#undef DLSORT_PREFIX
#undef DLSORT_TYPE_T




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const size_t BUFFERSIZE = 0x1000;
static const char COMMENT_CHARS[256] = {
  ['#']=1,
  ['%']=1,
  ['\'']=1,
  ['"']=1,
  ['/']=1
};


/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __quicksort_edges(vtx_t * const keys,
    wgt_t * const vals, const size_t n)
{
  vtx_t mid;
  size_t i,j,k;
  i = 1;
  j = n-1;
  k = n >> 1;
  mid = keys[k];
  keys[k] = keys[0];
  while (i < j) {
    if (keys[i] > mid) { /* keys[i] is on the wrong side */
      if (keys[j] <= mid) {
        dl_swap(keys[i],keys[j]);
        dl_swap(vals[i],vals[j]);
        ++i;
      }
      --j;
    } else {
      if (keys[j] > mid) { /* keys[j] is on the right side */
        --j;
      }
      ++i;
    }
  }
  if (keys[i] > mid) {
    --i;
  }
  keys[0] = keys[i];
  keys[i] = mid;
  if (i > 1) {
    __quicksort_edges(keys,vals,i);
  }
  ++i; /* skip the pivot element */
  if (n-i > 1) {
    __quicksort_edges(keys+i,vals+i,n-i);
  }
} 






/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int read_snap_graph(const char * const filename, vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, vtx_t ** const r_adjncy, wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  file_t * file;
  int rv, do_adjwgt;
  vtx_t nvtxs, i, j, u, v, maxv;
  adj_t nedges;
  ssize_t ll;
  size_t bufsize;
  char * line = NULL,*eptr, * sptr;
  wgt_t w;

  vtx_t * labels = NULL;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(filename,"r",&file)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }
  line = char_alloc(bufsize);

  do_adjwgt = 0;

  nedges = 0;

  /* first pass to count edges and vertices */
  maxv = 0;
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (vtx_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    j = (vtx_t)strtoull(sptr,&eptr,10);
    sptr = eptr;
    w = (wgt_t)strtod(sptr,&eptr);
    if (eptr != sptr) {
      /* if anyone is missing weight, we'll assign 1 */
      printf("Parsing edge weight from SNAP graph\n");
      do_adjwgt = 1;
    }
    if (i > maxv) {
      maxv = i;
    }
    if (j > maxv) {
      maxv = j;
    }
    /* count edges in both directions */
    nedges+=2;
  }
  nvtxs = maxv+1;

  labels = vtx_init_alloc(NULL_VTX,nvtxs);
  xadj = adj_calloc(nvtxs+1);
  adjncy = vtx_alloc(nedges);
  if (do_adjwgt) {
    adjwgt = wgt_alloc(nedges);
  }
  dl_reset_file(file);

  /* second pass to populate xadj */
  nvtxs = 0;
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = strtoull(sptr,&eptr,10);
    if (labels[i] == NULL_VTX) {
      u = labels[i] = nvtxs++;
    } else {
      u = labels[i];
    }
    sptr = eptr;
    j = strtoull(sptr,&eptr,10);
    if (labels[j] == NULL_VTX) {
      v = labels[j] = nvtxs++;
    } else {
      v = labels[j];
    }
    ++xadj[u+1];
    ++xadj[v+1];
  }
  xadj[0] = 0;
  adj_prefixsum_exc(xadj+1,nvtxs);
  dl_reset_file(file);

  DL_ASSERT(xadj[1] == 0,"Broken xadj for loading\n");

  /* final pass */
  w = 0; /* stops uninitialized complaints */
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    u = labels[strtoull(sptr,&eptr,10)];
    sptr = eptr;
    v = labels[strtoull(sptr,&eptr,10)];
    sptr = eptr;
    if (do_adjwgt) {
      w = (wgt_t)strtod(sptr,&eptr);
      if (eptr == sptr) {
        /* if anyone is missing weight, we'll assign 1 */
        w = 1;
      }
    }
    adjncy[xadj[u+1]] = v;
    adjncy[xadj[v+1]] = u;
    if (do_adjwgt) {
      adjwgt[xadj[u+1]] = w;
      adjwgt[xadj[v+1]] = w;
    }
    ++xadj[u+1];
    ++xadj[v+1];
  }
  dl_free(labels);
  dl_free(line);
  dl_close_file(file);

  /* sort the edges */
  if (do_adjwgt) {
    for (i=0;i<nvtxs;++i) {
      __quicksort_edges(adjncy+xadj[i],adjwgt+xadj[i],xadj[i+1]-xadj[i]); 
    }
  } else {
    for (i=0;i<nvtxs;++i) {
      vtx_quicksort(adjncy+xadj[i],xadj[i+1]-xadj[i]); 
    }
  }

  /* remove duplicates */
  nedges = 0;
  for (i=0;i<nvtxs;++i) {
    j=xadj[i];
    xadj[i] = nedges;
    if (j<xadj[i+1]) {
      adjncy[nedges] = adjncy[j];
      if (do_adjwgt) {
        adjwgt[nedges] = adjwgt[j];
      }
      ++nedges;
      ++j;
    }
    for (;j<xadj[i+1];++j) {
      if (adjncy[j] != adjncy[j-1]) {
        adjncy[nedges] = adjncy[j];
        if (do_adjwgt) {
          adjwgt[nedges] = adjwgt[j];
        }
        ++nedges;
      }
    }
  }
  xadj[nvtxs] = nedges;

  if (r_nvtxs) {
    *r_nvtxs = nvtxs;
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
  /* doesn't support vwgt */
  if (r_vwgt) {
    *r_vwgt = NULL;
  }
  if (r_adjwgt) {
    *r_adjwgt = adjwgt;
  } else if (adjwgt) {
    dl_free(adjwgt);
  }

  return BOWSTRING_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }
  if (labels) {
    dl_free(labels);
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


int write_snap_graph(const char * const filename, const vtx_t nvtxs, 
    const adj_t * const xadj, const vtx_t * const adjncy, 
    const wgt_t * const vwgt, const wgt_t * const adjwgt)
{
  return BOWSTRING_ERROR_UNIMPLEMENTED;
}




#endif

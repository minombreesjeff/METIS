/**
 * @file cloud9.c
 * @brief Functions for reading and writing graphs in the Cloud9 format
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2014-02-13
 */




#ifndef BOWSTRING_IO_CLOUD9_C
#define BOWSTRING_IO_CLOUD9_C




#include "cloud9.h"




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
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int read_cloud9_graph(const char * const filename, vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, vtx_t ** const r_adjncy, wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  return BOWSTRING_ERROR_UNIMPLEMENTED;
}


int write_cloud9_graph(const char * const filename, const vtx_t nvtxs, 
    const adj_t * const xadj, const vtx_t * const adjncy, 
    const wgt_t * const vwgt, const wgt_t * const adjwgt)
{
  file_t * file;
  size_t i, j;
  int do_vwgt, do_adjwgt, rv;

  int wgt_flags = 0;
  do_vwgt = do_adjwgt = 0;

  /* check if we should write the weights */
  if (vwgt) {
    for (i=0;i<nvtxs;++i) {
      if (vwgt[i] != 1.0) {
        break;
      }
    }
    if (i < nvtxs) {
      wgt_flags |= VWGT_FLAG;
      do_vwgt = 1;
    }
  }
  if (adjwgt) {
    for (i=0;i<nvtxs;++i) {
      if (adjwgt[i] != 1.0) {
        break;
      }
    }
    if (i < nvtxs) {
      wgt_flags |= ADJWGT_FLAG;
      do_adjwgt = 1;
    }
  }

  if ((rv = __open_file(filename,"w",&file)) != BOWSTRING_SUCCESS) {
    return rv;
  }

  /* print out the graph */
  for (i=0;i<nvtxs;++i) {
    dl_fprintf(file,"%zu ",i);
    if (do_vwgt) {
      dl_fprintf(file,"%lf ",(double)vwgt[i]);
    }
    for (j=(size_t)xadj[i];j<(size_t)xadj[i+1];++j) {
      dl_fprintf(file,"%zu ",(size_t)(adjncy[j]+1));
      if (do_adjwgt) {
        dl_fprintf(file,"%lf ",(double)adjwgt[j]);
      }
    }
    dl_fprintf(file,"\n");
  }
  dl_close_file(file);

  return BOWSTRING_SUCCESS;
}


#endif

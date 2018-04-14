/**
 * @file io.h
 * @brief IO File operations for graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-04-09
 */




#ifndef BOWSTRING_IO_H
#define BOWSTRING_IO_H




#include "../base.h"




/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/

/* top level reads */

int io_read_graph(
    char const * filename, 
    int graphtype, 
    vtx_t * nvtxs, 
    adj_t ** xadj, 
    vtx_t ** adjncy, 
    wgt_t ** vwgt, 
    wgt_t ** adjwgt);


int io_write_graph(
    char const * filename, 
    int graphtype, 
    vtx_t nvtxs, 
    const adj_t * xadj, 
    const vtx_t * adjncy, 
    const wgt_t * vwgt, 
    const wgt_t * adjwgt);


/* labels */

int read_vertex_labels(
    char const * filename, 
    vtx_t * nvtxs, 
    vlbl_t ** labels);


int write_vertex_labels(
    char const * filename, 
    vtx_t nvtxs, 
    const vlbl_t * labels);


int read_edge_labels(
    char const * filename, 
    adj_t * nedges, 
    elbl_t ** labels);


int write_edge_labels(
    char const * filename, 
    adj_t nedges, 
    const elbl_t * labels);


int write_weights(
    char const * filename,
    wgt_t const * wgts,
    size_t nwgts);


/* translations */

int translate_metis_to_snap(
    char const * infile, 
    char const * outfile);


int translate_metis_to_cloud9(
    char const * infile, 
    char const * outfile);


int translate_metis_to_dimacs(
    char const * infile, 
    char const * outfile);


int translate_csr_to_snap(
    char const * infile, 
    char const * outfile);


int translate_csr_to_cloud9(
    char const * infile, 
    char const * outfile);


int translate_csr_to_pegasus(
    char const * infile, 
    char const * outfile);


int translate_csr_to_matrixmarket(
    char const * infile, 
    char const * outfile);


int translate_edgelist_to_snap(
    char const * infile, 
    char const * outfile);


/* misc */

int tile_csr(
    char const * infile, 
    char const * outfile, 
    vtx_t nx, 
    vtx_t ny);




/******************************************************************************
* INLINE FUNCTIONS ************************************************************
******************************************************************************/


/**
 * @brief Convert bitwise weight flags to base 10 format:
 *    1 = 1, 2 = 10, 4 = 100, etc.
 *
 * @param flags The bitwised flags
 *
 * @return The base 10 representation suitable for printing
 */
static inline int __flag_to_int(
    int flags) 
{
  int iflags = 0;

  int p10 = 1;
  int p2 = 1;
  while (flags > 0) {
    if (flags & p2) {
      iflags += p10;
      flags ^= p2;
    }
    /* obnoxious way to multiply by 10 */
    p10 = (p10 << 3) + (p10 << 1);
    p2 <<= 1;
  }
  
  return iflags;
}


/**
 * @brief Wrapper around the dl_file routines to handle error codes
 *
 * @param filename The file to open
 * @param mode The mode to open the file in
 * @param r_file The reference to the file pointer
 *
 * @return BOWSTRING_SUCCESS, or an error code
 */
static inline int __open_file(
    char const * const filename, 
    char const * const mode, 
    file_t ** const r_file)
{
  int rv;
  if ((rv = dl_open_file(filename,mode,r_file)) != DL_FILE_SUCCESS) {
    switch (rv) {
      case DL_FILE_BAD_PARAMETERS:
      case DL_FILE_PATH_PARSE_FAILURE:
        eprintf("Bad filename '%s'\n",filename);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        break;
      case DL_FILE_PATH_BAD:
      case DL_FILE_FILENOTFOUND:
        eprintf("File not found '%s'\n",filename);
        rv = BOWSTRING_ERROR_FILENOTFOUND;
        break;
      case DL_FILE_PATH_ACCESS_DENIED:
      case DL_FILE_READ_ACCESS_DENIED:
      case DL_FILE_WRITE_ACCESS_DENIED:
        eprintf("Permission denied '%s'\n",filename);
        rv = BOWSTRING_ERROR_FILEPERMISSIONDENIED;
        break;
      default:
        eprintf("Unknown failure: %d opening '%s'\n",rv,filename);
        rv = BOWSTRING_ERROR_UNKNOWN;
        break;
    }
  } else {
    rv = BOWSTRING_SUCCESS;
  }
  return rv;
}




#endif

/**
 * @file io.c
 * @brief IO file operations for graphs
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2012
 * @version 1
 * @date 2013-04-09
 */




#ifndef BOWSTRING_IO_C
#define BOWSTRING_IO_C




#include "io.h"
#include "metis.h"
#include "nerstrand.h"
#include "stp.h"
#include "dimacs.h"
#include "snap.h"
#include "csr.h"
#include "matrix.h"
#include "cloud9.h"
#include <math.h>




/******************************************************************************
* DOMLIB MACROS ***************************************************************
******************************************************************************/


#define DLTREE_STATIC
#define DLTREE_PREFIX cv
#define DLTREE_KEY_T char const *
#define DLTREE_VAL_T vtx_t 
#include "dltree_headers.h"
#undef DLTREE_PREFIX
#undef DLTREE_KEY_T
#undef DLTREE_VAL_T
#undef DLTREE_STATIC




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define __READ_INT64(var,buffer,file,msg,fn) \
  do { \
    if (fread((buffer),8,1,(file)) != 1) { \
      dl_error("Failed reading "msg" from '%s'\n",(fn)); \
    } \
    dl_from_bytes(&(var),buffer,8); \
  } while (0)

#define __READ_INT(var,buffer,file,msg,fn) \
  do { \
    if (fread((buffer),4,1,(file)) != 1) { \
      dl_error("Failed reading "msg" from '%s'\n",(fn)); \
    } \
    dl_from_bytes(&(var),buffer,4); \
  } while (0)

#define __READ_BYTE(var,buffer,file,msg,fn) \
  do { \
    if (fread((buffer),1,1,(file)) != 1) { \
      dl_error("Failed reading "msg" from '%s'\n",(fn)); \
    } \
    dl_from_bytes(&(var),buffer,1); \
  } while (0)

#define __WRITE_INT64(var,buffer,file,msg,fn) \
  do { \
    dl_to_bytes((buffer),&(var),8); \
    if (fwrite((buffer),8,1,(file)) != 1) { \
      dl_error("Failed writing "msg" to '%s'\n",(fn)); \
    } \
  } while (0)

#define __WRITE_INT(var,buffer,file,msg,fn) \
  do { \
    dl_to_bytes((buffer),&(var),4); \
    if (fwrite((buffer),4,1,(file)) != 1) { \
      dl_error("Failed writing "msg" to '%s'\n",(fn)); \
    } \
  } while (0)

#define __WRITE_BYTE(var,buffer,file,msg,fn) \
  do { \
    dl_to_bytes((buffer),&(var),1); \
    if (fwrite((buffer),1,1,(file)) != 1) { \
      dl_error("Failed writing "msg" to '%s'\n",(fn)); \
    } \
  } while (0)

#define __CONVERT_ARRAY(var,buffer,bsize,width,arr,i,type) \
  do { \
    size_t j; \
    for ((j)=0;(j)<(bsize);++(j),++(i)) { \
      dl_from_bytes(&(var),(buffer)+((j)*(width)),width); \
      (arr)[i] = (type)(var); \
    } \
  } while (0)

#define __SETSTATS(stat,stat_max,stat_avg,level,level_max) \
  do { \
    if ((stat) > (stat_max)) { \
      (stat_max) = (stat); \
      (level_max) = (level); \
    } \
    (stat_avg) += (stat); \
  } while (0)




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
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int io_read_graph(
    char const * const filename, 
    int const graphtype, 
    vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  switch (graphtype) {
    case BOWSTRING_FORMAT_METIS:
      return read_metis_graph(filename,r_nvtxs,r_xadj,r_adjncy,r_vwgt,
          r_adjwgt);
    case BOWSTRING_FORMAT_CSR:
      return read_csr_graph(filename,r_nvtxs,r_xadj,r_adjncy,r_vwgt,r_adjwgt);
    case BOWSTRING_FORMAT_CLOUD9:
      return read_cloud9_graph(filename,r_nvtxs,r_xadj,r_adjncy,r_vwgt,
          r_adjwgt);
    case BOWSTRING_FORMAT_SNAP:
      return read_snap_graph(filename,r_nvtxs,r_xadj,r_adjncy,r_vwgt,r_adjwgt);
    case BOWSTRING_FORMAT_DIMACS:
      return read_dimacs_graph(filename,r_nvtxs,r_xadj,r_adjncy,r_vwgt,
          r_adjwgt);
    case BOWSTRING_FORMAT_STP:
      return read_stp_graph(filename,r_nvtxs,r_xadj,r_adjncy,r_vwgt,r_adjwgt);
    case BOWSTRING_FORMAT_NERSTRAND:
      return read_nerstrand_graph(filename,r_nvtxs,r_xadj,r_adjncy,r_vwgt,
          r_adjwgt,0);
    case BOWSTRING_FORMAT_NBG:
      return read_nerstrand_graph(filename,r_nvtxs,r_xadj,r_adjncy,r_vwgt,
          r_adjwgt,1);
    default:
      eprintf("Unknown graph type '%d'\n",graphtype);
      return BOWSTRING_ERROR_INVALIDINPUT;
  }
}


int io_write_graph(
    char const * const filename, 
    int graphtype, 
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    wgt_t const * const adjwgt)
{
  switch (graphtype) {
    case BOWSTRING_FORMAT_METIS:
      return write_metis_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt);
    case BOWSTRING_FORMAT_CSR:
      return write_csr_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt);
    case BOWSTRING_FORMAT_CLOUD9:
      return write_cloud9_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt);
    case BOWSTRING_FORMAT_SNAP:
      return write_snap_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt);
    case BOWSTRING_FORMAT_DIMACS:
      return write_dimacs_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt);
    case BOWSTRING_FORMAT_STP:
      return write_stp_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt);
    case BOWSTRING_FORMAT_NBG:
      return write_nerstrand_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt,1);
    #ifdef XXX
    case BOWSTRING_FORMAT_PEGASUS:
      return write_pegasus_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt);
    #endif
    default:
      eprintf("Unknown graph type '%d'\n",graphtype);
      return BOWSTRING_ERROR_INVALIDINPUT;
  }
}


/* META FILES ****************************************************************/

int read_vertex_labels(
    char const * const filename, 
    vtx_t * const r_nvtxs, 
    vlbl_t ** const r_labels)
{
  file_t * file = NULL;
  vtx_t nl, nvtxs;
  vlbl_t l;
  int rv = BOWSTRING_SUCCESS;
  ssize_t ll;
  size_t bufsize;
  char * line = NULL;
  vlbl_t * labels = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(filename,"r",&file)) != BOWSTRING_SUCCESS) {
    goto CLEANUP;
  }
  line = char_alloc(bufsize);
  ll = dl_get_next_line(file,&line,&bufsize);
  /* skip comments */
  while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
    ll = dl_get_next_line(file,&line,&bufsize);
  }
  if (!r_nvtxs || *r_nvtxs == 0) {
    /* determine the size of the file */
    nl = 0;
    while (ll > 0 && sscanf(line,PF_VLBL_T,&l) == 1) {
      ++nl; 
      ll = dl_get_next_line(file,&line,&bufsize);
    }
    dl_reset_file(file);
    ll = dl_get_next_line(file,&line,&bufsize);
    /* skip comments */
    while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
      ll = dl_get_next_line(file,&line,&bufsize);
    }
  } else {
    /* the size of the file was specified */
    nl = *r_nvtxs;
  }
  nvtxs = nl;
  labels = vlbl_alloc(nvtxs);
  nl = 0;
  /* read in the file */
  while (ll > 0 && sscanf(line,PF_VLBL_T,&l) == 1) {
    if (nl == nvtxs) {
      eprintf("More vertices in partition file than specified.\n");
      rv = BOWSTRING_ERROR_INVALIDINPUT;
      goto CLEANUP;
    }
    labels[nl++] = l;
    ll = dl_get_next_line(file,&line,&bufsize);
  }

  dl_close_file(file);
  file = NULL;

  if (r_nvtxs) {
    *r_nvtxs = nl;
  }

  if (r_labels) {
    *r_labels = labels;
    labels = NULL;
  } else {
    dl_free(labels);
  }

  CLEANUP:

  if (file) {
    dl_close_file(file);
  }
  if (line) {
    dl_free(line);
  } 
  if (labels) {
    dl_free(labels);
  }

  return rv;
}


int write_vertex_labels(
    char const * const filename, 
    vtx_t const nvtxs, 
    vlbl_t const * const labels)
{
  int rv;
  vtx_t i;
  file_t * file = NULL;

  if ((rv = __open_file(filename,"w",&file)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  for (i=0;i<nvtxs;++i) {
    dl_fprintf(file,PF_VLBL_T"\n",labels[i]);
  }

  dl_close_file(file);
  file = NULL;

  return BOWSTRING_SUCCESS;

  ERROR:

  if (file) {
    dl_close_file(file);
  }

  return rv;
}


int read_edge_labels(
    char const * const filename, 
    adj_t * const r_nedges, 
    elbl_t ** const r_labels)
{
  file_t * file;
  adj_t nl;
  elbl_t l;
  double d;
  int rv;
  ssize_t ll;
  size_t bufsize;
  char * line, * eptr, * sptr;

  elbl_t * labels = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(filename,"r",&file)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }
  line = char_alloc(bufsize);
  ll = dl_get_next_line(file,&line,&bufsize);
  /* skip comments */
  while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
    ll = dl_get_next_line(file,&line,&bufsize);
  }
  if (*r_nedges == 0) {
    /* determine the size of the file */
    nl = 0;
    while (ll > 0) {
      eptr = sptr = line;
      d = strtod(sptr,&eptr);
      if (isfinite(d)) {
        rv = BOWSTRING_ERROR_INVALIDVALUE;
        goto ERROR;
      }
      while (eptr != sptr) {
        ++nl; 
        d = strtod(sptr,&eptr);
        if (isfinite(d)) {
          rv = BOWSTRING_ERROR_INVALIDVALUE;
          goto ERROR;
        }
      }
      ll = dl_get_next_line(file,&line,&bufsize);
    }
    dl_reset_file(file);
    ll = dl_get_next_line(file,&line,&bufsize);
    /* skip comments */
    while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
      ll = dl_get_next_line(file,&line,&bufsize);
    }
  } else {
    /* the size of the file was specified */
    nl = *r_nedges;
  }
  labels = elbl_alloc(nl);
  nl = 0;
  /* read in the file */
  while (ll > 0) {
    eptr = sptr = line;
    l = (elbl_t)strtod(sptr,&eptr);
    while (eptr != sptr) {
      labels[nl++] = l;
      l = (elbl_t)strtod(sptr,&eptr);
    }
    ll = dl_get_next_line(file,&line,&bufsize);
  }

  *r_nedges = nl;

  dl_free(line);

  if (r_labels) {
    *r_labels = labels;
  } else {
    dl_free(labels);
  }

  return BOWSTRING_SUCCESS; 

  ERROR:

  if (line) {
    dl_free(line);
  } 
  if (labels) {
    dl_free(labels);
  }

  return rv;
}


int write_edge_labels(
    char const * const filename, 
    adj_t const nedges, 
    elbl_t const * const labels)
{
  int rv;
  adj_t i;
  file_t * file = NULL;

  if ((rv = __open_file(filename,"w",&file)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  for (i=0;i<nedges;++i) {
    dl_fprintf(file,PF_ELBL_T"\n",labels[i]);
  }

  dl_close_file(file);
  file = NULL;

  return BOWSTRING_SUCCESS;

  ERROR:

  if (file) {
    dl_close_file(file);
  }

  return rv;
}


int write_weights(
    char const * const filename,
    wgt_t const * const wgts,
    size_t const nwgts)
{
  int rv;
  size_t i;
  file_t * file = NULL;

  if ((rv = __open_file(filename,"w",&file)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  for (i=0;i<nwgts;++i) {
    dl_fprintf(file,PF_WGT_T"\n",wgts[i]);
  }

  dl_close_file(file);
  file = NULL;

  return BOWSTRING_SUCCESS;

  ERROR:

  if (file) {
    dl_close_file(file);
  }

  return rv;
}



/* TRANSLATION ***************************************************************/

int translate_metis_to_snap(
    char const * const infile, 
    char const * const outfile)
{
  file_t * fin = NULL, *fout = NULL;
  int w, rv;
  size_t n,m,i, j, ln, ne, num_header;
  ssize_t ll;
  size_t bufsize;
  char * line,*eptr, * sptr;
  wgt_t weight;

  bufsize = BUFFERSIZE;

  /* prep input file */
  if ((rv = __open_file(infile,"r",&fin)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  line = char_alloc(bufsize);
  ll = dl_get_next_line(fin,&line,&bufsize);
  /* skip comments */
  while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
    ll = dl_get_next_line(fin,&line,&bufsize);
  }
  num_header = sscanf(line,"%zu %zu %d",&n,&m,&w);

  int do_vwgt, do_adjwgt;

  if (num_header == 3) {
    do_vwgt = ((w & VWGT_FLAG) ? 1 : 0);
    do_adjwgt = ((w & ADJWGT_FLAG) ? 1 : 0);
  } else if (num_header == 2) {
    do_vwgt = 0;
    do_adjwgt = 0;
  } else {
    eprintf("Invalid header line '%s'\n",line);
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  /* prep output file */
  if ((rv = __open_file(outfile,"w",&fout)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  /* make edges count twice */
  m *= 2;
  if (do_vwgt) {
    printf("Ignoring vertex weights when converting to SNAP\n");
  }

  ln = j = 0;
  /* scan through the file */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      if (ln == n) {
        /* if we've read all vertices, ignore trailing blank lines */
        continue;
      }
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    ne = 0;
    i = strtoull(sptr,&eptr,10)-1;
    while (eptr != sptr) {
      sptr = eptr;
      if (i >= n) {
        eprintf("Invalid vertex %zu at line %zu\n",i+1,ln);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        goto ERROR;
      }
      if (j >= m) {
        eprintf("Edge %zu found at vertex %zu is greater "
            "than specified number of edges (%zu)\n", j,ln,m);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        goto ERROR;
      }
      dl_fprintf(fout,"%zu %zu",ln,i);
      /* edge weight */ 
      if (do_adjwgt) {
        weight = (wgt_t)strtod(sptr,&eptr);
        if (eptr == sptr) {
          eprintf("Missing edge weight for vertex %zu, edge %zu\n",ln,ne);
          rv = BOWSTRING_ERROR_INVALIDINPUT;
          goto ERROR;
        }
        sptr = eptr;
        dl_fprintf(fout," "PF_WGT_T,weight);
      }
      dl_fprintf(fout,"\n");
      ++j;
      ++ne;
      i = strtoull(sptr,&eptr,10)-1;
    }
    ++ln;
  }
  dl_close_file(fin);
  fin = NULL;
  dl_close_file(fout);
  fout = NULL;
  if (ln < n) {
    eprintf("Only found %zu out of %zu vertices in graph file\n",ln,n);
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  } else if (j != m) {
    eprintf("Only found %zu out of %zu edges in graph file\n",j,m);
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  dl_free(line);

  return BOWSTRING_SUCCESS;

  ERROR:
  if (fin) {
    dl_close_file(fin);
    fin = NULL;
  }
  if (fout) {
    dl_close_file(fout);
    fout = NULL;
  }
  dl_free(line);
  return rv;
}


int translate_metis_to_cloud9(
    char const * infile, 
    char const * outfile)
{ 
  file_t * fin = NULL, *fout = NULL;
  int w, rv;
  size_t n,m,i, j, ln, ne, num_header;
  ssize_t ll;
  size_t bufsize;
  char * line = NULL;
  char * eptr, * sptr;
  wgt_t weight;

  bufsize = BUFFERSIZE;

  /* prep input file */
  if ((rv = __open_file(infile,"r",&fin)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }
  line = char_alloc(bufsize);
  ll = dl_get_next_line(fin,&line,&bufsize);
  /* skip comments */
  while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
    ll = dl_get_next_line(fin,&line,&bufsize);
  }
  num_header = sscanf(line,"%zu %zu %d",&n,&m,&w);

  int do_vwgt, do_adjwgt;

  if (num_header == 3) {
    do_vwgt = ((w & VWGT_FLAG) ? 1 : 0);
    do_adjwgt = ((w & ADJWGT_FLAG) ? 1 : 0);
  } else if (num_header == 2) {
    do_vwgt = 0;
    do_adjwgt = 0;
  } else {
    eprintf("Invalid header line '%s'\n",line);
    dl_close_file(fin);
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  /* prep output file */
  if ((rv = __open_file(outfile,"w",&fout)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  /* make edges count twice */
  m *= 2;
  if (do_vwgt) {
    wprintf("Ignoring vertex weights when converting to Cloud9\n");
  }

  ln = j = 0;
  /* scan through the file */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      if (ln == n) {
        /* if we've read all vertices, ignore trailing blank lines */
        continue;
      }
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    ne = 0;
    dl_fprintf(fout,"%zu",ln);
    i = strtoull(sptr,&eptr,10)-1;
    while (eptr != sptr) {
      sptr = eptr;
      if (i >= n) {
        eprintf("Invalid vertex %zu at line %zu\n",i+1,ln);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        goto ERROR;
      }
      if (j >= m) {
        eprintf("Edge %zu found at vertex %zu is greater "
            "than specified number of edges (%zu)\n", j,ln,m);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        goto ERROR;
      }
      dl_fprintf(fout," %zu",i);
      /* edge weight */ 
      if (do_adjwgt) {
        weight = (wgt_t)strtod(sptr,&eptr);
        if (eptr == sptr) {
          eprintf("Missing edge weight for vertex %zu, edge %zu\n",ln,ne);
          rv = BOWSTRING_ERROR_INVALIDINPUT;
          goto ERROR;
        }
        sptr = eptr;
        dl_fprintf(fout," "PF_WGT_T,weight);
      }
      ++j;
      ++ne;
      i = strtoull(sptr,&eptr,10)-1;
    }
    dl_fprintf(fout,"\n");
    ++ln;
  }
  dl_close_file(fin);
  fin = NULL;
  dl_close_file(fout);
  fout = NULL;
  if (ln < n) {
    eprintf("Only found %zu out of %zu vertices in graph file\n",ln,n);
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  } else if (j != m) {
    eprintf("Only found %zu out of %zu edges in graph file\n",j,m);
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  dl_free(line);

  return BOWSTRING_SUCCESS;

  ERROR:
  if (fin) {
    dl_close_file(fin);
    fin = NULL;
  }
  if (fout) {
    dl_close_file(fout);
    fout = NULL;
  }
  if (line) {
    dl_free(line);
  }
  return rv;
}


int translate_metis_to_dimacs(
    char const * infile, 
    char const * outfile)
{
  return BOWSTRING_ERROR_UNIMPLEMENTED;
}



int translate_csr_to_snap(
    char const * infile, 
    char const * outfile)
{
  file_t * fin = NULL, *fout = NULL;
  int rv;
  size_t i, j, ln, ne;
  ssize_t ll;
  size_t bufsize;
  char * line,*eptr, * sptr;

  bufsize = BUFFERSIZE;
  line = char_alloc(bufsize);

  /* prep input file */
  if ((rv = __open_file(infile,"r",&fin)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  ll = dl_get_next_line(fin,&line,&bufsize);
  /* skip comments */
  while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
    ll = dl_get_next_line(fin,&line,&bufsize);
  }

  /* prep output file */
  if ((rv =__open_file(outfile,"w",&fout)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  /* make edges count twice */
  ln = j = 0;
  /* scan through the file */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* no way to tell if its a trailing blank line or island vertice */
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    ne = 0;
    i = strtoull(sptr,&eptr,10);
    while (eptr != sptr) {
      sptr = eptr;
      dl_fprintf(fout,"%zu %zu",ln,i);
      dl_fprintf(fout,"\n");
      ++j;
      ++ne;
      i = strtoull(sptr,&eptr,10);
    }
    ++ln;
  }
  dl_close_file(fin);
  fin = NULL;
  dl_close_file(fout);
  fout = NULL;

  dl_free(line);

  return BOWSTRING_SUCCESS;

  ERROR:
  if (fin) {
    dl_close_file(fin);
    fin = NULL;
  }
  if (fout) {
    dl_close_file(fout);
    fout = NULL;
  }
  dl_free(line);
  return rv;
}


int translate_csr_to_cloud9(
    char const * infile, 
    char const * outfile)
{ 
  file_t * fin = NULL, *fout = NULL;
  int rv;
  size_t i, j, ln, ne;
  ssize_t ll;
  size_t bufsize;
  char * line = NULL, *eptr, * sptr;

  bufsize = BUFFERSIZE;
  line = char_alloc(bufsize);

  /* prep input file */
  if ((rv = __open_file(infile,"r",&fin)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }
  ll = dl_get_next_line(fin,&line,&bufsize);

  /* prep output file */
  if ((rv = __open_file(outfile,"w",&fout)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  /* make edges count twice */
  ln = j = 0;
  /* scan through the file */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    ne = 0;
    dl_fprintf(fout,"%zu",ln);
    i = strtoull(sptr,&eptr,10);
    while (eptr != sptr) {
      sptr = eptr;
      dl_fprintf(fout," %zu",i);
      ++j;
      ++ne;
      i = strtoull(sptr,&eptr,10);
    }
    dl_fprintf(fout,"\n");
    ++ln;
  }
  dl_close_file(fin);
  fin = NULL;
  dl_close_file(fout);
  fout = NULL;

  dl_free(line);

  return BOWSTRING_SUCCESS;

  ERROR:
  if (fin) {
    dl_close_file(fin);
    fin = NULL;
  }
  if (fout) {
    dl_close_file(fout);
    fout = NULL;
  }
  dl_free(line);
  return rv;
}


int translate_csr_to_pegasus(
    char const * infile, 
    char const * outfile)
{
  file_t * fin = NULL, *fout = NULL;
  int rv;
  size_t i, j, ln, ne;
  ssize_t ll;
  size_t bufsize;
  char * line,*eptr, * sptr;

  bufsize = BUFFERSIZE;
  line = char_alloc(bufsize);

  /* prep input file */
  if ((rv = __open_file(infile,"r",&fin)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  ll = dl_get_next_line(fin,&line,&bufsize);
  /* skip comments */
  while (ll > 0 && COMMENT_CHARS[(unsigned int)line[0]]) {
    ll = dl_get_next_line(fin,&line,&bufsize);
  }

  /* prep output file */
  if ((rv =__open_file(outfile,"w",&fout)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  /* make edges count twice */
  ln = j = 0;
  /* scan through the file */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* no way to tell if its a trailing blank line or island vertice */
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    ne = 0;
    i = strtoull(sptr,&eptr,10);
    while (eptr != sptr) {
      sptr = eptr;
      dl_fprintf(fout,"%zu\t%zu",ln,i);
      dl_fprintf(fout,"\n");
      ++j;
      ++ne;
      i = strtoull(sptr,&eptr,10);
    }
    ++ln;
  }
  dl_close_file(fin);
  fin = NULL;
  dl_close_file(fout);
  fout = NULL;

  dl_free(line);

  return BOWSTRING_SUCCESS;

  ERROR:
  if (fin) {
    dl_close_file(fin);
    fin = NULL;
  }
  if (fout) {
    dl_close_file(fout);
    fout = NULL;
  }
  dl_free(line);
  return rv;

}


int translate_csr_to_matrixmarket(
    char const * infile, 
    char const * outfile)
{
  file_t * fin = NULL, *fout = NULL;
  int rv;
  adj_t ne;
  vtx_t minidx,i,ln,lm;
  ssize_t ll;
  size_t bufsize;
  char * line,*eptr, * sptr;
  wgt_t w;

  bufsize = BUFFERSIZE;
  line = char_alloc(bufsize);

  /* prep input file */
  if ((rv = __open_file(infile,"r",&fin)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  /* make edges count twice */
  minidx = 1;
  lm = ln = 0;
  ne = 0;
  /* scan through the file */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* no way to tell if its a trailing blank line or island vertice */
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    eptr = sptr = line;
    i = strtoull(sptr,&eptr,10);
    sptr = eptr;
    w = (wgt_t)strtod(sptr,&eptr);
    while (eptr != sptr) {
      if (i < minidx) {
        minidx = i;
      }
      if (i > lm) {
        lm = i;
      }
      ++ne;
      sptr = eptr;
      i = strtoull(sptr,&eptr,10);
      sptr = eptr;
      w = (wgt_t)strtod(sptr,&eptr);
    }
    ++ln;
  }
  if (minidx == 0) {
    ++lm;
  }
  dl_reset_file(fin);

  /* prep output file */
  if ((rv =__open_file(outfile,"w",&fout)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  dl_fprintf(fout,"%%%%MatrixMarket matrix coordinate real general\n");
  dl_fprintf(fout,PF_VTX_T" "PF_VTX_T" "PF_ADJ_T"\n",ln,lm,ne);

  /* make edges count twice */
  ln = 0;
  /* scan through the file */
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* no way to tell if its a trailing blank line or island vertice */
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = (vtx_t)(strtoull(sptr,&eptr,10)-minidx);
    sptr = eptr;
    w = (wgt_t)strtod(sptr,&eptr);
    while (eptr != sptr) {
      dl_fprintf(fout,PF_VTX_T" "PF_VTX_T" "PF_WGT_T"\n",ln,i,w);
      sptr = eptr;
      i = (vtx_t)(strtoull(sptr,&eptr,10)-minidx);
      sptr = eptr;
      w = (wgt_t)strtod(sptr,&eptr);
    }
    ++ln;
  }
  dl_close_file(fin);
  fin = NULL;
  dl_close_file(fout);
  fout = NULL;

  dl_free(line);

  return BOWSTRING_SUCCESS;

  ERROR:
  if (fin) {
    dl_close_file(fin);
    fin = NULL;
  }
  if (fout) {
    dl_close_file(fout);
    fout = NULL;
  }
  dl_free(line);
  return rv;
}


int tile_csr(
    char const * infile, 
    char const * outfile, 
    vtx_t const nx, 
    vtx_t const ny)
{
  file_t * fin = NULL, *fout = NULL;
  vtx_t ncol, nrow, j, i, x, y;
  wgt_t w;
  adj_t nnz;
  int rv;
  ssize_t ll;
  size_t bufsize;
  char * line,*eptr, * sptr;

  bufsize = BUFFERSIZE;
  line = char_alloc(bufsize);

  /* prep input file */
  if ((rv = __open_file(infile,"r",&fin)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  /* prep output file */
  if ((rv =__open_file(outfile,"w",&fout)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  dl_fprintf(fout,"%%%%MatrixMarket matrix coordinate real general\n");

  /* scan through the file to find ncol */
  ncol = 0;
  nrow = 0;
  nnz = 0;
  while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* no way to tell if its a trailing blank line or island vertice */
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    sptr = line;
    i = strtoull(sptr,&eptr,10);
    sptr = eptr;
    w = strtod(sptr,&eptr);
    while (eptr != sptr) {
      if (i > ncol) {
        ncol = i;
      }
      sptr = eptr;
      i = strtoull(sptr,&eptr,10);
      sptr = eptr;
      w = strtod(sptr,&eptr);
      ++nnz;
    }
    ++nrow;
  }

  ncol += 1;

  dl_fprintf(fout,PF_VTX_T" "PF_VTX_T" "PF_ADJ_T"\n",nrow*ny,ncol*nx,nnz*nx*ny);

  j= 0;
  for (y=0;y<ny;++y) {
    dl_reset_file(fin);
    /* scan through the file */
    while((ll = dl_get_next_line(fin,&line,&bufsize)) >= 0) {
      if (ll == 0) {
        /* no way to tell if its a trailing blank line or island vertice */
      } else {
        if (COMMENT_CHARS[(unsigned int)line[0]]) {
          /* skip comments */
          continue;
        }
      }
      for (x=0;x<nx;++x) {
        sptr = line;
        i = strtoull(sptr,&eptr,10);
        sptr = eptr;
        w = strtod(sptr,&eptr);
        while (eptr != sptr) {
          dl_fprintf(fout,PF_VTX_T" "PF_VTX_T" "PF_WGT_T"\n",j,i+(x*ncol),w);
          sptr = eptr;
          i = strtoull(sptr,&eptr,10);
          sptr = eptr;
          w = strtod(sptr,&eptr);
        }
      }
      ++j;
    }
  }
  dl_close_file(fin);
  fin = NULL;
  dl_close_file(fout);
  fout = NULL;

  dl_free(line);

  return BOWSTRING_SUCCESS;

  ERROR:
  if (fin) {
    dl_close_file(fin);
    fin = NULL;
  }
  if (fout) {
    dl_close_file(fout);
    fout = NULL;
  }
  dl_free(line);
  return rv;
}


int translate_edgelist_to_snap(
    char const * const infile, 
    char const * const outfile)
{
  int rv;
  vtx_t i, j, nvtxs;
  ssize_t ll;
  char * line = NULL;
  char * v;
  char * key;
  size_t bufsize;
  file_t * ifile = NULL, * ofile = NULL;
  cv_tree_t * tree = NULL;

  bufsize = BUFFERSIZE;

  if ((rv = __open_file(infile,"r",&ifile)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }
  line = char_alloc(bufsize);

  tree = cv_tree_create(strcmp);

  nvtxs = 0;
  /* do a pass over the vertices to count nedges and O(nvtxs) */
  while((ll = dl_get_next_line(ifile,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    v = strtok(line,"\t ");
    key = (char*)malloc(strlen(v)+1);
    memcpy(key,v,strlen(v)+1);
    if (cv_tree_add(key,nvtxs,tree)) {
      ++nvtxs;
    } else {
      dl_free(key);
    }
    key = (char*)malloc(strlen(v)+1);
    memcpy(key,v,strlen(v)+1);
    v = strtok(NULL,"\t ");
    key = (char*)malloc(strlen(v)+1);
    memcpy(key,v,strlen(v)+1);
    if (cv_tree_add(key,nvtxs,tree)) {
      ++nvtxs;
    } else {
      dl_free(key);
    }
  }

  dl_reset_file(ifile);

  if ((rv = __open_file(outfile,"w",&ofile)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  while((ll = dl_get_next_line(ifile,&line,&bufsize)) >= 0) {
    if (ll == 0) {
      /* ignore blank lines */
      continue;
    } else {
      if (COMMENT_CHARS[(unsigned int)line[0]]) {
        /* skip comments */
        continue;
      }
    }
    v = strtok(line,"\t ");
    i = cv_tree_get(v,tree);
    v = strtok(NULL,"\t ");
    j = cv_tree_get(v,tree);
    dl_fprintf(ofile,PF_VTX_T" "PF_VTX_T"\n",i,j);
  }

  dl_close_file(ifile);
  ifile = NULL;
  dl_close_file(ofile);
  ofile = NULL;

  dl_free(line);
  cv_tree_free(tree);

  return BOWSTRING_SUCCESS;

  ERROR:

  if (line) {
    dl_free(line);
  }

  if (ifile) {
    dl_close_file(ifile);
  }

  if (ofile) {
    dl_close_file(ofile);
  }

  if (tree) {
    cv_tree_free(tree);
  }

  return rv;
}




#endif

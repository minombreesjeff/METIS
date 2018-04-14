/**
 * @file nerstrand.c
 * @brief Functions for reading and writing graphs in the Nerstrand graph
 * format.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2014-02-13
 */




#ifndef BOWSTRING_IO_NERSTRAND_C
#define BOWSTRING_IO_NERSTRAND_C




#include "nerstrand.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/

static const uint32_t __FORMAT_VERSION = 1;
static const uint32_t __MAGIC_NUMBER = 0x646f6d3;
static size_t const BUFFERSIZE = 0x8000;
static char const COMMENT_CHARS[256] = {
  ['#']=1,
  ['%']=1,
  ['\'']=1,
  ['"']=1,
  ['/']=1
};




/******************************************************************************
* BSCOMP MACROS ***************************************************************
******************************************************************************/


#define BSCOMP_TYPE_T adj_t
#define BSCOMP_PREFIX adj
#include "bscomp_funcs.h"
#undef BSCOMP_TYPE_T
#undef BSCOMP_PREFIX


#define BSCOMP_TYPE_T vtx_t
#define BSCOMP_PREFIX vtx
#include "bscomp_funcs.h"
#undef BSCOMP_TYPE_T
#undef BSCOMP_PREFIX


#define BSCOMP_TYPE_T wgt_t
#define BSCOMP_PREFIX wgt
#define BSCOMP_TYPE_FLOAT
#include "bscomp_funcs.h"
#undef BSCOMP_TYPE_FLOAT
#undef BSCOMP_TYPE_T
#undef BSCOMP_PREFIX




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/

static int __read_text_graph(
    char const * const filename, 
    vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  file_t * file = NULL;
  int w, rv;
  size_t n,m,i, j, ln, ne, num_header;
  ssize_t ll;
  size_t bufsize;
  char * line,*eptr, * sptr;
  wgt_t weight;

  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * vwgt = NULL;
  wgt_t * adjwgt = NULL;

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
    dl_close_file(file);
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  /* make edges count twice */
  m *= 2;
  xadj = adj_alloc(n+1);
  adjncy = vtx_alloc(m);
  if (do_vwgt) {
    vwgt = wgt_alloc(n);
  }
  if (do_adjwgt) {
    adjwgt = wgt_alloc(m);
  }

  xadj[0] = 0;
  ln = j = 0;
  /* scan through the file */
  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
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
    /* vertex weight */
    if (do_vwgt) {
      weight = (wgt_t)strtod(sptr,&eptr);
      if (eptr == sptr) {
        eprintf("Missing vertex weight for vertex %zu\n",ln);
        dl_close_file(file);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        goto ERROR;
      }
      sptr = eptr;
      vwgt[ln] = weight;
    }
    i = strtoull(sptr,&eptr,10)-1;
    while (eptr != sptr) {
      sptr = eptr;
      if (i >= n) {
        eprintf("Invalid vertex %zu at line %zu\n",i+1,ln);
        dl_close_file(file);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        goto ERROR;
      }
      if (j >= m) {
        eprintf("Edge %zu found at vertex %zu is greater "
            "than specified number of edges (%zu)\n", j,ln,m);
        dl_close_file(file);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        goto ERROR;
      }
      adjncy[j] = (vtx_t)i;
      /* edge weight */ 
      if (do_adjwgt) {
        weight = (wgt_t)strtod(sptr,&eptr);
        if (eptr == sptr) {
          eprintf("Missing edge weight for vertex %zu, edge %zu\n",ln,ne);
          dl_close_file(file);
          rv = BOWSTRING_ERROR_INVALIDINPUT;
          goto ERROR;
        }
        sptr = eptr;
        adjwgt[j] = weight;
      }
      ++j;
      ++ne;
      i = strtoull(sptr,&eptr,10)-1;
    }
    xadj[++ln] = (adj_t)j;
  }
  dl_close_file(file);
  file = NULL;

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

  if (r_nvtxs) {
    *r_nvtxs = n;
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
  if (r_vwgt) {
    *r_vwgt = vwgt;
  } else if (vwgt) {
    dl_free(vwgt);
  }
  if (r_adjwgt) {
    *r_adjwgt = adjwgt;
  } else if (adjwgt) {
    dl_free(adjwgt);
  }

  return BOWSTRING_SUCCESS;

  ERROR:

  dl_free(line);
  if (xadj) {
    dl_free(xadj);
  }
  if (adjncy) {
    dl_free(adjncy);
  }
  if (vwgt) {
    dl_free(vwgt);
  }
  if (adjwgt) {
    dl_free(adjwgt);
  }
  if (file) {
    dl_close_file(file);
  }

  return rv;
}


/**
 * @brief Function for reading graph in the Nerstrand Binary Graph format. The
 * format is as follows:
 *
 * | bytes | description
 * | 0-3   | NBG Identifier (0x646f6d3)
 * | 4-7   | Format version
 * | 8     | Weight bytes (0 for not presnt)
 * | 9-16  | Number of vertices (determines size of vtx_t)
 * | 17-24 | Number of edges (determines size of adj_t)
 * | 25-28 | Compression used (0 for none)
 * | 29-36 | Size of compressed xadj
 * | N     | compressed xadj
 * | N+8   | Size of compressed adjncy
 * | M     | compressed adjncy 
 * | M+8   | Size of compressed vwgt (if present)
 * | V     | vwgt
 * | V+8   | Size of compressed adjwgt (if present) 
 * | W     | adjwgt
 *
 * Someday I need to make it so that missing r_* pointers results in the field
 * being skipped.
 *
 * @param filename
 * @param r_nvtxs
 * @param r_xadj
 * @param r_adjncy
 * @param r_vwgt
 * @param r_adjwgt
 *
 * @return 
 */
static int __read_binary_graph(
    char const * const filename, 
    vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  int rv, bigedges, bigvertices;
  uint64_t nedges, nvtxs;
  uint32_t b32, compression;
  uint8_t ww;

  file_t * file = NULL;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * vwgt = NULL;
  wgt_t * adjwgt = NULL;

  if ((rv = __open_file(filename,"r",&file)) != BOWSTRING_SUCCESS) {
    goto ERROR;
  }

  /* verify that its good graph file */
  if (fread(&b32,sizeof(uint32_t),1,file->fd) != 1) {
    rv = BOWSTRING_ERROR_FILEREAD;
    goto ERROR;
  }
  if (b32 != __MAGIC_NUMBER) {
    eprintf("Not a Nerstrand Binary Graph file (bad magic number).\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }
  if (fread(&b32,sizeof(uint32_t),1,file->fd) != 1) {
    rv = BOWSTRING_ERROR_FILEREAD;
    goto ERROR;
  }
  if (b32 != __FORMAT_VERSION) {
    eprintf("Bad format version in Nerstrand Binary Graph file.\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  /* load graph information */
  if (fread(&ww,sizeof(uint8_t),1,file->fd) != 1) {
    eprintf("Premature end of Nerstrand Binary Graph file (missing weight "
        "size).\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }
  dprintf("Nerstrand Binary Graph file with weight size of %zu\n",(size_t)ww);

  if (fread(&nvtxs,sizeof(uint64_t),1,file->fd) != 1) {
    eprintf("Premature end of Nerstrand Binary Graph file (missing number of "
        "vertices).\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }
  if (fread(&nedges,sizeof(uint64_t),1,file->fd) != 1) {
    eprintf("Premature end of Nerstrand Binary Graph file (missing number of "
        "edges).\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  dprintf("Nerstrand Binary Graph file with %"PRIu64" vertices and %"PRIu64
      " edges.\n",nvtxs,nedges);

  xadj = adj_alloc(nvtxs+1);
  adjncy = vtx_alloc(nedges);

  /* determine vertex and edge size */
  if (nvtxs >= UINT32_MAX /*0xFFFFFFFFULL*/) {
    bigvertices = 1;
  } else {
    bigvertices = 0;
  }
  dprintf("Nerstrand Binary Graph file with bigvertices = %d\n",bigvertices);
  if (nedges >= UINT32_MAX /*0xFFFFFFFFULL*/) {
    bigedges = 1;
  } else {
    bigedges = 0;
  }
  dprintf("Nerstrand Binary Graph file with bigedges = %d\n",bigedges);

  /* ensure we support the compression type */
  if (fread(&compression,sizeof(uint32_t),1,file->fd) != 1) {
    eprintf("Premature end of Nerstrand Binary Graph file (missing "
        "compression type).\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto ERROR;
  }

  /* read in xadj */
  if ((rv = __adj_read(xadj,nvtxs+1,file,bigedges,compression)) !=
      BOWSTRING_SUCCESS) {
    eprintf("Failed to read xadj from Nerstrand Binary Graph file.\n");
    goto ERROR;
  }

  /* read in adjncy */
  if ((rv = __vtx_read(adjncy,nedges,file,bigvertices,compression)) !=
      BOWSTRING_SUCCESS) {
    eprintf("Failed to read adjncy from Nerstrand Binary Graph file.\n");
    goto ERROR;
  }

  /* read in vwgt */
  if ((rv = __wgt_read(vwgt,nvtxs,file,sizeof(wgt_t)==sizeof(double),
        compression)) != BOWSTRING_SUCCESS) {
    eprintf("Failed to read vwgt from Nerstrand Binary Graph file.\n");
    goto ERROR;
  }

  /* read in adjwgt */
  if ((rv = __wgt_read(adjwgt,nedges,file,sizeof(wgt_t)==sizeof(double),
        compression)) != BOWSTRING_SUCCESS) {
    eprintf("Failed to read adjwgt from Nerstrand Binary Graph file.\n");
    goto ERROR;
  }

  dl_close_file(file);

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
  if (r_vwgt) {
    *r_vwgt = vwgt;
  } else if (vwgt) {
    dl_free(vwgt);
  }
  if (r_adjwgt) {
    *r_adjwgt = adjwgt;
  } else if (adjwgt) {
    dl_free(adjwgt);
  }

  return BOWSTRING_SUCCESS;

  ERROR:

  if (xadj) {
    dl_free(xadj);
  }
  if (adjncy) {
    dl_free(adjncy);
  }
  if (vwgt) {
    dl_free(vwgt);
  }
  if (adjwgt) {
    dl_free(adjwgt);
  }
  if (file) {
    dl_close_file(file);
  }

  return rv;
}


static int __write_text_graph(
    char const * const filename, 
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    wgt_t const * const adjwgt)
{
  int rv;
  file_t * file;
  size_t i, j;

  int do_vwgt, do_adjwgt, iflags;

  int wgt_flags = 0;
  do_vwgt = do_adjwgt = 0;

  adj_t const nedges = xadj[nvtxs];

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

  /* print the header line */
  dl_fprintf(file,"%zu %zu",(size_t)nvtxs,(size_t)(nedges/2));
  iflags = __flag_to_int(wgt_flags);
  if (iflags) {
    dl_fprintf(file," %d\n",iflags);
  } else {
    dl_fprintf(file,"\n");
  }

  /* print out the graph */
  for (i=0;i<nvtxs;++i) {
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


/**
 * @brief Function for writing nerstrand graphs in binary format. The format is
 * as follows:
 *
 * | bytes | description
 * | 0-3   | NBG Identifier (0x646f6d3)
 * | 4-7   | Format version
 * | 8     | Weight bytes (0 for not presnt)
 * | 9-16  | Number of vertices (determines size of vtx_t)
 * | 17-24 | Number of edges (determines size of adj_t)
 * | 25-28 | Compression used (0 for none)
 * | 29-36 | Size of compressed xadj
 * | N     | compressed xadj
 * | N+8   | Size of compressed adjncy
 * | M     | compressed adjncy 
 * | M+8   | Size of compressed vwgt (if present)
 * | V     | vwgt
 * | V+8   | Size of compressed adjwgt (if present) 
 * | W     | adjwgt
 *
 * @param filename
 * @param nvtxs
 * @param xadj
 * @param adjncy
 * @param vwgt
 * @param adjwgt
 * @param compression
 *
 * @return 
 */
static int __write_binary_graph(
    char const * const filename, 
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    wgt_t const * const adjwgt,
    const uint32_t compression)
{
  int rv;
  int bigedges, bigvertices;
  file_t * file = NULL;
  size_t i;
  uint64_t bnvtxs, bnedges;
  uint8_t size;

  int do_vwgt, do_adjwgt;

  do_vwgt = do_adjwgt = 0;

  adj_t const nedges = xadj[nvtxs];

  /* check if we should write the weights */
  if (vwgt) {
    for (i=0;i<nvtxs;++i) {
      if (vwgt[i] != 1.0) {
        break;
      }
    }
    if (i < nvtxs) {
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
      do_adjwgt = 1;
    }
  }

  bnvtxs = nvtxs;
  bnedges = xadj[nvtxs];

  /* determine vertex and edge size */
  if (bnvtxs >= 0xFFFFFFFFULL) {
    bigvertices = 1;
  } else {
    bigvertices = 0;
  }
  dprintf("Nerstrand Binary Graph file with bigvertices = %d\n",bigvertices);
  if (bnedges >= 0xFFFFFFFFULL) {
    bigedges = 1;
  } else {
    bigedges = 0;
  }
  dprintf("Nerstrand Binary Graph file with bigedges = %d\n",bigedges);

  if ((rv = __open_file(filename,"w",&file)) != BOWSTRING_SUCCESS) {
    return rv;
  }

  /* header */
  if (fwrite(&__MAGIC_NUMBER,sizeof(uint32_t),1,file->fd) != 1) {
    eprintf("Failed to write magic number to Nerstrand Binary Graph file.\n");
    rv = BOWSTRING_ERROR_FILEWRITE; 
    goto ERROR;
  }
  if (fwrite(&__FORMAT_VERSION,sizeof(uint32_t),1,file->fd) != 1) {
    eprintf("Failed to write format version to Nerstrand Binary Graph "
        "file.\n");
    rv = BOWSTRING_ERROR_FILEWRITE; 
    goto ERROR;
  }

  /* size of weights */
  size = sizeof(wgt_t);
  if (fwrite(&size,sizeof(uint8_t),1,file->fd) != 1) {
    eprintf("Failed to write weight size to Nerstrand Binary Graph file.\n");
    rv = BOWSTRING_ERROR_FILEWRITE; 
    goto ERROR;
  }
  dprintf("Nerstrand Binary Graph file with weight size of %zu\n",
      (size_t)size);

  /* number of vertices and edges */
  if (fwrite(&bnvtxs,sizeof(uint64_t),1,file->fd) != 1) {
    eprintf("Failed to write number of vertices to Nerstrand Binary Graph "
        "file.\n");
    rv = BOWSTRING_ERROR_FILEWRITE; 
    goto ERROR;
  }
  if (fwrite(&bnedges,sizeof(uint64_t),1,file->fd) != 1) {
    eprintf("Failed to write number of edges to Nerstrand Binary Graph "
        "file.\n");
    rv = BOWSTRING_ERROR_FILEWRITE; 
    goto ERROR;
  }
  dprintf("Nerstrand Binary Graph file with %"PRIu64" vertices and %"PRIu64
      " edges.\n",bnvtxs,bnedges);

  /* compression used */
  if (fwrite(&compression,sizeof(uint32_t),1,file->fd) != 1) {
    eprintf("Failed to write compression type to Nerstrand Binary Graph "
        "file.\n");
    rv = BOWSTRING_ERROR_FILEWRITE; 
    goto ERROR;
  }

  /* xadj */
  if ((rv = __adj_write(xadj,nvtxs+1,file,bigedges,compression)) != 
      BOWSTRING_SUCCESS) {
    eprintf("Failed to write xadj to Nerstrand Binary Graph file.\n");
    goto ERROR;
  }

  /* adjncy */
  if ((rv = __vtx_write(adjncy,nedges,file,bigvertices,compression)) !=
      BOWSTRING_SUCCESS) {
    eprintf("Failed to write adjncy to Nerstrand Binary Graph file.\n");
    goto ERROR;
  }

  /* vwgt */
  if ((rv = __wgt_write(vwgt,nvtxs*do_vwgt,file,sizeof(wgt_t)==sizeof(double),
          compression)) != BOWSTRING_SUCCESS) {
    eprintf("Failed to write vwgt to Nerstrand Binary Graph file.\n");
    goto ERROR;
  }

  /* adjwgt */
  if ((rv =__wgt_write(adjwgt,nedges*do_adjwgt,file,
        sizeof(wgt_t)==sizeof(double),compression)) != BOWSTRING_SUCCESS) {
    eprintf("Failed to write adjwgt to Nerstrand Binary Graph file.\n");
    goto ERROR;
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




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int read_nerstrand_graph(
    char const * const filename, 
    vtx_t * const r_nvtxs, 
    adj_t ** const xadj,
    vtx_t ** const adjncy, 
    wgt_t ** const vwgt, 
    wgt_t ** const adjwgt, 
    int const binary)
{
  if (binary) {
    return __read_binary_graph(filename,r_nvtxs,xadj,adjncy,vwgt,adjwgt);
  } else {
    return __read_text_graph(filename,r_nvtxs,xadj,adjncy,vwgt,adjwgt);
  }
}


int write_nerstrand_graph(
    char const * const filename, 
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    wgt_t const * const adjwgt, 
    int const binary)
{
  if (binary) {
    return __write_binary_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt,
        BOWSTRING_COMPRESSION_GZIP6);
  } else {
    return __write_text_graph(filename,nvtxs,xadj,adjncy,vwgt,adjwgt);
  }
}


int read_nerstrand_clustering(
    char const * const filename, 
    vtx_t * const nvtxs, 
    vlbl_t ** const labels,
    int const binary)
{
  if (binary) {
    return BOWSTRING_ERROR_UNIMPLEMENTED;
  } else {
    return read_vertex_labels(filename,nvtxs,labels);
  }
}


int write_nerstrand_clustering(
    char const * const filename, 
    vtx_t const nvtxs, 
    vlbl_t const * const labels,
    int const binary)
{
  if (binary) {
    return BOWSTRING_ERROR_UNIMPLEMENTED;
  } else {
    return write_vertex_labels(filename,nvtxs,labels);
  }
}




#endif

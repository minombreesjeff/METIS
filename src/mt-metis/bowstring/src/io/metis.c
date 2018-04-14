/**
 * @file metis.c
 * @brief Functions for reading and writing metis graphs and partition files
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2014-02-13
 */




#ifndef BOWSTRING_IO_METIS_C
#define BOWSTRING_IO_METIS_C




#include "metis.h"




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


int read_metis_graph(
    const char * const filename, 
    vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  file_t * file;
  int w, rv;
  size_t n,m,i, j, ln, ne, num_header;
  ssize_t ll;
  size_t bufsize;
  char * line = NULL, *eptr, * sptr;
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
    if (ln >= n) {
      eprintf("Greater than the specificied %zu vertices in header found.\n",
          n);
      rv = BOWSTRING_ERROR_INVALIDINPUT;
      goto ERROR;
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

  if (line) {
    dl_free(line);
  }
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
  return rv;
}


int write_metis_graph(
    const char * const filename, 
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const vwgt, 
    const wgt_t * const adjwgt)
{
  int rv;
  file_t * file;
  vtx_t i, k;
  adj_t j, nedges;

  int do_vwgt, do_adjwgt, iflags;

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

  /* filter self loops */
  nedges = xadj[nvtxs];
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k == i) {
        --nedges;
      }
    }
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
      dl_fprintf(file,"%zu ",(size_t)vwgt[i]);
    }
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k != i) {
        dl_fprintf(file,"%zu ",(size_t)(k+1));
        if (do_adjwgt) {
          dl_fprintf(file,"%zu ",(size_t)k);
        }
      }
    }
    dl_fprintf(file,"\n");
  }
  dl_close_file(file);

  return BOWSTRING_SUCCESS;
}


int read_metis_partition(
    char const * const filename, 
    vtx_t * const r_nvtxs, 
    vlbl_t ** const r_labels)
{
  file_t * file;
  vtx_t nl;
  vlbl_t l;
  int rv;
  ssize_t ll;
  size_t bufsize;
  char * line;

  vlbl_t * labels = NULL;

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
  if (!r_nvtxs || *r_nvtxs == 0) {
    /* determine the size of the file */
    nl = 0;
    while (ll > 0 && sscanf(line,RF_VLBL_T,&l) == 1) {
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
  labels = vlbl_alloc(nl);
  nl = 0;
  /* read in the file */
  while (ll > 0 && sscanf(line,RF_VLBL_T,&l) == 1) {
    labels[nl++] = l;
    ll = dl_get_next_line(file,&line,&bufsize);
  }
  if (r_nvtxs) {
    *r_nvtxs = nl;
  }

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


int write_metis_partition(const char * const filename, const vtx_t nvtxs, 
    const vlbl_t * const labels)
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


#endif

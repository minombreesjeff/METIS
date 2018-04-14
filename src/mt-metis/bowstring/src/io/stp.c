/**
 * @file stp.c
 * @brief Functions for reading and writing STP files
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2014-02-13
 */




#ifndef BOWSTRING_STP_C
#define BOWSTRING_STP_C




#include "stp.h"




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define __SECTION_STR "SECTION"
#define __STP_HEADER "33D32945 STP File, STP Format Version 1.0"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const BUFFERSIZE = 0x1000;




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


/**
 * @brief Move the current curser in the file to a specific section.
 *
 * @param file The file to move the cursor of.
 * @param section The name of the section to move to.
 * @param r_bufsize A reference to the current/resulting buffer size.
 *
 * @return 1 if the section is found.
 */
static int __find_section(
    file_t * file, 
    char const * const section, 
    size_t * const r_bufsize)
{
  int found = 0;
  ssize_t ll;
  char * line = NULL;
  size_t bufsize;
  char str[512];

  bufsize = *r_bufsize;

  while((ll = dl_get_next_line(file,&line,&bufsize)) >= 0) {
    dl_string_upperize(line);
    if (sscanf(line,__SECTION_STR"%*[ \t]%512s",str) != 1) {
      continue;
    } else {
      if (strcmp(section,str) == 0) {
        /* found the section */
        found = 1;
        break;
      } else {
        /* some other section */
      }
    }
  }
  dl_free(line);

  *r_bufsize = bufsize;

  return found;
}



/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int read_stp_graph(
    char const * const filename, 
    vtx_t * const r_nvtxs, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_vwgt, 
    wgt_t ** const r_adjwgt)
{
  int rv;
  vtx_t i, k, nvtxs;
  adj_t nedges, j, l;
  wgt_t w;
  ssize_t ll;
  char * line = NULL;
  size_t bufsize = BUFFERSIZE;
  file_t * ifile = NULL;

  vtx_t * x = NULL;
  vtx_t * y = NULL;
  wgt_t * z = NULL;

  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL;
  wgt_t * vwgt = NULL;


  if ((rv = __open_file(filename,"r",&ifile)) != BOWSTRING_SUCCESS) {
    goto END;
  }
  line = char_alloc(bufsize);

  /* find the graph section */
  if (!__find_section(ifile,"GRAPH",&bufsize)) {
    eprintf("Could not find graph section\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT; 
    goto END;
  }

  /* read the number of vertices/nodes */
  if ((ll = dl_get_next_line(ifile,&line,&bufsize)) < 0 || 
      sscanf(line,"Nodes%*[ \t]"RF_VTX_T,&nvtxs) != 1) {
    eprintf("Could not read number of nodes from graph section\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT; 
    goto END;
  }

  /* read the number of edges */
  if ((ll = dl_get_next_line(ifile,&line,&bufsize)) < 0 || 
      sscanf(line,"Edges%*[ \t]"RF_ADJ_T,&nedges) != 1) {
    eprintf("Could not read number of edges from graph section\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT; 
    goto END;
  }
  nedges *=2;

  /* allocate graph */
  xadj = adj_calloc(nvtxs+1);
  adjncy = vtx_alloc(nedges);
  adjwgt = wgt_alloc(nedges);


  /* allocate graph buffers */
  x = vtx_alloc(nedges/2);
  y = vtx_alloc(nedges/2);
  z = wgt_alloc(nedges/2);

  /* read in edges */
  for (j=0;j<nedges/2;++j) {
    while ((ll = dl_get_next_line(ifile,&line,&bufsize)) == 0) {
      /* blank line */
    }
    if (ll < 0) {
      eprintf("Could not find edge #"PF_ADJ_T": '%s'\n",j,line);
      rv = BOWSTRING_ERROR_INVALIDINPUT; 
      goto END;
    } else {
      if (dl_strncmp_nocase(line,"END",3) == 0) {
        wprintf("Only half the number of edges listed :"PF_ADJ_T"/"PF_ADJ_T \
            ".\n",j,nedges/2);
        nedges = j*2;
        goto NODEWEIGHTS;
      } else if (line[0] == 'E') { 
        if (sscanf(line,"E%*[ \t]"RF_VTX_T"%*[ \t]"RF_VTX_T"%*[ \t]"RF_WGT_T, \
            &i,&k,&w) != 3) {
          eprintf("Could not read edge #"PF_ADJ_T": '%s'\n",j,line);
          rv = BOWSTRING_ERROR_INVALIDINPUT; 
          goto END;
        }
      } else if (line[0] == 'A') {
        /* handle broken files with arcs listed instead of edges */
        if (sscanf(line,"A%*[ \t]"RF_VTX_T"%*[ \t]"RF_VTX_T"%*[ \t]"RF_WGT_T, \
            &i,&k,&w) != 3) {
          eprintf("Could not read arc #"PF_ADJ_T": '%s'\n",j,line);
          rv = BOWSTRING_ERROR_INVALIDINPUT; 
          goto END;
        }
      } else {
        eprintf("Bad label for edge #"PF_ADJ_T": '%s'\n",j,line);
        rv = BOWSTRING_ERROR_INVALIDINPUT; 
        goto END;
      }
      ++xadj[i];
      ++xadj[k];
      --i;
      --k;
      x[j] = i;
      y[j] = k;
      z[j] = w;
    }
  }

  /* skip blank lines */
  while ((ll = dl_get_next_line(ifile,&line,&bufsize)) == 0) {
    /* blank line */
  }

  /* check the end of the edge list */
  if (ll < 0 || dl_strncmp_nocase(line,"END",3) != 0) {
    eprintf("Invalid line at end of section: '%s'\n",line);
    rv = BOWSTRING_ERROR_INVALIDINPUT; 
    goto END;
  }

  NODEWEIGHTS:

  /* find the node weight section if it exists */
  if (r_vwgt && __find_section(ifile,"NODEWEIGHTS",&bufsize)) {
    vwgt = wgt_alloc(nvtxs);
    for (i=0;i<nvtxs;++i) {
      while ((ll = dl_get_next_line(ifile,&line,&bufsize)) == 0) {
        /* blank line */
      }
      if (ll < 0 || sscanf(line,"NW%*[ \t]"RF_WGT_T,&w) != 1) {
        eprintf("Could not read vertex weight for vertex "PF_VTX_T"\n",i);
        rv = BOWSTRING_ERROR_INVALIDINPUT;
        goto END;
      } else {
        vwgt[i] = w;
      }
    }

    /* skip blank lines */
    while ((ll = dl_get_next_line(ifile,&line,&bufsize)) == 0) {
      /* blank line */
    }

    /* check the end of the vertex weight list */
    if (ll < 0 || dl_strncmp_nocase(line,"END",3) != 0) {
      eprintf("Invalid line at end of section: '%s'\n",line);
      rv = BOWSTRING_ERROR_INVALIDINPUT; 
      goto END;
    }
  }

  /* close the file */
  dl_close_file(ifile);
  ifile = NULL;

  /* build edges */
  adj_prefixsum_exc(xadj+1,nvtxs);
  for (j=0;j<nedges/2;++j) {
    i = x[j];
    k = y[j];
    w = z[j];
    /* vertex a */
    l = xadj[i+1]++;
    adjncy[l] = k;
    adjwgt[l] = w;
    /* vertex b */
    l = xadj[k+1]++;
    adjncy[l] = i;
    adjwgt[l] = w;
  }

  *r_nvtxs = nvtxs;

  *r_xadj = xadj;
  xadj = NULL;

  *r_adjncy = adjncy;
  adjncy = NULL;

  if (r_adjwgt) {
    *r_adjwgt = adjwgt;
  } else { 
    dl_free(adjwgt);
  }
  adjwgt = NULL;

  if (r_vwgt) {
    *r_vwgt = vwgt;
  }
  vwgt = NULL;

  END:

  if (line) {
    dl_free(line);
  }
  if (ifile) {
    dl_close_file(ifile);
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
  if (vwgt) {
    dl_free(vwgt);
  }
  if (x) {
    dl_free(x);
  }
  if (y) {
    dl_free(y);
  }
  if (z) {
    dl_free(z);
  }

  return rv;
}


int write_stp_graph(
    char const * const filename, 
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    wgt_t const * const adjwgt)
{
  int rv;
  file_t * file;
  vtx_t i, k;
  adj_t j;


  if ((rv = __open_file(filename,"w",&file)) != BOWSTRING_SUCCESS) {
    return rv;
  }

  /* print the header line */
  dl_fprintf(file,"%s\n",__STP_HEADER);
  dl_fprintf(file,"\n");

  /* print comment section */
  dl_fprintf(file,"%s Comment\n",__SECTION_STR);
  dl_fprintf(file,"Name \"Saved Graph File\"\n");
  dl_fprintf(file,"Creator \"Bowstring v%d.%d.%d\"\n",BOWSTRING_VER_MAJOR, \
      BOWSTRING_VER_MINOR,BOWSTRING_VER_SUBMINOR);
  dl_fprintf(file,"END\n");
  dl_fprintf(file,"\n");


  /* print the graph numbers */
  dl_fprintf(file,"%s Graph\n",__SECTION_STR);
  dl_fprintf(file,"Nodes "PF_VTX_T"\n",nvtxs);
  dl_fprintf(file,"Edges "PF_ADJ_T"\n",xadj[nvtxs]/2);

  /* print out the edges */
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (i < k) {
        if (adjwgt) {
          dl_fprintf(file,"E "PF_VTX_T" "PF_VTX_T" "PF_WGT_T"\n",i+1,k+1, \
              adjwgt[j]);
        } else {
          dl_fprintf(file,"E "PF_VTX_T" "PF_VTX_T" 1\n",i+1,k+1);
        }
      }
    }
  }
  dl_fprintf(file,"END\n");
  dl_fprintf(file,"\n");

  /* print node weights */
  if (vwgt) {
    dl_fprintf(file,"%s NodeWeights\n",__SECTION_STR);
    for (i=0;i<nvtxs;++i) {
      dl_fprintf(file,"NW "PF_WGT_T"\n",vwgt[i]);
    }
    dl_fprintf(file,"END\n");
  }

  /* print footer */
  dl_fprintf(file,"\n");
  dl_fprintf(file,"EOF");

  dl_close_file(file);

  return BOWSTRING_SUCCESS;
}




#endif

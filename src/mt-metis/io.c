/**
 * @file io.c
 * @brief IO routines
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */


#include "includes.h"

idx_t ReadGraph(const char * filename, idx_t ** rxadj, idx_t ** radjncy,
    idx_t ** rvwgt, idx_t ** radjwgt)
{
  idx_t i, j,k, l, fmt, ncon, nfields, readew, readvw, readvs, edge,
      ewgt, nedges, nvtxs;
  char *line=NULL, fmtstr[256], *curstr, *newstr;
  size_t lnlen=0;
  FILE *fpin;

  if (!gk_fexists((char*)filename)) 
    errexit("File %s does not exist!\n", filename);

  fpin = gk_fopen((char*)filename, "r", "ReadGRaph: Graph");

  /* Skip comment lines until you get to the first valid line */
  do {
    if (gk_getline(&line, &lnlen, fpin) == -1) 
      errexit("Premature end of input file: file: %s\n", filename);
  } while (line[0] == '%');


  fmt = ncon = 0;
  nfields = sscanf(line, "%"SCIDX" %"SCIDX" %"SCIDX" %"SCIDX, 
                &(nvtxs), &(nedges), &fmt, &ncon);

  if (nfields < 2) 
    errexit("The input file does not specify the number of vertices and edges."
        "\n");

  if (nvtxs <= 0 || nedges <= 0) 
    errexit("The supplied nvtxs:%"PRIDX" and nedges:%"PRIDX" must be positive."
        "\n", nvtxs, nedges);
        
  if (fmt > 111) 
    errexit("Cannot read this type of file format [fmt=%"PRIDX"]!\n", fmt);

  sprintf(fmtstr, "%03"PRIDX, fmt%1000);
  readvs = (fmtstr[0] == '1');
  readvw = (fmtstr[1] == '1');
  readew = (fmtstr[2] == '1');
    
  if (ncon > 0 && !readvw) 
    errexit(
      "---------------------------------------------------------------------\n"
      "***  I detected an error in your input file  ***\n\n"
      "You specified ncon=%"PRIDX", but the fmt parameter does not specify "
      "vertex weights\n" 
      "Make sure that the fmt parameter is set to either 10 or 11.\n"
      "---------------------------------------------------------------------\n"
      , ncon);

  nedges *=2;

  idx_t * const xadj   = *rxadj = imalloc(nvtxs+1, "ReadGraph: xadj");
  idx_t * const adjncy = *radjncy = imalloc(nedges, "ReadGraph: adjncy");
  idx_t * const vwgt   = *rvwgt = ismalloc(nvtxs, 1, "ReadGraph: vwgt");
  idx_t * const adjwgt = *radjwgt = ismalloc(nedges, 1, "ReadGraph: adjwgt");
  idx_t * const vsize  = ismalloc(nvtxs, 1, "ReadGraph: vsize");

  xadj[0] = 0;


  /*----------------------------------------------------------------------
   * Read the sparse graph file
   *---------------------------------------------------------------------*/
  for (k=0, i=0; i<nvtxs; i++) {
    do {
      if (gk_getline(&line, &lnlen, fpin) == -1) 
        errexit("Premature end of input file while reading vertex %"PRIDX".\n",
            i+1);
    } while (line[0] == '%');

    curstr = line;
    newstr = NULL;

    /* Read vertex sizes */
    if (readvs) {
      vsize[i] = strtoidx(curstr, &newstr, 10);
      if (newstr == curstr)
        errexit("The line for vertex %"PRIDX" does not have vsize information"
            "\n", i+1);
      if (vsize[i] < 0)
        errexit("The size for vertex %"PRIDX" must be >= 0\n", i+1);
      curstr = newstr;
    }


    /* Read vertex weights */
    if (readvw) {
      for (l=0; l<ncon; l++) {
        vwgt[i*ncon+l] = strtoidx(curstr, &newstr, 10);
        if (newstr == curstr)
          errexit("The line for vertex %"PRIDX" does not have enough weights "
                  "for the %"PRIDX" constraints.\n", i+1, ncon);
        if (vwgt[i*ncon+l] < 0)
          errexit("The weight vertex %"PRIDX" and constraint %"PRIDX" must be "
            ">= 0\n", i+1, l);
        curstr = newstr;
      }
    }

    while (1) {
      edge = strtoidx(curstr, &newstr, 10);
      if (newstr == curstr)
        break; /* End of line */
      curstr = newstr;

      if (edge < 1 || edge > nvtxs)
        errexit("Edge %"PRIDX" for vertex %"PRIDX" is out of bounds\n", edge,
            i+1);

      ewgt = 1;
      if (readew) {
        ewgt = strtoidx(curstr, &newstr, 10);
        if (newstr == curstr)
          errexit("Premature end of line for vertex %"PRIDX"\n", i+1);
        if (ewgt <= 0)
          errexit("The weight (%"PRIDX") for edge (%"PRIDX", %"PRIDX") must be"
            "positive.\n", ewgt, i+1, edge);
        curstr = newstr;
      }

      if (k == nedges)
        errexit("There are more edges in the file than the %"PRIDX" specified."
            "\n", nedges/2); 
      adjncy[k] = edge-1;
      adjwgt[k] = ewgt;
      k++;
    } 
    xadj[i+1] = k;
  }
  gk_fclose(fpin);

  if (k != nedges) {
    printf(
"--------------------------------------------------------------------------\n"
"***  I detected an error in your input file  ***\n\n"
"In the first line of the file, you specified that the graph contained\n"
"%"PRIDX" edges. However, I only found %"PRIDX" edges in the file.\n", 
           nedges/2, k/2);

    if (2*k == nedges) {
      printf(
"\n *> I detected that you specified twice the number of edges that you have\n"
"    in the file. Remember that the number of edges specified in the first\n"
"    line counts each edge between vertices v and u only once.\n\n"
      );

    }
    printf(
"Please specify the correct number of edges in the first line of the file.\n"
"---------------------------------------------------------------------------\n"
    );
    exit(0);
  }

  ASSERT(xadj[i] == nedges);

  gk_free((void**)&vsize,LTERM);

  #ifdef HAVE_GETLINE_
  gk_free((void**)&line, LTERM);
  #else
  free(line);
  #endif

  return nvtxs;
}


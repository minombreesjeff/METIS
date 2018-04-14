/**
 * @file animikii.c
 * @brief Main driver function
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */


#include "includes.h"

static void usage(const char * name, FILE * out)
{
  fprintf(out,"USAGE:\n");
  fprintf(out,"%s <inputgraph> [outputfile | - [nparts [nthreads [seed "
      "[nruns]]]]:\n",name);
}

idx_t main(idx_t argc, char ** argv) {
  idx_t rv, i,j,n;
  idx_t seed = 1;
  idx_t ncuts = 1;
  idx_t nparts = 2;
  idx_t nthreads = omp_get_max_threads();
  idx_t * owhere = NULL;
  char * outfile = NULL;
  FILE * fout = NULL;

  if (argc < 2) {
    fprintf(stderr,"No graph specified to partition!\n");
    usage(argv[0],stderr);
    return 1;
  } else if (strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0) {
    usage(argv[0],stdout);
    return 0;
  } else if (argc > 7) {
    fprintf(stderr,"Invalid number of arguments!\n");
    usage(argv[0],stderr);
    return 2;
  }

  if (argc > 2) {
    if (strlen(argv[2]) == 1 && argv[2][0] == '-') {
      /* don't set the output file */
    } else {
      outfile = argv[2];
    }
  } 
  if (argc > 3) {
    nparts = atoi(argv[3]);
  } 
  if (argc > 4) {
    nthreads = atoi(argv[4]);
  }
  if (argc > 5) {
    seed = atoi(argv[5]);
  }
  if (argc > 6) {
    ncuts = atoi(argv[6]);
  }

  omp_set_nested(0);
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);

  real_t means[5] = {0.0};
  real_t stddevs[5] = {0.0};
  real_t max[5] = {0.0};
  real_t min[5] = {0.0};
  real_t timers[5];
  real_t delta;

  idx_t cut;
  
  printf("Reading graph from %s\n", argv[1]);
  idx_t * xadj, * vwgt, *adjncy, * adjwgt;
  idx_t nvtxs = ReadGraph(argv[1],&xadj,&adjncy,&vwgt,&adjwgt);

  if (outfile) {
    owhere = imalloc(nvtxs,"X");
  }
  
  for (i=0;i<ncuts;++i) {
    printf("**************************************************************\n");
    printf(" Using %"PRIDX" threads\n",nthreads);
    printf("**************************************************************\n");

    idx_t buffer[nthreads+1];

    dctrl_t * dctrl = AllocateDCtrl(nvtxs,nparts,nthreads); 
    dctrl->ctrl->seed = seed+i;
    dctrl->ctrl->nparts = nparts;

    dgraph_t * graph;

    #pragma omp parallel shared(cut,owhere) num_threads(nthreads)
    {
      idx_t x;
      ParSetupGraph(&graph,buffer,nvtxs,xadj,adjncy,adjwgt,vwgt);

      x = ParPartGraphKway(dctrl,graph,owhere);

      ParFreeGraph(&graph);

      #pragma omp master
      {
        cut = x;
      }
    }

    /* do stats */
    n = i + 1;
    timers[0] = dctrl->ctrl->TotalTmr;
    timers[1] = dctrl->ctrl->CoarsenTmr;
    timers[2] = dctrl->ctrl->InitPartTmr;
    timers[3] = dctrl->ctrl->UncoarsenTmr;
    timers[4] = cut;
    for (j=0;j<5;++j) {
      if (i==0) {
        max[j] = timers[j];
        min[j] = timers[j];
      } else {
        max[j] = gk_max(max[j],timers[j]);
        min[j] = gk_min(min[j],timers[j]);
      }
      delta = timers[j] - means[j];
      means[j] += (delta/n);
      stddevs[j] += delta*(timers[j]-means[j]);
    }

    FreeDCtrl(&dctrl,nthreads);
  }

  if (outfile) {
    fout = fopen(outfile,"w");
    for (i=0;i<nvtxs;++i) {
      fprintf(fout,"%"PRIDX"\n",owhere[i]);
    }
    fclose(fout);
    gk_free((void**)&owhere,LTERM);
  }

  if (ncuts > 1) {
    for (j=0;j<5;++j) {
      stddevs[j] = sqrt(stddevs[j]/ncuts);
    }
    /* print stats */
    printf("SUMMARY STATS\n");
    printf("Number of cuts = %"PRIDX"\n",ncuts);
    printf("edgecut avg = %7.5f, stdev = %7.5f, min = %7.5f, max = %7.5f\n",
        means[4],stddevs[4],min[4],max[4]);
    printf("total avg = %7.5f, stdev = %7.5f, min = %7.5f, "
        "max = %7.5f\n",means[0],stddevs[0],min[0],max[0]);
    printf("coarsening avg = %7.5f, stdev = %7.5f, min = %7.5f, "
        "max = %7.5f\n",means[1], stddevs[1],min[1],max[1]);
    printf("initpart avg = %7.5f, stdev = %7.5f, "
        "min = %7.5f, max = %7.5f\n", means[2],stddevs[2],min[2],max[2]);
    printf("uncoarsening avg = %7.5f, stdev = %7.5f, min = %7.5f, "
        "max = %7.5f\n",means[3], stddevs[3], min[3],max[3]);
  }

  gk_free((void**)&xadj,&adjncy,&vwgt,&adjwgt,LTERM);
  return 0;
}

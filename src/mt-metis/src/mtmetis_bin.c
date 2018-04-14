/**
 * @file mtmetis_bin.c
 * @brief Main driver function
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */





#ifndef MTMETIS_BIN_C
#define MTMETIS_BIN_C




#include "base.h"
#include "partition.h"
#include "order.h"
#include "graph.h"
#include "ctrl.h"
#include "internal.h"
#include <bowstring.h>




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define __ARRAY_SIZE(a) \
  (sizeof(a) > 0 ? (sizeof(a) / sizeof((a)[0])) : 0)




/******************************************************************************
* OPTIONS *********************************************************************
******************************************************************************/


static const cmd_opt_pair_t CTYPE_CHOICES[] = {
  {MTMETIS_STR_CTYPE_RM,"Random Matching",MTMETIS_CTYPE_RM},
  {MTMETIS_STR_CTYPE_SHEM,"Sorted Heavy Edge Matching",MTMETIS_CTYPE_SHEM},
  {MTMETIS_STR_CTYPE_FC,"FirstChoice Grouping",MTMETIS_CTYPE_FC}
};


static const cmd_opt_pair_t CONTYPE_CHOICES[] = {
  {MTMETIS_STR_CONTYPE_CLS,"Hash-table with linear scanning", \
      MTMETIS_CONTYPE_CLS},
  {MTMETIS_STR_CONTYPE_DENSE,"Dense vector",MTMETIS_CONTYPE_DENSE},
  {MTMETIS_STR_CONTYPE_SORT,"Sort and merge",MTMETIS_CONTYPE_SORT}
};


static const cmd_opt_pair_t RTYPE_CHOICES[] = {
  {MTMETIS_STR_RTYPE_GREEDY,"Greedy refinement",MTMETIS_RTYPE_GREEDY},
  {MTMETIS_STR_RTYPE_FM,"FM based serial refinement (pairwise for kway)", \
      MTMETIS_RTYPE_FM},
  {MTMETIS_STR_RTYPE_SFM,"Segmented FM based parallel refinement", \
      MTMETIS_RTYPE_SFM},
  {MTMETIS_STR_RTYPE_SFG,"Segmented FM plus Greedy parallel refinement", \
      MTMETIS_RTYPE_SFG},
};


static const cmd_opt_pair_t PTYPE_CHOICES[] = {
  {MTMETIS_STR_PTYPE_KWAY,"K-Way Edgecut",MTMETIS_PTYPE_KWAY},
  {MTMETIS_STR_PTYPE_ESEP,"Edge Separator",MTMETIS_PTYPE_ESEP},
  {MTMETIS_STR_PTYPE_RB,"RB Edgecut",MTMETIS_PTYPE_RB},
  {MTMETIS_STR_PTYPE_VSEP,"Vertex Separator",MTMETIS_PTYPE_VSEP},
  {MTMETIS_STR_PTYPE_ND,"Nested Dissection",MTMETIS_PTYPE_ND}
};


static const cmd_opt_pair_t VERBOSITY_CHOICES[] = {
  {MTMETIS_STR_VERBOSITY_NONE,"Do not print any runtime information.", \
      MTMETIS_VERBOSITY_NONE},
  {MTMETIS_STR_VERBOSITY_LOW,"Print only metric summary.", \
      MTMETIS_VERBOSITY_LOW},
  {MTMETIS_STR_VERBOSITY_MEDIUM,"Print summary run information.", \
      MTMETIS_VERBOSITY_MEDIUM},
  {MTMETIS_STR_VERBOSITY_HIGH,"Print verbose run information.", \
      MTMETIS_VERBOSITY_HIGH},
  {MTMETIS_STR_VERBOSITY_MAXIMUM,"Print everything.", \
      MTMETIS_VERBOSITY_MAXIMUM}
};


static const cmd_opt_pair_t DISTRIBUTION_CHOICES[] = {
  {MTMETIS_STR_DISTRIBUTION_BLOCK,"Distribute the vertices in continous " \
      "from the initial ordering.",MTMETIS_DISTRIBUTION_BLOCK},
  {MTMETIS_STR_DISTRIBUTION_CYCLIC,"Distribute the vertices in a cyclic " \
      "fashion.",MTMETIS_DISTRIBUTION_CYCLIC},
  {MTMETIS_STR_DISTRIBUTION_BLOCKCYCLIC,"Distribute the vertices in a " \
      "blockcyclic fashion.",MTMETIS_DISTRIBUTION_BLOCKCYCLIC}
};


static const cmd_opt_pair_t IGNOREWEIGHTS_CHOICES[] = {
  {MTMETIS_STR_IGNORE_NONE,"Use all weights normally", \
      MTMETIS_IGNORE_NONE},
  {MTMETIS_STR_IGNORE_VERTEXWEIGHTS,"Force all vertex weights to be one", \
      MTMETIS_IGNORE_VERTEXWEIGHTS},
  {MTMETIS_STR_IGNORE_EDGEWEIGHTS,"Force all edge weights to be one", \
      MTMETIS_IGNORE_EDGEWEIGHTS},
  {MTMETIS_STR_IGNORE_BOTH,"Force all weights to be one", \
      MTMETIS_IGNORE_VERTEXWEIGHTS | MTMETIS_IGNORE_EDGEWEIGHTS}
};


static const cmd_opt_t OPTS[] = {
  {MTMETIS_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {MTMETIS_OPTION_CTYPE,'c',"ctype","The type of coarsening.", \
      CMD_OPT_CHOICE,CTYPE_CHOICES,__ARRAY_SIZE(CTYPE_CHOICES)},
  {MTMETIS_OPTION_CONTYPE,'d',"contype","How to merge adjacency lists " \
      "during contraction.",CMD_OPT_CHOICE,CONTYPE_CHOICES, \
        __ARRAY_SIZE(CONTYPE_CHOICES)},
  {MTMETIS_OPTION_RTYPE,'r',"rtype","The type of refinement.", \
      CMD_OPT_CHOICE,RTYPE_CHOICES,__ARRAY_SIZE(RTYPE_CHOICES)},
  {MTMETIS_OPTION_SEED,'s',"seed","The random seed to use.",CMD_OPT_INT,NULL, \
      0},
  {MTMETIS_OPTION_NCUTS,'N',"cuts","The number of cuts to " \
      "generate using succesive random seeds.",CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_NRUNS,'n',"runs","The number of partitionings to " \
      "generate using succesive random seeds.",CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_NINITSOLUTIONS,'i',"initialcuts","The number of " \
      "initial cuts to generate at the coarsest level.",CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_NITER,'R',"nrefpass","The maximum number of refinement " \
      "passes.",CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_TIME,'t',"times","Print timing information",CMD_OPT_FLAG, \
      NULL,0},
  {MTMETIS_OPTION_NTHREADS,'T',"threads","The number of threads to use.", \
      CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_VERBOSITY,'v',"verbosity","The amount of information to " \
      "print during partitioning.",CMD_OPT_CHOICE,VERBOSITY_CHOICES, \
      __ARRAY_SIZE(VERBOSITY_CHOICES)},
  {MTMETIS_OPTION_DISTRIBUTION,'D',"distribution","The distribution to use " \
      "for assigning vertices to threads.",CMD_OPT_CHOICE, \
      DISTRIBUTION_CHOICES,__ARRAY_SIZE(DISTRIBUTION_CHOICES)},
  {MTMETIS_OPTION_UBFACTOR,'b',"balance","The balance constraint (1.03 " \
      "means allowing for a 3% imbalance).",CMD_OPT_FLOAT,NULL,0},
  {MTMETIS_OPTION_PTYPE,'p',"ptype","The type of partition to compute", \
      CMD_OPT_CHOICE,PTYPE_CHOICES,__ARRAY_SIZE(PTYPE_CHOICES)},
  {MTMETIS_OPTION_RUNSTATS,'C',"partstats","Statics on quality of " \
      "partitions",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_METIS,'M',"metis","When run with one thread, call Metis " \
      "directly.",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_LEAFMATCH,'L',"leafmatch","Match leaf vertices together " \
      "if there are too many unmatched vertices.",CMD_OPT_BOOL,NULL,0},
  {MTMETIS_OPTION_REMOVEISLANDS,'I',"removeislands","Remove island vertices " \
      "before partitioning.",CMD_OPT_BOOL,NULL,0},
  {MTMETIS_OPTION_VWGTDEGREE,'V',"vwgtdegree","Use the degree of each " \
      "vertex as its weight.",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_IGNORE,'W',"ignoreweights","Ignore input weights " \
      "on a graph file.",CMD_OPT_CHOICE,IGNOREWEIGHTS_CHOICES, \
      __ARRAY_SIZE(IGNOREWEIGHTS_CHOICES)}
};


static const size_t NOPTS = __ARRAY_SIZE(OPTS);


#undef __ARRAY_SIZE




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static char const * __bool2str(
    int const b)
{
  if (b) {
    return "Yes";
  } else {
    return "No";
  }
}


static int __usage(
    char const * const name)
{
  fprintf(stderr,"USAGE:\n");
  fprintf(stderr,"%s [options] <graphfile> <nparts> [ <partfile> | - ]\n", \
      name);
  fprintf(stderr,"\n");
  fprintf(stderr,"Options:\n");
  fprint_cmd_opts(stderr,OPTS,NOPTS);
  return 1;
}


static double * __parse_args(
    cmd_arg_t * args, 
    size_t nargs,
    char const ** r_input, 
    char const ** r_output)
{
  size_t i, xarg;
  double * options = NULL;
  const char * input_file = NULL, * output_file = NULL;

  options = mtmetis_init_options();

  /* set default verbosity to low */
  options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_LOW;
  options[MTMETIS_OPTION_NPARTS] = 2.0;
  options[MTMETIS_OPTION_PTYPE] = MTMETIS_PTYPE_KWAY;

  for (i=0;i<nargs;++i) {
    switch (args[i].type) {
      case CMD_OPT_CHOICE:
        options[args[i].id] = (double)args[i].val.o;
        break;
      case CMD_OPT_BOOL:
        options[args[i].id] = (double)args[i].val.b;
        break;
      case CMD_OPT_INT:
        options[args[i].id] = (double)args[i].val.i;
        break;
      case CMD_OPT_FLOAT:
        options[args[i].id] = (double)args[i].val.f;
        break;
      case CMD_OPT_FLAG:
        options[args[i].id] = 1.0;
        break;
      default:
        break;
    }
  }

  xarg = 0;
  for (i=0;i<nargs;++i) {
    /* check for help */
    if (args[i].id == MTMETIS_OPTION_HELP) {
      goto CLEANUP;
    }
    if (args[i].type == CMD_OPT_XARG) {
      if (xarg == 0) {
        input_file = args[i].val.s;
      } else {
        if (options[MTMETIS_OPTION_PTYPE] == MTMETIS_PTYPE_KWAY || \
            options[MTMETIS_OPTION_PTYPE] == MTMETIS_PTYPE_RB) {
          if (xarg == 1) {
            options[MTMETIS_OPTION_NPARTS] = (pid_t)atoll(args[i].val.s);
          } else if (xarg == 2) {
            output_file = args[i].val.s;
            if (strcmp(output_file,"-") == 0) {
              /* if we are going to print to stdout, don't print anything else */
              options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_NONE;
            }
          } else {
            eprintf("Unknown extra argument '%s'\n",args[i].val.s);
            goto CLEANUP;
          }
        } else {
          if (xarg == 1) {
            output_file = args[i].val.s;
            if (strcmp(output_file,"-") == 0) {
              /* if we are going to print to stdout, don't print anything else */
              options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_NONE;
            }
          } else {
            eprintf("Unknown extra argument '%s'\n",args[i].val.s);
            goto CLEANUP;
          }
        }
      }
      ++xarg;
    }
  }

  if (input_file == NULL) {
    eprintf("Must supply at least an input graph to partition\n");
    goto CLEANUP;
  }

  *r_output = output_file;
  *r_input = input_file;

  return options;

  CLEANUP:
  dl_free(options);
  *r_output = NULL;
  *r_input = NULL;

  return NULL;
}




/******************************************************************************
* MAIN ************************************************************************
******************************************************************************/


int main(
    int argc, 
    char ** argv) 
{
  int rv;
  size_t nargs;
  vtx_t nvtxs, i;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * vwgt = NULL, * adjwgt = NULL;
  ctrl_t * ctrl = NULL;
  double * options = NULL;
  cmd_arg_t * args = NULL;
  pid_t * owhere = NULL;
  char const * output_file = NULL, * input_file = NULL;
  graph_t * graph = NULL;
  timers_t * timers;

  /* parse user specified options */
  rv = cmd_parse_args(argc-1,argv+1,OPTS,NOPTS,&args,&nargs);
  if (rv != DL_CMDLINE_SUCCESS) {
    __usage(argv[0]);
    rv = 1;
    goto CLEANUP;
  }
  options = __parse_args(args,nargs,&input_file,&output_file);
  if (options == NULL) {
    __usage(argv[0]);
    rv = 2;
    goto CLEANUP;
  }
  if (ctrl_parse(options,&ctrl) != MTMETIS_SUCCESS) {
    __usage(argv[0]);
    rv = 3;
    goto CLEANUP;
  }
  dl_free(options);
  options = NULL;

  timers = &(ctrl->timers);

  /* start timers */
  dl_start_timer(&timers->total);
  dl_start_timer(&timers->io);

  if (ctrl->verbosity >= MTMETIS_VERBOSITY_LOW) {
    dl_print_header("PARAMETERS",'%');
    printf("Number of Threads: %"PF_TID_T" | Verbosity: %s\n",ctrl->nthreads, \
        trans_verbosity_string(ctrl->verbosity));
    printf("Number of Runs: %zu | Random Seed: %u\n",ctrl->nruns, \
        ctrl->seed);
    printf("Number of Partitions: %"PF_PID_T" | Partition Type: %s\n", \
        ctrl->nparts,trans_ptype_string(ctrl->ptype));
    printf("Coarsening Type: %s | Contraction Type: %s\n", \
        trans_ctype_string(ctrl->ctype),trans_contype_string(ctrl->contype));
    printf("Refinement Type: %s | Number of Refinement Passes: %zu\n",
        trans_rtype_string(ctrl->rtype),ctrl->nrefpass);
    printf("Balance: %0.2lf\n",ctrl->ubfactor);
    printf("Leaf-Matching: %s | Remove Islands: %s\n", \
        __bool2str(ctrl->leafmatch),__bool2str(ctrl->removeislands));
    dl_print_footer('%');
  }

  /* read the input graph */
  vprintf(ctrl->verbosity,MTMETIS_VERBOSITY_LOW,"Reading '%s'\n", \
      input_file);
  rv = bowstring_read_graph(input_file,BOWSTRING_FORMAT_AUTO,&nvtxs,&xadj,
      &adjncy,(bowstring_wgt_t**)&vwgt,(bowstring_wgt_t**)&adjwgt);
  if (rv != BOWSTRING_SUCCESS) {
    eprintf("Error reading from input file '%s'\n",input_file);
    rv = 4;
    goto CLEANUP;
  }

  ctrl_setup(ctrl,NULL,nvtxs);

  dl_stop_timer(&timers->io);
  
  if (output_file) {
    owhere = pid_alloc(nvtxs);
  }

  if (ctrl->ignore & MTMETIS_IGNORE_EDGEWEIGHTS) {
    if (adjwgt) {
      dl_free(adjwgt);
      adjwgt = NULL;
    }
  }
  if (ctrl->ignore & MTMETIS_IGNORE_VERTEXWEIGHTS) {
    if (vwgt) {
      dl_free(vwgt);
      vwgt = NULL;
    }
  }

  if (ctrl->vwgtdegree) {
    if (vwgt) {
      /* discard file vertex weights */
      dl_free(vwgt);
    }
    vwgt = wgt_alloc(nvtxs);
    for (i=0;i<nvtxs;++i) {
      vwgt[i] = xadj[i+1] - xadj[i];
    }
  }

  graph = graph_distribute(MTMETIS_DISTRIBUTION_BLOCKCYCLIC,nvtxs,xadj, \
      adjncy,vwgt,adjwgt,ctrl->nthreads);

  switch (ctrl->ptype) {
    case MTMETIS_PTYPE_RB:
      mtmetis_partition_rb_int(ctrl,graph,owhere);
      break;
    case MTMETIS_PTYPE_KWAY:
      mtmetis_partition_kway_int(ctrl,graph,owhere);
      break;
    case MTMETIS_PTYPE_ESEP:
      mtmetis_partition_esep_int(ctrl,graph,owhere);
      break;
    case MTMETIS_PTYPE_VSEP:
      mtmetis_partition_vsep_int(ctrl,graph,owhere);
      break;
    case MTMETIS_PTYPE_ND:
      order_nd(ctrl,graph,owhere);
      break;
    default:
      dl_error("Unknown ptype '%d'",ctrl->ptype);
  }

  dl_start_timer(&timers->io);

  if (output_file) {
    if (strcmp(output_file,"-") == 0) {
      /* write to stdout */
      for (i=0;i<graph->nvtxs;++i) {
        printf("%"PF_PID_T"\n",owhere[i]);
      }
    } else {
      /* save to file */
      bowstring_write_vertex_labels(output_file,nvtxs,owhere);
    }
  }

  dl_stop_timer(&timers->io);

  if (ctrl->runstats) {
    dl_print_header("STATISTICS",'&');
    printf("Best Objective: %"PF_WGT_T"\n",wgt_min_value(ctrl->runs, \
          ctrl->nruns));
    printf("Worst Objective: %"PF_WGT_T"\n",wgt_max_value(ctrl->runs, \
          ctrl->nruns));
    printf("Median Objective: %"PF_WGT_T"\n",wgt_median(ctrl->runs, \
          ctrl->nruns));
    printf("Mean Objective - Geo.: %0.2lf - Ari.: %.2lf\n", \
        wgt_geometric_mean(ctrl->runs,ctrl->nruns), \
        wgt_arithmetic_mean(ctrl->runs,ctrl->nruns));
    dl_print_footer('&');
  }

  if (ctrl->time && ctrl->verbosity > MTMETIS_VERBOSITY_NONE) {
    dl_print_header("TIMING",'$');
    printf("Total Time: %.05fs\n",dl_poll_timer(&(timers->total)));
    printf("\tIO: %.05fs\n",dl_poll_timer(&(timers->io)));
    printf("\tPreprocessing: %.05fs\n",dl_poll_timer(&(timers->preprocess)));
    if (ctrl->ptype == MTMETIS_PTYPE_ND) {
      printf("\tOrdering: %.05fs\n",dl_poll_timer(&(timers->ordering)));
    }
    printf("\tPartitioning: %.05fs\n",dl_poll_timer(&(timers->partitioning)));
    printf("\t\tCoarsening: %.05fs\n",dl_poll_timer(&(timers->coarsening)));
    printf("\t\t\tMatching: %.05fs\n",dl_poll_timer(&(timers->matching)));
    printf("\t\t\tContraction: %.05fs\n", \
        dl_poll_timer(&(timers->contraction)));
    printf("\t\tInitial Partitioning: %.05fs\n",
        dl_poll_timer(&(timers->initpart))); \
    printf("\t\tUncoarsening: %.05fs\n", \
        dl_poll_timer(&(timers->uncoarsening)));
    printf("\t\t\tProjection: %.05fs\n",dl_poll_timer(&(timers->projection)));
    printf("\t\t\tRefinement: %.05fs\n",dl_poll_timer(&(timers->refinement)));
    if (ctrl->ptype == MTMETIS_PTYPE_ND || ctrl->ptype == MTMETIS_PTYPE_RB) {
      printf("\t\tRecursion Overhead: %.05fs\n", \
          dl_poll_timer(&(timers->recursion)));
    }
    if (ctrl->ptype == MTMETIS_PTYPE_ND || ctrl->ptype == MTMETIS_PTYPE_RB || \
        ctrl->metis_serial) {
      printf("\tMetis: %.05fs\n",dl_poll_timer(&(timers->metis)));
    }
    dl_print_footer('$');
  }

  CLEANUP:

  if (ctrl) {
    ctrl_free(ctrl);
  }
  if (graph) {
    graph_free(graph);
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
  if (owhere) {
    dl_free(owhere);
  }
  if (args) {
    dl_free(args);
  }

  return 0;
}




#endif

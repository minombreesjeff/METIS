/**
 * @file bowstring_bin.c
 * @brief Command line interface for bowstring
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-04-09
 */




#ifndef BOWSTRING_C
#define BOWSTRING_C




#include "base.h"
#include "order.h"
#include "tree.h"
#include "graph.h"
#include "sparsen.h"
#include "analyze.h"
#include "coordinates.h"
#include "cut.h"
#include "generate.h"
#include "flow.h"
#include "io/io.h"





/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define ARRAY_SIZE(a) \
  (sizeof(a) > 0 ?  (sizeof(a) / sizeof((a)[0])) : 0)




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


/* COMMANDS ******************************************************************/

typedef enum command_t {
  COMMAND_HELP,
  COMMAND_CONVERT,
  COMMAND_ORDER,
  COMMAND_ANALYSIS,
  COMMAND_SPARSEN,
  COMMAND_COORDINATES,
  COMMAND_FLOW,
  COMMAND_GENERATE
} command_t;


/* CONVERT ****************************************************************/

typedef enum convert_option_t {
  CONVERT_OPTION_HELP,
  CONVERT_OPTION_INFILE,
  CONVERT_OPTION_OUTFILE,
  CONVERT_OPTION_INFORMAT,
  CONVERT_OPTION_OUTFORMAT,
  CONVERT_OPTION_DROPWEIGHTS
} convert_option_t;


/* ORDER *********************************************************************/

typedef enum order_option_t {
  ORDER_OPTION_HELP,
  ORDER_OPTION_INFILE,
  ORDER_OPTION_OUTFILE,
  ORDER_OPTION_INFORMAT,
  ORDER_OPTION_OUTFORMAT,
  ORDER_OPTION_SEED,
  ORDER_OPTION_TYPE
} order_option_t;


/* ANALYSIS ******************************************************************/

typedef enum analysis_t {
  ANALYSIS_GRAPHSTATS,
  ANALYSIS_PARTSTATS,
  ANALYSIS_BNDSIM,
  ANALYSIS_STATS,
  ANALYSIS_DEGDIST,
  ANALYSIS_NHDEGDIST,
  ANALYSIS_NTRIANGLES,
  ANALYSIS_ATOMICCYCLES,
  ANALYSIS_STARS,
  ANALYSIS_SUPERSTARS
} analysis_t;


typedef enum analysis_option_t {
  ANALYSIS_OPTION_HELP,
  ANALYSIS_OPTION_INFILE,
  ANALYSIS_OPTION_INFORMAT,
  ANALYSIS_OPTION_PARTFILE,
  ANALYSIS_OPTION_TYPE
} analysis_option_t;


/* SPARSEN *******************************************************************/

typedef enum sparsen_t {
  SPARSEN_DFS,
  SPARSEN_BFS,
  SPARSEN_MST,
  SPARSEN_NIR,
  SPARSEN_MSR,
  SPARSEN_ASR
} sparsen_t;


typedef enum sparsen_option_t {
  SPARSEN_OPTION_HELP,
  SPARSEN_OPTION_TYPE,
  SPARSEN_OPTION_REWEIGHT,
  SPARSEN_OPTION_SEED,
  SPARSEN_OPTION_INFILE,
  SPARSEN_OPTION_INFORMAT,
  SPARSEN_OPTION_OUTFILE,
  SPARSEN_OPTION_OUTFORMAT,
  SPARSEN_OPTION_FRACTION
} sparsen_option_t;


/* COORDINATES ***************************************************************/

typedef enum coordinates_t {
  COORDINATES_BFS
} coordinates_t;


typedef enum coordinates_option_t {
  COORDINATES_OPTION_HELP,
  COORDINATES_OPTION_TYPE,
  COORDINATES_OPTION_SEED,
  COORDINATES_OPTION_INFILE,
  COORDINATES_OPTION_INFORMAT,
  COORDINATES_OPTION_OUTFILE,
  COORDINATES_OPTION_OUTFORMAT,
  COORDINATES_OPTION_NDIM
} coordinates_option_t;


/* GENERATE ******************************************************************/

typedef enum generate_t {
  GENERATE_COMPLETE,
  GENERATE_HYPERCUBE,
  GENERATE_GRID,
  GENERATE_DUMBBELL,
  GENERATE_CAVEMAN
} generate_t;


typedef enum generate_option_t {
  GENERATE_OPTION_HELP,
  GENERATE_OPTION_TYPE,
  GENERATE_OPTION_OUTFILE,
  GENERATE_OPTION_OUTFORMAT,
  GENERATE_OPTION_PARAMS
} generate_option_t;


/* FLOW **********************************************************************/


typedef enum flow_t {
  FLOW_EDGE,
  FLOW_VERTEX
} flow_t;


typedef enum flow_option_t {
  FLOW_OPTION_HELP,
  FLOW_OPTION_TYPE,
  FLOW_OPTION_INFILE,
  FLOW_OPTION_INFORMAT,
  FLOW_OPTION_OUTFILE,
  FLOW_OPTION_OUTFORMAT,
  FLOW_OPTION_SRC,
  FLOW_OPTION_DST,
  FLOW_OPTION_TIME
} flow_option_t;




/******************************************************************************
* OPTION ARRAYS ***************************************************************
******************************************************************************/


/* COMMANDS ******************************************************************/

static const cmd_opt_pair_t COMMANDS[] = {
  [COMMAND_HELP] = {"help","Display list of available commands",COMMAND_HELP},
  [COMMAND_CONVERT] = {"convert","Convert a graph from one format to another",
      COMMAND_CONVERT},
  [COMMAND_ORDER] = {"order","Order a graph by a given method.",COMMAND_ORDER},
  [COMMAND_ANALYSIS] = {"analyze","Analyze a graph and/or partitions of a "
      "graph",COMMAND_ANALYSIS},
  [COMMAND_SPARSEN] = {"sparsen","Remove edges from a graph",COMMAND_SPARSEN},
  [COMMAND_COORDINATES] = {"coordinates","Embbed the vertices in a coordinate "
      "space",COMMAND_COORDINATES},
  [COMMAND_GENERATE] = {"generate","Generate a synthetic graph",
      COMMAND_GENERATE},
  [COMMAND_FLOW] = {"flow","Generate a flow on a graph given a source and " \
      "destination vertex.",COMMAND_FLOW}
};
static const size_t NCOMMANDS = ARRAY_SIZE(COMMANDS);


/* CONVERT ****************************************************************/

static const cmd_opt_pair_t FORMAT[] = {
  {"auto","Determine the format based on the extension.",
      BOWSTRING_FORMAT_AUTO},
  {"metis","The metis graph format.",BOWSTRING_FORMAT_METIS},
  {"dimacs","The dimacs graph format.",BOWSTRING_FORMAT_DIMACS},
  {"cloud9","The cloud9 graph format.",BOWSTRING_FORMAT_CLOUD9},
  {"snap","The Stanford Network Collection format (snap).",
      BOWSTRING_FORMAT_SNAP},
  {"csr","The compressed sparse row graph format.",BOWSTRING_FORMAT_CSR},
  {"pegasus","The pegusus graph format.",BOWSTRING_FORMAT_PEGASUS},
  {"matrixmarket","The Matrix-Market graph format.",
      BOWSTRING_FORMAT_MATRIXMARKET},
  {"edgelist","String based edge list format.",BOWSTRING_FORMAT_EDGELIST},
  {"nerstrand","Identical to metis graph format.",BOWSTRING_FORMAT_NERSTRAND},
  {"nbg","Nerstrand binary graph format.",BOWSTRING_FORMAT_NBG}
};


static const cmd_opt_t CONVERT_OPTS[] = {
  {CONVERT_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG,
    NULL,0},
  {CONVERT_OPTION_INFILE,'i',"infile","The input graph file.",
    CMD_OPT_STRING,NULL,0},
  {CONVERT_OPTION_OUTFILE,'o',"outfile","The output graph file.",
      CMD_OPT_STRING,NULL,0},
  {CONVERT_OPTION_INFORMAT,'I',"informat","The input graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {CONVERT_OPTION_OUTFORMAT,'O',"outformat","The output graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {CONVERT_OPTION_DROPWEIGHTS,'d',"dropweights","Drop weights (vertex and " \
    "edge from the graph).",CMD_OPT_FLAG,NULL,0}
};
static const size_t NCONVERT_OPTS = ARRAY_SIZE(CONVERT_OPTS);


/* ORDER *********************************************************************/

static const cmd_opt_pair_t ORDER[] = {
  {"random","Order a graph randomly",BOWSTRING_ORDER_RANDOM},
  {"bfs","Order a graph using a BFS search",BOWSTRING_ORDER_BFS},
  {"dfs","Order a graph using a DFS search",BOWSTRING_ORDER_DFS},
  {"post","Order a graph/tree using a post-order DFS search", \
      BOWSTRING_ORDER_POST},
  {"rcm","Order a graph using Reverse CuthillMckee",BOWSTRING_ORDER_RCM}
};


static const cmd_opt_t ORDER_OPTS[] = {
  {ORDER_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG,
      NULL,0},
  {ORDER_OPTION_INFILE,'i',"infile","The input graph.",CMD_OPT_STRING,
      NULL,0},
  {ORDER_OPTION_OUTFILE,'o',"outfile","The output graph.",CMD_OPT_STRING,
      NULL,0},
  {ORDER_OPTION_INFORMAT,'I',"informat","The input graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {ORDER_OPTION_INFORMAT,'O',"outformat","The output graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {ORDER_OPTION_SEED,'s',"seed","The seed to use for randomness.",CMD_OPT_INT,
      NULL,0},
  {ORDER_OPTION_TYPE,'t',"type","The ordering type to use.",CMD_OPT_CHOICE,
      ORDER,ARRAY_SIZE(ORDER)}
};
static const size_t NORDER_OPTS = ARRAY_SIZE(ORDER_OPTS);


/* ANALYSIS ******************************************************************/

static const cmd_opt_pair_t ANALYSIS[] = {
  {"graphstats","Calculate statistics of a graph",
      ANALYSIS_GRAPHSTATS},
  {"partstats","Calculate statistics of a graph partitioning",
      ANALYSIS_PARTSTATS},
  {"boundary","Calculate the boundary similarity of two partitions",
      ANALYSIS_BNDSIM},
  {"degdist","Calculate the degree distribution of the graph",
      ANALYSIS_DEGDIST},
  {"nhopdegdist","Calculate the n-hop degree distribution of the graph",
      ANALYSIS_NHDEGDIST},
  {"ntriangles","Calculate the number of triangles in the of the graph",
      ANALYSIS_NTRIANGLES},
  {"atomiccycles","Histogram the atomic cycles in the graph",
      ANALYSIS_ATOMICCYCLES},
  {"stars","Histogram of the size of star-subgraphs",ANALYSIS_STARS},
/*  {"superstars","Histogram of the size of superstar-subgraphs",
      ANALYSIS_SUPERSTARS}, */
};


static const cmd_opt_t ANALYSIS_OPTS[] = {
  {ANALYSIS_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG,NULL,
      0},
  {ANALYSIS_OPTION_INFILE,'i',"infile","The input graph.",CMD_OPT_STRING,
      NULL,0},
  {ANALYSIS_OPTION_INFORMAT,'I',"informat","The input graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {ANALYSIS_OPTION_PARTFILE,'p',"partfile","A file containing the partition.",
      CMD_OPT_STRING,NULL,0},
  {ANALYSIS_OPTION_TYPE,'t',"type","The type of analysis to perform.",
      CMD_OPT_CHOICE,ANALYSIS,ARRAY_SIZE(ANALYSIS)}
};
static const size_t NANALYSIS_OPTS = ARRAY_SIZE(ANALYSIS_OPTS);


/* SPARSEN *******************************************************************/

static const cmd_opt_pair_t SPARSEN[] = {
  {"dfs","Create a DFS from a graph",SPARSEN_DFS},
  {"bfs","Create a BFS from a graph",SPARSEN_BFS},
  {"mst","Create an MST from a graph",SPARSEN_MST},
  {"nir","Create an NIRank from a graph",SPARSEN_NIR},
  {"msr","Create an MSRank from a graph",SPARSEN_MSR},
  {"asr","Create an ASRank from a graph",SPARSEN_ASR}
};


static const cmd_opt_pair_t REWEIGHT[] = {
  {"none","Do not re-distribute the removed edge-weight.", \
      BOWSTRING_REWEIGHT_NONE},
  {"exact","Reweight using a BFS to find shortest paths.", \
      BOWSTRING_REWEIGHT_EXACT},
  {"approx","Reweight using DFS to approximate the shortest paths.", \
      BOWSTRING_REWEIGHT_APPROX}
};


static const cmd_opt_t SPARSEN_OPTS[] = {
  {SPARSEN_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG,NULL,
      0},
  {SPARSEN_OPTION_INFILE,'i',"infile","The input graph file.",CMD_OPT_STRING,
      NULL,0},
  {SPARSEN_OPTION_OUTFILE,'o',"outfile","The output graph file.",
      CMD_OPT_STRING,NULL,0},
  {SPARSEN_OPTION_INFORMAT,'I',"informat","The input graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {SPARSEN_OPTION_OUTFORMAT,'O',"outformat","The output graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {SPARSEN_OPTION_TYPE,'t',"type","The strategy to use to sparsen the graph.",
      CMD_OPT_CHOICE,SPARSEN,ARRAY_SIZE(SPARSEN)},
  {SPARSEN_OPTION_REWEIGHT,'r',"reweight","The weight distribution method to "
      "use.",CMD_OPT_CHOICE,REWEIGHT,ARRAY_SIZE(REWEIGHT)},
  {SPARSEN_OPTION_SEED,'s',"seed","The seed to use when making 'random' "
      "decisions.",CMD_OPT_INT,NULL,0},
  {SPARSEN_OPTION_FRACTION,'f',"fraction","The fraction of edges to keep in "
      "graph.",CMD_OPT_FLOAT,NULL,0}
};
static const size_t NSPARSEN_OPTS = ARRAY_SIZE(SPARSEN_OPTS);


/* COORDINATES ***************************************************************/

static const cmd_opt_pair_t COORDINATES[] = {
  {"bfs","Assign coordinates based on BFS depth.",COORDINATES_BFS}
};


static const cmd_opt_t COORDINATES_OPTS[] = {
  {COORDINATES_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG,
      NULL,0},
  {COORDINATES_OPTION_INFILE,'i',"infile","The input graph file.",
      CMD_OPT_STRING,NULL,0},
  {COORDINATES_OPTION_OUTFILE,'o',"outfile","The output graph file.",
      CMD_OPT_STRING,NULL,0},
  {COORDINATES_OPTION_INFORMAT,'I',"informat","The input graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {COORDINATES_OPTION_NDIM,'d',"dimensions","The number of dimensions to "
      "embbed the vertices in.",CMD_OPT_INT,NULL,0},
  {COORDINATES_OPTION_SEED,'s',"seed","The seed to use when making 'random' "
      "decisions.",CMD_OPT_INT,NULL,0},
  {COORDINATES_OPTION_TYPE,'t',"type","The method for assigning coordinates.",
      CMD_OPT_CHOICE,COORDINATES,ARRAY_SIZE(COORDINATES)}
};
static const size_t NCOORDINATES_OPTS = ARRAY_SIZE(COORDINATES_OPTS);


/* GENERATE ******************************************************************/

static const cmd_opt_pair_t GENERATE[] = {
  {"complete","Create a complete graph",GENERATE_COMPLETE},
  {"hypercube","Create a hyper cube",GENERATE_HYPERCUBE},
  {"grid","Create a 3D grid",GENERATE_GRID},
  {"dumbbell","Create a dumbbell graph",GENERATE_DUMBBELL},
  {"caveman","Create a caveman graph",GENERATE_CAVEMAN}
};


static const cmd_opt_t GENERATE_OPTS[] = {
  {GENERATE_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG,NULL,
      0},
  {GENERATE_OPTION_TYPE,'t',"type","The type of graph to generate.",
      CMD_OPT_CHOICE,GENERATE,ARRAY_SIZE(GENERATE)},
  {GENERATE_OPTION_OUTFILE,'o',"outfile","The output graph file.",
      CMD_OPT_STRING,NULL,0},
  {GENERATE_OPTION_OUTFORMAT,'O',"outformat","The output graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {GENERATE_OPTION_PARAMS,'p',"param","A parameter for the graph generation.",
      CMD_OPT_FLOAT,NULL,0}
};
static const size_t NGENERATE_OPTS = ARRAY_SIZE(GENERATE_OPTS);


/* FLOW **********************************************************************/

static const cmd_opt_pair_t FLOW[] = {
  {"edge","Create a flow on the edges of a graph.",FLOW_EDGE},
  {"vertex","Create a flow on the vertices of a graph",FLOW_VERTEX}
};


static const cmd_opt_t FLOW_OPTS[] = {
  {FLOW_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG,NULL,
      0},
  {FLOW_OPTION_TYPE,'f',"flow","The type of flow to generate.",
      CMD_OPT_CHOICE,FLOW,ARRAY_SIZE(FLOW)},
  {FLOW_OPTION_INFILE,'i',"infile","The input graph file.",
      CMD_OPT_STRING,NULL,0},
  {FLOW_OPTION_INFORMAT,'I',"informat","The input graph format.",
      CMD_OPT_CHOICE,FORMAT,ARRAY_SIZE(FORMAT)},
  {FLOW_OPTION_OUTFILE,'o',"outfile","The output flow file.",
      CMD_OPT_STRING,NULL,0},
  {FLOW_OPTION_SRC,'s',"src","The source vertex for the flow.",
      CMD_OPT_INT,NULL,0},
  {FLOW_OPTION_DST,'d',"dst","The destination vertex for the flow.",
      CMD_OPT_INT,NULL,0},
  {FLOW_OPTION_TIME,'t',"time","Time the flow operation.",CMD_OPT_FLAG,NULL,0}
};
static const size_t NFLOW_OPTS = ARRAY_SIZE(FLOW_OPTS);




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static int __usage(
    const char * const name, 
    FILE * fout)
{
  size_t i;

  fprintf(fout,"USAGE:\n");
  fprintf(fout,"%s <command> [options]\n",name);
  fprintf(fout,"\n");
  fprintf(fout,"Commands:\n");

  for (i=0;i<NCOMMANDS;++i) {
    fprintf(fout,"\t%s : %s\n",COMMANDS[i].str,COMMANDS[i].desc);
  }

  return 1;
}


static int __command_usage(
    const char * const name, 
    const char * const cmd,
    const cmd_opt_t * const opts, 
    const size_t nopts, 
    FILE * fout)
{
  fprintf(stdout,"USAGE:\n");
  fprintf(stdout,"%s %s [options]\n",name,cmd);
  fprintf(stdout,"\n");
  fprint_cmd_opts(fout,opts,nopts);

  return 1;
}


/* COMMAND FUNCTIONS *********************************************************/

static int __help(
    int argc, 
    char ** argv)
{
  __usage(argv[0],stdout);
  return BOWSTRING_SUCCESS;
}


static int __convert(
    int argc, 
    char ** argv)
{
  size_t nargs, i;
  int err, ondisk, dropweights = 0;
  bowstring_graph_type_t informat = BOWSTRING_FORMAT_AUTO, 
      outformat = BOWSTRING_FORMAT_AUTO;
  cmd_arg_t * args = NULL;
  const char * infile = NULL, * outfile = NULL;

  vtx_t nvtxs;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL, * vwgt = NULL;

  err = cmd_parse_args(argc-2,argv+2,CONVERT_OPTS,NCONVERT_OPTS,&args,
      &nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return BOWSTRING_ERROR_INVALIDINPUT;
  }

  err = BOWSTRING_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],CONVERT_OPTS,NCONVERT_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case CONVERT_OPTION_HELP:
        __command_usage(argv[0],argv[1],CONVERT_OPTS,NCONVERT_OPTS,stdout);
        goto END;
        break;
      case CONVERT_OPTION_INFILE:
        infile = args[i].val.s;
        break;
      case CONVERT_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case CONVERT_OPTION_INFORMAT:
        informat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case CONVERT_OPTION_OUTFORMAT:
        outformat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case CONVERT_OPTION_DROPWEIGHTS:
        dropweights = 1;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = BOWSTRING_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  /* check for on-disk convert */
  ondisk = 0;
  if (informat == BOWSTRING_FORMAT_METIS) {
    if (outformat == BOWSTRING_FORMAT_SNAP) {
      err = translate_metis_to_snap(infile,outfile);
      ondisk = 1;
    } else if (outformat == BOWSTRING_FORMAT_CLOUD9) {
      err = translate_metis_to_cloud9(infile,outfile);
      ondisk = 1;
    }
  } else if (informat == BOWSTRING_FORMAT_CSR) {
    if (outformat == BOWSTRING_FORMAT_SNAP) {
      err = translate_csr_to_snap(infile,outfile);
      ondisk = 1;
    } else if (outformat == BOWSTRING_FORMAT_CLOUD9) {
      err = translate_csr_to_cloud9(infile,outfile);
      ondisk = 1;
    } else if (outformat == BOWSTRING_FORMAT_PEGASUS) {
      err = translate_csr_to_pegasus(infile,outfile);
      ondisk = 1;
    } else if (outformat == BOWSTRING_FORMAT_MATRIXMARKET) {
      err = translate_csr_to_matrixmarket(infile,outfile);
      ondisk = 1;
    }
  } else if (informat == BOWSTRING_FORMAT_EDGELIST) {
    if (outformat == BOWSTRING_FORMAT_SNAP) {
      err = translate_edgelist_to_snap(infile,outfile);
      ondisk = 1;
    }
  }

  if (!ondisk) {
    /* otherwise load it into memory */
    err = bowstring_read_graph(infile,informat,&nvtxs,&xadj,&adjncy,&vwgt,
        &adjwgt);
    if (err != BOWSTRING_SUCCESS) {
      goto END;
    }

    if (dropweights) {
      if (vwgt) {
        dl_free(vwgt);
        vwgt = NULL;
      }
      if (adjwgt) {
        dl_free(adjwgt);
        adjwgt = NULL;
      }
    }

    err = bowstring_write_graph(outfile,outformat,nvtxs,xadj,adjncy,vwgt,
        adjwgt);
    if (err != BOWSTRING_SUCCESS) {
      goto END;
    }
  }

  END:

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

  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __order(
    int argc, 
    char ** argv) 
{
  size_t nargs, i;
  int err;
  vtx_t nvtxs;
  int informat = BOWSTRING_FORMAT_AUTO, outformat = BOWSTRING_FORMAT_AUTO;
  int type = BOWSTRING_ORDER_RANDOM;
  cmd_arg_t * args = NULL;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL, * perm = NULL;
  wgt_t * adjwgt = NULL, * vwgt = NULL;
  const char * infile = NULL, * outfile = NULL;

  err = cmd_parse_args(argc-2,argv+2,ORDER_OPTS,NORDER_OPTS,&args,&nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return BOWSTRING_ERROR_INVALIDINPUT;
  }

  err = BOWSTRING_SUCCESS;

  if (nargs < 1) {
    __command_usage(argv[0],argv[1],ORDER_OPTS,NORDER_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case ORDER_OPTION_HELP:
        __command_usage(argv[0],argv[1],ORDER_OPTS,NORDER_OPTS,stdout);
        goto END;
        break;
      case ORDER_OPTION_INFILE:
        infile = args[i].val.s;
        break;
      case ORDER_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case ORDER_OPTION_INFORMAT:
        informat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case ORDER_OPTION_OUTFORMAT:
        outformat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case ORDER_OPTION_SEED:
        dl_set_rand((unsigned int)args[i].val.i);
        break;
      case ORDER_OPTION_TYPE:
        type = (int)args[i].val.o;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = BOWSTRING_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  err = bowstring_read_graph(infile,informat,&nvtxs,&xadj,&adjncy,&vwgt,
      &adjwgt);
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  perm = vtx_alloc(nvtxs);
  bowstring_permutation(type,nvtxs,xadj,adjncy,vwgt,adjwgt,perm);
  bowstring_order_graph(nvtxs,xadj,adjncy,vwgt,adjwgt,perm);
  dl_free(perm);

  err = bowstring_write_graph(outfile,outformat,nvtxs,xadj,adjncy,vwgt,adjwgt);
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  END:

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

  if (args) {
    dl_free(args); 
  }

  return err;
}


static void __write_degree_distribution(
    FILE * out, 
    const vtx_t nvtxs, 
    const adj_t * const xadj)
{
  vtx_t maxdeg, i;
  vtx_t * degs;

  degree_distribution(nvtxs,xadj,&degs,&maxdeg);

  for (i=0;i<=maxdeg;++i) {
    fprintf(out,PF_VTX_T" "PF_VTX_T"\n",i,degs[i]);
  }

  dl_free(degs);
}


static int __analyze(
    int argc, 
    char ** argv)
{
  size_t i, nargs, npf = 0, ntriangles;
  size_t * cycles;
  vtx_t * stars;
  adj_t ncut;
  vlbl_t nparts;
  vtx_t maxcycle, maxstar, maxdeg, nislands, ncomponents, nleafs, deg;
  double tadjwgt, tvwgt, wcut;
  int * dc;
  int err;
  bowstring_graph_type_t informat = BOWSTRING_FORMAT_AUTO;
  cmd_arg_t * args = NULL;
  const char * infile = NULL;
  const char * partfiles[1024];
  analysis_t type = ANALYSIS_STATS;

  vtx_t nvtxs;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL, * vwgt = NULL;
  vlbl_t * wheres[1024]; 

  err = cmd_parse_args(argc-2,argv+2,ANALYSIS_OPTS,NANALYSIS_OPTS,&args,
      &nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return BOWSTRING_ERROR_INVALIDINPUT;
  }

  err = BOWSTRING_SUCCESS;

  if (nargs < 1) {
    __command_usage(argv[0],argv[1],ANALYSIS_OPTS,NANALYSIS_OPTS,stdout);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case ANALYSIS_OPTION_HELP:
        __command_usage(argv[0],argv[1],ANALYSIS_OPTS,NANALYSIS_OPTS,stdout);
        goto END;
        break;
      case ANALYSIS_OPTION_INFILE:
        infile = args[i].val.s;
        break;
      case ANALYSIS_OPTION_INFORMAT:
        informat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case ANALYSIS_OPTION_PARTFILE:
        partfiles[npf++] = args[i].val.s;
        break;
      case ANALYSIS_OPTION_TYPE:
        type = (analysis_t)args[i].val.o;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = BOWSTRING_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  err = bowstring_read_graph(infile,informat,&nvtxs,&xadj,&adjncy,&vwgt,
      &adjwgt);
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  for (i=0;i<npf;++i) {
    err = read_vertex_labels(partfiles[i],NULL,wheres+i);
    if (err != BOWSTRING_SUCCESS) {
      goto END;
    }
  }

  switch (type) {
    case ANALYSIS_GRAPHSTATS:
      maxdeg = 0;
      nleafs = 0;
      nislands = 0;
      for (i=0;i<nvtxs;++i) {
        deg = (vtx_t)(xadj[i+1] - xadj[i]);
        if (deg == 0) {
          ++nislands;
        } else {
          if (deg == 1) {
            ++nleafs;
          }
          if (maxdeg < xadj[i+1] - xadj[i]) {
            maxdeg = xadj[i+1] - xadj[i];
          }
        }
      }
      if (vwgt) { 
        tvwgt = wgt_fa_sum(vwgt,nvtxs); 
      } else {
        tvwgt = (double)nvtxs;
      }
      if (adjwgt) {
        tadjwgt = wgt_fa_sum(adjwgt,xadj[nvtxs])/2.0; 
      } else {
        tadjwgt = (double)xadj[nvtxs]/2.0;
      }
      label_components(nvtxs,xadj,adjncy,NULL,&ncomponents);
      printf("Number of Vertices   = %16zu\n",(size_t)nvtxs);
      printf("Number of Edges      = %16zu\n",(size_t)(xadj[nvtxs]/2));
      printf("Total Vertex Weight  = %16.3lf\n",tvwgt);
      printf("Total Edge Weight    = %16.3lf\n",tadjwgt);
      printf("Maximum Degree       = %16zu\n",(size_t)maxdeg); 
      printf("Connected Components = %16zu\n",(size_t)ncomponents);
      printf("Island Vertices      = %16zu\n",(size_t)nislands);
      printf("Leaf Vertices        = %16zu\n",(size_t)nleafs);
      break;
    case ANALYSIS_PARTSTATS:
      for (i=0;i<npf;++i) {
        nparts = vlbl_max_value(wheres[i],nvtxs) - \
                 vlbl_min_value(wheres[i],nvtxs)+1;
        if (npf > 1) {
          printf("Partition %zu\n",i);
        }
        dc = int_alloc(nparts);
        calc_domainconn(nvtxs,xadj,adjncy,nparts,wheres[i],dc);
        ncut = calc_edgecut(nvtxs,xadj,adjncy,NULL,wheres[i]);
        wcut = calc_edgecut(nvtxs,xadj,adjncy,adjwgt,wheres[i]);
        printf("Number of Partitions  = %15zu\n",(size_t)nparts);
        printf("Edgecut               = %15zu\n",(size_t)ncut);
        printf("Weighted Edgecut      = %19.3lf\n",wcut);
        printf("Internal edges        = %15zu\n",(size_t)(xadj[nvtxs]/2)-ncut);
        printf("Internal edges %%      = %19.3lf%%\n", \
            (100.0*((xadj[nvtxs]/2)-ncut))/(xadj[nvtxs]/2));
        printf("Max Domain Degree     = %19.3lf\n", \
            calc_max_domaindegree(nvtxs,xadj,adjncy,adjwgt,nparts,wheres[i]));
        printf("Domain Degree Balance = %19.3lf\n", \
            calc_degree_balance(nvtxs,xadj,adjncy,adjwgt,nparts,wheres[i]));
        printf("Max Domain Conn.      = %15zu\n", \
            (size_t)int_max_value(dc,nparts));
        printf("Avg. Domain Conn.     = %19.3lf\n", \
            int_sum(dc,nparts)/(double)nparts);
        printf("Communication Volume  = %19.3lf\n", \
            calc_communicationvolume(nvtxs,xadj,adjncy,adjwgt,nparts, \
              wheres[i]));
        printf("Max domain ComVol     = %19.3lf\n", \
            calc_max_domaincomvol(nvtxs,xadj,adjncy,adjwgt,nparts,wheres[i]));
        printf("Domain ComVol Balance = %19.3lf\n", \
            calc_comvol_balance(nvtxs,xadj,adjncy,adjwgt,nparts,wheres[i]));
        printf("Vertex Balance        = %19.3lf\n", \
            calc_vertex_balance(nvtxs,xadj,adjncy,vwgt,nparts,wheres[i]));
        printf("Edge Balance          = %19.3lf\n", \
            calc_edge_balance(nvtxs,xadj,adjncy,adjwgt,nparts,wheres[i]));
        printf("Modularity            = %23.7lf\n", \
            calc_modularity(nvtxs,xadj,adjncy,adjwgt,nparts,wheres[i]));
        printf("Partition Components  = %19.3lf\n", \
            calc_partition_components(nvtxs,xadj,adjncy,nparts,wheres[i]));
        dl_free(dc);
      }
      break;
    #ifdef XXX
    case ANALYSIS_BNDSIM:
      printf("Calculating bndsim of graph with "PF_VTX_T" vertices\n",
          nvtxs);
      read_vertex_labels(oper.outfile,&nvtxs,&where);
      anbnd = find_boundary_vertices(nvtxs,xadj,adjncy,where,&abnd);
      avgdeg = 0;
      for (i=0;i<anbnd;++i) {
        v = abnd[i];
        avgdeg += xadj[v+1]-xadj[v];
      }
      printf("A has average bnd degree of %f and total of %f\n",avgdeg/anbnd,
          avgdeg);
      dl_free(where);
      read_vertex_labels(oper.auxfile,&nvtxs,&where);
      bnbnd = find_boundary_vertices(nvtxs,xadj,adjncy,where,&bbnd);
      avgdeg = 0;
      for (i=0;i<bnbnd;++i) {
        v = bbnd[i];
        avgdeg += xadj[v+1]-xadj[v];
      }
      printf("B has average bnd degree of %f and total of %f\n",avgdeg/bnbnd,
          avgdeg);
      printf("Cut A has "PF_VTX_T" boundary vertices, cut B has "PF_VTX_T
          " vertices, and they share "PF_SIZE_T" vertices\n",anbnd,bnbnd,
          vtx_intersection_size(abnd,anbnd,bbnd,bnbnd));
      dl_free(where);
      break;
    #endif
    case ANALYSIS_DEGDIST:
      __write_degree_distribution(stdout,nvtxs,xadj);
      break;
    case ANALYSIS_NTRIANGLES:
      ntriangles = count_triangles(nvtxs,xadj,adjncy,NULL);
      printf("Number of Triangles: %zu\n",ntriangles);
      printf("Triangles Per Vertex: %0.03lf\n",ntriangles/(double)nvtxs);
      printf("Triangles Per Edge: %0.03lf\n",ntriangles/(double)xadj[nvtxs]);
      break;
    case ANALYSIS_ATOMICCYCLES:
      atomic_cycle_distribution(nvtxs,xadj,adjncy,NULL,&cycles,&maxcycle); 
      printf("# Atomic Cycle Histogram of '%s'\n",infile);
      for (i=0;i<=(size_t)maxcycle;++i) {
        printf(PF_SIZE_T" "PF_SIZE_T"\n",i,cycles[i]);
      }
      dl_free(cycles);
      break;
    case ANALYSIS_STARS:
      star_distribution(nvtxs,xadj,adjncy,&stars,&maxstar);
      printf("# Star Histogram of : %s\n",infile);
      for (i=0;i<=(size_t)maxstar;++i) {
        printf(PF_SIZE_T" "PF_VTX_T"\n",i,stars[i]);
      }
      dl_free(stars);
      break;
    default:
      eprintf("Unsupported/Unimplemented analysis '%d'\n",type);
      err = BOWSTRING_ERROR_UNIMPLEMENTED;
      break;
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  END:

  for (i=0;i<npf;++i) {
    if (wheres[i]) {
      dl_free(wheres[i]);
    }
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

  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __sparsen(
    int argc, 
    char ** argv)
{
  size_t nargs, i;
  int err;
  dl_timer_t tmr;

  elbl_t maxrank = 0;
  double frac = 0.5;
  bowstring_graph_type_t informat = BOWSTRING_FORMAT_AUTO, 
    outformat = BOWSTRING_FORMAT_AUTO;
  sparsen_t type = SPARSEN_MSR;
  bowstring_reweight_t reweight = BOWSTRING_REWEIGHT_EXACT;

  int * adjmap = NULL;
  cmd_arg_t * args = NULL;
  elbl_t * rank = NULL;
  const char * infile = NULL, * outfile = NULL;

  vtx_t nvtxs;
  adj_t * xadj = NULL, * sxadj = NULL;
  vtx_t * adjncy = NULL, * sadjncy = NULL;
  wgt_t * adjwgt = NULL, * vwgt = NULL, * sadjwgt = NULL;

  dl_init_timer(&tmr);

  err = cmd_parse_args(argc-2,argv+2,SPARSEN_OPTS,NSPARSEN_OPTS,&args,&nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return BOWSTRING_ERROR_INVALIDINPUT;
  }

  err = BOWSTRING_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],SPARSEN_OPTS,NSPARSEN_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case SPARSEN_OPTION_HELP:
        __command_usage(argv[0],argv[1],SPARSEN_OPTS,NSPARSEN_OPTS,stdout);
        goto END;
        break;
      case SPARSEN_OPTION_INFILE:
        infile = args[i].val.s;
        break;
      case SPARSEN_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case SPARSEN_OPTION_INFORMAT:
        informat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case SPARSEN_OPTION_OUTFORMAT:
        outformat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case SPARSEN_OPTION_SEED:
        dl_set_rand((unsigned int)args[i].val.i);
        break;
      case SPARSEN_OPTION_FRACTION:
        frac = args[i].val.f;
        break;
      case SPARSEN_OPTION_REWEIGHT:
        reweight = (bowstring_reweight_t)args[i].val.o;
        break;
      case SPARSEN_OPTION_TYPE:
        type = (int)args[i].val.o;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = BOWSTRING_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  err = bowstring_read_graph(infile,informat,&nvtxs,&xadj,&adjncy,&vwgt,
      &adjwgt);
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  dl_start_timer(&tmr);
  switch (type) {
    case SPARSEN_NIR:
      rank = elbl_calloc(xadj[nvtxs]);
      maxrank = build_nirank(nvtxs,xadj,adjncy,adjwgt,rank);
      break;
    case SPARSEN_MSR:
      rank = elbl_alloc(xadj[nvtxs]);
      maxrank = build_mstrank(nvtxs,xadj,adjncy,adjwgt,rank);
      break;
    case SPARSEN_ASR:
      rank = elbl_alloc(xadj[nvtxs]);
      maxrank = build_astrank(nvtxs,xadj,adjncy,adjwgt,rank);
      break;
    case SPARSEN_DFS:
      adjmap = int_calloc(xadj[nvtxs]);
      build_dfs_tree(nvtxs,xadj,adjncy,vtx_rand(0,nvtxs),NULL,NULL,NULL,
          adjmap);
      break;
    case SPARSEN_BFS:
      adjmap = int_calloc(xadj[nvtxs]);
      build_bfs_tree(nvtxs,xadj,adjncy,vtx_rand(0,nvtxs),NULL,NULL,NULL,
          adjmap);
      break;
    default:
      eprintf("Unsupported/Unimplemented sparsening '%d'\n",type);
      err = BOWSTRING_ERROR_UNIMPLEMENTED;
      break;
  }
  dl_stop_timer(&tmr);
  printf("Edge selection took %lfs (%lf edges / second)\n",dl_poll_timer(&tmr),
      xadj[nvtxs]/dl_poll_timer(&tmr));

  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  dl_reset_timer(&tmr);
  dl_start_timer(&tmr);
  switch (type) {
    case SPARSEN_NIR:
    case SPARSEN_MSR:
    case SPARSEN_ASR:
      prune_ranked_edges(nvtxs,xadj,adjncy,adjwgt,rank,maxrank,frac,&sxadj,
          &sadjncy,&sadjwgt,reweight);
      break;
    case SPARSEN_DFS:
    case SPARSEN_BFS:
      apply_edge_mask(nvtxs,xadj,adjncy,adjwgt,adjmap,&sxadj,&sadjncy,
          &sadjwgt);
      break;
    default:
      eprintf("Unsupported/Unimplemented sparsening '%d'\n",type);
      err = BOWSTRING_ERROR_UNIMPLEMENTED;
      break;
  }
  dl_stop_timer(&tmr);
  printf("Edge removal took %lfs (%lf edges / second)\n",dl_poll_timer(&tmr),
      xadj[nvtxs]/dl_poll_timer(&tmr));

  err = bowstring_write_graph(outfile,outformat,nvtxs,sxadj,sadjncy,vwgt,
      sadjwgt);
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  END:

  if (xadj) {
    dl_free(xadj);
  }

  if (sxadj) {
    dl_free(sxadj);
  }

  if (adjncy) {
    dl_free(adjncy);
  }

  if (sadjncy) {
    dl_free(sadjncy);
  }

  if (adjwgt) {
    dl_free(adjwgt);
  }

  if (sadjwgt) {
    dl_free(sadjwgt);
  }

  if (vwgt) {
    dl_free(vwgt);
  }

  if (rank) {
    dl_free(rank);
  }

  if (adjmap) {
    dl_free(adjmap);
  }

  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __coordinates(
    int argc, 
    char ** argv) 
{
  size_t nargs, i, ndim;
  int err;
  bowstring_graph_type_t informat = BOWSTRING_FORMAT_AUTO;
  coordinates_t type = COORDINATES_BFS;
  coord_t ** coords = NULL;
  cmd_arg_t * args = NULL;
  const char * infile = NULL, * outfile = NULL;

  vtx_t nvtxs;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL, * vwgt = NULL;

  ndim = 2;

  err = cmd_parse_args(argc-2,argv+2,COORDINATES_OPTS,NCOORDINATES_OPTS,&args,
      &nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return BOWSTRING_ERROR_INVALIDINPUT;
  }

  err = BOWSTRING_SUCCESS;

  if (nargs < 1) {
    __command_usage(argv[0],argv[1],COORDINATES_OPTS,NCOORDINATES_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case COORDINATES_OPTION_HELP:
        __command_usage(argv[0],argv[1],COORDINATES_OPTS,NCOORDINATES_OPTS,
            stdout);
        goto END;
        break;
      case COORDINATES_OPTION_INFILE:
        infile = args[i].val.s;
        break;
      case COORDINATES_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case COORDINATES_OPTION_INFORMAT:
        informat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case COORDINATES_OPTION_SEED:
        dl_set_rand((unsigned int)args[i].val.i);
        break;
      case COORDINATES_OPTION_NDIM:
        ndim = args[i].val.i;
        break;
      case COORDINATES_OPTION_TYPE:
        type = (coordinates_t)args[i].val.o;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = BOWSTRING_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  err = bowstring_read_graph(infile,informat,&nvtxs,&xadj,&adjncy,&vwgt,
      &adjwgt);
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  coords = r_coord_sym_alloc(nvtxs,ndim);

  switch (type) {
    case COORDINATES_BFS:
      err = bfs_coordinates(nvtxs,xadj,adjncy,adjwgt,ndim,coords);
      break;
    default:
      eprintf("Unsupported/Unimplemented ordering '%d'\n",type);
      err = BOWSTRING_ERROR_UNIMPLEMENTED;
      break;
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }


  /* delete me */
  vlbl_t * where = vlbl_alloc(nvtxs);

  for (i=0;i<ndim;++i) {
    coordinate_bisection(nvtxs,NULL,coords[i],nvtxs/2,where);
    printf("%zu: cut = "PF_WGT_T"\n",i,calc_edgecut(nvtxs,xadj,adjncy,adjwgt,
        where));
  }

  err = write_vertex_labels(outfile,nvtxs,where);

  //err = write_matrix(outfile,(const coord_t * const *)coords,nvtxs,ndim);
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }


  END:

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

  if (coords) {
    for (i=0;i<ndim;++i) {
      dl_free(coords[i]);
    }
    dl_free(coords);
  }

  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __generate(
    int argc, 
    char ** argv)
{
  size_t nargs, i, nparams;
  int err;
  double params[32];
  cmd_arg_t * args = NULL;

  vtx_t nvtxs;
  bowstring_graph_type_t outformat = BOWSTRING_FORMAT_AUTO; 
  generate_t type = GENERATE_COMPLETE;
  const char * outfile = NULL;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;

  err = cmd_parse_args(argc-2,argv+2,GENERATE_OPTS,NGENERATE_OPTS,&args,
      &nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return BOWSTRING_ERROR_INVALIDINPUT;
  }

  err = BOWSTRING_SUCCESS;

  nparams = 0;

  if (nargs < 1) {
    __command_usage(argv[0],argv[1],GENERATE_OPTS,NGENERATE_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case GENERATE_OPTION_HELP:
        __command_usage(argv[0],argv[1],GENERATE_OPTS,NGENERATE_OPTS,stdout);
        goto END;
        break;
      case GENERATE_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case GENERATE_OPTION_OUTFORMAT:
        outformat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case GENERATE_OPTION_TYPE:
        type = (generate_t)args[i].val.o;
        break;
      case GENERATE_OPTION_PARAMS:
        params[nparams++] = args[i].val.f;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = BOWSTRING_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  switch(type) {
    case GENERATE_COMPLETE:
      if (nparams != 1) {
        eprintf("Invalid number of parameters for generating complete "
            "graph\n");
        err = BOWSTRING_ERROR_INVALIDINPUT;
      } else {
        nvtxs = (vtx_t)params[0];
        generate_complete_graph(nvtxs,&xadj,&adjncy,NULL,NULL);
      }
      break;
    case GENERATE_GRID:
      if (nparams > 3) {
        eprintf("Invalid number of parameters for generating 3D grid\n");
        err = BOWSTRING_ERROR_INVALIDINPUT;
        goto END;
      } else {
        if (nparams < 3) {
          params[2] = 1.0;
        }
        if (nparams < 2) {
          params[1] = 1.0;
        }
        nvtxs = (vtx_t)(params[0]*params[1]*params[2]);
        generate_grid_graph(params[0],params[1],params[2],&xadj,&adjncy,
            NULL,NULL);
      }
      break;
    default:
      eprintf("Unknown argument '%s'\n",args[i].val.s);
      err = BOWSTRING_ERROR_INVALIDINPUT;
      break;
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  bowstring_write_graph(outfile,outformat,nvtxs,xadj,adjncy,NULL,NULL);

  END:

  if (xadj) {
    dl_free(xadj);
  }

  if (adjncy) {
    dl_free(adjncy);
  }

  if (args) {
    dl_free(args); 
  }

  return err;
}


static int __flow(
    int argc, 
    char ** argv)
{
  size_t nargs, i, nflow;
  int err, time = 0;
  dl_timer_t flowtmr, iotmr;
  vtx_t src = 0, dst = 1;
  wgt_t maxflow = 0;

  bowstring_graph_type_t informat = BOWSTRING_FORMAT_AUTO;
  flow_t type = FLOW_EDGE;

  cmd_arg_t * args = NULL;
  const char * infile = NULL, * outfile = NULL;

  vtx_t nvtxs;
  adj_t * xadj = NULL;
  vtx_t * adjncy = NULL;
  wgt_t * adjwgt = NULL, * vwgt = NULL, * flow = NULL;

  dl_init_timer(&flowtmr);
  dl_init_timer(&iotmr);

  err = cmd_parse_args(argc-2,argv+2,FLOW_OPTS,NFLOW_OPTS,&args,&nargs);
  if (err != DL_CMDLINE_SUCCESS) {
    return BOWSTRING_ERROR_INVALIDINPUT;
  }

  err = BOWSTRING_SUCCESS;

  if (nargs < 2) {
    __command_usage(argv[0],argv[1],FLOW_OPTS,NFLOW_OPTS,stderr);
    goto END;
  }
  for (i=0;i<nargs;++i) {
    switch (args[i].id) {
      case FLOW_OPTION_HELP:
        __command_usage(argv[0],argv[1],FLOW_OPTS,NFLOW_OPTS,stdout);
        goto END;
        break;
      case FLOW_OPTION_INFILE:
        infile = args[i].val.s;
        break;
      case FLOW_OPTION_OUTFILE:
        outfile = args[i].val.s;
        break;
      case FLOW_OPTION_INFORMAT:
        informat = (bowstring_graph_type_t)args[i].val.o;
        break;
      case FLOW_OPTION_TYPE:
        type = (int)args[i].val.o;
        break;
      case FLOW_OPTION_SRC:
        src = (vtx_t)args[i].val.i;
        break;
      case FLOW_OPTION_DST:
        dst = (vtx_t)args[i].val.i;
        break;
      case FLOW_OPTION_TIME:
        time = 1;
        break;
      default:
        eprintf("Unknown argument '%s'\n",args[i].val.s);
        err = BOWSTRING_ERROR_INVALIDINPUT;
        break;
    }
  }
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  dl_start_timer(&iotmr);

  err = bowstring_read_graph(infile,informat,&nvtxs,&xadj,&adjncy,&vwgt,
      &adjwgt);
  if (err != BOWSTRING_SUCCESS) {
    goto END;
  }

  if (src > nvtxs || src < 1) {
    eprintf("Invalid source vertex "PF_VTX_T" for graph of "PF_VTX_T \
        " vertices\n",src,nvtxs);
    goto END;
  }
  if (dst > nvtxs || dst < 1) {
    eprintf("Invalid destination vertex "PF_VTX_T" for graph of "PF_VTX_T \
        " vertices\n",dst,nvtxs);
    goto END;
  }
  if (dst == src) {
    eprintf("Source and destination vertices are both "PF_VTX_T"\n",dst);
    goto END;
  }

  dl_stop_timer(&iotmr);

  /* adjust vertex labeling */
  --src;
  --dst;

  dl_start_timer(&flowtmr);
  nflow = xadj[nvtxs];
  flow = wgt_alloc(nflow);
  switch (type) {
    case FLOW_EDGE:
      dprintf("Performing edge based flow from "PF_VTX_T" to "PF_VTX_T"\n", \
          src,dst);
      if (adjwgt == NULL) {
        adjwgt = wgt_init_alloc(1,xadj[nvtxs]);
      }
      maxflow = maxflow_edge(src,dst,nvtxs,xadj,adjncy,adjwgt,flow);
      break;
    case FLOW_VERTEX:
      dprintf("Performing vertex based flow from "PF_VTX_T" to "PF_VTX_T"\n", \
          src,dst);
      if (vwgt == NULL) {
        vwgt = wgt_init_alloc(1,nvtxs);
      }
      maxflow = maxflow_vertex(src,dst,nvtxs,xadj,adjncy,vwgt,flow);
      break;
    default:
      eprintf("Unsupported/Unimplemented flow '%d'\n",type);
      err = BOWSTRING_ERROR_UNIMPLEMENTED;
      goto END;
  }
  dl_stop_timer(&flowtmr);

  dl_start_timer(&iotmr);
  if (outfile) {
    err = write_weights(outfile,flow,nflow);
    if (err != BOWSTRING_SUCCESS) {
      goto END;
    }
  }
  printf("Maximum flow of "PF_WGT_T"\n",maxflow);
  dl_stop_timer(&iotmr);

  if (time) {
    printf("Flow Time: %0.05fs | IO Time: %0.05fs\n",dl_poll_timer(&flowtmr), \
        dl_poll_timer(&iotmr));
  }

  END:

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

  if (flow) {
    dl_free(flow);
  }

  if (args) {
    dl_free(args); 
  }

  return err;
}


typedef int (*__cmdfuncptr_t)(int,char**); 
static const __cmdfuncptr_t COMMAND_FUNCS[] = {
  [COMMAND_HELP] = __help,
  [COMMAND_CONVERT] = __convert,
  [COMMAND_ORDER] = __order,
  [COMMAND_ANALYSIS] = __analyze,
  [COMMAND_SPARSEN] = __sparsen,
  [COMMAND_COORDINATES] = __coordinates,
  [COMMAND_GENERATE] = __generate,
  [COMMAND_FLOW] = __flow
};


int main(
    int argc, 
    char ** argv)
{
  int err;
  char * cmdstr;
  size_t i;

  dl_init_rand();

  if (argc < 2) {
    eprintf("Must supply a command.\n");
    __usage(argv[0],stderr);
    return 1;
  }
  cmdstr = argv[1];

  for (i=0;i<NCOMMANDS;++i) {
    if (strcmp(cmdstr,COMMANDS[i].str) == 0) {
      err = COMMAND_FUNCS[i](argc,argv);
      break;
    }
  }
  if (i == NCOMMANDS) {
    eprintf("Unrecognized command '%s'.\n",cmdstr);
    __usage(argv[0],stderr);
    return 1;
  }

  if (err == BOWSTRING_SUCCESS) {
    return 0;
  } else {
    eprintf("Operation failed.\n");
    return 2;
  }
}

  


#endif

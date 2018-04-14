/*
 * cmdline_pmetis.c
 *
 * This file parses the command line arguments
 *
 * Started 8/9/02
 * George
 *
 * $Id: cmdline_pmetis.c,v 1.4 2002/08/12 15:20:51 karypis Exp $
 *
 */


#include <metisbin.h>


/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct option long_options[] = {
  {"mtype",          1,      0,      CMD_MTYPE},
  {"itype",          1,      0,      CMD_ITYPE},
  {"rtype",          1,      0,      CMD_RTYPE},

  {"balanced",       0,      0,      CMD_BALANCE},

  {"niter",          1,      0,      CMD_NITER},

  {"tpwgts",         1,      0,      CMD_TPWGTS},

  {"seed",           1,      0,      CMD_SEED},

  {"dbglvl",         1,      0,      CMD_DBGLVL},

  {"help",           0,      0,      CMD_HELP},
  {0,                0,      0,      0}
};



/*-------------------------------------------------------------------
 * Mappings for the various parameter values
 *-------------------------------------------------------------------*/
static StringMapType mtype_options[] = {
 {"rm",                 MTYPE_RM},
 {"hem",                MTYPE_HEM},
 {"shem",               MTYPE_SHEM},
 {"shebm",              MTYPE_SHEBM_ONENORM},
 {"sbhem",              MTYPE_SBHEM_ONENORM},
 {NULL,                 0}
};


static StringMapType itype_options[] = {
 {"greedy",              ITYPE_GGPKL},
 {"random",              ITYPE_RANDOM},
 {NULL,                 0}
};


static StringMapType rtype_options[] = {
 {"fm",            RTYPE_FM},
 {NULL,                 0}
};



/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: pmetis [options] <filename> <nparts>",
" ",
" Required parameters",
"    filename    Stores the graph to be partitioned.",
"    nparts      The number of partitions to split the graph.",
" ",
" Optional parameters",
"  -mtyep=string",
"     Specifies the scheme to be used to match the vertices of the graph",
"     during the coarsening.",
"     The possible values are:",
"        rm       - Random matching",
"        hem      - Heavy-edge matching",
"        shem     - Sorted heavy-edge matching [default]",
"        shebm    - Combination of shem and balanced matching for",
"                   multi-constraint.",
"        sbhem    - Similar as shebm but priority is given to balance",
" ",
"  -itype=string",
"     Specifies the scheme to be used to compute the initial partitioning",
"     of the graph.",
"     The possible values are:",
"        greedy   - Grow a bisection using a greedy strategy [default]",
"        random   - Compute a bisection at random",
" ",
"  -rtype=string",
"     Specifies the scheme to be used for refinement",
"     The possible values are:",
"        fm       - FM refinement",
" ",
"  -balance",
"     Specifies that the final partitioning should contain nparts-1 equal",
"     size partitions with the last partition having upto nparts-1 fewer",
"     vertices.",
" ",
"  -seed=int      ",
"     Selects the seed of the random number generator.  ",
" ",
"  -dbglvl=int      ",
"     Selects the dbglvl.  ",
" ",
"  -help",
"     Prints this message.",
""
};

static char shorthelpstr[][100] = {
" ",
"   Usage: pmetis [options] <filename> <nparts>",
"          use 'pmetis -help' for a summary of the options.",
""
};
 


/*************************************************************************
* This is the entry point of the command-line argument parser
**************************************************************************/
void parse_cmdline(ParamType *params, idxtype argc, char *argv[])
{
  idxtype i, j, k;
  idxtype c, option_index;

  /* initialize the params data structure */
  params->mtype         = PMETIS_CTYPE;
  params->itype         = PMETIS_ITYPE;
  params->rtype         = PMETIS_RTYPE;
  params->dbglvl        = PMETIS_DBGLVL;

  params->balance       = 0;
  params->seed            = -1;
  params->dbglvl          = 0;
  params->filename        = NULL;
  params->nparts          = 1;


  /* Parse the command line arguments  */
  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
      case CMD_MTYPE:
        if (optarg)
          if ((params->mtype = GetStringID(mtype_options, optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, optarg);
        break;
      case CMD_ITYPE:
        if (optarg)
          if ((params->itype = GetStringID(itype_options, optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, optarg);
        break;
      case CMD_RTYPE:
        if (optarg)
          if ((params->rtype = GetStringID(rtype_options, optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, optarg);
        break;

      case CMD_BALANCE:
        params->balance = 1;
        break;


      case CMD_SEED:
        if (optarg) params->seed = atoi(optarg);
        break;

      case CMD_DBGLVL:
        if (optarg) params->dbglvl = atoi(optarg);
        break;

      case CMD_HELP:
        for (i=0; strlen(helpstr[i]) > 0; i++)
          printf("%s\n", helpstr[i]);
        exit(0);
        break;
      case '?':
      default:
        printf("Illegal command-line option(s)\nUse %s -help for a summary of the options.\n", argv[0]);
        exit(0);
    }
  }

  if (argc-optind != 2) {
    printf("Missing parameters.");
    for (i=0; strlen(shorthelpstr[i]) > 0; i++)
      printf("%s\n", shorthelpstr[i]);
    exit(0);
  }

  params->filename = strdup(argv[optind++]);
  params->nparts   = atoi(argv[optind++]);
    
}



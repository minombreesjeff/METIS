#include "test.h"

#define NOPTS 7

sint_t test(void) 
{
  sint_t rv = 0;
  int success;
  size_t i, nargs,x;
  char * argv_b[] = {"--opt_str=hello","--opt_int=45","-y","arg1","-F","-Byes",
    "--opt_bool=yes","--opt_float=3.14159","arg2","-f","2342"};
  char * xargs[] = {"arg1","arg2"};
  char ** argv = (char**)argv_b;
  size_t argc = 9;

  cmd_opt_t * opts = cmd_opt_alloc(NOPTS); 

  opts[0].id = 1;
  opts[0].type = CMD_OPT_INT;
  opts[0].sflag = 'i';
  opts[0].lflag = "opt_int";
  opts[0].desc = "Some imteger option\n";

  opts[1].id = 0;
  opts[1].type = CMD_OPT_STRING;
  opts[1].sflag = 's';
  opts[1].lflag = "opt_str";
  opts[1].desc = "Some string option\n";

  opts[2].id = 2;
  opts[2].type = CMD_OPT_BOOL;
  opts[2].sflag = 'b';
  opts[2].lflag = "opt_bool";
  opts[2].desc = "Some bool option\n";

  opts[3].id = 3;
  opts[3].type = CMD_OPT_FLAG;
  opts[3].sflag = 'y';
  opts[3].lflag = "opt_flag";
  opts[3].desc = "Some flag option\n";

  opts[4].id = 4;
  opts[4].type = CMD_OPT_FLOAT;
  opts[4].sflag = 'f';
  opts[4].lflag = "opt_float";
  opts[4].desc = "Some flag option\n";

  opts[5].id = 5;
  opts[5].type = CMD_OPT_FLAG;
  opts[5].sflag = 'F';
  opts[5].lflag = "other_opt_flag";
  opts[5].desc = "Some other flag option\n";

  opts[6].id = 6;
  opts[6].type = CMD_OPT_BOOL;
  opts[6].sflag = 'B';
  opts[6].lflag = "other_opt_BOOL";
  opts[6].desc = "Some other bool option\n";

  cmd_arg_t * args;
  success = cmd_parse_args(argc, argv, opts, NOPTS, &args, &nargs);
  TESTEQUALS(success,DL_CMDLINE_SUCCESS,"%d");
  TESTEQUALS(nargs,(size_t)9,"%zu");
  TESTTRUE(args != NULL);

  dl_free(opts);

  x = 0;
  argc = 7;
  for (i=0;i<nargs;++i) { 
    switch (args[i].id) {
      case -1:
        printf("args[%d] = '%s'\n",(int)i,args[i].val.s);
        TESTEQUALS(args[i].type,CMD_OPT_XARG,"%d");
        TESTEQUALS(strcmp(xargs[x],args[i].val.s),0,"%d");
        ++x;
        break;
      case 0:
        TESTEQUALS(args[i].type,CMD_OPT_STRING,"%d");
        TESTEQUALS(strcmp(args[i].val.s,"hello"),0,"%d");
        --argc;
        break;
      case 1:
        TESTEQUALS(args[i].type,CMD_OPT_INT,"%d");
        TESTEQUALS(args[i].val.i,45LL,"%lld");
        --argc;
        break;
      case 2:
        TESTEQUALS(args[i].type,CMD_OPT_BOOL,"%d");
        TESTEQUALS((args[i].val.b ? 1 : 0),1,"%d");
        --argc;
        break;
      case 3:
        /* not sure how to test this yet */
        TESTEQUALS(args[i].type,CMD_OPT_FLAG,"%d");
        --argc;
        break;
      case 4:
        TESTEQUALS(args[i].type,CMD_OPT_FLOAT,"%d");
        TESTEQUALS(args[i].val.f,3.14159,"%lf");
        --argc;
        break;
      case 5:
        /* not sure how to test this yet */
        TESTEQUALS(args[i].type,CMD_OPT_FLAG,"%d");
        --argc;
        break;
      case 6:
        TESTEQUALS(args[i].type,CMD_OPT_BOOL,"%d");
        TESTEQUALS((args[i].val.b ? 1 : 0),1,"%d");
        --argc;
        break;
      default:
        /* do nothing */
        break;
    }
  }

  TESTEQUALS(argc,(size_t)0,"%zu");

  dl_free(args);

  return rv;
}

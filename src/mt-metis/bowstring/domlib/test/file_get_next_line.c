#include "test.h"

#define N 10L
#define M 265010L
#define FILENAMESMALL "data/10.txt"
#define FILENAMEBIG "data/big.txt"

sint_t test(void) 
{
  file_t * file;

  int rv = dl_open_file(FILENAMESMALL,"r",&file);
    
  TESTEQUALS(rv,DL_FILE_SUCCESS,"%d");

  ssize_t ll,i;
  size_t bs = 5;
  char * buffer = char_alloc(bs);

  i = 0;
  while ((ll = dl_get_next_line(file,&buffer,&bs)) > 0) {
    TESTEQUALS(ll,1L,PF_SSIZE_T);
    TESTEQUALS((ssize_t)atoi(buffer),i,PF_SSIZE_T);
    ++i;
  }

  dl_close_file(file);

  rv = dl_open_file(FILENAMEBIG,"r",&file);
  TESTEQUALS(rv,DL_FILE_SUCCESS,"%d");

  i = 0;
  while ((ll = dl_get_next_line(file,&buffer,&bs)) >= 0) {
    TESTEQUALS(buffer[ll],'\0',"%c");
    ++i;
  }

  TESTEQUALS(i,M,PF_SSIZE_T);

  dl_close_file(file);

  return 0;
}

#include "test.h"

#define B 'c'
#define I 0x5230ff23
#define D 153.5
#define N 1024*1024

sint_t test(void) 
{
  sint_t rv = 0;
  char * mem = (char*)malloc(sizeof(char)+ sizeof(int) + sizeof(double));

  char c = B, cb;
  int i = I, ib;
  double d = D, db;

  size_t j;

  for (j=0;j<N;++j) {
    cb = (char)(c + j);
    ib = (int)(i + j);
    db = (double)(d+j); 
    dl_to_bytes(mem,&cb,sizeof(char));
    dl_to_bytes(mem+sizeof(char),&ib,sizeof(int));
    dl_to_bytes(mem+sizeof(char)+sizeof(int),&db,sizeof(double));

    dl_from_bytes(&c,mem,sizeof(char));
    dl_from_bytes(&i,mem+sizeof(char),sizeof(int));
    dl_from_bytes(&d,mem+sizeof(char)+sizeof(int),sizeof(double));

    TESTEQUALS(c,cb,"%c"); 
    TESTEQUALS(i,ib,"%d"); 
    TESTEQUALS(d,db,"%lf"); 
  }

  free(mem);
  return rv;
}

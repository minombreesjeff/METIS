#include "test.h"


sint_t test(void)
{
  TESTTRUE(dl_string_endswith("hello_world.txt",".txt"));
  TESTTRUE(dl_string_endswith("hello_world.txt","hello_world.txt"));
  TESTTRUE(dl_string_endswith("hello_world.txt",""));
  TESTTRUE(dl_string_endswith("foo bar","oo bar"));
  TESTTRUE(dl_string_endswith("",""));
  TESTTRUE(!dl_string_endswith("def","abcdef"));
  TESTTRUE(!dl_string_endswith("hello_world.txt","world"));
  TESTTRUE(!dl_string_endswith("hello_world.txt","hello"));
  TESTTRUE(!dl_string_endswith("hello_world.txt","sxt"));
  TESTTRUE(!dl_string_endswith("hello_world.txt","txts"));
  TESTTRUE(!dl_string_endswith("hello_world.txt",NULL));
  TESTTRUE(!dl_string_endswith(NULL,"txt"));
  TESTTRUE(!dl_string_endswith(NULL,NULL));

  return 0;
}

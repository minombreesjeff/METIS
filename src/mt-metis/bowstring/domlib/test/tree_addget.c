#include "test.h"


static int compare(const char * str1, const char * str2)
{
  return strcmp(str1,str2);
}

sint_t test(void) 
{
  int x;
  sint_t y;
  sint_tree_t * tree = sint_tree_create(compare);

  sint_tree_add("hello",1,tree);
  sint_tree_add("world",2,tree);
  sint_tree_add("foo",0,tree);
  sint_tree_add("bar",3,tree);
  sint_tree_add("dog",8,tree);
  sint_tree_add("cat",80,tree);
  sint_tree_add("mouse",4,tree);
  sint_tree_add("rat",22,tree);
  sint_tree_add("bat",22,tree);

  x = sint_tree_add("foo",1,tree);

  TESTEQUALS(x,0,"%d");

  x = sint_tree_add("bats",22,tree);

  TESTTRUE(x != 0);

  y = sint_tree_get("rat",tree);

  TESTEQUALS(y,22,PF_SINT_T);

  y = sint_tree_get("bar",tree);

  TESTEQUALS(y,3,PF_SINT_T);

  y = sint_tree_get("foo",tree);

  TESTEQUALS(y,0,PF_SINT_T);

  y = sint_tree_get("hello",tree);

  TESTEQUALS(y,1,PF_SINT_T);

  y = sint_tree_get("dog",tree);

  TESTEQUALS(y,8,PF_SINT_T);

  sint_tree_free(tree);

  return 0;
}

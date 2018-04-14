#include "test.h"

#define N 10

sint_t test(void)
{
  const sint_t sorted[] = {-4,0,3,4,5,5,7,10,256,1104};

  TESTEQUALS(sint_binarysearch(sorted,-4,N),(ssize_t)0,PF_SSIZE_T);
  TESTEQUALS(sint_binarysearch(sorted,-3,N),(ssize_t)0,PF_SSIZE_T);
  TESTEQUALS(sint_binarysearch(sorted,10,N),(ssize_t)7,PF_SSIZE_T);
  TESTEQUALS(sint_binarysearch(sorted,1104,N),(ssize_t)9,PF_SSIZE_T);
  TESTEQUALS(sint_binarysearch(sorted,2048,N),(ssize_t)9,PF_SSIZE_T);

  return 0;
}

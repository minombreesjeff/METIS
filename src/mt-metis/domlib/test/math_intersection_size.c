#include "test.h"

#define N 1000
#define M 197

sint_t test(void) 
{
  sint_t an = 10; 
  sint_t a[] = {1,3,5,7,9,11,13,15,17,19};
  sint_t bn = 7;
  sint_t b[] = {0,1,2,3,4,5,6};
  sint_t cn = 12;
  sint_t c[] = {11,12,13,14,15,16,17,18,19,20,21,22};

  TESTEQUALS(sint_intersection_size(a,an,b,bn),(size_t)3,PF_SIZE_T);
  TESTEQUALS(sint_intersection_size(a,an,c,cn),(size_t)5,PF_SIZE_T);
  TESTEQUALS(sint_intersection_size(c,cn,b,bn),(size_t)0,PF_SIZE_T);

  return 0;
}

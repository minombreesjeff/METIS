#include "test.h"

#define N 12

sint_t test(void) 
{
  ss_ht_t * map = ss_ht_create(N,N);

  TESTEQUALS(0UL,map->size,PF_SIZE_T);
  TESTEQUALS(16UL,map->maxsize,PF_SIZE_T);
  TESTEQUALS(15UL,map->hashmask,PF_SIZE_T);
  TESTEQUALS(16UL,map->hashsize,PF_SIZE_T);
  TESTEQUALS(12,map->csize,"%d");
  TESTEQUALS(0,map->cidx,"%d");

  ss_ht_free(map);

  return 0;
}

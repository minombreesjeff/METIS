#include "test.h"

#define N 12

sint_t test(void) 
{
  sint_t i;

  ss_ht_t * map = ss_ht_create(N,N);

  ss_ht_put(0,10,map);
  ss_ht_put(18,20,map);
  ss_ht_put(19,30,map);
  ss_ht_put(122,40,map);
  i = ss_ht_put(10,50,map);

  TESTEQUALS(-1,i,PF_SINT_T);
  TESTEQUALS(5UL,map->size,PF_SIZE_T);

  i = ss_ht_put(10,60,map);

  TESTEQUALS(5UL,map->size,PF_SIZE_T);

  TESTEQUALS(50,i,PF_SINT_T);

  TESTEQUALS(0,map->elements[0].key,PF_SINT_T);
  TESTEQUALS(10,map->elements[0].val,PF_SINT_T);
  TESTEQUALS(-1,map->elements[0].next,"%d");

  TESTEQUALS(18,map->elements[2].key,PF_SINT_T);
  TESTEQUALS(20,map->elements[2].val,PF_SINT_T);
  TESTEQUALS(-1,map->elements[2].next,"%d");

  TESTEQUALS(19,map->elements[3].key,PF_SINT_T);
  TESTEQUALS(30,map->elements[3].val,PF_SINT_T);
  TESTEQUALS(-1,map->elements[3].next,"%d");

  TESTEQUALS(122,map->elements[10].key,PF_SINT_T);
  TESTEQUALS(40,map->elements[10].val,PF_SINT_T);
  TESTEQUALS(0,map->elements[10].next,"%d");

  TESTEQUALS(10,map->chain[0].key,PF_SINT_T);
  TESTEQUALS(60,map->chain[0].val,PF_SINT_T);
  TESTEQUALS(-1,map->chain[0].next,"%d");
  TESTEQUALS(1,map->cidx,"%d"); 

  ss_ht_free(map);

  return 0;
}

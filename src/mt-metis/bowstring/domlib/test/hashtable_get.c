#include "test.h"

#define N 12

sint_t test(void) 
{
  ss_ht_t * map = ss_ht_create(N,N);

  map->elements[0].key = 0;
  map->elements[0].val = 10;
  map->elements[0].next = -1;

  map->elements[2].key = 18;
  map->elements[2].val = 20;
  map->elements[0].next = -1;

  map->elements[3].key = 19;
  map->elements[3].val = 30;
  map->elements[0].next = -1;

  map->elements[10].key = 122;
  map->elements[10].val = 40;
  map->elements[10].next = 0;

  map->chain[0].key = 10;
  map->chain[0].val = 50;
  map->chain[0].next = -1;
  map->cidx = 1;

  map->size = 5;

  TESTEQUALS(10,ss_ht_get(0,map),PF_SINT_T);
  TESTEQUALS(20,ss_ht_get(18,map),PF_SINT_T);
  TESTEQUALS(30,ss_ht_get(19,map),PF_SINT_T);
  TESTEQUALS(40,ss_ht_get(122,map),PF_SINT_T);
  TESTEQUALS(50,ss_ht_get(10,map),PF_SINT_T);

  TESTEQUALS(-1,ss_ht_get(9,map),PF_SINT_T);
  TESTEQUALS(-1,ss_ht_get(16,map),PF_SINT_T);
  TESTEQUALS(-1,ss_ht_get(7,map),PF_SINT_T);
  TESTEQUALS(-1,ss_ht_get(57,map),PF_SINT_T);

  ss_ht_free(map);

  return 0;
}

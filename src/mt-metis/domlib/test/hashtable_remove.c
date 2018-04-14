#include "test.h"

#define N 12

sint_t test(void) 
{
  sint_t i;

  ss_ht_t * map = ss_ht_create(N,N);

  map->elements[0].key = 0;
  map->elements[0].val = 10;
  map->elements[0].next = -1;

  map->elements[2].key = 18;
  map->elements[2].val = 20;
  map->elements[2].next = -1;

  map->elements[3].key = 19;
  map->elements[3].val = 30;
  map->elements[3].next = -1;

  map->elements[10].key = 122;
  map->elements[10].val = 40;
  map->elements[10].next = 0;

  map->chain[0].key = 10;
  map->chain[0].val = 50;
  map->chain[0].next = 1;

  map->chain[1].key = 26;
  map->chain[1].val = 60;
  map->chain[1].next = -1;

  map->cidx = 2;

  map->size = 6;

  i = ss_ht_remove(18,map);

  TESTEQUALS(20,i,PF_SINT_T);
  TESTEQUALS(-1,map->elements[2].key,PF_SINT_T);
  TESTEQUALS(5UL,map->size,PF_SIZE_T);

  i = ss_ht_remove(0,map);

  TESTEQUALS(10,i,PF_SINT_T);
  TESTEQUALS(-1,map->elements[0].key,PF_SINT_T);
  TESTEQUALS(4UL,map->size,PF_SIZE_T);

  i = ss_ht_remove(19,map);

  TESTEQUALS(30,i,PF_SINT_T);
  TESTEQUALS(-1,map->elements[3].key,PF_SINT_T);
  TESTEQUALS(3UL,map->size,PF_SIZE_T);

  i = ss_ht_remove(10,map);

  TESTEQUALS(50,i,PF_SINT_T);
  TESTEQUALS(1,map->cidx,PF_SINT_T);
  TESTEQUALS(26,map->chain[0].key,PF_SINT_T);
  TESTEQUALS(60,map->chain[0].val,PF_SINT_T);
  TESTEQUALS(2UL,map->size,PF_SIZE_T);

  i = ss_ht_remove(122,map);

  TESTEQUALS(40,i,PF_SINT_T);
  TESTEQUALS(0,map->cidx,PF_SINT_T);
  TESTEQUALS(1UL,map->size,PF_SIZE_T);

  i = ss_ht_remove(26,map);

  TESTEQUALS(60,i,PF_SINT_T);
  TESTEQUALS(-1,map->elements[10].key,PF_SINT_T);
  TESTEQUALS(0UL,map->size,PF_SIZE_T);

  ss_ht_free(map);

  return 0;
}

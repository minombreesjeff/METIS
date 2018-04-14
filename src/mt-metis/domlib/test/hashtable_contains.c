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

  TESTTRUE(ss_ht_contains(0,map));
  TESTTRUE(ss_ht_contains(18,map));
  TESTTRUE(ss_ht_contains(19,map));
  TESTTRUE(ss_ht_contains(122,map));
  TESTTRUE(ss_ht_contains(10,map));

  TESTEQUALS(0,ss_ht_contains(9,map),"%d");
  TESTEQUALS(0,ss_ht_contains(16,map),"%d");
  TESTEQUALS(0,ss_ht_contains(7,map),"%d");
  TESTEQUALS(0,ss_ht_contains(57,map),"%d");

  ss_ht_free(map);

  return 0;
}

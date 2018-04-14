#include "test.h"


sint_t test(void)
{
  sint_list_t * list = sint_list_create();

  TESTEQUALS(list->head,NULL,"%p");
  TESTEQUALS(list->tail,NULL,"%p");
  TESTEQUALS(list->size,0UL,"%zu");

  sint_list_free(list);

  return 0;
}

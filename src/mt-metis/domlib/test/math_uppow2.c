#include "test.h"

sint_t test(void) 
{
  TESTEQUALS(1,sint_uppow2(0),PF_SINT_T);
  TESTEQUALS(1,sint_uppow2(1),PF_SINT_T);
  TESTEQUALS(2,sint_uppow2(2),PF_SINT_T);
  TESTEQUALS(4,sint_uppow2(3),PF_SINT_T);
  TESTEQUALS(4,sint_uppow2(4),PF_SINT_T);
  TESTEQUALS(8,sint_uppow2(5),PF_SINT_T);
  TESTEQUALS(16,sint_uppow2(12),PF_SINT_T);
  TESTEQUALS(128,sint_uppow2(123),PF_SINT_T);
  TESTEQUALS(2048,sint_uppow2(1230),PF_SINT_T);
  TESTEQUALS(4096,sint_uppow2(4096),PF_SINT_T);
  return 0;
}

#include "test.h"

sint_t test(void) 
{
  TESTEQUALS(1,sint_downpow2(0),PF_SINT_T);
  TESTEQUALS(1,sint_downpow2(1),PF_SINT_T);
  TESTEQUALS(2,sint_downpow2(2),PF_SINT_T);
  TESTEQUALS(2,sint_downpow2(3),PF_SINT_T);
  TESTEQUALS(4,sint_downpow2(4),PF_SINT_T);
  TESTEQUALS(4,sint_downpow2(5),PF_SINT_T);
  TESTEQUALS(8,sint_downpow2(12),PF_SINT_T);
  TESTEQUALS(64,sint_downpow2(123),PF_SINT_T);
  TESTEQUALS(1024,sint_downpow2(1230),PF_SINT_T);
  TESTEQUALS(4096,sint_downpow2(4096),PF_SINT_T);
  return 0;
}

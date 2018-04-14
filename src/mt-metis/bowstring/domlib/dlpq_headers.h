/**
 * @file dlpq_headers.h
 * @brief Priority queue function prototypes
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-10-04
 */


/* prefixing ugliness */
#define DLPQ_PRE2(prefix,suffix) prefix ## _ ## suffix
#define DLPQ_PRE1(prefix,suffix) DLPQ_PRE2(prefix,suffix)
#define DLPQ_PUB(name) DLPQ_PRE1(DLPQ_PREFIX,name)
#define DLPQ_PRI(name) DLPQ_PRE1(_,DLPQ_PRE1(DLPQ_PREFIX,name))




/******************************************************************************
* TYPES ***********************************************************************
******************************************************************************/


typedef struct DLPQ_PRI(pq_kv_t) {
  DLPQ_KEY_T key;
  DLPQ_VAL_T val;
} DLPQ_PRI(pq_kv_t);


typedef struct DLPQ_PUB(priority_queue_t) {
  DLPQ_PRI(pq_kv_t) * elements;
  size_t maxsize;
  size_t size;
  size_t * index;
  DLPQ_VAL_T min;
  DLPQ_VAL_T max;
} DLPQ_PUB(priority_queue_t);





#ifndef DLPQ_STATIC


/******************************************************************************
* FUNCTION PROTOTYPES *********************************************************
******************************************************************************/


DLPQ_PUB(priority_queue_t) * DLPQ_PUB(priority_queue_create)(
    DLPQ_VAL_T min,
    DLPQ_VAL_T max);


int DLPQ_PUB(maxpq_push)(
    DLPQ_KEY_T key, DLPQ_VAL_T val,
    DLPQ_PUB(priority_queue_t) * q);


int DLPQ_PUB(minpq_push)(
    DLPQ_KEY_T key, 
    DLPQ_VAL_T val,
    DLPQ_PUB(priority_queue_t) * q);


DLPQ_VAL_T DLPQ_PUB(maxpq_pop)(
    DLPQ_PUB(priority_queue_t) *q);


DLPQ_VAL_T DLPQ_PUB(minpq_pop)(
    DLPQ_PUB(priority_queue_t) *q);


DLPQ_VAL_T DLPQ_PUB(maxpq_peek)(
    DLPQ_PUB(priority_queue_t) const * q);


DLPQ_VAL_T DLPQ_PUB(minpq_peek)(
    DLPQ_PUB(priority_queue_t) const * q);


DLPQ_KEY_T DLPQ_PUB(maxpq_max)(
    DLPQ_PUB(priority_queue_t) const * q);


DLPQ_KEY_T DLPQ_PUB(minpq_min)(
    DLPQ_PUB(priority_queue_t) const * q);


int DLPQ_PUB(maxpq_update)(
    DLPQ_KEY_T p, 
    DLPQ_VAL_T v,
    DLPQ_PUB(priority_queue_t) * q);


int DLPQ_PUB(minpq_update)(
    DLPQ_KEY_T p, 
    DLPQ_VAL_T v,
    DLPQ_PUB(priority_queue_t) *q);


int DLPQ_PUB(maxpq_delete)(
    DLPQ_VAL_T v, 
    DLPQ_PUB(priority_queue_t) * q);


int DLPQ_PUB(minpq_delete)(
    DLPQ_VAL_T v, 
    DLPQ_PUB(priority_queue_t) * q);


int DLPQ_PUB(maxpq_clear)(
    DLPQ_PUB(priority_queue_t) * q);


int DLPQ_PUB(minpq_clear)(
    DLPQ_PUB(priority_queue_t) * q);


int DLPQ_PUB(maxpq_contains)(
    DLPQ_VAL_T v, 
    DLPQ_PUB(priority_queue_t) * q);


int DLPQ_PUB(minpq_contains)(
    DLPQ_VAL_T v, 
    DLPQ_PUB(priority_queue_t) * q);


void DLPQ_PUB(priority_queue_free)(
    DLPQ_PUB(priority_queue_t) * q);


#undef DLPQ_PRE2
#undef DLPQ_PRE1
#undef DLPQ_PUB
#undef DLPQ_PRI


#else


#undef DLPQ_PRE2
#undef DLPQ_PRE1
#undef DLPQ_PUB
#undef DLPQ_PRI


#define DLPQ_VISIBILITY static
#include "dlpq_funcs.h"
#undef DLPQ_VISIBILITY


#endif

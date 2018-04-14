/**
 * @file dlpq_funcs.h
 * @brief Functions for priority queues
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


#ifndef DLPQ_VISIBILITY
  #define DLPQ_VISIBILITY
  #define DLPQ_DEFVIS
#endif




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define __DL_HEAP_EMPTY_SLOT (-1)
#define __LEFTCHILDINDEX(i) ((2*i)+1)
#define __RIGHTCHILDINDEX(i) ((2*i)+2)
#define __PARENTINDEX(i) ((i-1)>>1)
#define __ISLESS(a,b) (a < b)
#define __ISMORE(a,b) (a > b)
#define __ISLESSKV(a,b) (a.key < b.key)
#define __ISMOREKV(a,b) (a.key > b.key)
#define __ADJINDEX(i,q) \
  (((ssize_t)i) - ((ssize_t)((q)->min)))


#define __MOVEINDEXLESS(i,j,h) \
  do { \
    h->elements[i] = h->elements[j]; \
    i = j; \
  } while (0)


#define __MOVEINDEXFUL(i,j,h) \
  do { \
    h->elements[i] = h->elements[j]; \
    h->index[__ADJINDEX(h->elements[j].val,h)] = i; \
    i = j; \
  } while (0)


#define __SETINDEXLESS(i,v,h) \
  do { \
    h->elements[i] = v; \
  } while (0)


#define __SETINDEXFUL(i,v,h) \
  do { \
    h->elements[i] = v; \
    h->index[__ADJINDEX(v.val,h)] = i; \
  } while (0)




/******************************************************************************
* META MACROS *****************************************************************
******************************************************************************/


#define __HEAP_FIX_DOWN(i,val,heap,isorder,move,set) \
  do { \
    size_t k,j; \
    while (1) { \
      set(i,val,heap); \
      j = __LEFTCHILDINDEX(i); \
      k = __RIGHTCHILDINDEX(i); \
      if (j < heap->size) { \
        if (k < heap->size && !isorder(val,heap->elements[k]) &&  \
            !isorder(heap->elements[j],heap->elements[k])) { \
          move(i,k,heap); \
        } else if (!isorder(val,heap->elements[j])) { \
          move(i,j,heap); \
        } else { \
          break; \
        } \
      } else { \
        break; \
      } \
    } \
  } while(0)


#define __HEAP_FIX_UP(i,val,heap,isorder,move,set) \
  do { \
    size_t j =__PARENTINDEX(i); \
    while (i > 0 && !isorder(heap->elements[j],val)) { \
      move(i,j,heap); \
      j=__PARENTINDEX(i); \
    } \
    set(i,val,heap); \
  } while (0)




/******************************************************************************
* MEMORY FUNCTIONS ************************************************************
******************************************************************************/


#define DLMEM_PREFIX DLPQ_PRI(pq_kv)
#define DLMEM_TYPE_T DLPQ_PRI(pq_kv_t)
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX


#define DLMEM_PREFIX DLPQ_PUB(priority_queue)
#define DLMEM_TYPE_T DLPQ_PUB(priority_queue_t)
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_TYPE_T
#undef DLMEM_PREFIX




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static size_t const DLPQ_PRI(pq_empty_slot) = (size_t)-1;




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


DLPQ_VISIBILITY DLPQ_PUB(priority_queue_t) * DLPQ_PUB(priority_queue_create)(
    DLPQ_VAL_T const min, 
    DLPQ_VAL_T const max)
{
  DLPQ_PUB(priority_queue_t) * q;

  size_t const n = (size_t)(((ssize_t)max) - ((ssize_t)min));

  q = DLPQ_PUB(priority_queue_alloc)(1);

  q->min = min;
  q->max = max;
  q->maxsize = n;
  q->size = 0;
  q->elements = DLPQ_PRI(pq_kv_alloc)(q->maxsize);
  q->index = size_init_alloc(DLPQ_PRI(pq_empty_slot),q->maxsize);

  return q;
}


DLPQ_VISIBILITY int DLPQ_PUB(maxpq_push)(
    DLPQ_KEY_T const key, 
    DLPQ_VAL_T const val, 
    DLPQ_PUB(priority_queue_t) * const q)
{
  DL_ASSERT(val >= q->min,"Attempting to add value below minimum to "
      "priority queue\n");
  DL_ASSERT(val < q->max,"Attempting to add value above maximum to "
      "priority queue\n");
  DL_ASSERT(q->index[__ADJINDEX(val,q)] == DLPQ_PRI(pq_empty_slot),
      "Value already exists in priority queue\n");

  DLPQ_PRI(pq_kv_t) kv;

  kv.key = key;
  kv.val = val;
  size_t i = q->size++;

  __HEAP_FIX_UP(i,kv,q,__ISMOREKV,__MOVEINDEXFUL,__SETINDEXFUL);

  return 1;
}


DLPQ_VISIBILITY int DLPQ_PUB(minpq_push)(
    DLPQ_KEY_T const key, 
    DLPQ_VAL_T const val, 
    DLPQ_PUB(priority_queue_t) * const q)
{
  DL_ASSERT(val >= q->min,"Attempting to add value below minimum to "
      "priority queue\n");
  DL_ASSERT(val < q->max,"Attempting to add value above maximum to "
      "priority queue\n");
  DL_ASSERT(q->index[__ADJINDEX(val,q)] == DLPQ_PRI(pq_empty_slot),
      "Value already exists in priority queue\n");

  DLPQ_PRI(pq_kv_t) kv;

  kv.key = key;
  kv.val = val;
  size_t i = q->size++;

  __HEAP_FIX_UP(i,kv,q,__ISLESSKV,__MOVEINDEXFUL,__SETINDEXFUL);

  return 1;
}


DLPQ_VISIBILITY DLPQ_VAL_T DLPQ_PUB(maxpq_pop)(
    DLPQ_PUB(priority_queue_t) * const q)
{
  DL_ASSERT(q->size > 0,"Trying to pop() from an empty queue\n");

  size_t i;
  
  i = 0;
  DLPQ_PRI(pq_kv_t) kv = q->elements[i];
  DLPQ_PRI(pq_kv_t) val = q->elements[--q->size];

  __HEAP_FIX_DOWN(i,val,q,__ISMOREKV,__MOVEINDEXFUL,__SETINDEXFUL);
  q->index[__ADJINDEX(kv.val,q)] = DLPQ_PRI(pq_empty_slot);

  return kv.val;
}


DLPQ_VISIBILITY DLPQ_VAL_T DLPQ_PUB(minpq_pop)(
    DLPQ_PUB(priority_queue_t) * const q)
{
  DL_ASSERT(q->size > 0,"Trying to pop() from an empty queue\n");

  size_t i;

  i = 0;
  DLPQ_PRI(pq_kv_t) kv = q->elements[i];
  DLPQ_PRI(pq_kv_t) val = q->elements[--q->size];

  __HEAP_FIX_DOWN(i,val,q,__ISLESSKV,__MOVEINDEXFUL,__SETINDEXFUL);
  q->index[__ADJINDEX(kv.val,q)] = DLPQ_PRI(pq_empty_slot);
  return kv.val;
}


DLPQ_VISIBILITY DLPQ_VAL_T DLPQ_PUB(maxpq_peek)(
    DLPQ_PUB(priority_queue_t) const * const q)
{
  DL_ASSERT(q->size > 0,"Attempting to look at first element of an empty " \
      "priority queue\n");

  return q->elements[0].val;
}


DLPQ_VISIBILITY DLPQ_VAL_T DLPQ_PUB(minpq_peek)(
    DLPQ_PUB(priority_queue_t) const * const q)
{
  DL_ASSERT(q->size > 0,"Attempting to look at first element of an empty " \
      "priority queue\n");

  return q->elements[0].val;
}


DLPQ_VISIBILITY DLPQ_KEY_T DLPQ_PUB(maxpq_max)(
    DLPQ_PUB(priority_queue_t) const * const q)
{
  DL_ASSERT(q->size > 0,"Attempting to look at first element of an empty " \
      "priority queue\n");

  return q->elements[0].key;
}


DLPQ_VISIBILITY DLPQ_KEY_T DLPQ_PUB(minpq_min)(
    const DLPQ_PUB(priority_queue_t) * q)
{
  DL_ASSERT(q->size > 0,"Attempting to look at first element of an empty " \
      "priority queue\n");
  return q->elements[0].key;
}


DLPQ_VISIBILITY int DLPQ_PUB(maxpq_update)(
    DLPQ_KEY_T const key, 
    DLPQ_VAL_T const val, 
    DLPQ_PUB(priority_queue_t) * const q)
{
  DL_ASSERT(q->index[__ADJINDEX(val,q)] != DLPQ_PRI(pq_empty_slot), \
      "Can't update value not in priority queue\n");

  size_t i; 
  DLPQ_PRI(pq_kv_t) kv;

  i = q->index[__ADJINDEX(val,q)];
  kv.key = key;
  kv.val = val;

  if (__ISMOREKV(kv,q->elements[i])) { /* goes up */
    __HEAP_FIX_UP(i,kv,q,__ISMOREKV,__MOVEINDEXFUL,__SETINDEXFUL);
  } else {
    __HEAP_FIX_DOWN(i,kv,q,__ISMOREKV,__MOVEINDEXFUL,__SETINDEXFUL);
  } /* otherwise do nothing */

  return 1;
}


DLPQ_VISIBILITY int DLPQ_PUB(minpq_update)(
    DLPQ_KEY_T const key, 
    DLPQ_VAL_T const val, 
    DLPQ_PUB(priority_queue_t) * const q)
{
  DL_ASSERT(q->index[__ADJINDEX(val,q)] != DLPQ_PRI(pq_empty_slot), \
      "Can't update value not in priority queue\n");

  size_t i;
  DLPQ_PRI(pq_kv_t) kv;

  i = q->index[__ADJINDEX(val,q)];
  kv.key = key;
  kv.val = val;

  if (__ISLESSKV(kv,q->elements[i])) { /* goes up */
    __HEAP_FIX_UP(i,kv,q,__ISLESSKV,__MOVEINDEXFUL,__SETINDEXFUL);
  } else {
    __HEAP_FIX_DOWN(i,kv,q,__ISLESSKV,__MOVEINDEXFUL,__SETINDEXFUL);
  } /* otherwise do nothing */

  return 1;
}


DLPQ_VISIBILITY int DLPQ_PUB(maxpq_delete)(
    DLPQ_VAL_T const v,
    DLPQ_PUB(priority_queue_t) * const q)
{
  DL_ASSERT(q->index[__ADJINDEX(v,q)] !=  DLPQ_PRI(pq_empty_slot), \
      "Can't delete value not in priority queue\n");

  size_t i;
  DLPQ_PRI(pq_kv_t) val;

  i = q->index[__ADJINDEX(v,q)];
  val = q->elements[--q->size];

  if (val.key > q->elements[i].key) {
    __HEAP_FIX_UP(i,val,q,__ISMOREKV,__MOVEINDEXFUL,__SETINDEXFUL);
  } else {
    __HEAP_FIX_DOWN(i,val,q,__ISMOREKV,__MOVEINDEXFUL,__SETINDEXFUL);
  }
  q->index[__ADJINDEX(v,q)] = DLPQ_PRI(pq_empty_slot);

  return 1;
}


DLPQ_VISIBILITY int DLPQ_PUB(minpq_delete)(
    DLPQ_VAL_T const v,
    DLPQ_PUB(priority_queue_t) * const q)
{
  DL_ASSERT(q->index[__ADJINDEX(v,q)] !=  DLPQ_PRI(pq_empty_slot), \
      "Can't delete value not in priority queue\n");

  size_t i;
  DLPQ_PRI(pq_kv_t) val;

  i = q->index[__ADJINDEX(v,q)];
  val = q->elements[--q->size];
  
  if (val.key < q->elements[i].key) {
    __HEAP_FIX_UP(i,val,q,__ISLESSKV,__MOVEINDEXFUL,__SETINDEXFUL);
  } else {
    __HEAP_FIX_DOWN(i,val,q,__ISLESSKV,__MOVEINDEXFUL,__SETINDEXFUL);
  }
  q->index[__ADJINDEX(v,q)] = DLPQ_PRI(pq_empty_slot);

  return 1;
}


DLPQ_VISIBILITY int DLPQ_PUB(maxpq_clear)(
    DLPQ_PUB(priority_queue_t) * const q)
{
  DLPQ_VAL_T v;

  if (q->size > 0.125 * q->maxsize) {
    size_set(q->index,DLPQ_PRI(pq_empty_slot),q->maxsize);
    q->size = 0;
  } else {
    while (q->size>0) {
      v = q->elements[--q->size].val;
      q->index[__ADJINDEX(v,q)] = DLPQ_PRI(pq_empty_slot);
    }
  }

  return 1;
}


DLPQ_VISIBILITY int DLPQ_PUB(minpq_clear)(
    DLPQ_PUB(priority_queue_t) * const q) 
{
  DLPQ_VAL_T v;

  if (q->size > 0.125 * q->maxsize) {
    size_set(q->index,DLPQ_PRI(pq_empty_slot),q->maxsize);
    q->size = 0;
  } else {
    while (q->size>0) {
      v = q->elements[--q->size].val;
      q->index[__ADJINDEX(v,q)] = DLPQ_PRI(pq_empty_slot);
    }
  }

  return 1;
}


DLPQ_VISIBILITY int DLPQ_PUB(maxpq_contains)(
    DLPQ_VAL_T const v, 
    DLPQ_PUB(priority_queue_t) * const q)
{
  return q->index[__ADJINDEX(v,q)] != DLPQ_PRI(pq_empty_slot);
}


DLPQ_VISIBILITY int DLPQ_PUB(minpq_contains)(
    DLPQ_VAL_T const v, 
    DLPQ_PUB(priority_queue_t) * const q)
{
  return q->index[__ADJINDEX(v,q)] != DLPQ_PRI(pq_empty_slot);
}


DLPQ_VISIBILITY void DLPQ_PUB(priority_queue_free)(
    DLPQ_PUB(priority_queue_t) * q)
{
  dl_free(q->elements);
  dl_free(q->index);
  dl_free(q);
}


#ifdef DLPQ_DEFVIS
  #undef DLPQ_VISIBILITY
  #undef DLPQ_DEFVIS
#endif


#undef DLPQ_PRE2
#undef DLPQ_PRE1
#undef DLPQ_PUB
#undef DLPQ_PRI




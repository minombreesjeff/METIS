/**
 * @file dlht_funcs.h
 * @brief Hash table functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-10-04
 */




/* prefixing ugliness */
#define DLHT_PRE2(prefix,suffix) prefix ## _ ## suffix
#define DLHT_PRE1(prefix,suffix) DLHT_PRE2(prefix,suffix)
#define DLHT_PUB(name) DLHT_PRE1(DLHT_PREFIX,name)
#define DLHT_PRI(name) DLHT_PRE1(_,DLHT_PRE1(DLHT_PREFIX,name))


#ifndef DLHT_VISIBILITY
  #define DLHT_DEFVIS
  #define DLHT_VISIBILITY
#endif



/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/

static const DLHT_KEY_T DLHT_PUB(ht_null_key) = (DLHT_KEY_T)-1;
static const DLHT_VAL_T DLHT_PUB(ht_null_value) = (DLHT_VAL_T)-1;
static const DLHT_VAL_T DLHT_PUB(ht_fail_value) = (DLHT_VAL_T)-2;



/******************************************************************************
* MEMORY FUNCTIONS ************************************************************
******************************************************************************/


#define DLMEM_PREFIX DLHT_PRI(ht_kv)
#define DLMEM_TYPE_T DLHT_PRI(ht_kv_t)
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T


#define DLMEM_PREFIX DLHT_PUB(ht)
#define DLMEM_TYPE_T DLHT_PUB(ht_t)
#define DLMEM_STATIC
#include "dlmem_headers.h"
#undef DLMEM_STATIC
#undef DLMEM_PREFIX
#undef DLMEM_TYPE_T




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void DLHT_PRI(ht_fixchain)(const int cidx, 
    DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(cidx >= 0, "Bad chain index passed to fixchain\n");
  size_t hk;
  DLHT_PRI(ht_kv_t) * oldkv;
  --map->cidx;
  if (cidx != map->cidx) {
    map->chain[cidx] = map->chain[map->cidx];
    hk = ((size_t)map->chain[cidx].key) & map->hashmask;
    /* I know that hk has a chain */
    oldkv = map->elements+hk;
    while (oldkv->next != map->cidx) {
      DL_ASSERT(oldkv->next >= 0,"Broken chain encountered when "
          "performing fix\n");
      DL_ASSERT(oldkv->next < map->csize,"Out of range chain index when "
          "fixing chains\n");
      oldkv = map->chain+oldkv->next;
    }
    oldkv->next = cidx;
  }
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/

static const int __NULL_NEXT = -1;



DLHT_VISIBILITY DLHT_PUB(ht_t) * DLHT_PUB(ht_create)(const size_t size,
    const int csize)
{
  size_t i;
  DLHT_PUB(ht_t) * const map = DLHT_PUB(ht_alloc)(1);
  map->size = 0;
  map->hashsize = size_uppow2(size);
  map->maxsize = map->hashsize;
  map->hashmask = map->hashsize - 1;
  map->elements = DLHT_PRI(ht_kv_alloc)(map->maxsize);
  map->csize = csize;
  map->cidx = 0;
  map->chain = DLHT_PRI(ht_kv_alloc)(map->csize);
  /* empty my hash table */
  for (i=0;i<map->maxsize;++i) {
    map->elements[i].key = DLHT_PUB(ht_null_key);
  }
  return map;
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_get)(const DLHT_KEY_T key,
    const DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to search a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to search a hashmap "
      "with a null key\n");
  int cidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* life is good */
    return map->elements[hk].val;
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* we don't have the key */
    return DLHT_PUB(ht_null_value);
  } else {
    /* search the stupid chain */
    cidx = map->elements[hk].next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        return map->chain[cidx].val;
      }
      cidx = map->chain[cidx].next;
    }
    return DLHT_PUB(ht_null_value);
  }
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_put)(
    DLHT_KEY_T const key, 
    DLHT_VAL_T const val, 
    DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to search a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to set a null key "
      "in the hasmap\n");
  DL_ASSERT(val != DLHT_PUB(ht_null_value),"Attempt to set a null "
      "value in the hashmap\n");

  int cidx;
  DLHT_PRI(ht_kv_t) * oldkv;
  DLHT_VAL_T oldval;
  size_t hk;

  hk = ((size_t)key) & map->hashmask;

  if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* new kv in empty slot */
    map->elements[hk].key = key;
    map->elements[hk].val = val;
    map->elements[hk].next = __NULL_NEXT;
    ++map->size;
    return DLHT_PUB(ht_null_value);
  } else if (map->elements[hk].key == key) {
    /* replace old kv */
    oldval = map->elements[hk].val;
    map->elements[hk].val = val;
    return oldval;
  } else {
    /* search the stupid chain */
    oldkv = map->elements+hk;
    cidx = oldkv->next;

    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        /* replace old kv */
        oldval = map->chain[cidx].val;
        map->chain[cidx].val = val;
        return oldval;
      }
      oldkv = map->chain+cidx;
      cidx = map->chain[cidx].next;
    }

    if (map->cidx < map->csize) {
      /* add new kv to chain */
      cidx = (oldkv->next = map->cidx++);
      map->chain[cidx].key = key;
      map->chain[cidx].val = val;
      map->chain[cidx].next = __NULL_NEXT;
      ++map->size;
      return DLHT_PUB(ht_null_value);
    } else {
      /* our chain is full return an error */
      return DLHT_PUB(ht_fail_value);
    }
  }
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_min)(const DLHT_KEY_T key, 
    const DLHT_VAL_T av, DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to search a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to search a hashmap "
      "with a null key\n");
  int cidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* life is good */
    if (map->chain[hk].val > av) {
      map->chain[hk].val = av;
    }
    return map->chain[hk].val;
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* we don't have the key */
    DLHT_PUB(ht_put)(key,av,map);
    return av;
  } else {
    /* search the stupid chain */
    cidx = map->elements[hk].next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        if (map->chain[cidx].val > av) {
          map->chain[cidx].val = av;
        }
        return map->chain[cidx].val;
      }
      cidx = map->chain[cidx].next;
    }
    return DLHT_PUB(ht_null_value);
  }
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_max)(const DLHT_KEY_T key, 
    const DLHT_VAL_T av, DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to search a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to search a hashmap "
      "with a null key\n");
  int cidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* life is good */
    if (map->chain[hk].val < av) {
      map->chain[hk].val = av;
    }
    return map->chain[hk].val;
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* we don't have the key */
    DLHT_PUB(ht_put)(key,av,map);
    return av;
  } else {
    /* search the stupid chain */
    cidx = map->elements[hk].next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        if (map->chain[cidx].val < av) {
          map->chain[cidx].val = av;
        }
        return map->chain[cidx].val;
      }
      cidx = map->chain[cidx].next;
    }
    return DLHT_PUB(ht_null_value);
  }
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_add)(const DLHT_KEY_T key, 
    const DLHT_VAL_T av, DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to search a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to search a hashmap "
      "with a null key\n");
  int cidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* life is good */
    return (map->elements[hk].val += av);
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* we don't have the key */
    DLHT_PUB(ht_put)(key,av,map);
    return av;
  } else {
    /* search the stupid chain */
    cidx = map->elements[hk].next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        return (map->chain[cidx].val += av);
      }
      cidx = map->chain[cidx].next;
    }
    return DLHT_PUB(ht_null_value);
  }
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_multiply)(const DLHT_KEY_T key, 
    const DLHT_VAL_T av, DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to search a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to search a hashmap "
      "with a null key\n");
  int cidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* life is good */
    return (map->elements[hk].val *= av);
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* we don't have the key */
    DLHT_PUB(ht_put)(key,av,map);
    return av;
  } else {
    /* search the stupid chain */
    cidx = map->elements[hk].next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        return (map->chain[cidx].val *= av);
      }
      cidx = map->chain[cidx].next;
    }
    return DLHT_PUB(ht_null_value);
  }
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_and)(const DLHT_KEY_T key, 
    const DLHT_VAL_T av, DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to search a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to search a hashmap "
      "with a null key\n");
  int cidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* life is good */
    return (map->elements[hk].val = 
        (((int)map->elements[hk].val) && ((int)av)));
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* we don't have the key */
    DLHT_PUB(ht_put)(key,av,map);
    return av;
  } else {
    /* search the stupid chain */
    cidx = map->elements[hk].next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        return (map->chain[cidx].val =
            (((int)map->chain[cidx].val) && ((int)av)));
      }
      cidx = map->chain[cidx].next;
    }
    return DLHT_PUB(ht_null_value);
  }
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_or)(const DLHT_KEY_T key, 
    const DLHT_VAL_T av, DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to search a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to search a hashmap "
      "with a null key\n");
  int cidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* life is good */
    return (map->elements[hk].val = 
        ((int)map->elements[hk].val) || ((int)map->elements[hk].val));
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* we don't have the key */
    DLHT_PUB(ht_put)(key,av,map);
    return av;
  } else {
    /* search the stupid chain */
    cidx = map->elements[hk].next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        return (map->chain[cidx].val =
            (((int)map->chain[cidx].val) || ((int)av)));
      }
      cidx = map->chain[cidx].next;
    }
    return DLHT_PUB(ht_null_value);
  }
}


DLHT_VISIBILITY DLHT_VAL_T DLHT_PUB(ht_remove)(const DLHT_KEY_T key,
    DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to remove from a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to remove a null "
      "key in the hasmap\n");
  DLHT_VAL_T oldval;
  DLHT_PRI(ht_kv_t) * oldkv;
  int cidx,nidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* key in the table -- remove it */
    oldval = map->elements[hk].val;
    if ((cidx = map->elements[hk].next) >= 0) { 
      map->elements[hk] = map->chain[cidx];
      DLHT_PRI(ht_fixchain)(cidx,map);
    } else {
      map->elements[hk].key = DLHT_PUB(ht_null_key);
    }
    --map->size;
    return oldval;
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* key not in the table */
    return DLHT_PUB(ht_null_value);
  } else {
    /* search the stupid chain */
    oldkv = map->elements+hk;
    cidx = oldkv->next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        /* delete the kv */
        oldval = map->chain[cidx].val;
        nidx = map->chain[cidx].next;
        if (nidx >= 0) {
          map->chain[cidx] = map->chain[nidx];
          DLHT_PRI(ht_fixchain)(nidx,map);
        } else {
          oldkv->next = __NULL_NEXT;
        }
        --map->size;
        return oldval;
      }
      oldkv = map->chain+cidx;
      cidx = map->chain[cidx].next;
    }
    return DLHT_PUB(ht_null_value);
  }
}


DLHT_VISIBILITY int DLHT_PUB(ht_contains)(const DLHT_KEY_T key,
    const DLHT_PUB(ht_t) * const map)
{
  DL_ASSERT(map != NULL,"Attempt to remove from a NULL hashmap\n");
  DL_ASSERT(key != DLHT_PUB(ht_null_key),"Attempt to remove a null "
      "key in the hasmap\n");
  int cidx;
  size_t hk = ((size_t)key) & map->hashmask;
  if (map->elements[hk].key == key) {
    /* life is good */
    return 1;
  } else if (map->elements[hk].key == DLHT_PUB(ht_null_key)) {
    /* we don't have the key */
    return 0;
  } else {
    /* search the stupid chain */
    cidx = map->elements[hk].next;
    while (cidx >= 0) {
      if (map->chain[cidx].key == key) {
        return 1;
      }
      cidx = map->chain[cidx].next;
    }
    return 0;
  }
}


DLHT_VISIBILITY size_t DLHT_PUB(ht_clear)(DLHT_PUB(ht_t) * const map)
{
  size_t i, oldsize;
  oldsize = map->size;
  for (i=0;i<map->hashsize;++i) {
    map->elements[i].key = DLHT_PUB(ht_null_key);
  }
  map->size = 0;
  map->cidx = 0;
  return oldsize;
}


DLHT_VISIBILITY void DLHT_PUB(ht_clear_chains)(DLHT_PUB(ht_t) * map)
{
  map->cidx = 0;
}


DLHT_VISIBILITY void DLHT_PUB(ht_clear_slot)(DLHT_KEY_T key, 
    DLHT_PUB(ht_t) * map)
{
  size_t hk = ((size_t)key) & map->hashmask;
  map->elements[hk].key = DLHT_PUB(ht_null_key);
  map->elements[hk].next = __NULL_NEXT;
}


DLHT_VISIBILITY int DLHT_PUB(ht_adjust_size)(const size_t newsize,
    DLHT_PUB(ht_t) * const map)
{
  if (map->size != 0 || map->maxsize < newsize) {
    return 0;
  } else {
    map->hashsize = size_uppow2(newsize);
    map->hashmask = map->hashsize - 1;
    return 1;
  }
}


DLHT_VISIBILITY int DLHT_PUB(ht_free)(DLHT_PUB(ht_t) * map)
{
  dl_free(map->elements);
  dl_free(map->chain);
  dl_free(map);
  return 1;
}



#ifdef DLHT_DEFVIS
  #undef DLHT_DEFVIS
  #undef DLHT_VISIBILITY
#endif



#undef DLHT_PRE1
#undef DLHT_PRE2
#undef DLHT_PUB
#undef DLHT_PRI



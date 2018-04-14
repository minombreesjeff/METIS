/**
 * @file util.c
 * @brief Utility functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2012, Regents of the University of Minnesota
 * @version 1
 * @date 2012-07-10
 */

#include "includes.h"

idx_t irandInRange_r(const idx_t max, idx_t * seed) 
{
  if (max <= 0) {
    return 0;
  } else {
    return rand_r(seed)%max;
  }
}

void irandArrayPermute_r(const idx_t n, idx_t * p, const idx_t nshuffles,
    const idx_t flag, idx_t*seed) {
  idx_t i,u,v;
  idx_t tmp;
  if (flag) {
    for (i=0;i<n;++i) {
      p[i] = i;
    }
  }
  if (n > 4) {
    for (i=0;i<nshuffles; ++i) {
      v = irandInRange_r(n-3,seed);
      u = irandInRange_r(n-3,seed);
      gk_SWAP(p[v+0],p[u+2],tmp);
      gk_SWAP(p[v+1],p[u+3],tmp);
      gk_SWAP(p[v+2],p[u+0],tmp);
      gk_SWAP(p[v+3],p[u+1],tmp);
    }
  }
}

void BucketSortKeysInc_t(idx_t * const counts, idx_t n, idx_t max, idx_t *keys, 
         idx_t *tperm, idx_t *perm)
{
  idx_t i,ii;

  for (i=0; i<n; i++) {
    counts[keys[i]]++;
  }
  MAKECSR(i, max+1, counts);

  for (ii=0; ii<n; ii++) {
    i = tperm[ii];
    perm[counts[keys[ii]]++] = i;
  }
}




/**
* @file macros.h
* @brief Macros and some micro-Macros...
* @author Dominique LaSalle <lasalle@cs.umn.edu>
* Copyright 2012, Regents of the University of Minnesota
* @version 1
* @date 2012-06-13
*/


#ifndef MACROS_H
#define MACROS_H

#define ONES 0xFFFFFFFF

/***************************************************************************
* VERTEX DISTRIBUTION MACROS ***********************************************
***************************************************************************/
#define GBID_2_THRID(gbid,dmask) \
  ( (gbid) & (dmask) )

#define GBID_2_LBID(gbid,dshift) \
  ( (gbid) >> (dshift) )

#define GBID_2_THRID(gbid,dmask) \
  ( (gbid) & (dmask) )

#define LBID_2_GBID(lbid,thrid,dshift) \
  ( ((lbid) << (dshift)) + (thrid) )

#define GVTX_2_GBID(gvtx) \
  ( (gvtx)>>(BLOCKSHIFT) )

#define GVTX_2_THRID(gvtx,dmask) \
  ( GBID_2_THRID(GVTX_2_GBID(gvtx),dmask) )

#define GVTX_2_LBID(gvtx,dshift) \
  ( GBID_2_LBID(GVTX_2_GBID(gvtx),dshift) )

#define GVTX_2_LVTX(gvtx,dshift) \
  ( ((gvtx)&BLOCKMASK) + (GVTX_2_LBID(gvtx,dshift)<<BLOCKSHIFT) )

#define LVTX_2_LBID(lvtx) \
  ( (lvtx)>>BLOCKSHIFT )

#define LVTX_2_GBID(lvtx,thrid,dshift) \
  ( LBID_2_GBID(LVTX_2_LBID(lvtx),thrid,dshift) )

#define LVTX_2_GVTX(lvtx,thrid,dshift) \
  ( ((lvtx) & BLOCKMASK) + \
    (LVTX_2_GBID(lvtx,thrid,dshift)<<BLOCKSHIFT) )

/******************************************************************************
* COARSENING MACROS ***********************************************************
******************************************************************************/
#define MY_CVTX(v,u,dv,du) \
  ( ((v) == (u)) || ((dv) > (du)) || \
    ( ( ((dv) == (du)) && \
     ( ((((v)+(u))%2==0)&&((v)>(u))) || \
       ((((v)+(u))%2==1)&&((v)<(u))) ) ) ) ) 

#define skipEdge(v,u) \
  if ((v) == (u)) continue

#define cEdgeMask(start,j,jj,k,kk,m,nedges,mask,htable,adjncy,adjwgt, \
    cadjncy,cadjwgt) \
  do { \
    kk = k&mask; \
    if ((m = htable[kk]) == -1) { \
      cadjncy[nedges] = k; \
      cadjwgt[nedges] = adjwgt[j]; \
      htable[kk]      = nedges++; \
    } else if (cadjncy[m] == k) { \
      cadjwgt[m] += adjwgt[j]; \
    } else { \
      for (jj=start; jj<nedges; ++jj) { \
        if (cadjncy[jj] == k) { \
          cadjwgt[jj] += adjwgt[j]; \
          break; \
        } \
      } \
      if (jj == nedges) { \
        cadjncy[nedges]   = k; \
        cadjwgt[nedges++] = adjwgt[j]; \
      } \
    } \
  } while(0)

#define cEdgeNoMask(j,k,m,nedges,htable,adjncy,adjwgt,cadjncy,cadjwgt) \
  do { \
    if ((m = htable[k]) == -1) { \
      cadjncy[nedges] = k; \
      cadjwgt[nedges] = adjwgt[j]; \
      htable[k] = nedges++; \
    } else { \
      cadjwgt[m] += adjwgt[j]; \
    } \
  } while(0)

/******************************************************************************
* OMP/PARALLEL MACROS *********************************************************
******************************************************************************/
#define INIT_PARALLEL() \
  const idx_t myid = omp_get_thread_num(); \
  const idx_t nthreads = omp_get_num_threads() \
  /* semi colon dropped here */

#define startwctimer(tmr) \
  do { \
    _Pragma ("omp master") \
    { \
      gk_startwctimer(tmr); \
    } \
  } while (0)

#define stopwctimer(tmr) \
  do { \
    _Pragma ("omp master") \
    { \
      gk_stopwctimer(tmr); \
    } \
  } while (0)

#define startcputimer(tmr) \
  do { \
    _Pragma ("omp master") \
    { \
      gk_startcputimer(tmr); \
    } \
  } while (0)

#define stopcputimer(tmr) \
  do { \
    _Pragma ("omp master") \
    { \
      gk_stopcputimer(tmr); \
    } \
  } while (0)

#define pprintf(fmt, ...) \
  do { \
    _Pragma ("omp master") \
    { \
      printf(fmt, ##__VA_ARGS__); \
    } \
  } while (0)

#ifdef DEBUG
  #define dprintf(fmt, ...) \
  do { \
    _Pragma("omp master") \
    { \
      printf(fmt, ##__VA_ARGS__); \
    } \
  } while (0)
#else
  #define dprintf(fmt, ...)
#endif

/* omp junk */
/* return the id of the thread with the minimum value */
#define dl_omp_minreduce_index(myid,idx,val,buffer,nthreads) \
  do { \
    _Pragma ("omp barrier") \
    buffer[myid] = val; \
    idx = myid; \
    _Pragma ("omp barrier") \
    _Pragma ("omp master") \
    { \
      idx_t _i; \
      for (_i=0;_i<nthreads;++_i) { \
        if (buffer[_i] < val || (buffer[_i] == val && _i < idx)) { \
          idx = _i; \
          val = buffer[_i]; \
        } \
      } \
      buffer[0] = idx; \
    } \
    _Pragma("omp barrier") \
    (idx) = buffer[0]; \
  } while(0) 

/* return the id of the thread with the maximum value */ \
#define dl_omp_maxreduce_index(values, maxidx, nthreads) \
  do { \
    _Pragma ("omp barrier") \
    idx_t _nbrid,_i,_j,_nidx,_midx; \
    const idx_t _n = uplog(nthreads); \
    const idx_t _myid = omp_get_thread_num(); \
    _midx = maxidx[_myid] = _myid; \
    for (_j=0;_j<(1<<_n == nthreads ? 1 : 2);++_j) { \
      for (_i=0;_i<_n;++_i) { \
        _nbrid = _myid ^ (1<<_i); \
        _Pragma ("omp barrier") \
        if (_nbrid < nthreads) { \
          _nidx = maxidx[_nbrid]; \
          if (values[_midx] < values[_nidx] || \
              (values[_midx] == values[_nidx] && _midx > _nidx)) { \
            _midx = maxidx[_myid] = _nidx; \
          } \
        } \
      } \
    } \
    _Pragma ("omp barrier") \
  } while(0) 

#define dl_omp_maxreduce(myid,val,buffer,nthreads) \
  do { \
    _Pragma("omp barrier") \
    idx_t _i; \
    (buffer)[myid] = (val); \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      for (_i=0;_i<(nthreads);++_i) { \
        if ((val) < (buffer)[_i]) { \
          (val) = (buffer)[_i]; \
        } \
      } \
      (buffer)[0] = (val); \
    } \
    _Pragma("omp barrier") \
    (val) = (buffer)[0]; \
  } while(0) 

#define dl_omp_minreduce(myid,val,buffer,nthreads) \
  do { \
    _Pragma("omp barrier") \
    idx_t _i; \
    (buffer)[myid] = (val); \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    _Pragma("omp master") \
    { \
      for (_i=0;_i<(nthreads);++_i) { \
        if ((val) < (buffer)[_i]) { \
          (val) = (buffer)[_i]; \
        } \
      } \
      (buffer)[0] = (val); \
    } \
    _Pragma("omp barrier") \
    (val) = (buffer)[0]; \
  } while(0) 

#define dl_omp_maxreduce_fs(myid,val,buffer,nthreads) \
  do { \
    _Pragma("omp barrier") \
    idx_t _i; \
    falseShareSet(buffer,myid,val); \
    _Pragma("omp barrier") \
    /* works for small number of threads */ \
    for (_i=0;_i<(nthreads);++_i) { \
      if ((val) < falseShareGet(buffer,_i)) { \
        (val) = falseShareGet(buffer,_i); \
      } \
    } \
  } while(0) 

#define dl_omp_sumreduce(myid,val,buffer,nthreads) \
    do { \
      _Pragma("omp barrier") \
      idx_t _i; \
      (buffer)[myid] = (val); \
      _Pragma("omp barrier") \
      /* works for small number of threads */ \
      for ((val)=_i=0;_i<(nthreads);++_i) { \
        (val) += (buffer)[_i]; \
      } \
    } while (0)

#define dl_omp_prefixsum(myid,val,buffer,nthreads) \
do { \
  _Pragma("omp barrier") \
  buffer[myid+1] = val; \
  _Pragma("omp barrier") \
  _Pragma("omp single") \
  { \
    idx_t _i; \
    buffer[0] = 0; \
    for (_i=1;_i<nthreads+1;++_i) { \
      buffer[_i] += buffer[_i-1];  \
    } \
  } \
  val = buffer[myid]; \
} while(0) 

/* prefix sum */
#define dl_omp_pfreduce_off(myid,values,nvalues,partial,buffer,voff,noff, \
    nthreads) \
  do { \
    _Pragma ("omp barrier") \
    idx_t _i; \
    idx_t _mystart = (((nvalues)/(nthreads))*(myid)) + \
      gk_min((myid),(nvalues)%(nthreads)); \
    idx_t _mysize = ((nvalues)/(nthreads)) + \
      ((myid) < (nvalues)%(nthreads) ? 1 : 0); \
    for (_i=_mystart+1;_i<_mystart+_mysize;++_i) { \
      (values)[_i*(voff)] += (values)[(_i-1)*(voff)]; \
    } \
    (partial)[(myid)*(noff)] = (values)[(_i-1)*(voff)]; \
    _Pragma("omp barrier") \
    for (_i=1;_i<nthreads;_i<<=1) { \
      if ((myid)>=_i) { \
        (buffer)[(myid)*(noff)] = (partial)[(myid)*(noff)] + \
            (partial)[(myid-_i)*(noff)]; \
      } \
      _Pragma("omp barrier") \
      if (myid > 0) { \
        partial[((myid)*(noff))] = (buffer)[((myid)*(noff))]; \
      } \
      _Pragma("omp barrier") \
    } \
    for (_i=_mystart;_i<_mystart+_mysize;++_i) { \
      (values)[_i*(voff)] += (partial)[(myid)*(noff)] - \
          (values)[(_mystart+_mysize-1)*(voff)]; \
    } \
    _Pragma ("omp barrier") \
  } while(0)

#define dl_omp_pfreduce(myid,values,nvalues,partial,buffer,nthreads) \
  dl_omp_pfreduce_off(myid,values,nvalues,partial,buffer,1,1,nthreads)

#define dl_omp_pfreduce_fs(myid,values,nvalues,partial,buffer,nthreads) \
  dl_omp_pfreduce_off(myid,values,nvalues,partial,buffer, \
    falseShareSize(*(values)),falseShareSize(*(buffer)), \
    nthreads)

#define dl_omp_balanced_set(ptr,val,n) \
  do { \
    idx_t _i; \
    _Pragma("omp for schedule(guided)") \
    for (_i=0;_i<(n);++_i) { \
      (ptr)[_i] = (val); \
    } \
  } while (0) \

#define dl_omp_balanced_copy(dst,src,n) \
  do { \
    idx_t _i; \
    _Pragma("omp for schedule(guided)") \
    for (_i=0;_i<(n);++_i) { \
      (dst)[_i] = (src)[_i]; \
    } \
  } while (0) \



#define dl_omp_pfreduce_struct(myid,values,nvalues,partial,buffer,mem, \
    nthreads) \
  do { \
    _Pragma ("omp barrier") \
    idx_t _i; \
    idx_t _mystart = (((nvalues)/(nthreads))*(myid)) + \
      gk_min((myid),(nvalues)%(nthreads)); \
    idx_t _mysize = ((nvalues)/(nthreads)) + \
      ((myid) < (nvalues)%(nthreads) ? 1 : 0); \
    for (_i=_mystart+1;_i<_mystart+_mysize;++_i) { \
      (values)[_i].mem += (values)[_i-1].mem; \
    } \
    (partial)[myid] = (values)[_i-1].mem; \
    _Pragma("omp barrier") \
    for (_i=1;_i<nthreads;_i<<=1) { \
      if ((myid)>=_i) { \
        (buffer)[myid] = (partial)[myid] + \
            (partial)[(myid)-_i]; \
      } \
      _Pragma("omp barrier") \
      if ((myid) > 0) { \
        partial[myid] = (buffer)[myid]; \
      } \
      _Pragma("omp barrier") \
    } \
    for (_i=_mystart;_i<_mystart+_mysize;++_i) { \
      (values)[_i].mem += (partial)[myid] - \
          (values)[_mystart+_mysize-1].mem; \
    } \
    _Pragma ("omp barrier") \
  } while(0)

#define dl_omp_mergelists(myid,lst,blst,lstsize,nlst,partial,buffer, \
    nthreads) \
  do { \
    _Pragma("omp barrier") \
    idx_t _mysize = (lstsize)[(myid)+1]; \
    dl_omp_pfreduce((myid),((lstsize)+1),(nthreads),(partial),(buffer), \
        (nthreads)); \
    memcpy((blst)+(lstsize)[myid],(lst),_mysize*sizeof(*(lst))); \
    ASSERT(lstsize[nthreads] <= nlst); \
    _Pragma("omp barrier") \
  } while(0)

#define dl_omp_mergebnd(myid,lst,plst,blst,lstsize,nlst,partial,buffer, \
    nthreads) \
  do { \
    idx_t _j, _k; \
    dl_omp_mergelists(myid,lst,blst,lstsize,nlst,partial,buffer, \
      nthreads); \
    _Pragma("omp for schedule(guided,CACHE_LINE_INTS)") \
    for(_j=0;_j<(lstsize)[nthreads];++_j) { \
      _k = (blst)[_j]; \
      (plst)[_k] = _j; \
    } \
  } while(0)



/******************************************************************************
* JUNK ************************************************************************
******************************************************************************/
#define RIGHT_SIDE(c,to,from) \
  ((c == 0 && to > from) || (c == 1 && to < from))

#define MYSTART(myid,n,m) \
  ((((n)/(m))*(myid)) + gk_min((n)%(m),(myid)))

#define MYSIZE(myid,n,m) \
  (((n)/(m))+((myid)<(n)%(m)?1:0))

#define falseShareSize(t) \
  (CACHE_LINE_SIZE/sizeof(t))

#define falseShareArraySize(n,type) \
  ((n)*falseShareSize(type))

#define falseShareInc(ptr,i) \
  (++(ptr)[(i)*falseShareSize(*(ptr))])

#define falseShareDec(ptr,i) \
  (--(ptr)[(i)*falseShareSize(*(ptr))])

#define falseShareGet(ptr,i) \
  ((ptr)[(i)*falseShareSize(*(ptr))])

#define falseShareSet(ptr,i,v) \
  ((ptr)[(i)*falseShareSize(*(ptr))] = v)

#define falseShareOffset(ptr,o) \
  ((ptr)+((o)*falseShareSize(*(ptr))))

#define falseShareMemset(ptr,v,n) \
  memset((ptr),v,(n)*CACHE_LINE_SIZE);

#define invertGet(ptr,i) \
  ((ptr)[i] ^ ONES)

#define invertExists(ptr,i) \
  ((ptr)[i] < 0)

#define invertSet(ptr,i,v) \
  ((ptr)[i] = (v) ^ ONES)

#define ProperSide(c, from, other) \
  (((c) == 0 && ((from)-(other) < 0)) || \
    ((c) == 1 && ((from)-(other) > 0)))

#define uplog(i) \
  ((idx_t)ceil(log2(i)))

#define PerturbedMatch(ctrl,s,maxidx) \
  (maxidx == -1 || s == 0 || rand_r(&ctrl->seed) % s)

/* SHEM and RM */
#define MATCH(j,maxidx,cnvtx,cmap,match) \
  do { \
    /* order here actualyl matters for race conditions */ \
    match[j] = maxidx; \
    match[maxidx] = j; \
    cmap[j] = cmap[maxidx] = cnvtx++; \
  } while (0)

#define SELFMATCH(i,cnvtx,cmap,match) \
  do { \
    match[i] = i; \
    cmap[i] = cnvtx++; \
  } while (0)



#endif

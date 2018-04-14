/**
 * @file dlmacros.h
 * @brief Misc macros
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2013-10-08
 */



#ifndef DL_MACRO_H
#define DL_MACRO_H

/* check whether or not we can use C99 stuf */
#if (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L)
  #define __DL_RESTRICT restrict
#else
  #define __DL_RESTRICT
#endif

/******************************************************************************
* Useful Macros ***************************************************************
******************************************************************************/
#define dl_bitsize(a) (sizeof(a)*8u)

#define dl_maxshift(a,b) (dl_min(b,(dl_bitsize(a)-1)))

#define dl_min(a,b) ((a) > (b) ? (b) : (a))

#define dl_max(a,b) ((a) < (b) ? (b) : (a))

#define dl_swap(a,b) \
  do { \
    __typeof__(b) _swap_var = (b); \
    (b) = (a); \
    (a) = _swap_var; \
  } while(0)

#define dl_safe_diff(a,b) \
  (((a) > (b)) ? ((a) - (b)) : ((b) - (a)))

#define dl_storemax(max,v) (((max) < (v)) ? ((max) = (v)) : (max))

#define dl_storemin(min,v) (((min) > (v)) ? ((min) = (v)) : (min))

/* count leading zeros of arbitrary size number */
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
  #define dl_clz(a) \
    ((size_t)( (a) == 0 ? dl_bitsize(a) : \
      (sizeof(a) == sizeof(unsigned long long) ? \
        (size_t)__builtin_clzll((unsigned long long)a) \
        : (size_t)(__builtin_clz((unsigned int)a) - \
          (dl_bitsize(unsigned int) - dl_bitsize(a))) \
      ) \
    ))
#else
  /* cause trouble */
#endif

#define dl_near_equal(a,b) \
  ((a) == (b) ? 1 : ( \
      (fabs(a) > fabs(b) ? (fabs(((a)-(b))/(fabs(a)+1.0)) < 0.001) : \
        (fabs(((a)-(b))/(fabs(b)+1.0)) < 0.001))))

   
#endif

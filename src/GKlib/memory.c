/*!
\file  memory.c
\brief This file contains various allocation routines 

The allocation routines included are for 1D and 2D arrays of the 
most datatypes that GKlib support. Many of these routines are 
defined with the help of the macros in gk_memory.h. These macros 
can be used to define other memory allocation routines.

\date   Started 4/3/2007
\author George
\version\verbatim $Id: memory.c 7750 2009-12-22 15:02:25Z karypis $ \endverbatim
*/


#include <GKlib.h>

#ifdef USE_DLMALLOC
#ifdef GKMSPACE
/* This is the mspace for all the gk_malloc() calls. This is a thread local allocation */
static __thread mspace gk_mspace = 0;

/* This function is mostly for debugging */
void gk_printmspaceaddr() { printf("mspace: %p\n", (void *)gk_mspace); }
#endif
#endif


/*************************************************************************/
/*! Define the set of memory allocation routines for each data type */
/**************************************************************************/
GK_MKALLOC(gk_c,   char)
GK_MKALLOC(gk_i,   int)
GK_MKALLOC(gk_i32, int32_t)
GK_MKALLOC(gk_i64, int64_t)
GK_MKALLOC(gk_f,   float)
GK_MKALLOC(gk_d,   double)
GK_MKALLOC(gk_idx, gk_idx_t)

GK_MKALLOC(gk_ckv,   gk_ckv_t)
GK_MKALLOC(gk_ikv,   gk_ikv_t)
GK_MKALLOC(gk_i32kv, gk_i32kv_t)
GK_MKALLOC(gk_i64kv, gk_i64kv_t)
GK_MKALLOC(gk_fkv,   gk_fkv_t)
GK_MKALLOC(gk_dkv,   gk_dkv_t)
GK_MKALLOC(gk_skv,   gk_skv_t)
GK_MKALLOC(gk_idxkv, gk_idxkv_t)






/*************************************************************************
* This function allocates a two-dimensional matrix 
**************************************************************************/
void gk_AllocMatrix(void ***r_matrix, size_t elmlen, size_t ndim1, size_t ndim2)
{
  gk_idx_t i, j;
  void **matrix;

  *r_matrix = NULL;

  if ((matrix = (void **)gk_malloc(ndim1*sizeof(void *), "gk_AllocMatrix: matrix")) == NULL)
    return;

  for (i=0; i<ndim1; i++) {
    if ((matrix[i] = (void *)gk_malloc(ndim2*elmlen, "gk_AllocMatrix: matrix[i]")) == NULL) {
      for (j=0; j<i; j++) 
        gk_free((void **)&matrix[j], LTERM);
      return;
    }
  }

  *r_matrix = matrix;
}


/*************************************************************************
* This function frees a two-dimensional matrix 
**************************************************************************/
void gk_FreeMatrix(void ***r_matrix, size_t ndim1, size_t ndim2)
{
  gk_idx_t i;
  void **matrix;

  if ((matrix = *r_matrix) == NULL)
    return;

  for (i=0; i<ndim1; i++) 
    gk_free((void **)&matrix[i], LTERM);

  gk_free((void **)r_matrix, LTERM); 

}


/*************************************************************************/
/*! This function is my wrapper around malloc that provides the following
    enhancements over malloc:
    * It always allocates one byte of memory, even if 0 bytes are requested.
      This is to ensure that checks of returned values do not lead to NULL
      due to 0 bytes requested.
    * It zeros-out the memory that is allocated. This is for a quick init
      of the underlying datastructures.
*/
/**************************************************************************/
void *gk_malloc(size_t nbytes, char *msg)
{
  void *ptr=NULL;

  if (nbytes == 0)
    nbytes++;  /* This was added to make all the mallocs to actually allocate some memory */

#ifdef USE_DLMALLOC
#ifdef GKMSPACE
  if (gk_mspace == 0)
    gk_mspace = create_mspace(0, 0);

  if (gk_mspace == NULL) {
    gk_errexit(SIGMEM, "***Memory allocation failed for creating gk_mspace.");
    return NULL;
  }

  ptr = (void *)mspace_malloc(gk_mspace, nbytes);
#else
  ptr = (void *)dlmalloc(nbytes);
#endif
#else
  ptr = (void *)malloc(nbytes);
#endif

  if (ptr == NULL) {
    fprintf(stderr, "   Current memory used:  %10"PRId64" bytes\n", (int64_t)gk_GetCurMemoryUsed());
    fprintf(stderr, "   Maximum memory used:  %10"PRId64" bytes\n", (int64_t)gk_GetMaxMemoryUsed());

    gk_errexit(SIGMEM, "***Memory allocation failed for %s. Requested size: %"PRId64" bytes", 
        msg, (int64_t)nbytes);
    return NULL;
  }

  /* zero-out the allocated space */
  memset(ptr, 0, nbytes);

  return ptr;
}


/*************************************************************************
* This function is my wrapper around realloc
**************************************************************************/
void *gk_realloc(void *oldptr, size_t nbytes, char *msg)
{
  void *ptr=NULL;

  nbytes++;  /* This was added to make all the mallocs to actually allocate some memory */

  if (nbytes == 0) {
    gk_free((void **)&oldptr, LTERM);
    return NULL;
  }

#ifdef USE_DLMALLOC
#ifdef GKMSPACE
  if (gk_mspace == 0)
    gk_mspace = create_mspace(0, 0);

  if (gk_mspace == NULL) {
    gk_errexit(SIGMEM, "***Memory allocation failed for creating gk_mspace.");
    return NULL;
  }

  ptr = (void *)mspace_realloc(gk_mspace, oldptr, nbytes);
#else
  ptr = (void *)dlrealloc(oldptr, nbytes);
#endif
#else
  ptr = (void *)realloc(oldptr, nbytes);
#endif

  if (ptr == NULL) {
    fprintf(stderr, "   Maximum memory used:              %10ju bytes\n", gk_GetMaxMemoryUsed());
    fprintf(stderr, "   Current memory used:              %10ju bytes\n", gk_GetCurMemoryUsed());

    gk_errexit(SIGMEM, "***Memory re-allocation failed for %s. Requested size: %zd bytes", msg, nbytes);
    return NULL;
  }

  return ptr;
}


/*************************************************************************
* This function is my wrapper around free, allows multiple pointers    
**************************************************************************/
void gk_free(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
#ifdef USE_DLMALLOC
#ifdef GKMSPACE
    mspace_free(gk_mspace, *ptr1);
#else
    dlfree(*ptr1);
#endif
#else
    free(*ptr1);
#endif
  *ptr1 = NULL;

  va_start(plist, ptr1);

  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
#ifdef USE_DLMALLOC
#ifdef GKMSPACE
      mspace_free(gk_mspace, *ptr);
#else
      dlfree(*ptr);
#endif
#else
      free(*ptr);
#endif
    *ptr = NULL;
  }

  va_end(plist);
}            


/*************************************************************************
* This function cleans up the memory that has been allocated thus far.
* This work only if code has been compiled with GKMSPACE 
**************************************************************************/
void gk_malloc_cleanup()
{
#ifdef USE_DLMALLOC
#ifdef GKMSPACE
  if (gk_mspace != 0)
    destroy_mspace(gk_mspace);
  gk_mspace = 0;
#endif
#endif
}


/*************************************************************************
* This function returns the current ammount of dynamically allocated
* memory that is used by the system
**************************************************************************/
uintmax_t gk_GetCurMemoryUsed()
{
  size_t cused=0;
#ifdef USE_DLMALLOC
  struct mallinfo meminfo;

#ifdef GKMSPACE
  if (gk_mspace != 0) {
    meminfo = mspace_mallinfo(gk_mspace);
    cused = meminfo.uordblks;
  }
#else
  meminfo = dlmallinfo();
  cused = meminfo.uordblks;
#endif
#endif

  return (uintmax_t) cused;
}


/*************************************************************************
* This function returns the maximum ammount of dynamically allocated 
* memory that was used by the system
**************************************************************************/
uintmax_t gk_GetMaxMemoryUsed()
{
  size_t mused=0;
#ifdef USE_DLMALLOC
  struct mallinfo meminfo;

#ifdef GKMSPACE
  if (gk_mspace != 0) {
    meminfo = mspace_mallinfo(gk_mspace);
    mused = meminfo.usmblks;
  }
#else
  meminfo = dlmallinfo();
  mused = meminfo.usmblks;
#endif
#endif

  return (uintmax_t) mused;
}



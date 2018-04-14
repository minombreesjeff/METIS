/*!
\file 
\brief Functions dealing with creating and allocating mcores

\date Started 5/30/11
\author George
\author Copyright 1997-2011, Regents of the University of Minnesota 
\version $Id: mcore.c 10020 2011-05-30 17:27:20Z karypis $
*/

#include <GKlib.h>


/*************************************************************************/
/*! This function creates an mcore */
/*************************************************************************/
gk_mcore_t *gk_mcoreCreate(size_t coresize)
{
  gk_mcore_t *mcore;

  mcore = (gk_mcore_t *) gk_malloc(sizeof(gk_mcore_t), "gk_mcoreCreate: mcore");
  memset(mcore, 0, sizeof(gk_mcore_t));

  mcore->coresize = coresize;
  mcore->corecpos = 0;

  mcore->core = gk_malloc(mcore->coresize, "gk_mcoreCreate: core");

  /* allocate the memory for keeping track of malloc ops */
  mcore->nmops = 2048;
  mcore->cmop  = 0;
  mcore->mops  = (gk_mop_t *)gk_malloc(mcore->nmops*sizeof(gk_mop_t), 
                             "gk_mcoreCreate: mops");

  return mcore;
}


/*************************************************************************/
/*! This function destroyes an mcore */
/*************************************************************************/
void gk_mcoreDestroy(gk_mcore_t **r_mcore, int showstats)
{
  gk_mcore_t *mcore = *r_mcore;

  if (showstats)
    printf("\n gk_mcore statistics\n" 
           "           coresize: %12zu         nmops: %12zu  cmop: %6zu\n"
           "        num_callocs: %12zu   num_hallocs: %12zu\n"
           "       size_callocs: %12zu  size_hallocs: %12zu\n"
           "        cur_callocs: %12zu   cur_hallocs: %12zu\n"
           "        max_callocs: %12zu   max_hallocs: %12zu\n",
           mcore->coresize, mcore->nmops, mcore->cmop,
           mcore->num_callocs,  mcore->num_hallocs,
           mcore->size_callocs, mcore->size_hallocs,
           mcore->cur_callocs,  mcore->cur_hallocs,
           mcore->max_callocs,  mcore->max_hallocs);

  if (mcore->cur_callocs != 0 || mcore->cur_hallocs != 0 || mcore->cmop != 0) {
    printf("***Warning: mcore memory was not fully freed when destroyed.\n"
           " cur_callocs: %6zu  cur_hallocs: %6zu cmop: %6zu\n",
           mcore->cur_callocs,  mcore->cur_hallocs, mcore->cmop);
  }

  gk_free((void **)&mcore->core, &mcore->mops, &mcore, LTERM);

  *r_mcore = NULL;
}


/*************************************************************************/
/*! This function allocate space from the core/heap */
/*************************************************************************/
void *gk_mcoreMalloc(gk_mcore_t *mcore, size_t nbytes)
{
  void *ptr;

  /* pad to make pointers 8-byte aligned */
  nbytes += (nbytes%8 == 0 ? 0 : 8 - nbytes%8);

  if (mcore->cmop == mcore->nmops) {
    mcore->nmops *= 2;
    mcore->mops = 
        gk_realloc(mcore->mops, mcore->nmops*sizeof(gk_mop_t), "gk_mcoremalloc: mops");
  }

  if (mcore->corecpos + nbytes < mcore->coresize) {
    /* service this request from the core */
    ptr = ((char *)mcore->core)+mcore->corecpos;
    mcore->corecpos += nbytes;
    mcore->mops[mcore->cmop].flag   = 1;
    mcore->mops[mcore->cmop].nbytes = nbytes;
    mcore->mops[mcore->cmop].ptr    = ptr;
    mcore->cmop++;

    mcore->num_callocs++;
    mcore->size_callocs += nbytes;
    mcore->cur_callocs  += nbytes;
    if (mcore->max_callocs < mcore->cur_callocs)
      mcore->max_callocs = mcore->cur_callocs;
  }
  else {
    /* service this request from the heap */
    ptr = gk_malloc(nbytes, "gk_mcoremalloc: ptr");
    mcore->mops[mcore->cmop].flag   = 2;
    mcore->mops[mcore->cmop].nbytes = nbytes;
    mcore->mops[mcore->cmop].ptr    = ptr;
    mcore->cmop++;

    mcore->num_hallocs++;
    mcore->size_hallocs += nbytes;
    mcore->cur_hallocs  += nbytes;
    if (mcore->max_hallocs < mcore->cur_hallocs)
      mcore->max_hallocs = mcore->cur_hallocs;
  }

  /*
  printf("MCMALLOC: %zu %d %8zu\n", mcore->cmop-1, 
      mcore->mops[mcore->cmop-1].flag, mcore->mops[mcore->cmop-1].nbytes);
  */

  return ptr;
}


/*************************************************************************/
/*! This function sets a marker in the stack of malloc ops to be used
    subsequently for freeing purposes */
/*************************************************************************/
void gk_mcorePush(gk_mcore_t *mcore)
{
  if (mcore->cmop == mcore->nmops) {
    mcore->nmops *= 2;
    mcore->mops = gk_realloc(mcore->mops, mcore->nmops*sizeof(gk_mop_t), 
                      "gk_mcorepush: mops");
  }

  mcore->mops[mcore->cmop].flag   = 0;
  mcore->mops[mcore->cmop].nbytes = 0;
  mcore->mops[mcore->cmop].ptr    = NULL;
  mcore->cmop++;

  /* printf("MCPPUSH:   %zu\n", mcore->cmop-1); */
}


/*************************************************************************/
/*! This function frees all mops since the last push */
/*************************************************************************/
void gk_mcorePop(gk_mcore_t *mcore)
{
  while (mcore->cmop > 0) {
    mcore->cmop--;
    switch (mcore->mops[mcore->cmop].flag) {
      case 0: /* push marker */
        goto DONE;
        break; 
      case 1: /* core free */
        mcore->corecpos    -= mcore->mops[mcore->cmop].nbytes;
        mcore->cur_callocs -= mcore->mops[mcore->cmop].nbytes;
        if (mcore->corecpos < 0)
          errexit("Internal Error: wspace's core over-freed [%zu, %zu]\n",
              mcore->coresize, mcore->corecpos);
        break;
      case 2: /* heap free */
        gk_free((void **)&mcore->mops[mcore->cmop].ptr, LTERM);
        mcore->cur_hallocs -= mcore->mops[mcore->cmop].nbytes;
        break;
      default:
        errexit("Unknown mop flag of %d\n", mcore->mops[mcore->cmop].flag);
    }
  }

DONE:
  ;
  /*printf("MCPPOP:    %zu\n", mcore->cmop); */
}


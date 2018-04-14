/*!
\file  
\brief Various routines for providing portable 32 and 64 bit random number
       generators.

\date   Started 5/17/2007
\author George
\version\verbatim $Id: random.c 1902 2007-05-28 14:21:18Z karypis $ \endverbatim
*/


#include <GKlib.h>


/*************************************************************************/
/*! Create the various random number functions */
/*************************************************************************/
GK_MKRANDOM(gk_c,   size_t, char)
GK_MKRANDOM(gk_i,   size_t, int)
GK_MKRANDOM(gk_f,   size_t, float)
GK_MKRANDOM(gk_d,   size_t, double)
GK_MKRANDOM(gk_idx, size_t, gk_idx_t)

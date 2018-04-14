/**
\file
\brief This file contains various helper API routines for using METIS.

\date   Started 5/12/2011
\author George  
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version\verbatim $Id: auxapi.c 10103 2011-06-07 18:28:27Z karypis $ \endverbatim
*/


#include "metislib.h"


/*************************************************************************/
/*! This function free memory that was allocated by METIS and retuned
    to the application.
    
    \param ptr points to the memory that was previously allocated by
           METIS.
*/
/*************************************************************************/
int METIS_Free(void *ptr)
{
  gk_free((void **)&ptr, LTERM);

  return METIS_OK;
}


/*************************************************************************/
/*! This function sets the default values for the options.
    
    \param options points to an array of size at least METIS_NOPTIONS.
*/
/*************************************************************************/
int METIS_SetDefaultOptions(idx_t *options)
{
  iset(METIS_NOPTIONS, -1, options);

  return METIS_OK;
}



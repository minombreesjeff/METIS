/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * par_metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: parmetis.h,v 1.1 1998/05/25 17:55:05 karypis Exp $
 */

#define PARMETIS_MAJOR_VERSION        3
#define PARMETIS_MINOR_VERSION        0 

/*
#define DEBUG			1
#define DMALLOC			1
*/

#include <stdheaders.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include <rename.h>
#include <defs.h>
#include <struct.h>
#include <macros.h>
#include <proto.h>


/**
 * @file defs.h
 * @brief Definitions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */


#ifndef DEFS_H
#define DEFS_H

#define MAXDEGREEFOR2HOP 20

/* sys stuff */
#define CACHE_LINE_SIZE 64
#define CACHE_LINE_INTS (CACHE_LINE_SIZE / sizeof(idx_t)) 
#define CACHE_LINE_REALS (CACHE_LINE_SIZE / sizeof(real_t))
#define CACHE_SIZE 12000000
#define CACHE_DATA_SIZE (CACHE_SIZE*0.75f) /* assume we only get to use 3/4 */
#define PROC_SIZE 4
#define CACHE_PER_THREAD (CACHE_DATA_SIZE / PROC_SIZE)
#define PAGESIZE 4096

/* performance parameters */
#define BLOCKSHIFT 12 
#define BLOCKSIZE (1<<BLOCKSHIFT) 
#define BLOCKMASK (BLOCKSIZE-1)
#define NBRPOOL_EXP_RATE 1.5

/* parallel stuff */
#define PAR_COARSEN
#define PAR_PARTITION
#define PAR_REFINE

/* coarsening */
#define ROOT -2
#define LWGT 2
/* #define NOMASK */
/* #define NOSHEM */

/* initial partitioning */
#define NSOLUTIONS 16
#define PAR_COARSEN_FRACTION 0.75
#define PAR_COARSEN_FACTOR 128
#define TPPRATIO 4
/* #define LGCOPY */

/* refinement stuff */
#define MOVEBUFFERSIZE 512 
#define MINUPDATES 1024

/* other stuff */
#define LISTUNINDEXED -2
#define LISTEMPTY -1
#define UBFACTOR 1.030

#endif

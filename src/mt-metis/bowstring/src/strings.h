/**
 * @file strings.h
 * @brief String literals
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Dominique LaSalle
 * @version 1
 * @date 2014-09-24
 */




#ifndef BOWSTRING_STRINGS_H
#define BOWSTRING_STRINGS_H




/******************************************************************************
* TYPES FORMATTING ************************************************************
******************************************************************************/


#define PF_FLOAT_T "%0.05f"
#define RF_FLOAT_T "%f"
#define PF_DOUBLE_T "%0.05lf"
#define RF_DOUBLE_T "%lf"


#ifdef BOWSTRING_64BIT_VERTICES
#define PF_VTX_T "%"PRIu64
#define RF_VTX_T "%"PRIu64
#else
#define PF_VTX_T "%"PRIu32
#define RF_VTX_T "%"PRIu32
#endif
#ifdef BOWSTRING_64BIT_EDGES
#define PF_ADJ_T "%"PRIu64
#define RF_ADJ_T "%"PRIu64
#else
#define PF_ADJ_T "%"PRIu32
#define RF_ADJ_T "%"PRIu32
#endif
#ifdef BOWSTRING_64BIT_THREADS
#define PF_TID_T "%"PRIu64
#define RF_TID_T "%"PRIu64
#else
#define PF_TID_T "%"PRIu32
#define RF_TID_T "%"PRIu32
#endif
#ifdef BOWSTRING_INT_WEIGHTS
#ifdef BOWSTRING_64BIT_WEIGHTS
#define PF_WGT_T "%"PRId64
#define RF_WGT_T "%"PRId64
#else
#define PF_WGT_T "%"PRId32
#define RF_WGT_T "%"PRId32
#endif
#else
#ifdef BOWSTRING_64BIT_WEIGHTS
#define PF_WGT_T PF_DOUBLE_T
#define RF_WGT_T RF_DOUBLE_T
#else
#define PF_WGT_T PF_FLOAT_T
#define RF_WGT_T RF_FLOAT_T
#endif
#endif


#ifdef BOWSTRING_64BIT_VLABELS
#define PF_VLBL_T "%"PRIu64
#define RF_VLBL_T "%"PRIu64
#else
#define PF_VLBL_T "%"PRIu32
#define RF_VLBL_T "%"PRIu32
#endif /* BOWSTRING_64BIT_VLABELS */
#ifdef BOWSTRING_64BIT_ELABELS
#define PF_ELBL_T "%"PRIu64
#define RF_ELBL_T "%"PRIu64
#else
#define PF_ELBL_T "%"PRIu32
#define RF_ELBL_T "%"PRIu32
#endif /* BOWSTRING_64BIT_ELABELS */


#ifdef BOWSTRING_DOUBLE_COORDS
#define PF_COORD_T PF_DOUBLE_T
#define RF_COORD_T RF_DOUBLE_T
#else /* BOWSTRING_DOUBLE_COORDS */
#define PF_COORD_T PF_FLOAT_T
#define RF_COORD_T RF_FLOAT_T
#endif /* BOWSTRING_DOUBLE_COORDS */





#endif

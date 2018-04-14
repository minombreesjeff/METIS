/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * Started 8/27/94
 * George
 *
 * $Id: defs.h 10188 2011-06-13 14:30:26Z karypis $
 *
 */

#ifndef _LIBMETIS_DEFS_H_
#define _LIBMETIS_DEFS_H_

#define METISTITLE              "METIS 5.0 Copyright 1998-11, Regents of the University of Minnesota\n"
#define MAXLINE			1280000

#define LTERM			(void **) 0	/* List terminator for gk_free() */

#define MAXNCON			16		/* The maximum number of constrains */
#define MAXNOBJ			16		/* The maximum number of objectives */

#define HTLENGTH		((1<<11)-1)

#define INIT_MAXNAD             200     /* Initial number of maximum number of 
                                           adjacent domains. This number will be
                                           adjusted as required. */

/* Types of boundaries */
#define BNDTYPE_REFINE          1       /* Used for k-way refinement-purposes */
#define BNDTYPE_BALANCE         2       /* Used for k-way balancing purposes */

/* Mode of optimization */
#define OMODE_REFINE            1       /* Optimize the objective function */
#define OMODE_BALANCE           2       /* Balance the subdomains */

/* Types of vertex statues in the priority queue */
#define VPQSTATUS_PRESENT      1       /* The vertex is in the queue */
#define VPQSTATUS_EXTRACTED    2       /* The vertex has been extracted from the queue */
#define VPQSTATUS_NOTPRESENT   3       /* The vertex is not present in the queue and
                                          has not been extracted before */

#define UNMATCHED		-1

#define HTABLE_EMPTY    	-1

#define NGR_PASSES		4	/* Number of greedy refinement passes */
#define NLGR_PASSES		5	/* Number of GR refinement during IPartition */

#define LARGENIPARTS		8	/* Number of random initial partitions */
#define SMALLNIPARTS		3	/* Number of random initial partitions */

#define COARSEN_FRACTION	0.75	/* Node reduction between succesive coarsening levels */
#define COARSEN_FRACTION2	0.99	/* Node reduction between succesive coarsening levels */


#define UNBALANCE_FRACTION		1.03

#define COMPRESSION_FRACTION		0.85

#define ORDER_UNBALANCE_FRACTION	1.10

#define MMDSWITCH		        120

/* Default ufactors for the various operational modes */
#define PMETIS_DEFAULT_UFACTOR          1
#define MCPMETIS_DEFAULT_UFACTOR        10
#define KMETIS_DEFAULT_UFACTOR          30
#define OMETIS_DEFAULT_UFACTOR          30

#endif

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains header files
 *
 * Started 10/2/97
 * George
 *
 * $Id: rename.h,v 1.4 1997/11/05 19:44:24 karypis Exp $
 *
 */


/* coarsen.c */
#define Coarsen2Way			__Coarsen2Way

/* fm.c */
#define FM_2WayEdgeRefine		__FM_2WayEdgeRefine

/* balance.c */
#define Balance2Way			__Balance2Way
#define Bnd2WayBalance			__Bnd2WayBalance
#define General2WayBalance		__General2WayBalance

/* initpart.c */
#define Init2WayPartition		__Init2WayPartition
#define InitSeparator			__InitSeparator
#define GrowBisection			__GrowBisection
#define GrowBisectionNode		__GrowBisectionNode
#define RandomBisection			__RandomBisection

/* match.c */
#define Match_RM			__Match_RM
#define Match_HEM			__Match_HEM
#define Match_SHEM			__Match_SHEM

/* ccgraph.c */
#define CreateCoarseGraph		__CreateCoarseGraph
#define CreateCoarseGraphNoMask		__CreateCoarseGraphNoMask

/* memory.c */
#define AllocateWorkSpace		__AllocateWorkSpace
#define FreeWorkSpace			__FreeWorkSpace
#define idxwspacemalloc			__idxwspacemalloc
#define idxwspacefree			__idxwspacefree
#define CreateGraph			__CreateGraph
#define InitGraph			__InitGraph
#define FreeGraph			__FreeGraph

/* pmetis.c */
#define MlevelRecursiveBisection	__MlevelRecursiveBisection
#define MlevelEdgeBisection		__MlevelEdgeBisection
#define SplitGraphPart			__SplitGraphPart

/* pqueue.c */
#define PQueueInit			__PQueueInit
#define PQueueReset			__PQueueReset
#define PQueueFree			__PQueueFree
#define PQueueInsert			__PQueueInsert
#define PQueueDelete			__PQueueDelete
#define PQueueUpdate			__PQueueUpdate
#define PQueueUpdateUp			__PQueueUpdateUp
#define PQueueGetMax			__PQueueGetMax
#define PQueueSeeMax			__PQueueSeeMax
#define CheckHeap			__CheckHeap


/* refine.c */
#define Refine2Way			__Refine2Way
#define Allocate2WayPartitionMemory	__Allocate2WayPartitionMemory
#define Compute2WayPartitionParams	__Compute2WayPartitionParams
#define Project2WayPartition		__Project2WayPartition



/* util.c */
#define errexit				__errexit
#define GKfree				__GKfree
#ifndef DMALLOC
#define imalloc				__imalloc
#define idxmalloc			__idxmalloc
#define fmalloc				__fmalloc
#define ismalloc			__ismalloc
#define idxsmalloc			__idxsmalloc
#define GKmalloc			__GKmalloc
#endif
#define iset				__iset
#define idxset				__idxset
#define idxamax				__idxamax
#define idxamin				__idxamin
#define idxasum				__idxasum
#define idxsum				__idxsum
#define charsum				__charsum
#define isum				__isum
#define ssum				__ssum
#define sscale				__sscale
#define snorm2				__snorm2
#define sdot				__sdot
#define saxpy				__saxpy
#define RandomPermute			__RandomPermute
#define ispow2				__ispow2
#define InitRandom			__InitRandom


/* timing.c */
#define InitTimers			__InitTimers
#define PrintTimers			__PrintTimers
#define seconds				__seconds

/* debug.c */
#define ComputeCut			__ComputeCut
#define CheckBnd			__CheckBnd
#define CheckBnd2			__CheckBnd2
#define CheckNodeBnd			__CheckNodeBnd
#define CheckRInfo			__CheckRInfo
#define CheckNodePartitionParams	__CheckNodePartitionParams
#define IsSeparable			__IsSeparable

/* bucketsort.c */
#define BucketSortKeysInc		__BucketSortKeysInc


/* graph.c */
#define SetUpGraph			__SetUpGraph
#define RandomizeGraph			__RandomizeGraph
#define IsConnected			__IsConnected
#define IsConnectedSubdomain		__IsConnectedSubdomain
#define IsConnected2			__IsConnected2
#define FindComponents			__FindComponents
#define AnalyzeGraph			__AnalyzeGraph


/* stat.c */
#define ComputePartitionInfo		__ComputePartitionInfo
#define ComputePartitionBalance		__ComputePartitionBalance
#define ComputeElementBalance		__ComputeElementBalance

/* kmetis.c */
#define MlevelKWayPartitioning		__MlevelKWayPartitioning

/* kwayrefine.c */
#define RefineKWay			__RefineKWay
#define AllocateKWayPartitionMemory	__AllocateKWayPartitionMemory
#define ComputeKWayPartitionParams	__ComputeKWayPartitionParams
#define ProjectKWayPartition		__ProjectKWayPartition
#define IsBalanced			__IsBalanced
#define ComputeKWayBoundary		__ComputeKWayBoundary
#define ComputeKWayBalanceBoundary	__ComputeKWayBalanceBoundary

/* kwayfm.c */
#define Random_KWayEdgeRefine		__Random_KWayEdgeRefine
#define Random_KWayEdgeBalance		__Random_KWayEdgeBalance
#define Greedy_KWayEdgeRefine		__Greedy_KWayEdgeRefine
#define Greedy_KWayEdgeBalance		__Greedy_KWayEdgeBalance

/* ometis.c */
#define MlevelNestedDissection		__MlevelNestedDissection
#define MlevelNestedDissectionCC	__MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple	__MlevelNodeBisectionMultiple
#define MlevelNodeBisection		__MlevelNodeBisection
#define SplitGraphOrder			__SplitGraphOrder
#define MMDOrder			__MMDOrder
#define SplitGraphOrderCC		__SplitGraphOrderCC

/* srefine.c */
#define Refine2WayNode			__Refine2WayNode
#define Allocate2WayNodePartitionMemory	__Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams	__Compute2WayNodePartitionParams
#define Project2WayNodePartition	__Project2WayNodePartition

/* sfm.c */
#define FM_2WayNodeRefine		__FM_2WayNodeRefine
#define FM_2WayNodeRefineEqWgt		__FM_2WayNodeRefineEqWgt
#define FM_2WayNodeRefine_OneSided	__FM_2WayNodeRefine_OneSided
#define FM_2WayNodeBalance		__FM_2WayNodeBalance
#define ComputeMaxNodeGain		__ComputeMaxNodeGain


/* separator.c */
#define ConstructSeparator		__ConstructSeparator
#define ConstructMinCoverSeparator0	__ConstructMinCoverSeparator0
#define ConstructMinCoverSeparator	__ConstructMinCoverSeparator


/* mincover.o */
#define MinCover			__MinCover
#define MinCover_Augment		__MinCover_Augment
#define MinCover_Decompose		__MinCover_Decompose
#define MinCover_ColDFS			__MinCover_ColDFS
#define MinCover_RowDFS			__MinCover_RowDFS

/* mmd.c */
#define genmmd				__genmmd
#define mmdelm				__mmdelm
#define mmdint				__mmdint
#define mmdnum				__mmdnum
#define mmdupd				__mmdupd


/* compress.c */
#define CompressGraph			__CompressGraph
#define PruneGraph			__PruneGraph

/* svrefine.c */
#define VSepRefine			__VSepRefine

/* mesh.c */
#define TRIDUALMETIS			__TRIDUALMETIS
#define TETDUALMETIS			__TETDUALMETIS
#define HEXDUALMETIS			__HEXDUALMETIS
#define TRINODALMETIS			__TRINODALMETIS
#define TETNODALMETIS			__TETNODALMETIS
#define HEXNODALMETIS			__HEXNODALMETIS


/* fortran.c */
#define Change2CNumbering		__Change2CNumbering
#define Change2FNumbering		__Change2FNumbering
#define Change2FNumbering2		__Change2FNumbering2
#define Change2FNumberingOrder		__Change2FNumberingOrder
#define ChangeMesh2CNumbering		__ChangeMesh2CNumbering
#define ChangeMesh2FNumbering		__ChangeMesh2FNumbering
#define ChangeMesh2FNumbering2		__ChangeMesh2FNumbering2


/* myqsort.c */
#define iidxsort			__iidxsort
#define iintsort			__iintsort
#define ikeysort			__ikeysort
#define ikeyvalsort			__ikeyvalsort

/* ometispar.c */
#define MlevelNestedDissectionP		__MlevelNestedDissectionP

/* estmem.c */
#define EstimateCFraction		__EstimateCFraction
#define ComputeCoarseGraphSize		__ComputeCoarseGraphSize

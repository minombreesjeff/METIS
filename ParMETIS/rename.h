/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains renaming #defines to remove any conflicts that the
 * library may have with the users functions.
 *
 * Started 7/17/97
 * George
 *
 * $Id: rename.h,v 1.2 1997/07/18 00:32:12 karypis Exp $
 */

/* debug.c */
#define PrintVector PrintVector__
#define PrintGraph PrintGraph__
#define PrintGraph2 PrintGraph2__
#define PrintSetUpInfo PrintSetUpInfo__

/* comm.c */
#define CommInterfaceData CommInterfaceData__
#define CommChangedInterfaceData CommChangedInterfaceData__
#define GlobalSEMaxDouble GlobalSEMaxDouble__
#define GlobalSEMin GlobalSEMin__
#define GlobalSESum GlobalSESum__


/* util.c */
#define errexit errexit__
#define myprintf myprintf__
#define rprintf rprintf__
#ifndef DMALLOC
#define imalloc imalloc__
#define idxmalloc idxmalloc__
#define fmalloc fmalloc__
#define ismalloc ismalloc__
#define idxsmalloc idxsmalloc__
#define GKmalloc GKmalloc__
#endif
#define GKfree GKfree__
#define iset iset__
#define idxset idxset__
#define idxamax idxamax__
#define idxamin idxamin__
#define idxasum idxasum__
#define snorm2 snorm2__
#define sdot sdot__
#define saxpy saxpy__
#define ikeyvalsort_org ikeyvalsort_org__
#define IncKeyValueCmp IncKeyValueCmp__
#define dkeyvalsort dkeyvalsort__
#define DecKeyValueCmp DecKeyValueCmp__
#define BSearch BSearch__
#define RandomPermute RandomPermute__
#define FastRandomPermute FastRandomPermute__
#define ispow2 ispow2__


/* qsort_special.c */
#define iidxsort iidxsort__
#define iintsort iintsort__
#define ikeysort ikeysort__
#define ikeyvalsort ikeyvalsort__

/* memory.c */
#define PreAllocateMemory PreAllocateMemory__
#define FreeWSpace FreeWSpace__
#define CreateGraph CreateGraph__
#define InitGraph InitGraph__
#define FreeGraph FreeGraph__
#define FreeInitialGraph FreeInitialGraph__

/* kmetis.c */
#define SetUpGraph SetUpGraph__
#define SetUpCtrl SetUpCtrl__
#define ChangeNumbering ChangeNumbering__
#define InitTimers InitTimers__
#define PrintTimingInfo PrintTimingInfo__
#define PrintTimer PrintTimer__
#define GraphRandomPermute GraphRandomPermute__
#define ComputeMoveStatistics ComputeMoveStatistics__

/* setup.c */
#define SetUp SetUp__
#define Home_PE Home_PE__

/* coarsen.c */
#define GlobalMatch_HEM GlobalMatch_HEM__
#define Global_CreateCoarseGraph Global_CreateCoarseGraph__
#define LocalMatch_HEM LocalMatch_HEM__
#define Local_CreateCoarseGraph Local_CreateCoarseGraph__


/* edge_refine.c */
#define ProjectPartition ProjectPartition__
#define ComputePartitionParams ComputePartitionParams__
#define ComputeCut ComputeCut__
#define KWayRefine KWayRefine__
#define KWayRefineClean KWayRefineClean__
#define KWayAdaptiveRefineClean KWayAdaptiveRefineClean__

/* node_refine.c */
#define ComputeNodePartitionParams0 ComputeNodePartitionParams0__
#define ComputeNodePartitionParams ComputeNodePartitionParams__
#define KWayNodeRefine0 KWayNodeRefine0__
#define KWayNodeRefine KWayNodeRefine__
#define KWayNodeRefine2 KWayNodeRefine2__
#define PrintNodeBalanceInfo PrintNodeBalanceInfo__

/* initpart.c */
#define InitPartition InitPartition__
#define InitPartition_RB InitPartition_RB__
#define AssembleGraph AssembleGraph__
#define KeepPart KeepPart__

/* initmsection.c */
#define InitMultisection InitMultisection__
#define AssembleMultisectedGraph AssembleMultisectedGraph__
#define ConstructSeparator ConstructSeparator__
#define ConstructSeparator1 ConstructSeparator1__


/* serkmetis.c */
#define Ser_KMetis Ser_KMetis__
#define Ser_OMetis Ser_OMetis__ 
#define Ser_Coarsen Ser_Coarsen__
#define Ser_Match_HEM Ser_Match_HEM__
#define Ser_CreateCoarseGraph Ser_CreateCoarseGraph__
#define Ser_InitPartition Ser_InitPartition__
#define Ser_RecursiveBisection Ser_RecursiveBisection__
#define Ser_EdgeFM Ser_EdgeFM__
#define Ser_MlevelBisectGraph Ser_MlevelBisectGraph__
#define Ser_BisectGraph Ser_BisectGraph__
#define Ser_SplitGraphPart Ser_SplitGraphPart__
#define Ser_Refine Ser_Refine__
#define Ser_ProjectPartition Ser_ProjectPartition__
#define Ser_ComputePartitionParams Ser_ComputePartitionParams__
#define Ser_KWayRefine Ser_KWayRefine__


/* drivers.c */
#define Global_Partition Global_Partition__
#define FldGlobal_Partition FldGlobal_Partition__
#define Refine_Partition Refine_Partition__
#define AdaptiveUndirected_Partition AdaptiveUndirected_Partition__
#define AdaptiveDirected_Partition AdaptiveDirected_Partition__
#define SmallerSubGraph SmallerSubGraph__

/* fold.c */
#define EnoughMemory EnoughMemory__
#define FoldGraph FoldGraph__
#define UnFoldGraph UnFoldGraph__


/* move.c */
#define MoveGraph MoveGraph__
#define ProjectInfoBack ProjectInfoBack__
#define FindVtxPerm FindVtxPerm__

/* initdiff.c */
#define InitDiffusion InitDiffusion__
#define AssembleAdaptiveGraph AssembleAdaptiveGraph__
#define KWay_InitialDiffuser KWay_InitialDiffuser__
#define KWay_InitialDiffuser2 KWay_InitialDiffuser2__
#define ConjGrad ConjGrad__
#define mvMult mvMult__
#define setupLaplace setupLaplace__
#define Ser_KWayUpdateDegrees Ser_KWayUpdateDegrees__


/* order.c */
#define MultilevelOrder MultilevelOrder__
#define LabelSeparators LabelSeparators__
#define CompactGraph CompactGraph__
#define LocalOrder LocalOrder__


/* pqueue.c */
#define PQueueInit PQueueInit__
#define PQueueReset PQueueReset__
#define PQueueFree PQueueFree__
#define PQueueInsert PQueueInsert__
#define PQueueUpdate PQueueUpdate__
#define PQueueDelete PQueueDelete__
#define PQueueGetMax PQueueGetMax__
#define PQueueSeeMax PQueueSeeMax__
#define PQueueCheck PQueueCheck__

/* serometis.c */
#define Ser_NodeOMetis Ser_NodeOMetis__
#define Ser_NodeBisection Ser_NodeBisection__
#define Ser_NodeRefine Ser_NodeRefine__
#define Ser_NodeComputePartitionParams Ser_NodeComputePartitionParams__
#define Ser_NodeFM Ser_NodeFM__
#define Ser_ConstructSeparator Ser_ConstructSeparator__
#define CheckPartitionParams CheckPartitionParams__

/* mmd.c */
#define genmmd genmmd__
#define mmdelm mmdelm__
#define mmdint mmdint__
#define mmdnum mmdnum__
#define mmdupd mmdupd__


/* xyzpart.c */
#define Coordinate_Partition Coordinate_Partition__
#define PartSort PartSort__

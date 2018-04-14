/*
 * Copyright 1997-2015, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains header files
 *
 * Started 10/2/97
 * George
 *
 *
 * - Dominique LaSalle 2015-02-28
 * Modified to obfuscate names using the '__mtmetis_' prefix.
 *
 */


#ifndef _LIBMETIS_RENAME_H_
#define _LIBMETIS_RENAME_H_


/* balance.c */
#define Balance2Way			__mtmetis_Balance2Way
#define Bnd2WayBalance			__mtmetis_Bnd2WayBalance
#define General2WayBalance		__mtmetis_General2WayBalance
#define McGeneral2WayBalance            __mtmetis_McGeneral2WayBalance

/* bucketsort.c */
#define BucketSortKeysInc		__mtmetis_BucketSortKeysInc

/* checkgraph.c */
#define CheckGraph                      __mtmetis_CheckGraph
#define CheckInputGraphWeights          __mtmetis_CheckInputGraphWeights
#define FixGraph                        __mtmetis_FixGraph

/* coarsen.c */
#define CoarsenGraph			__mtmetis_CoarsenGraph
#define Match_RM                        __mtmetis_Match_RM
#define Match_SHEM                      __mtmetis_Match_SHEM
#define Match_2Hop                      __mtmetis_Match_2Hop
#define Match_2HopAny                   __mtmetis_Match_2HopAny
#define Match_2HopAll                   __mtmetis_Match_2HopAll
#define PrintCGraphStats                __mtmetis_PrintCGraphStats
#define CreateCoarseGraph		__mtmetis_CreateCoarseGraph
#define CreateCoarseGraphNoMask		__mtmetis_CreateCoarseGraphNoMask
#define CreateCoarseGraphPerm		__mtmetis_CreateCoarseGraphPerm
#define SetupCoarseGraph		__mtmetis_SetupCoarseGraph
#define ReAdjustMemory			__mtmetis_ReAdjustMemory

/* compress.c */
#define CompressGraph			__mtmetis_CompressGraph
#define PruneGraph			__mtmetis_PruneGraph

/* contig.c */
#define FindPartitionInducedComponents  __mtmetis_FindPartitionInducedComponents   
#define IsConnected                     __mtmetis_IsConnected
#define IsConnectedSubdomain            __mtmetis_IsConnectedSubdomain
#define FindSepInducedComponents        __mtmetis_FindSepInducedComponents
#define EliminateComponents             __mtmetis_EliminateComponents
#define MoveGroupContigForCut           __mtmetis_MoveGroupContigForCut
#define MoveGroupContigForVol           __mtmetis_MoveGroupContigForVol

/* debug.c */
#define ComputeCut			__mtmetis_ComputeCut
#define ComputeVolume			__mtmetis_ComputeVolume
#define ComputeMaxCut			__mtmetis_ComputeMaxCut
#define CheckBnd			__mtmetis_CheckBnd
#define CheckBnd2			__mtmetis_CheckBnd2
#define CheckNodeBnd			__mtmetis_CheckNodeBnd
#define CheckRInfo			__mtmetis_CheckRInfo
#define CheckNodePartitionParams	__mtmetis_CheckNodePartitionParams
#define IsSeparable			__mtmetis_IsSeparable
#define CheckKWayVolPartitionParams     __mtmetis_CheckKWayVolPartitionParams

/* fm.c */
#define FM_2WayRefine                   __mtmetis_FM_2WayRefine
#define FM_2WayCutRefine                __mtmetis_FM_2WayCutRefine
#define FM_Mc2WayCutRefine              __mtmetis_FM_Mc2WayCutRefine
#define SelectQueue                     __mtmetis_SelectQueue
#define Print2WayRefineStats            __mtmetis_Print2WayRefineStats

/* fortran.c */
#define Change2CNumbering		__mtmetis_Change2CNumbering
#define Change2FNumbering		__mtmetis_Change2FNumbering
#define Change2FNumbering2		__mtmetis_Change2FNumbering2
#define Change2FNumberingOrder		__mtmetis_Change2FNumberingOrder
#define ChangeMesh2CNumbering		__mtmetis_ChangeMesh2CNumbering
#define ChangeMesh2FNumbering		__mtmetis_ChangeMesh2FNumbering
#define ChangeMesh2FNumbering2		__mtmetis_ChangeMesh2FNumbering2

/* graph.c */
#define SetupGraph			__mtmetis_SetupGraph
#define SetupGraph_adjrsum              __mtmetis_SetupGraph_adjrsum
#define SetupGraph_tvwgt                __mtmetis_SetupGraph_tvwgt
#define SetupGraph_label                __mtmetis_SetupGraph_label
#define SetupSplitGraph                 __mtmetis_SetupSplitGraph
#define CreateGraph                     __mtmetis_CreateGraph
#define InitGraph                       __mtmetis_InitGraph
#define FreeRData                       __mtmetis_FreeRData
#define FreeGraph                       __mtmetis_FreeGraph
#define graph_WriteToDisk               __mtmetis_graph_WriteToDisk
#define graph_ReadFromDisk              __mtmetis_graph_ReadFromDisk

/* initpart.c */
#define Init2WayPartition		__mtmetis_Init2WayPartition
#define InitSeparator			__mtmetis_InitSeparator
#define RandomBisection			__mtmetis_RandomBisection
#define GrowBisection			__mtmetis_GrowBisection
#define McRandomBisection               __mtmetis_McRandomBisection
#define McGrowBisection                 __mtmetis_McGrowBisection
#define GrowBisectionNode		__mtmetis_GrowBisectionNode

/* kmetis.c */
#define MlevelKWayPartitioning		__mtmetis_MlevelKWayPartitioning
#define InitKWayPartitioning            __mtmetis_InitKWayPartitioning
#define InitKWayPartitioningRB          __mtmetis_InitKWayPartitioningRB
#define InitKWayPartitioningGrow        __mtmetis_InitKWayPartitioningGrow

/* kwayfm.c */
#define Greedy_KWayOptimize		__mtmetis_Greedy_KWayOptimize
#define Greedy_KWayCutOptimize		__mtmetis_Greedy_KWayCutOptimize
#define Greedy_KWayVolOptimize          __mtmetis_Greedy_KWayVolOptimize
#define Greedy_McKWayCutOptimize        __mtmetis_Greedy_McKWayCutOptimize
#define Greedy_McKWayVolOptimize        __mtmetis_Greedy_McKWayVolOptimize
#define IsArticulationNode              __mtmetis_IsArticulationNode
#define KWayVolUpdate                   __mtmetis_KWayVolUpdate

/* kwayrefine.c */
#define RefineKWay			__mtmetis_RefineKWay
#define AllocateKWayPartitionMemory	__mtmetis_AllocateKWayPartitionMemory
#define ComputeKWayPartitionParams	__mtmetis_ComputeKWayPartitionParams
#define ProjectKWayPartition		__mtmetis_ProjectKWayPartition
#define ComputeKWayBoundary		__mtmetis_ComputeKWayBoundary
#define ComputeKWayVolGains             __mtmetis_ComputeKWayVolGains
#define IsBalanced			__mtmetis_IsBalanced

/* mcutil */
#define rvecle                          __mtmetis_rvecle
#define rvecge                          __mtmetis_rvecge
#define rvecsumle                       __mtmetis_rvecsumle
#define rvecmaxdiff                     __mtmetis_rvecmaxdiff
#define ivecle                          __mtmetis_ivecle
#define ivecge                          __mtmetis_ivecge
#define ivecaxpylez                     __mtmetis_ivecaxpylez
#define ivecaxpygez                     __mtmetis_ivecaxpygez
#define BetterVBalance                  __mtmetis_BetterVBalance
#define BetterBalance2Way               __mtmetis_BetterBalance2Way
#define BetterBalanceKWay               __mtmetis_BetterBalanceKWay
#define ComputeLoadImbalance            __mtmetis_ComputeLoadImbalance
#define ComputeLoadImbalanceDiff        __mtmetis_ComputeLoadImbalanceDiff
#define ComputeLoadImbalanceDiffVec     __mtmetis_ComputeLoadImbalanceDiffVec
#define ComputeLoadImbalanceVec         __mtmetis_ComputeLoadImbalanceVec

/* mesh.c */
#define CreateGraphDual                 __mtmetis_CreateGraphDual
#define FindCommonElements              __mtmetis_FindCommonElements
#define CreateGraphNodal                __mtmetis_CreateGraphNodal
#define FindCommonNodes                 __mtmetis_FindCommonNodes
#define CreateMesh                      __mtmetis_CreateMesh
#define InitMesh                        __mtmetis_InitMesh
#define FreeMesh                        __mtmetis_FreeMesh

/* meshpart.c */
#define InduceRowPartFromColumnPart     __mtmetis_InduceRowPartFromColumnPart

/* minconn.c */
#define ComputeSubDomainGraph           __mtmetis_ComputeSubDomainGraph
#define UpdateEdgeSubDomainGraph        __mtmetis_UpdateEdgeSubDomainGraph
#define PrintSubDomainGraph             __mtmetis_PrintSubDomainGraph
#define EliminateSubDomainEdges         __mtmetis_EliminateSubDomainEdges
#define MoveGroupMinConnForCut          __mtmetis_MoveGroupMinConnForCut
#define MoveGroupMinConnForVol          __mtmetis_MoveGroupMinConnForVol

/* mincover.c */
#define MinCover			__mtmetis_MinCover
#define MinCover_Augment		__mtmetis_MinCover_Augment
#define MinCover_Decompose		__mtmetis_MinCover_Decompose
#define MinCover_ColDFS			__mtmetis_MinCover_ColDFS
#define MinCover_RowDFS			__mtmetis_MinCover_RowDFS

/* mmd.c */
#define genmmd				__mtmetis_genmmd
#define mmdelm				__mtmetis_mmdelm
#define mmdint				__mtmetis_mmdint
#define mmdnum				__mtmetis_mmdnum
#define mmdupd				__mtmetis_mmdupd


/* ometis.c */
#define MlevelNestedDissection		__mtmetis_MlevelNestedDissection
#define MlevelNestedDissectionCC	__mtmetis_MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple	__mtmetis_MlevelNodeBisectionMultiple
#define MlevelNodeBisectionL2		__mtmetis_MlevelNodeBisectionL2
#define MlevelNodeBisectionL1		__mtmetis_MlevelNodeBisectionL1
#define SplitGraphOrder			__mtmetis_SplitGraphOrder
#define SplitGraphOrderCC		__mtmetis_SplitGraphOrderCC
#define MMDOrder			__mtmetis_MMDOrder

/* options.c */
#define SetupCtrl                       __mtmetis_SetupCtrl
#define SetupKWayBalMultipliers         __mtmetis_SetupKWayBalMultipliers
#define Setup2WayBalMultipliers         __mtmetis_Setup2WayBalMultipliers
#define PrintCtrl                       __mtmetis_PrintCtrl
#define FreeCtrl                        __mtmetis_FreeCtrl
#define CheckParams                     __mtmetis_CheckParams

/* parmetis.c */
#define MlevelNestedDissectionP		__mtmetis_MlevelNestedDissectionP
#define FM_2WayNodeRefine1SidedP        __mtmetis_FM_2WayNodeRefine1SidedP
#define FM_2WayNodeRefine2SidedP        __mtmetis_FM_2WayNodeRefine2SidedP

/* pmetis.c */
#define MlevelRecursiveBisection	__mtmetis_MlevelRecursiveBisection
#define MultilevelBisect		__mtmetis_MultilevelBisect
#define SplitGraphPart			__mtmetis_SplitGraphPart

/* refine.c */
#define Refine2Way			__mtmetis_Refine2Way
#define Allocate2WayPartitionMemory	__mtmetis_Allocate2WayPartitionMemory
#define Compute2WayPartitionParams	__mtmetis_Compute2WayPartitionParams
#define Project2WayPartition		__mtmetis_Project2WayPartition

/* separator.c */
#define ConstructSeparator		__mtmetis_ConstructSeparator
#define ConstructMinCoverSeparator	__mtmetis_ConstructMinCoverSeparator

/* sfm.c */
#define FM_2WayNodeRefine2Sided         __mtmetis_FM_2WayNodeRefine2Sided 
#define FM_2WayNodeRefine1Sided         __mtmetis_FM_2WayNodeRefine1Sided
#define FM_2WayNodeBalance              __mtmetis_FM_2WayNodeBalance

/* srefine.c */
#define Refine2WayNode			__mtmetis_Refine2WayNode
#define Allocate2WayNodePartitionMemory	__mtmetis_Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams	__mtmetis_Compute2WayNodePartitionParams
#define Project2WayNodePartition	__mtmetis_Project2WayNodePartition

/* stat.c */
#define ComputePartitionInfoBipartite   __mtmetis_ComputePartitionInfoBipartite
#define ComputePartitionBalance		__mtmetis_ComputePartitionBalance
#define ComputeElementBalance		__mtmetis_ComputeElementBalance

/* timing.c */
#define InitTimers			__mtmetis_InitTimers
#define PrintTimers			__mtmetis_PrintTimers

/* util.c */
#define irandInRange                    __mtmetis_irandInRange
#define irandArrayPermute               __mtmetis_irandArrayPermute
#define iargmax_strd                    __mtmetis_iargmax_strd 
#define iargmax_nrm                     __mtmetis_iargmax_nrm
#define iargmax2_nrm                    __mtmetis_iargmax2_nrm
#define rargmax2                        __mtmetis_rargmax2
#define InitRandom                      __mtmetis_InitRandom
#define metis_rcode                     __mtmetis_metis_rcode

/* wspace.c */
#define AllocateWorkSpace               __mtmetis_AllocateWorkSpace                  
#define AllocateRefinementWorkSpace     __mtmetis_AllocateRefinementWorkSpace
#define FreeWorkSpace                   __mtmetis_FreeWorkSpace
#define wspacemalloc                    __mtmetis_wspacemalloc
#define wspacepush                      __mtmetis_wspacepush
#define wspacepop                       __mtmetis_wspacepop
#define iwspacemalloc                   __mtmetis_iwspacemalloc
#define rwspacemalloc                   __mtmetis_rwspacemalloc
#define ikvwspacemalloc                 __mtmetis_ikvwspacemalloc
#define cnbrpoolReset                   __mtmetis_cnbrpoolReset
#define cnbrpoolGetNext                 __mtmetis_cnbrpoolGetNext
#define vnbrpoolReset                   __mtmetis_vnbrpoolReset
#define vnbrpoolGetNext                 __mtmetis_vnbrpoolGetNext


/* Hide toplevel names */
#define METIS_PartGraphRecursive __METIS_PartGraphRecursive 
#define METIS_PartGraphKway __METIS_PartGraphKway
#define METIS_NodeND __METIS_NodeND
#define METIS_Free __METIS_Free
#define METIS_SetDefaultOptions __METIS_SetDefaultOptions
#define METIS_ComputeVertexSeparator __METIS_ComputeVertexSeparator
#define METIS_NodeRefine __METIS_NodeRefine
 

#endif



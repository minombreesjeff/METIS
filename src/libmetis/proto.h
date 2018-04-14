/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h,v 1.27 2003/05/03 16:10:48 karypis Exp $
 *
 */

/* balance.c */
void Balance2Way(CtrlType *, GraphType *, idxtype *, float);
void Bnd2WayBalance(CtrlType *, GraphType *, idxtype *);
void General2WayBalance(CtrlType *, GraphType *, idxtype *);

/* bucketsort.c */
void BucketSortKeysInc(int, int, idxtype *, idxtype *, idxtype *);

/* ccgraph.c */
void CreateCoarseGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraphNoMask(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraph_NVW(CtrlType *, GraphType *, int, idxtype *, idxtype *);
GraphType *SetUpCoarseGraph(GraphType *, int, int);
void ReAdjustMemory(GraphType *, GraphType *, int);

/* coarsen.c */
GraphType *Coarsen2Way(CtrlType *, GraphType *);

/* compress.c */
void CompressGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *, idxtype *, idxtype *);
void PruneGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *, idxtype *, float);

/* debug.c */
idxtype ComputeCut(GraphType *, idxtype *);
idxtype ComputeMaxCut(GraphType *graph, idxtype nparts, idxtype *where);
idxtype CheckBnd(GraphType *);
idxtype CheckBnd2(GraphType *);
idxtype CheckNodeBnd(GraphType *, int);
idxtype CheckRInfo(RInfoType *);
idxtype CheckNodePartitionParams(GraphType *);
idxtype IsSeparable(GraphType *);

/* estmem.c */
void EstimateCFraction(int, idxtype *, idxtype *, float *, float *);
idxtype ComputeCoarseGraphSize(int, idxtype *, idxtype *, int, idxtype *, idxtype *, idxtype *);

/* fm.c */
void FM_2WayEdgeRefine(CtrlType *, GraphType *, idxtype *, int);

/* fortran.c */
void Change2CNumbering(int, idxtype *, idxtype *);
void Change2FNumbering(int, idxtype *, idxtype *, idxtype *);
void Change2FNumbering2(int, idxtype *, idxtype *);
void Change2FNumberingOrder(int, idxtype *, idxtype *, idxtype *, idxtype *);
void ChangeMesh2CNumbering(int, idxtype *);
void ChangeMesh2FNumbering(int, idxtype *, int, idxtype *, idxtype *);
void ChangeMesh2FNumbering2(int, idxtype *, int, int, idxtype *, idxtype *);
void ChangeMesh2FNumbering3(int, idxtype *);

/* frename.c */
void METIS_PARTGRAPHRECURSIVE(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphrecursive(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphrecursive_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphrecursive__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_WPARTGRAPHRECURSIVE(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphrecursive(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphrecursive_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphrecursive__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void METIS_PARTGRAPHKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_WPARTGRAPHKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void METIS_EDGEND(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_edgend(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_edgend_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_edgend__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_NODEND(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodend(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodend_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodend__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_NODEWND(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodewnd(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodewnd_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodewnd__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_PARTMESHNODAL(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partmeshnodal(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partmeshnodal_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partmeshnodal__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_PARTMESHDUAL(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, idxtype *);
void metis_partmeshdual(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, idxtype *);
void metis_partmeshdual_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, idxtype *);
void metis_partmeshdual__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, idxtype *);
void METIS_MESHTONODAL(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtonodal(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtonodal_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtonodal__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_MESHTODUAL(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtodual(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtodual_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtodual__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_ESTIMATEMEMORY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_estimatememory(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_estimatememory_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_estimatememory__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_MCPARTGRAPHRECURSIVE(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphrecursive(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphrecursive_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphrecursive__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_MCPARTGRAPHKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void METIS_PARTGRAPHVKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partgraphvkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partgraphvkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partgraphvkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_WPARTGRAPHVKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_wpartgraphvkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_wpartgraphvkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_wpartgraphvkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);

/* graph.c */
void SetUpGraph(GraphType *, int, int, int, idxtype *, idxtype *, idxtype *, idxtype *, int);
void SetUpGraphKway(GraphType *, int, idxtype *, idxtype *);
void SetUpGraph2(GraphType *, int, int, idxtype *, idxtype *, float *, idxtype *);
void VolSetUpGraph(GraphType *, int, int, int, idxtype *, idxtype *, idxtype *, idxtype *, int);
void RandomizeGraph(GraphType *);
idxtype IsConnectedSubdomain(CtrlType *, GraphType *, int, int);
idxtype IsConnected(CtrlType *, GraphType *, int);
idxtype IsConnected2(GraphType *, int);
idxtype FindComponents(CtrlType *, GraphType *, idxtype *, idxtype *);

/* initpart.c */
void Init2WayPartition(CtrlType *, GraphType *, idxtype *, float);
void InitSeparator(CtrlType *, GraphType *, float);
void GrowBisection(CtrlType *, GraphType *, idxtype *, float);
void GrowBisectionNode(CtrlType *, GraphType *, float);
void RandomBisection(CtrlType *, GraphType *, idxtype *, float);

/* kmetis.c */
idxtype MlevelKWayPartitioning(CtrlType *, GraphType *, int, idxtype *, float *, float);

/* kvmetis.c */
idxtype MlevelVolKWayPartitioning(CtrlType *, GraphType *, int, idxtype *, float *, float);

/* kwayfm.c */
void Random_KWayEdgeRefine(CtrlType *, GraphType *, int, float *, float, int, int);
void Greedy_KWayEdgeRefine(CtrlType *, GraphType *, int, float *, float, int);
void Greedy_KWayEdgeBalance(CtrlType *, GraphType *, int, float *, float, int);

/* kwayrefine.c */
void RefineKWay(CtrlType *, GraphType *, GraphType *, int, float *, float);
void RefineKWayRefinement(CtrlType *, GraphType *, GraphType *, int, float *, float);
void AllocateKWayPartitionMemory(CtrlType *, GraphType *, int);
void ComputeKWayPartitionParams(CtrlType *, GraphType *, int);
void ProjectKWayPartition(CtrlType *, GraphType *, int);
idxtype IsBalanced(idxtype *, int, float *, float);
void ComputeKWayBoundary(CtrlType *, GraphType *, int);
void ComputeKWayBalanceBoundary(CtrlType *, GraphType *, int);

/* kwayvolfm.c */
void Random_KWayVolRefine(CtrlType *, GraphType *, int, float *, float, int, int);
void Random_KWayVolRefineMConn(CtrlType *, GraphType *, int, float *, float, int, int);
void Greedy_KWayVolBalance(CtrlType *, GraphType *, int, float *, float, int);
void Greedy_KWayVolBalanceMConn(CtrlType *, GraphType *, int, float *, float, int);
void KWayVolUpdate(CtrlType *, GraphType *, int, int, int, idxtype *, idxtype *, idxtype *);
void ComputeKWayVolume(GraphType *, int, idxtype *, idxtype *, idxtype *);
idxtype ComputeVolume(GraphType *, idxtype *);
void CheckVolKWayPartitionParams(CtrlType *, GraphType *, int);
void ComputeVolSubDomainGraph(GraphType *, int, idxtype *, idxtype *);
void EliminateVolSubDomainEdges(CtrlType *, GraphType *, int, float *);
void EliminateVolComponents(CtrlType *, GraphType *, int, float *, float);

/* kwayvolrefine.c */
void RefineVolKWay(CtrlType *, GraphType *, GraphType *, int, float *, float);
void AllocateVolKWayPartitionMemory(CtrlType *, GraphType *, int);
void ComputeVolKWayPartitionParams(CtrlType *, GraphType *, int);
void ComputeKWayVolGains(CtrlType *, GraphType *, int);
void ProjectVolKWayPartition(CtrlType *, GraphType *, int);
void ComputeVolKWayBoundary(CtrlType *, GraphType *, int);
void ComputeVolKWayBalanceBoundary(CtrlType *, GraphType *, int);

/* match.c */
void Match_RM(CtrlType *, GraphType *);
void Match_RM_NVW(CtrlType *, GraphType *);
void Match_HEM(CtrlType *, GraphType *);
void Match_SHEM(CtrlType *, GraphType *);

/* mbalance.c */
void MocBalance2Way(CtrlType *, GraphType *, float *, float);
void MocGeneral2WayBalance(CtrlType *, GraphType *, float *, float);

/* mbalance2.c */
void MocBalance2Way2(CtrlType *, GraphType *, float *, float *);
void MocGeneral2WayBalance2(CtrlType *, GraphType *, float *, float *);
void SelectQueue3(int, float *, float *, idxtype *, idxtype *, PQueueType [MAXNCON][2], float *);

/* mcoarsen.c */
GraphType *MCCoarsen2Way(CtrlType *, GraphType *);

/* memory.c */
void AllocateWorkSpace(CtrlType *, GraphType *, int);
void FreeWorkSpace(CtrlType *, GraphType *);
idxtype WspaceAvail(CtrlType *);
idxtype *idxwspacemalloc(CtrlType *, int);
void idxwspacefree(CtrlType *, int);
float *fwspacemalloc(CtrlType *, int);
void fwspacefree(CtrlType *, int);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);

/* mesh.c */
idxtype GENDUALMETIS_COUNT(idxtype , idxtype , idxtype , idxtype *, idxtype *);
void GENDUALMETIS(int, int, int, idxtype *, idxtype *, idxtype *, idxtype *adjncy);
void TRINODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void TETNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void HEXNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void QUADNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void LINENODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);

/* meshpart.c */

/* mfm.c */
void MocFM_2WayEdgeRefine(CtrlType *, GraphType *, float *, int);
void SelectQueue(int, float *, float *, idxtype *, idxtype *, PQueueType [MAXNCON][2]);
idxtype BetterBalance(int, float *, float *, float *);
float Compute2WayHLoadImbalance(int, float *, float *);
void Compute2WayHLoadImbalanceVec(int, float *, float *, float *);

/* mfm2.c */
void MocFM_2WayEdgeRefine2(CtrlType *, GraphType *, float *, float *, int);
void SelectQueue2(int, float *, float *, idxtype *, idxtype *, PQueueType [MAXNCON][2], float *);
idxtype IsBetter2wayBalance(int, float *, float *, float *);

/* mincover.o */
void MinCover(idxtype *, idxtype *, int, int, idxtype *, idxtype *);
idxtype MinCover_Augment(idxtype *, idxtype *, int, idxtype *, idxtype *, idxtype *, int);
void MinCover_Decompose(idxtype *, idxtype *, int, int, idxtype *, idxtype *, idxtype *);
void MinCover_ColDFS(idxtype *, idxtype *, int, idxtype *, idxtype *, int);
void MinCover_RowDFS(idxtype *, idxtype *, int, idxtype *, idxtype *, int);

/* minitpart.c */
void MocInit2WayPartition(CtrlType *, GraphType *, float *, float);
void MocGrowBisection(CtrlType *, GraphType *, float *, float);
void MocRandomBisection(CtrlType *, GraphType *, float *, float);
void MocInit2WayBalance(CtrlType *, GraphType *, float *);
idxtype SelectQueueOneWay(idxtype ncon, float *npwgts, float *tpwgts, idxtype from, PQueueType queues[MAXNCON][2]);

/* minitpart2.c */
void MocInit2WayPartition2(CtrlType *, GraphType *, float *, float *);
void MocGrowBisection2(CtrlType *, GraphType *, float *, float *);
void MocGrowBisectionNew2(CtrlType *, GraphType *, float *, float *);
void MocInit2WayBalance2(CtrlType *, GraphType *, float *, float *);
idxtype SelectQueueOneWay2(idxtype ncon, float *pto, PQueueType queues[MAXNCON][2], float *ubvec);



/* mkmetis.c */
idxtype MCMlevelKWayPartitioning(CtrlType *, GraphType *, int, idxtype *, float *);

/* mkwayfmh.c */
void MCRandom_KWayEdgeRefineHorizontal(CtrlType *, GraphType *, int, float *, int);
void MCGreedy_KWayEdgeBalanceHorizontal(CtrlType *, GraphType *, int, float *, int);
idxtype AreAllHVwgtsBelow(int, float, float *, float, float *, float *);
idxtype AreAllHVwgtsAbove(int, float, float *, float, float *, float *);
void ComputeHKWayLoadImbalance(int, int, float *, float *);
idxtype MocIsHBalanced(int, int, float *, float *);
idxtype IsHBalanceBetterFT(int, int, float *, float *, float *, float *);
idxtype IsHBalanceBetterTT(int, int, float *, float *, float *, float *);

/* mkwayrefine.c */
void MocRefineKWayHorizontal(CtrlType *, GraphType *, GraphType *, int, float *);
void MocAllocateKWayPartitionMemory(CtrlType *, GraphType *, int);
void MocComputeKWayPartitionParams(CtrlType *, GraphType *, int);
void MocProjectKWayPartition(CtrlType *, GraphType *, int);
void MocComputeKWayBalanceBoundary(CtrlType *, GraphType *, int);

/* mmatch.c */
void MCMatch_RM(CtrlType *, GraphType *);
void MCMatch_HEM(CtrlType *, GraphType *);
void MCMatch_SHEM(CtrlType *, GraphType *);
void MCMatch_SHEBM(CtrlType *, GraphType *, int);
void MCMatch_SBHEM(CtrlType *, GraphType *, int);
float BetterVBalance(int, int, float *, float *, float *);
idxtype AreAllVwgtsBelowFast(int, float *, float *, float);

/* mmd.c */
void genmmd(int, idxtype *, idxtype *, idxtype *, idxtype *, idxtype , idxtype *, idxtype *, idxtype *, idxtype *, int, idxtype *);
void mmdelm(int, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, int);
idxtype  mmdint(int, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void mmdnum(int, idxtype *, idxtype *, idxtype *);
void mmdupd(int, int, idxtype *, idxtype *, int, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, idxtype *tag);

/* mpmetis.c */
idxtype MCMlevelRecursiveBisection(CtrlType *, GraphType *, int, idxtype *, float, int);
idxtype MCHMlevelRecursiveBisection(CtrlType *, GraphType *, int, idxtype *, float *, int);
void MCMlevelEdgeBisection(CtrlType *, GraphType *, float *, float);
void MCHMlevelEdgeBisection(CtrlType *, GraphType *, float *, float *);

/* mrefine.c */
void MocRefine2Way(CtrlType *, GraphType *, GraphType *, float *, float);
void MocAllocate2WayPartitionMemory(CtrlType *, GraphType *);
void MocCompute2WayPartitionParams(CtrlType *, GraphType *);
void MocProject2WayPartition(CtrlType *, GraphType *);

/* mrefine2.c */
void MocRefine2Way2(CtrlType *, GraphType *, GraphType *, float *, float *);

/* mutil.c */
idxtype AreAllVwgtsBelow(int, float, float *, float, float *, float);
idxtype AreAnyVwgtsBelow(int, float, float *, float, float *, float);
idxtype AreAllVwgtsAbove(int, float, float *, float, float *, float);
float ComputeLoadImbalance(int, int, float *, float *);
idxtype AreAllBelow(int, float *, float *);

/* myqsort.c */
void iidxsort(int, idxtype *);
void iintsort(int, idxtype *);
void ikeysort(int, KeyValueType *);
void ikeyvalsort(int, KeyValueType *);

/* ometis.c */
void MlevelNestedDissection(CtrlType *, GraphType *, idxtype *, float, int);
void MlevelNestedDissectionCC(CtrlType *, GraphType *, idxtype *, float, int);
void MlevelNodeBisectionMultiple(CtrlType *, GraphType *, idxtype *, float);
void MlevelNodeBisection(CtrlType *, GraphType *, idxtype *, float);
void SplitGraphOrder(CtrlType *, GraphType *, GraphType *, GraphType *);
void MMDOrder(CtrlType *, GraphType *, idxtype *, int);
idxtype SplitGraphOrderCC(CtrlType *, GraphType *, GraphType *, int, idxtype *, idxtype *);

/* parmetis.c */
void MlevelNestedDissectionP(CtrlType *, GraphType *, idxtype *, int, int, int, idxtype *);

/* pmetis.c */
idxtype MlevelRecursiveBisection(CtrlType *, GraphType *, int, idxtype *, float *, float, int);
void MlevelEdgeBisection(CtrlType *, GraphType *, idxtype *, float);
void SplitGraphPart(CtrlType *, GraphType *, GraphType *, GraphType *);
void SetUpSplitGraph(GraphType *, GraphType *, int, int);

/* pqueue.c */
void PQueueInit(CtrlType *ctrl, PQueueType *, int, int);
void PQueueReset(PQueueType *);
void PQueueFree(CtrlType *ctrl, PQueueType *);
idxtype PQueueGetSize(PQueueType *);
idxtype PQueueInsert(PQueueType *, int, int);
idxtype PQueueDelete(PQueueType *, int, int);
idxtype PQueueUpdate(PQueueType *, int, int, int);
void PQueueUpdateUp(PQueueType *, int, int, int);
idxtype PQueueGetMax(PQueueType *);
idxtype PQueueSeeMax(PQueueType *);
idxtype PQueueGetKey(PQueueType *);
idxtype CheckHeap(PQueueType *);

/* refine.c */
void Refine2Way(CtrlType *, GraphType *, GraphType *, idxtype *, float ubfactor);
void Allocate2WayPartitionMemory(CtrlType *, GraphType *);
void Compute2WayPartitionParams(CtrlType *, GraphType *);
void Project2WayPartition(CtrlType *, GraphType *);

/* separator.c */
void ConstructSeparator(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator0(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator(CtrlType *, GraphType *, float);

/* sfm.c */
void FM_2WayNodeRefine(CtrlType *, GraphType *, float, int);
void FM_2WayNodeRefineEqWgt(CtrlType *, GraphType *, int);
void FM_2WayNodeRefine_OneSided(CtrlType *, GraphType *, float, int);
void FM_2WayNodeBalance(CtrlType *, GraphType *, float);
idxtype ComputeMaxNodeGain(int, idxtype *, idxtype *, idxtype *);

/* srefine.c */
void Refine2WayNode(CtrlType *, GraphType *, GraphType *, float);
void Allocate2WayNodePartitionMemory(CtrlType *, GraphType *);
void Compute2WayNodePartitionParams(CtrlType *, GraphType *);
void Project2WayNodePartition(CtrlType *, GraphType *);

/* stat.c */
void ComputePartitionInfo(GraphType *, int, idxtype *);
void ComputePartitionInfoBipartite(GraphType *, int, idxtype *);
void ComputePartitionBalance(GraphType *, int, idxtype *, float *);
float ComputeElementBalance(int, int, idxtype *);

/* subdomains.c */
void Random_KWayEdgeRefineMConn(CtrlType *, GraphType *, int, float *, float, int, int);
void Greedy_KWayEdgeBalanceMConn(CtrlType *, GraphType *, int, float *, float, int);
void PrintSubDomainGraph(GraphType *, int, idxtype *);
void ComputeSubDomainGraph(GraphType *, int, idxtype *, idxtype *);
void EliminateSubDomainEdges(CtrlType *, GraphType *, int, float *);
void MoveGroupMConn(CtrlType *, GraphType *, idxtype *, idxtype *, int, int, int, idxtype *);
void EliminateComponents(CtrlType *, GraphType *, int, float *, float);
void MoveGroup(CtrlType *, GraphType *, int, int, int, idxtype *, idxtype *);

/* timing.c */
void InitTimers(CtrlType *);
void PrintTimers(CtrlType *);
double seconds(void);

/* util.c */
void errexit(char *,...);
idxtype *imalloc(size_t, char *);
idxtype *idxmalloc(size_t, char *);
float *fmalloc(size_t, char *);
double *dmalloc(size_t, char *);
idxtype *ismalloc(size_t, int, char *);
idxtype *idxsmalloc(size_t, idxtype, char *);
double *dsmalloc(size_t, double, char *);
void *GKmalloc(size_t, char *);
void GKfree(void **,...); 
idxtype *iset(size_t n, idxtype val, idxtype *x);
idxtype *idxset(size_t n, idxtype val, idxtype *x);
double *dset(size_t n, double val, double *x);
float *sset(size_t n, float val, float *x);
idxtype iamax(size_t, idxtype *);
idxtype idxamax(size_t, idxtype *);
idxtype idxamax_strd(size_t, idxtype *, int);
idxtype samax(size_t, float *);
idxtype samax2(size_t, float *);
idxtype idxamin(size_t, idxtype *);
idxtype samin(size_t, float *);
idxtype idxsum(size_t, idxtype *);
idxtype idxsum_strd(size_t, idxtype *, int);
void idxadd(size_t, idxtype *, idxtype *);
idxtype charsum(size_t, char *);
idxtype isum(size_t, idxtype *);
float ssum(size_t, float *);
float ssum_strd(size_t n, float *x, int);
void sscale(size_t n, float, float *x);
float snorm2(size_t, float *);
float sdot(size_t n, float *, float *);
void saxpy(size_t, float, float *, int, float *, int);
void RandomPermute(size_t, idxtype *, int);
double drand48();
void srand48(long);
idxtype ispow2(int);
void InitRandom(int);
idxtype log2i(int);
FILE *GKfopen(char *fname, char *mode, char *msg);
void GKfclose(FILE *fp);











/***************************************************************
* Programs Directory
****************************************************************/

/* io.c */
void ReadGraph(GraphType *, char *, idxtype *);
void WritePartition(char *, idxtype *, int, int);
void WriteMeshPartition(char *, int, int, idxtype *, int, idxtype *);
void WritePermutation(char *, idxtype *, int);
idxtype CheckGraph(GraphType *);
idxtype MeshType(char *);
idxtype *ReadWgt(char *, idxtype *, idxtype *, idxtype *);
idxtype *ReadMesh(char *, idxtype *, idxtype *, idxtype *);
idxtype *ReadMeshWgt(char *, idxtype *, idxtype *, idxtype *, idxtype *);
idxtype *ReadMixedMesh(char *, idxtype *, idxtype *, idxtype *);
idxtype *ReadMixedMeshWgt(char *, idxtype *, idxtype *, idxtype *, idxtype *);
void WriteGraph(char *, int, idxtype *, idxtype *);
idxtype MixedElements(char *);
idxtype *ReadMgcnums(char *);
void WriteWgtGraph(char *, idxtype , idxtype *, idxtype *, idxtype *);


/* smbfactor.c */
void ComputeFillIn(GraphType *, idxtype *);
idxtype ComputeFillIn2(GraphType *, idxtype *);
idxtype smbfct(int, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);

/* kfmetis.c */
void BalanceFillIn(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part);
GraphType *ExtractPartitionGraph(GraphType *graph, idxtype *part, idxtype pid, idxtype *vmap, idxtype *vimap);
void ComputePartitionFillIn(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part, idxtype *spart, idxtype pid, idxtype *r_fill, idxtype *r_subfill);
void RefineTopLevelSeparators(CtrlType *ctrl, GraphType *graph, idxtype nparts,  idxtype *part, idxtype *spart, idxtype *, idxtype *, idxtype *, idxtype *);


/* rkmetis.c */
void METIS_RefineGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part);
void METIS_WRefineGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part);
idxtype MlevelKWayRefinement(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part, float *tpwgts, float ubfactor);


/* cmetis.c */
void *METIS_PartGraphForContact(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy,
                double *xyzcoords, idxtype *sflag, idxtype *numflag, idxtype *nparts,
                idxtype *options, idxtype *edgecut, idxtype *part);
void METIS_UpdateContactInfo(void *raw_cinfo, idxtype *nvtxs, double *xyzcoords, idxtype *sflag);
void *METIS_SetupContact0(idxtype *nvtxs, double *xyzcoords, idxtype *sflag,
                idxtype *nparts, idxtype *part);
void *METIS_SetupContact(idxtype *nvtxs, double *xyzcoords, idxtype *sflag,
                idxtype *nparts, idxtype *part);
void METIS_FindContacts(void *raw_cinfo, idxtype *nboxes, double *boxcoords, idxtype *nparts, 
               idxtype **r_cntptr, idxtype **r_cntind);
void METIS_FreeContactInfo(void *raw_cinfo);
GraphType *CreatePartitionGraphForContact(idxtype nvtxs, idxtype *xadj, idxtype *adjncy, 
                idxtype *vwgt, idxtype *adjwgt, idxtype cnvtxs, idxtype *part);
idxtype InduceDecissionTree(idxtype nvtxs, DKeyValueType **xyzcand, idxtype *sflag, idxtype nparts, 
          idxtype *part, idxtype maxnvtxs, idxtype minnvtxs, float minfrac, idxtype *r_nnodes, idxtype *r_nlnodes,
          DTreeNodeType *dtree, idxtype *dtpart, idxtype *dtipart, idxtype *r_nclean,
          idxtype *r_naclean, idxtype *r_ndirty, idxtype *r_maxdepth, idxtype *marker);
idxtype FindBoxContacts(ContactInfoType *cinfo, double *coords, idxtype *stack, idxtype *cntind, idxtype *marker);
void BuildDTLeafContents(ContactInfoType *cinfo, idxtype *sflag);
void  CheckDTree(idxtype nvtxs, double *xyzcoords, idxtype *part, ContactInfoType *cinfo);
void  CheckDTreeSurface(idxtype nvtxs, double *xyzcoords, idxtype *part, ContactInfoType *cinfo, idxtype *sflag);
void idkeysort_qsort(idxtype n, DKeyValueType *cand);
idxtype cmpr_idkey(const void *v1, const void *v2);
void ikeyvalsort_qsort(idxtype n, KeyValueType *cand);
idxtype cmpr_ikeyval(const void *v1, const void *v2);
void *METIS_PartSurfForContactRCB(idxtype *nvtxs, double *xyzcoords, idxtype *sflag,
                idxtype *nparts, idxtype *part, idxtype *bestdims);
idxtype InduceRCBTree(idxtype nvtxs, DKeyValueType **xyzcand, idxtype firstPID, idxtype nparts,
          idxtype *r_nnodes, idxtype *r_nlnodes, DTreeNodeType *dtree, idxtype *leafpart,
          idxtype *part, idxtype *marker, idxtype *oldBestDims);



/* mrkmetis.c */
void METIS_mCRefineGraphKway(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy,
                          idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag,
                          idxtype *nparts, float *rubvec, idxtype *options, idxtype *edgecut,
                          idxtype *part);
idxtype MCMlevelKWayRefinement(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part,
      float *rubvec);



/* idkeysort.o */
void idkeysort(idxtype total_elems, DKeyValueType *pbase);



/***************************************************************
* Test Directory
****************************************************************/
void Test_PartGraph(int, idxtype *, idxtype *);
idxtype VerifyPart(int, idxtype *, idxtype *, idxtype *, idxtype *, int, int, idxtype *);
idxtype VerifyWPart(int, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, int, idxtype *);
void Test_PartGraphV(int, idxtype *, idxtype *);
idxtype VerifyPartV(int, idxtype *, idxtype *, idxtype *, idxtype *, int, int, idxtype *);
idxtype VerifyWPartV(int, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, int, idxtype *);
void Test_PartGraphmC(int, idxtype *, idxtype *);
idxtype VerifyPartmC(int, int, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, int, idxtype *);
void Test_ND(int, idxtype *, idxtype *);
idxtype VerifyND(int, idxtype *, idxtype *);


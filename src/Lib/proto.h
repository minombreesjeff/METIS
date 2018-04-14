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
 * $Id: proto.h,v 1.8 1998/01/21 16:31:50 karypis Exp $
 *
 */


/* coarsen.c */
GraphType *Coarsen2Way(CtrlType *, GraphType *);

/* fm.c */
void FM_2WayEdgeRefine(CtrlType *, GraphType *, int *, int);

/* balance.c */
void Balance2Way(CtrlType *, GraphType *, int *, float);
void Bnd2WayBalance(CtrlType *, GraphType *, int *);
void General2WayBalance(CtrlType *, GraphType *, int *);

/* initpart.c */
void Init2WayPartition(CtrlType *, GraphType *, int *, float);
void InitSeparator(CtrlType *, GraphType *, float);
void GrowBisection(CtrlType *, GraphType *, int *, float);
void GrowBisectionNode(CtrlType *, GraphType *, float);
void RandomBisection(CtrlType *, GraphType *, int *, float);

/* io.c */
void ReadGraph(GraphType *, char *);
void WritePartition(char *, idxtype *, int, int);
void WriteMeshPartition(char *, int, int, idxtype *, int, idxtype *);
void WritePermutation(char *, idxtype *, int);
int CheckGraph(GraphType *);
idxtype *ReadMesh(char *, int *, int *, int *);
void WriteGraph(char *, int, idxtype *, idxtype *);

/* match.c */
void Match_RM(CtrlType *, GraphType *);
void Match_HEM(CtrlType *, GraphType *);
void Match_SHEM(CtrlType *, GraphType *);

/* ccgraph.c */
void CreateCoarseGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraphNoMask(CtrlType *, GraphType *, int, idxtype *, idxtype *);

/* memory.c */
void AllocateWorkSpace(CtrlType *, GraphType *, int);
void FreeWorkSpace(CtrlType *, GraphType *);
idxtype *idxwspacemalloc(CtrlType *, int);
void idxwspacefree(CtrlType *, int);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);

/* pmetis.c */
void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
int MlevelRecursiveBisection(CtrlType *, GraphType *, int, idxtype *, float *, float, int);
void MlevelEdgeBisection(CtrlType *, GraphType *, int *, float);
void SplitGraphPart(CtrlType *, GraphType *, GraphType *, GraphType *);

/* pqueue.c */
void PQueueInit(CtrlType *ctrl, PQueueType *, int, int);
void PQueueReset(PQueueType *);
void PQueueFree(CtrlType *ctrl, PQueueType *);
int PQueueInsert(PQueueType *, int, int);
int PQueueDelete(PQueueType *, int, int);
int PQueueUpdate(PQueueType *, int, int, int);
void PQueueUpdateUp(PQueueType *, int, int, int);
int PQueueGetMax(PQueueType *);
int PQueueSeeMax(PQueueType *);
int CheckHeap(PQueueType *);


/* refine.c */
void Refine2Way(CtrlType *, GraphType *, GraphType *, int *, float ubfactor);
void Allocate2WayPartitionMemory(CtrlType *, GraphType *);
void Compute2WayPartitionParams(CtrlType *, GraphType *);
void Project2WayPartition(CtrlType *, GraphType *);



/* util.c */
void errexit(char *,...);
#ifndef DMALLOC
int *imalloc(int, char *);
idxtype *idxmalloc(int, char *);
float *fmalloc(int, char *);
int *ismalloc(int, int, char *);
idxtype *idxsmalloc(int, idxtype, char *);
void *GKmalloc(int, char *);
#endif
/*void GKfree(void **,...); */
int *iset(int n, int val, int *x);
idxtype * idxset(int n, idxtype val, idxtype *x);
int idxamax(int, idxtype *);
int idxamin(int, idxtype *);
int idxasum(int, idxtype *);
int idxmax(int, idxtype *);
int idxmin(int, idxtype *);
int idxsum(int, idxtype *);
/* idxtype *idxcopy(int, idxtype *, idxtype *); */
float ssum(int n, float *x);
void sscale(int n, float, float *x);
float snorm2(int, float *);
float sdot(int n, float *, float *);
void saxpy(int, float, float *, float *);
void RandomPermute(int, idxtype *, int);
double drand48();
void srand48(long);
int ispow2(int);


/* timing.c */
void InitTimers(CtrlType *);
void PrintTimers(CtrlType *);
double seconds(void);

/* debug.c */
int ComputeCut(GraphType *, idxtype *);
int CheckBnd(GraphType *);
int CheckBnd2(GraphType *);
int CheckNodeBnd(GraphType *, int);
int CheckRInfo(RInfoType *);
int CheckNodePartitionParams(GraphType *);
int IsSeparable(GraphType *);

/* bucketsort.c */
void BucketSortKeysInc(int, int, idxtype *, idxtype *, idxtype *);


/* graph.c */
void SetUpGraph(GraphType *, int, int, idxtype *, idxtype *, idxtype *, idxtype *, int);
void RandomizeGraph(GraphType *);
int IsConnected(CtrlType *, GraphType *, int);
int IsConnectedSubdomain(CtrlType *, GraphType *, int, int);
int IsConnected2(GraphType *, int);
int FindComponents(CtrlType *, GraphType *, idxtype *, idxtype *);
void AnalyzeGraph(GraphType *);


/* stat.c */
void ComputePartitionInfo(GraphType *, int, idxtype *);
float ComputePartitionBalance(GraphType *, int, idxtype *);
float ComputeElementBalance(int, int, idxtype *);

/* kmetis.c */
void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
int MlevelKWayPartitioning(CtrlType *, GraphType *, int, idxtype *, float *, float);

/* kwayrefine.c */
void RefineKWay(CtrlType *, GraphType *, GraphType *, int, float *, float);
void AllocateKWayPartitionMemory(CtrlType *, GraphType *, int);
void ComputeKWayPartitionParams(CtrlType *, GraphType *, int);
void ProjectKWayPartition(CtrlType *, GraphType *, int);
int IsBalanced(idxtype *, int, float *, float);
void ComputeKWayBoundary(CtrlType *, GraphType *, int);
void ComputeKWayBalanceBoundary(CtrlType *, GraphType *, int);

/* kwayfm.c */
void Random_KWayEdgeRefine(CtrlType *, GraphType *, int, float *, float, int, int);
void Random_KWayEdgeBalance(CtrlType *, GraphType *, int, float *, float, int);
void Greedy_KWayEdgeRefine(CtrlType *, GraphType *, int, float *, float, int);
void Greedy_KWayEdgeBalance(CtrlType *, GraphType *, int, float *, float, int);

/* ometis.c */
void METIS_EdgeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_NodeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_NodeWND(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void MlevelNestedDissection(CtrlType *, GraphType *, idxtype *, float, int);
void MlevelNestedDissectionCC(CtrlType *, GraphType *, idxtype *, float, int);
void MlevelNodeBisectionMultiple(CtrlType *, GraphType *, int *, float);
void MlevelNodeBisection(CtrlType *, GraphType *, int *, float);
void SplitGraphOrder(CtrlType *, GraphType *, GraphType *, GraphType *);
void MMDOrder(CtrlType *, GraphType *, idxtype *, int);
int SplitGraphOrderCC(CtrlType *, GraphType *, GraphType *, int, idxtype *, idxtype *);

/* srefine.c */
void Refine2WayNode(CtrlType *, GraphType *, GraphType *, float);
void Allocate2WayNodePartitionMemory(CtrlType *, GraphType *);
void Compute2WayNodePartitionParams(CtrlType *, GraphType *);
void Project2WayNodePartition(CtrlType *, GraphType *);

/* sfm.c */
void FM_2WayNodeRefine(CtrlType *, GraphType *, float, int);
void FM_2WayNodeRefineEqWgt(CtrlType *, GraphType *, int);
void FM_2WayNodeRefine_OneSided(CtrlType *, GraphType *, float, int);
void FM_2WayNodeBalance(CtrlType *, GraphType *, float);
int ComputeMaxNodeGain(int, idxtype *, idxtype *, idxtype *);

/* separator.c */
void ConstructSeparator(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator0(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator(CtrlType *, GraphType *, float);


/* mincover.o */
void MinCover(idxtype *, idxtype *, int, int, idxtype *, int *);
int MinCover_Augment(idxtype *, idxtype *, int, idxtype *, idxtype *, idxtype *, int);
void MinCover_Decompose(idxtype *, idxtype *, int, int, idxtype *, idxtype *, int *);
void MinCover_ColDFS(idxtype *, idxtype *, int, idxtype *, idxtype *, int);
void MinCover_RowDFS(idxtype *, idxtype *, int, idxtype *, idxtype *, int);

/* mmd.c */
void genmmd(int, idxtype *, idxtype *, idxtype *, idxtype *, int , idxtype *, idxtype *, idxtype *, idxtype *, int, int *);
void mmdelm(int, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, int);
int  mmdint(int, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void mmdnum(int, idxtype *, idxtype *, idxtype *);
void mmdupd(int, int, idxtype *, idxtype *, int, int *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, int *tag);

/* smbfactor.c */
void ComputeFillIn(GraphType *, idxtype *);
idxtype ComputeFillIn2(GraphType *, idxtype *);
int smbfct(int, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, idxtype *, idxtype *, int *);

/* compress.c */
void CompressGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *, idxtype *, idxtype *);
void PruneGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *, idxtype *, float);

/* svrefine.c */
void VSepRefine(CtrlType *, GraphType *, float);

/* mesh.c */
void METIS_MeshToDual(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void METIS_MeshToNodal(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void GENDUALMETIS(int, int, int, idxtype *, idxtype *, idxtype *adjncy);
void TRINODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void TETNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void HEXNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);
void QUADNODALMETIS(int, int, idxtype *, idxtype *, idxtype *adjncy);

/* meshpart.c */
void METIS_PartMeshNodal(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void METIS_PartMeshDual(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);

/* fortran.c */
void Change2CNumbering(int, idxtype *, idxtype *);
void Change2FNumbering(int, idxtype *, idxtype *, idxtype *);
void Change2FNumbering2(int, idxtype *, idxtype *);
void Change2FNumberingOrder(int, idxtype *, idxtype *, idxtype *, idxtype *);
void ChangeMesh2CNumbering(int, idxtype *);
void ChangeMesh2FNumbering(int, idxtype *, int, idxtype *, idxtype *);
void ChangeMesh2FNumbering2(int, idxtype *, int, int, idxtype *, idxtype *);

/* frename.c */
void METIS_PARTGRAPHRECURSIVE(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphrecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphrecursive_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphrecursive__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPARTGRAPHRECURSIVE(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphrecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphrecursive_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphrecursive__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void METIS_PARTGRAPHKWAY(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphkway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphkway_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void metis_partgraphkway__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *); 
void METIS_WPARTGRAPHKWAY(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphkway(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphkway_(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void metis_wpartgraphkway__(int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *); 
void METIS_EDGEND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_edgend(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_edgend_(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_edgend__(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_NODEND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodend(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodend_(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodend__(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_NODEWND(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodewnd(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodewnd_(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void metis_nodewnd__(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *); 
void METIS_PARTMESHNODAL(int *, int *, idxtype *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_partmeshnodal(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshnodal_(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshnodal__(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void METIS_PARTMESHDUAL(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshdual(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshdual_(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void metis_partmeshdual__(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
void METIS_MESHTONODAL(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtonodal(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtonodal_(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtonodal__(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void METIS_MESHTODUAL(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtodual(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtodual_(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void metis_meshtodual__(int *, int *, idxtype *, int *, int *, idxtype *, idxtype *);
void METIS_ESTIMATEMEMORY(int *, idxtype *, idxtype *, int *, int *, int *);
void metis_estimatememory(int *, idxtype *, idxtype *, int *, int *, int *);
void metis_estimatememory_(int *, idxtype *, idxtype *, int *, int *, int *);
void metis_estimatememory__(int *, idxtype *, idxtype *, int *, int *, int *);

/* ometispar.c */
void METIS_NodeNDP(int, idxtype *, idxtype *, int, int *, idxtype *, idxtype *, idxtype *);
void MlevelNestedDissectionP(CtrlType *, GraphType *, idxtype *, int, int, int, idxtype *);


/* estmem.c */
void METIS_EstimateMemory(int *, idxtype *, idxtype *, int *, int *, int *);
void EstimateCFraction(int, idxtype *, idxtype *, float *, float *);
int ComputeCoarseGraphSize(int, idxtype *, idxtype *, int, idxtype *, idxtype *, idxtype *);

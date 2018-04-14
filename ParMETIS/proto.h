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
 * $Id: proto.h,v 1.9 1997/07/18 00:32:12 karypis Exp $
 *
 */

/* debug.c */
void PrintVector(CtrlType *, int, int, idxtype *, char *);
void PrintGraph(CtrlType *, GraphType *);
void PrintGraph2(CtrlType *, GraphType *);
void PrintSetUpInfo(CtrlType *ctrl, GraphType *graph);

/* comm.c */
void CommInterfaceData(CtrlType *, GraphType *, idxtype *, idxtype *, idxtype *);
void CommChangedInterfaceData(CtrlType *, GraphType *, int, idxtype *, idxtype *, KeyValueType *, KeyValueType *, idxtype *);
int GlobalSEMax(CtrlType *, int);
double GlobalSEMaxDouble(CtrlType *, double);
int GlobalSEMin(CtrlType *, int);
int GlobalSESum(CtrlType *, int);

/* io.c */
void ReadGraph(GraphType *, char *, MPI_Comm);
void ReadPartitionedGraph(GraphType *, char *, MPI_Comm);
float *ReadCoordinates(GraphType *, char *, MPI_Comm);
void WritePVector(char *, idxtype *, idxtype *, MPI_Comm);


/* util.c */
void errexit(char *,...);
void myprintf(CtrlType *, char *f_str,...);
void rprintf(CtrlType *, char *f_str,...);
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
int idxamax(int n, idxtype *x);
int idxamin(int n, idxtype *x);
int idxasum(int n, idxtype *x);
float snorm2(int, float *);
float sdot(int n, float *, float *);
void saxpy(int, float, float *, float *);
void ikeyvalsort_org(int, KeyValueType *);
int IncKeyValueCmp(const void *, const void *);
void dkeyvalsort(int, KeyValueType *);
int DecKeyValueCmp(const void *, const void *);
int BSearch(int, idxtype *, int);
void RandomPermute(int, idxtype *, int);
void FastRandomPermute(int, idxtype *, int);
double drand48();
void srand48(long);
int ispow2(int);


/* qsort_special.c */
void iidxsort(int, idxtype *);
void iintsort(int, int *);
void ikeysort(int, KeyValueType *);
void ikeyvalsort(int, KeyValueType *);

/* memory.c */
void PreAllocateMemory(CtrlType *, GraphType *, WorkSpaceType *);
void FreeWSpace(WorkSpaceType *);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);
void FreeInitialGraph(GraphType *, int, int);

/* kmetis.c */
GraphType *SetUpGraph(CtrlType *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void SetUpCtrl(CtrlType *ctrl, int *, MPI_Comm);
void ChangeNumbering(idxtype *, idxtype *, idxtype *, idxtype *, int, int, int);
void InitTimers(CtrlType *);
void PrintTimingInfo(CtrlType *);
void PrintTimer(CtrlType *, timer, char *);
void GraphRandomPermute(GraphType *);
void ComputeMoveStatistics(CtrlType *, GraphType *, int *, int *, int *);

/* setup.c */
void SetUp(CtrlType *, GraphType *, WorkSpaceType *);
int Home_PE(int, int, idxtype *, int);

/* coarsen.c */
void GlobalMatch_HEM(CtrlType *, GraphType *, WorkSpaceType *);
void Global_CreateCoarseGraph(CtrlType *, GraphType *, WorkSpaceType *, int);
void LocalMatch_HEM(CtrlType *, GraphType *, WorkSpaceType *);
void Local_CreateCoarseGraph(CtrlType *, GraphType *, WorkSpaceType *, int);


/* edge_refine.c */
void ProjectPartition(CtrlType *, GraphType *, WorkSpaceType *);
void ComputePartitionParams(CtrlType *, GraphType *, WorkSpaceType *);
void ComputeCut(CtrlType *, GraphType *, WorkSpaceType *);
void KWayRefine(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayRefineClean(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayAdaptiveRefineClean(CtrlType *, GraphType *, WorkSpaceType *, int, float);

/* node_refine.c */
void ComputeNodePartitionParams0(CtrlType *, GraphType *, WorkSpaceType *);
void ComputeNodePartitionParams(CtrlType *, GraphType *, WorkSpaceType *);
void KWayNodeRefine0(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayNodeRefine(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayNodeRefine2(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void PrintNodeBalanceInfo(CtrlType *, int, idxtype *, idxtype *, idxtype *, int);

/* initpart.c */
void InitPartition(CtrlType *, GraphType *, WorkSpaceType *, int);
void InitPartition_RB(CtrlType *, GraphType *, WorkSpaceType *, int);
GraphType *AssembleGraph(CtrlType *, GraphType *, WorkSpaceType *, int);
void KeepPart(CtrlType *, GraphType *, WorkSpaceType *, idxtype *, int);

/* initmsection.c */
void InitMultisection(CtrlType *, GraphType *, WorkSpaceType *);
GraphType *AssembleMultisectedGraph(CtrlType *, GraphType *, WorkSpaceType *);
void ConstructSeparator(CtrlType *, GraphType *, int);
void ConstructSeparator1(CtrlType *, GraphType *, int);


/* serkmetis.c */
void Ser_KMetis(GraphType *, int, float);
void Ser_OMetis(GraphType *, float);
GraphType *Ser_Coarsen(GraphType *, int, int);
void Ser_Match_HEM(GraphType *, int); 
void Ser_CreateCoarseGraph(GraphType *, int, idxtype *);
void Ser_InitPartition(GraphType *, int, float);
void Ser_RecursiveBisection(GraphType *, idxtype *, int, float, int, idxtype *);
void Ser_EdgeFM(GraphType *, int *, int *, int);
void Ser_MlevelBisectGraph(GraphType *, float, idxtype *);
void Ser_BisectGraph(GraphType *, int, float, idxtype *);
void Ser_SplitGraphPart(GraphType *, idxtype *, GraphType *, GraphType *);
void Ser_Refine(GraphType *, GraphType *, int, float, EdgeType *degrees);
void Ser_ProjectPartition(GraphType *);
void Ser_ComputePartitionParams(GraphType *, int, EdgeType *);
void Ser_KWayRefine(GraphType *, int, float, int);


/* drivers.c */
void Global_Partition(CtrlType *, GraphType *, WorkSpaceType *);
void FldGlobal_Partition(CtrlType *, GraphType *, WorkSpaceType *, int);
void Refine_Partition(CtrlType *, GraphType *, WorkSpaceType *);
void AdaptiveUndirected_Partition(CtrlType *, GraphType *, WorkSpaceType *);
void AdaptiveDirected_Partition(CtrlType *, GraphType *, WorkSpaceType *);
int SmallerSubGraph(CtrlType *, GraphType *);

/* fold.c */
int EnoughMemory(CtrlType *, int, int);
GraphType *FoldGraph(CtrlType *, GraphType *, WorkSpaceType *);
void UnFoldGraph(CtrlType *, GraphType *, GraphType *, WorkSpaceType *);

/* balance.c */
void BalancePartition(CtrlType *, GraphType *, WorkSpaceType *);

/* move.c */
GraphType *MoveGraph(CtrlType *, GraphType *, WorkSpaceType *);
void ProjectInfoBack(CtrlType *, GraphType *, idxtype *, idxtype *, WorkSpaceType *);
void FindVtxPerm(CtrlType *, GraphType *, idxtype *, WorkSpaceType *);

/* initdiff.c */
void InitDiffusion(CtrlType *, GraphType *, WorkSpaceType *);
GraphType *AssembleAdaptiveGraph(CtrlType *, GraphType *, WorkSpaceType *);
void KWay_InitialDiffuser(CtrlType *, GraphType *, int, float);
void KWay_InitialDiffuser2(GraphType *, int, float);
int ConjGrad(int, idxtype *, idxtype *, float *, float *, float *, float);
void mvMult(int, idxtype *, idxtype *, float *, float *, float *);
void setupLaplace(GraphType *, int, idxtype *, idxtype *, float *, idxtype *);
void Ser_KWayUpdateDegrees(GraphType *, int, int);

/* test.c */
void TestParMetis(char *, MPI_Comm);
void TestAdaptiveMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int, MPI_Comm);
void AdaptGraph(GraphType *, int, MPI_Comm);
int ComputeRealCut(idxtype *, idxtype *, char *, MPI_Comm);
int ComputeRealCut2(idxtype *, idxtype *, idxtype *, idxtype *, char *, MPI_Comm);
void TestMoveGraph(GraphType *, GraphType *, idxtype *, MPI_Comm);


/* order.c */
void MultilevelOrder(CtrlType *, GraphType *, idxtype *, idxtype *, WorkSpaceType *);
void LabelSeparators(CtrlType *, GraphType *, idxtype *, idxtype *, idxtype *, idxtype *, WorkSpaceType *);
void CompactGraph(CtrlType *, GraphType *, idxtype *, WorkSpaceType *);
void LocalOrder(CtrlType *, GraphType *, idxtype *, int, WorkSpaceType *);


/* pqueue.c */
void PQueueInit(PQueueType *, int);
void PQueueReset(PQueueType *);
void PQueueFree(PQueueType *);
int PQueueInsert(PQueueType *, int, int);
int PQueueUpdate(PQueueType *, int, int, int);
int PQueueDelete(PQueueType *, int);
int PQueueGetMax(PQueueType *);
int PQueueSeeMax(PQueueType *);
int PQueueCheck(PQueueType *);

/* serometis.c */
void Ser_NodeOMetis(GraphType *, float);
void Ser_NodeBisection(GraphType *, float);
void Ser_NodeRefine(GraphType *, GraphType *, float);
void Ser_NodeComputePartitionParams(GraphType *);
void Ser_NodeFM(GraphType *, float, int);
void Ser_ConstructSeparator(GraphType *);
int CheckPartitionParams(GraphType *);

/* mmd.c */
void genmmd(int, idxtype *, idxtype *, idxtype *, idxtype *, int , idxtype *, idxtype *, idxtype *, idxtype *, int, int *);
void mmdelm(int, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, int);
int  mmdint(int, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void mmdnum(int, idxtype *, idxtype *, idxtype *);
void mmdupd(int, int, idxtype *, idxtype *, int, int *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, int *tag);


/* xyzpart.c */
void Coordinate_Partition(CtrlType *, GraphType *, int, float *, int, WorkSpaceType *);
void PartSort(CtrlType *, GraphType *, KeyValueType *, WorkSpaceType *);

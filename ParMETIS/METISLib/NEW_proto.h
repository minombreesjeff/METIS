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
 * $Id: proto.h,v 1.6 1998/09/20 17:36:13 karypis Exp $
 *
 */

/* NEW_parmetis.c */
void METIS_mCPartGraphRecursive2(int *, int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *);
int MCMlevelRecursiveBisection2(CtrlType *, GraphType *, int, float *, idxtype *, float, int); 

/* NEW_stats.c */
void Moc_ComputePartitionBalance(GraphType *, int, idxtype *, float *);

/* NEW_checkgraph.c */
int CheckGraph(GraphType *);


/**
* @file proto.h
* @brief Prototypes for functions
* @author Dominique LaSalle <lasalle@cs.umn.edu> 
* Copyright 2012, Regents of the University of Minnesota
* @version 1
* @date 2012-06-13
*/


#ifndef PROTO_H
#define PROTO_H

#include "includes.h"

/* check.c */
idx_t checkWhere(dctrl_t * dctrl, const dgraph_t * graph, 
    const idx_t *const * where);
idx_t checkInfo(dctrl_t * dctrl, const dgraph_t * graph, 
    const idx_t * const* where); 
idx_t checkNED(dctrl_t * dctrl, const ckrinfo_t * myrinfo, idx_t i);
idx_t checkGraph(const dgraph_t * graph);
idx_t checkBND(const dgraph_t * graph); 

/* coarsen.c */
dgraph_t * ParCoarsenGraph(dctrl_t * dctrl, dgraph_t * graph);
idx_t ParMatch_SHEM(dctrl_t * dctrl, dgraph_t *graph, idx_t * const * gmatch,
    idx_t * fcmap);
idx_t ParMatch_RM(dctrl_t * dctrl, dgraph_t * graph, idx_t * const * gmatch,
    idx_t * fcmap);
void ParCreateCoarseGraph(dctrl_t * dctrl, dgraph_t * graph, idx_t cnvtxs,
    const idx_t * const * gmatch, const idx_t * fcmap);
void ParCreateCoarseGraphNoMask(dctrl_t * dctrl, dgraph_t * graph, idx_t cnvtxs,
    const idx_t * const * gmatch, const idx_t * fcmap);

/* graph.c */
dgraph_t * SerCreateGraph(void);
void SerInitGraph(dgraph_t * graph);
dgraph_t *ParSetupGraph(dgraph_t ** r_graph, idx_t * buffer, idx_t nvtxs, 
    const idx_t * xadj, const idx_t *adjncy, const idx_t * adjwgt, 
    const idx_t *vwgt);
dgraph_t * ParSetupCoarseGraph(idx_t * buffer, dgraph_t ** r_graph, 
    dgraph_t * graph, idx_t cnvtxs);
void ParFreeGraph(dgraph_t ** r_graph);
void ParFreeRData(dgraph_t * graph);
void ParSetupGraph_tvwgt(idx_t * buffer, dgraph_t *graph);
void ParSetupGraph_label(dgraph_t *graph);
real_t ParComputeLoadImbalance(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t nparts, const real_t *pijbm);
real_t ParComputeLoadImbalanceDiff(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t nparts, const real_t *pijbm, const real_t *ubvec);
real_t SerComputeLoadImbalance(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t nparts, const real_t *pijbm);
real_t SerComputeLoadImbalanceDiff(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t nparts, const real_t *pijbm, const real_t *ubvec);
idx_t ParComputeCut(dctrl_t * dctrl, const dgraph_t *graph, 
    const idx_t * const * where);
idx_t SerComputeCut(const dgraph_t *graph, const idx_t * const * where);
idx_t ParComputeVolume(const dgraph_t *graph, const idx_t *const*where);
idx_t SerComputeVolume(const dgraph_t *graph, const idx_t *const*where);
idx_t ParIsBalanced(dctrl_t *ctrl, const dgraph_t *graph, const real_t ffactor);
graph_t * ParConvertDGraph(dctrl_t * dctrl, dgraph_t * ograph, graph_t ** rg);
idx_t computeEEW(dgraph_t * graph);
void nParReAdjustMemory(const idx_t myid, dctrl_t * dctrl, dgraph_t * graph,
    idx_t nthreads);

/* io.c */
idx_t ReadGraph(const char * filename, idx_t ** xadj, idx_t ** vwgt,
    idx_t ** adjncy, idx_t ** adjwgt);

/* initpart.c */
idx_t ParInitKWayPartitioning(dctrl_t * dctrl, dgraph_t * graph);

/* kpart.c */
idx_t ParPartGraphKway(dctrl_t * dctrl, dgraph_t * graph, idx_t * where);
idx_t ParMlevelKWayPartitioning(dctrl_t * dctrl, dgraph_t * graph, 
    idx_t *const * part);
void ParSetupKWayBalMultipliers(dctrl_t * dctrl, dgraph_t *graph);
void PrintPartInfo(ctrl_t * ctrl, const dgraph_t * graph, 
    const idx_t * const * where);
void PrintExtraTimers(dctrl_t * dctrl); 

/* kwayfm.c */
idx_t ParGreedy_KWayOptimize(dctrl_t * dctrl, dgraph_t *graph, idx_t niter,
    real_t ffactor, idx_t * const * const * refcon);
idx_t updateVertex(dctrl_t * dctrl, idx_t myid, 
    idx_t k, idx_t to, idx_t from, idx_t ewgt, dgraph_t * graph, 
    idx_t * rnbnd, idx_t * rnnbrpool, rpq_t * queue, idx_t done);

/* kwayrefine.c */
void ParRefineKWay(dctrl_t * dctrl, dgraph_t * orggraph, dgraph_t * graph); 
void ParComputeKWayPartitionParams(dctrl_t * dctrl, dgraph_t * graph,
    idx_t * const * gpwgts);
void ParProjectKWayPartition(dctrl_t * dctrl, dgraph_t * graph);

/* mem.c */
idx_t ** pimalloc(size_t n,const char * msg);
move_t * mset(const idx_t n, const idx_t x, move_t * ptr); 
void ParAllocateRefinementWorkSpace(dctrl_t * dctrl, const dgraph_t * graph);
void ParAllocateKWayPartitionMemory(ctrl_t *ctrl, dgraph_t *graph);
dctrl_t * AllocateDCtrl(idx_t nvtxs, idx_t nparts, idx_t nthreads);
void FreeDCtrl(dctrl_t ** dctrl, idx_t nthreads);
void ParCnbrpoolReset(dctrl_t * dctrl);

/* util.c */
idx_t irandInRange_r(const idx_t max, idx_t * seed);
void irandArrayPermute_r(const idx_t n, idx_t *p, const idx_t nshuffles, 
  const idx_t flag, idx_t * seed);
void BucketSortKeysInc_t(idx_t * const counts, idx_t n, idx_t max, idx_t *keys, 
         idx_t *tperm, idx_t *perm);

#endif

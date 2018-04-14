/**
 * @file struct.h
 * @brief Structures
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */


#ifndef STRUCT_H
#define STRUCT_H

typedef struct dctrl_t {
  double contractParPartTmr, refineParPartTmr,matchParPartTmr,matchLoopTmr,
    parPartTmr,initPartTmr,projectParPartTmr,convertTmr,ctrlDupTmr,cpuTmr,
    rmMatchTmr,shemMatchTmr,emptyQueueTmr,rpqTmr,bndTmr,projectWhereTmr,
    projectCKRTmr,refHdrTmr,updNeiTmr,makMovTmr,checkMovTmr,parBisTmr,
    renParTmr,cmapTmr,rmHeaderTmr,shemHeaderTmr,undoMovTmr,finMovTmr,
    pwgtsTmr,conLoopTmr,mmapTmr,ipTmr,projectEdTmr,projectIdTmr,lCopyTmr;
  ctrl_t * ctrl;
  idx_t * buffer1, * buffer2;
  real_t * buffer3;
  idx_t ** tptr_idx;
  real_t ** tptr_real;
  void ** tptr_void;
  idx_t freeCtrlParts;
  /* nbr stuff */
  cnbr_t ** nbrpool;
  idx_t * nnbrpool;
  idx_t * maxnnbrpool;
} dctrl_t;

typedef struct move_t {
  idx_t to,from,mvid;
} move_t;

typedef struct update_t {
  idx_t to,from,ewgt,adjid;
} update_t;

/*************************************************************************/
/*! This data structure holds a graph */
/*************************************************************************/
typedef struct dgraph_t {
  idx_t nvtxs, nedges, maxdeg, nbnd,gnvtxs;

  idx_t *tvwgt;         /* The sum of the vertex weights in the graph */
  real_t *invtvwgt;     /* The inverse of the sum of the vertex weights in the graph */

  idx_t ndist;
  idx_t dsize;
  idx_t dmask;
  idx_t dshift;

  idx_t * mynvtxs;
  idx_t * mymaxdeg;

  idx_t ** xadj;
  idx_t ** vwgt;
  idx_t ** adjncy;
  idx_t ** adjwgt;
  idx_t ** cmap;
  idx_t ** where;

  idx_t * mynbnd;
  idx_t ** bndind;
  idx_t ** bndptr;
  ckrinfo_t ** ckrinfo;

  idx_t ** rename;


  /* These are to keep track control if the corresponding fields correspond to
     application or library memory */
  int free_xadj, free_vwgt, free_vsize, free_adjncy, free_adjwgt;

  idx_t *label;

  /* Partition parameters */
  idx_t mincut, minvol;
  idx_t *pwgts;

  /* Bisection refinement parameters */
  idx_t *id, *ed;

  struct dgraph_t *coarser, *finer;
} dgraph_t;

#endif

/**
 * @file check.c
 * @brief This file contains various sanity checks intended for ASSERT
 * statements
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2012, Regents of the University of Minnesota
 * @version 1
 * @date 2012-07-02
 */




#ifndef MTMETIS_CHECK_C
#define MTMETIS_CHECK_C




#include "check.h"




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


int check_info(
    ucinfo_t const * const ucinfo,
    graph_t const * const graph,
    pid_t const * const * const where)
{
  vtx_t i,k,mynvtxs,gvtx,lvtx;
  adj_t j;
  wgt_t ted, tid, m;
  pid_t me, other, nnbrs, nbrid;
  adj_t const * xadj;
  vtx_t const * adjncy;
  wgt_t const * adjwgt;
  pid_t const * mywhere;
  nbrinfo_t const * myrinfo;
  adjinfo_t const * mynbrs;

  pid_t const nparts = ucinfo->nparts;

  tid_t const myid = omp_get_thread_num();

  wgt_t htable[nparts];

  mynvtxs = graph->mynvtxs[myid];
  xadj = graph->xadj[myid];
  adjncy = graph->adjncy[myid];
  adjwgt = graph->adjwgt[myid];
  mywhere = where[myid];
  for (i=0;i<mynvtxs;++i) {
    nnbrs = tid = ted = 0;
    wgt_set(htable,NULL_WGT,nparts);
    me = mywhere[i];
    gvtx = lvtx_to_gvtx(i,myid,graph->dist);
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k < mynvtxs) {
        lvtx = k;
        nbrid = myid;
      } else {
        lvtx = gvtx_to_lvtx(k,graph->dist);
        nbrid = gvtx_to_tid(k,graph->dist);
      }
      other = where[nbrid][lvtx];
      if (other != me) {
        if ((m = htable[other]) != NULL_WGT) {
          htable[other] += adjwgt[j];
        } else {
          htable[other] = adjwgt[j];
          ++nnbrs;
        }
        ted += adjwgt[j];
      } else {
        tid += adjwgt[j];
      }
    }
    myrinfo = ucinfo->nbrinfo + i;
    mynbrs = ucinfo->nbrpool + myrinfo->nbrstart;
    if (myrinfo->id != tid) {
      printf("[%"PF_TID_T"] Mismatch of tid_ternal degree of vertex %" \
          PF_VTX_T"(%"PF_VTX_T") in p%"PF_PID_T" (expected [%"PF_WGT_T \
          ":%"PF_WGT_T"], but got [%"PF_WGT_T":%"PF_WGT_T"])\n",myid,gvtx, \
          i,me,tid,ted,myrinfo->id,myrinfo->ed);
      printf("[%"PF_TID_T"] R-Neighbors = {",myid);
      for (i=0;i<nparts;++i) {
        if (htable[i] > 0) {
          printf("[p%"PF_PID_T",w%"PF_WGT_T"]",i,htable[i]);
        }
      }
      printf("}\n");
      printf("[%"PF_TID_T"] E-Neighbors = {",myid);
      for (i=0;i<myrinfo->nnbrs;++i) {
        printf("[p%"PF_PID_T",w%"PF_WGT_T"]",mynbrs[i].pid,mynbrs[i].ed);
      }
      printf("}\n");

      return 0;
    } else if (myrinfo->ed != ted) {
      printf("[%"PF_TID_T"] Mismatch of external degree of vertex %" \
          PF_VTX_T"(%"PF_VTX_T") in p%"PF_PID_T" (expected [%"PF_WGT_T":%" \
          PF_WGT_T"], but got [%"PF_WGT_T":%"PF_WGT_T"])\n",myid,gvtx,i,me, \
          tid,ted,myrinfo->id,myrinfo->ed);
      printf("[%"PF_TID_T"] R-Neighbors = {",myid);
      for (i=0;i<nparts;++i) {
        if (htable[i] > 0) {
          printf("[p%"PF_PID_T",w%"PF_WGT_T"]",i,htable[i]);
        }
      }
      printf("}\n");
      printf("[%"PF_TID_T"] E-Neighbors = {",myid);
      for (i=0;i<myrinfo->nnbrs;++i) {
        printf("[p%"PF_PID_T",w%"PF_WGT_T"]",mynbrs[i].pid,mynbrs[i].ed);
      }
      printf("}\n");

      return 0;
    } else if (myrinfo->nnbrs != nnbrs) {
      printf("[%"PF_TID_T"] Mismatch of number of neighbors of vertex %" \
          PF_VTX_T"(%"PF_VTX_T") in %"PF_PID_T" (expected %"PF_PID_T", " \
          "but got %"PF_PID_T")\n",myid,gvtx,i,me,nnbrs,myrinfo->nnbrs);
      printf("[%"PF_TID_T"] R-Neighbors = {",myid);
      for (i=0;i<nparts;++i) {
        if (htable[i] > 0) {
          printf("[p%"PF_PID_T",w%"PF_WGT_T"]",i,htable[i]);
        }
      }
      printf("}\n");
      printf("[%"PF_TID_T"] E-Neighbors = {",myid);
      for (i=0;i<myrinfo->nnbrs;++i) {
        printf("[p%"PF_PID_T",w%"PF_WGT_T"]",mynbrs[i].pid,mynbrs[i].ed);
      }
      printf("}\n");
      return 0;
    } else {
      for (j=0;j<nnbrs;++j) {
        other = mynbrs[j].pid;
        if (htable[other] == NULL_WGT) {
          printf("[%"PF_TID_T"] Vertex %"PF_VTX_T"(%"PF_VTX_T")[p%"PF_PID_T \
              "] thinks its connected to %"PF_PID_T"/%"PF_WGT_T"",myid, \
              gvtx,i,me,other,htable[other]);
          printf(" : {");
          for (k=0;k<nparts;++k) {
            if (htable[k] != NULL_WGT) {
              printf("[p%"PF_PID_T":w%"PF_WGT_T"]",k,htable[k]);
            }
          }
          printf("}\n");
          return 0;
        } else if (htable[other] != mynbrs[j].ed) {
          printf("[%"PF_TID_T"] Mismatch of neighbor p%"PF_PID_T" weight " \
              "of vertex %"PF_VTX_T"(%"PF_VTX_T") of p%"PF_PID_T" " \
              "(expected %"PF_WGT_T", but got %"PF_WGT_T")\n",myid,other, \
              gvtx,i,me,htable[other],mynbrs[j].ed);
          return 0;
        }
      }
    }
  }
  return 1;
}


int check_graph(
    graph_t const * const graph)
{
  adj_t * xadj, * xudj;
  vtx_t * adjncy, * udjncy;
  wgt_t * adjwgt, * udjwgt;

  vtx_t mynvtxs, i, v, u, k, m;
  adj_t j, l;
  tid_t nbrid;

  tid_t const myid = omp_get_thread_num();

  mynvtxs = graph->mynvtxs[myid];
  xadj = graph->xadj[myid];
  adjncy = graph->adjncy[myid];
  adjwgt = graph->adjwgt[myid];
  for (i=0;i<mynvtxs;++i) {
    v = lvtx_to_gvtx(i,myid,graph->dist);
    for (j=xadj[i];j<xadj[i+1];++j) {
      u = adjncy[j];
      if (u < mynvtxs) {
        k = u;
        nbrid = myid;
      } else {
        k = gvtx_to_lvtx(u,graph->dist);
        nbrid = gvtx_to_tid(u,graph->dist);
        if (nbrid == myid) {
          printf("Local vertex is stored as remote\n");
          return 0;
        }
      }
      xudj = graph->xadj[nbrid];
      udjncy = graph->adjncy[nbrid];
      udjwgt = graph->adjwgt[nbrid];
      /* find the reverse edge */
      for (l=xudj[k];l<xudj[k+1];++l) {
        m = udjncy[l];
        if (m < graph->mynvtxs[nbrid]) {
          m = lvtx_to_gvtx(m,nbrid,graph->dist);
        }
        if (m == v) {
          if (udjwgt[l] != adjwgt[j]) {
            printf("[%"PF_TID_T"] Adjwgt of edge {%"PF_VTX_T"/%"PF_VTX_T":%" \
                PF_TID_T",%"PF_VTX_T"/%"PF_VTX_T":%"PF_TID_T"} is uneven (%" \
                PF_WGT_T":%"PF_WGT_T")\n",myid,i,mynvtxs,myid,k, \
                graph->mynvtxs[nbrid],nbrid,adjwgt[j],udjwgt[l]);
            return 0;
          } else {
            break;
          }
        }
      }
      if (l == xudj[k+1]) {
        printf("[%"PF_TID_T"] Edge {%"PF_VTX_T"/%"PF_VTX_T":%"PF_TID_T",%" \
            PF_VTX_T"/%"PF_VTX_T":%"PF_TID_T"} is only in one direction\n", \
            myid,i,mynvtxs,myid,k,graph->mynvtxs[nbrid],nbrid);
        return 0;
      }
    }
  }

  return 1;
}


int check_bnd(
    vtx_iset_t const * const bnd,
    graph_t const * const graph)
{
  vtx_t mynvtxs, i, lvtx, gvtx, k;
  adj_t j;
  pid_t nbrid, me, other;
  wgt_t tid, ted;
  adj_t * xadj;
  vtx_t * adjncy;
  wgt_t * adjwgt;

  pid_t const * const * const gwhere = (pid_t const **)graph->where;
  tid_t const myid = omp_get_thread_num();
  
  xadj = graph->xadj[myid];
  adjncy = graph->adjncy[myid];
  adjwgt = graph->adjwgt[myid];
  mynvtxs = graph->mynvtxs[myid];
  for (i=0;i<mynvtxs;++i) {
    tid = ted =0;
    gvtx = lvtx_to_gvtx(i,myid,graph->dist);
    me = gwhere[myid][i];
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (k <mynvtxs) {
        lvtx = k;
        nbrid = myid;
      } else {
        lvtx = gvtx_to_lvtx(k,graph->dist);
        nbrid = gvtx_to_tid(k,graph->dist);
      }
      other = gwhere[nbrid][lvtx];
      if (me != other) {
        ted += adjwgt[j];
      } else {
        tid += adjwgt[j];
      }
    }
    if (ted >= tid) {
      if (!vtx_iset_contains(i,bnd)) {
        printf("[%"PF_TID_T"] vertex %"PF_VTX_T"(%"PF_VTX_T") should be " \
            "on the border [%"PF_WGT_T":%"PF_WGT_T"]\n",myid,gvtx,i,tid,ted);
        return 0;
      }
    } else if (vtx_iset_contains(i,bnd)) {
      printf("[%"PF_TID_T"] vertex %"PF_VTX_T"(%"PF_VTX_T") should not be " \
          "on the border [%"PF_WGT_T":%"PF_WGT_T"]\n",myid,gvtx,i,tid,ted);
      return 0;
    }
  }

  return 1;
}




#endif

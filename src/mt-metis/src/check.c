/**
 * @file check.c
 * @brief This file contains various sanity checks intended for ASSERT
 * statements
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2012, Regents of the University of Minnesota
 * @version 1
 * @date 2012-07-02
 */

#include "includes.h"

idx_t checkInfo(dctrl_t * dctrl, const dgraph_t * graph, 
    const idx_t *const* where) 
{
  idx_t myid,i,j,k,m,ted,tid,mynvtxs,gvtx,lvtx,nbrid,nnbrs,me,other;
  idx_t * xadj, *adjncy, *adjwgt;

  const idx_t ndist = graph->ndist;
  const idx_t nparts = dctrl->ctrl->nparts;
  const idx_t * mywhere;

  ckrinfo_t * myrinfo;
  cnbr_t * mynbrs;

  idx_t htable[nparts];
  for (myid=0;myid<ndist;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    xadj = graph->xadj[myid];
    adjncy = graph->adjncy[myid];
    adjwgt = graph->adjwgt[myid];
    mywhere = where[myid];
    for (i=0;i<mynvtxs;++i) {
      nnbrs = tid = ted = 0;
      iset(nparts,-1,htable);
      me = mywhere[i];
      gvtx = LVTX_2_GVTX(i,myid,graph->dshift);
      for (j=xadj[i];j<xadj[i+1];++j) {
        k = adjncy[j];
        lvtx = GVTX_2_LVTX(k,graph->dshift);
        nbrid = GVTX_2_THRID(k,graph->dmask);
        other = where[nbrid][lvtx];
        if (other != me) {
          if ((m = htable[other]) != -1) {
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
      myrinfo = graph->ckrinfo[myid] + i;
      mynbrs = dctrl->nbrpool[myid] + myrinfo->inbr;
      if (myrinfo->id != tid) {
        printf("[%"PRIDX"] Mismatch of internal degree of vertex %"PRIDX"(%"PRIDX") in p%"PRIDX" "
            "(expected [%"PRIDX":%"PRIDX"], but got [%"PRIDX":%"PRIDX"])\n",myid,gvtx,i,me,tid,ted,
            myrinfo->id,myrinfo->ed);
        printf("[%"PRIDX"] R-Neighbors = {",myid);
        for (i=0;i<nparts;++i) {
          if (htable[i] > 0) {
            printf("[p%"PRIDX",w%"PRIDX"]",i,htable[i]);
          }
        }
        printf("}\n");
        printf("[%"PRIDX"] E-Neighbors = {",myid);
        for (i=0;i<myrinfo->nnbrs;++i) {
          printf("[p%"PRIDX",w%"PRIDX"]",mynbrs[i].pid,mynbrs[i].ed);
        }
        printf("}\n");

        return 0;
      } else if (myrinfo->ed != ted) {
        printf("[%"PRIDX"] Mismatch of external degree of vertex %"PRIDX"(%"PRIDX") in p%"PRIDX" "
            "(expected [%"PRIDX":%"PRIDX"], but got [%"PRIDX":%"PRIDX"])\n",myid,gvtx,i,me,tid,ted,
            myrinfo->id,myrinfo->ed);
        printf("[%"PRIDX"] R-Neighbors = {",myid);
        for (i=0;i<nparts;++i) {
          if (htable[i] > 0) {
            printf("[p%"PRIDX",w%"PRIDX"]",i,htable[i]);
          }
        }
        printf("}\n");
        printf("[%"PRIDX"] E-Neighbors = {",myid);
        for (i=0;i<myrinfo->nnbrs;++i) {
          printf("[p%"PRIDX",w%"PRIDX"]",mynbrs[i].pid,mynbrs[i].ed);
        }
        printf("}\n");

        return 0;
      } else if (myrinfo->nnbrs != nnbrs) {
        printf("[%"PRIDX"] Mismatch of number of neighbors of vertex %"PRIDX"(%"PRIDX") in %"PRIDX""
            "(expected %"PRIDX", but got %"PRIDX")\n",myid,gvtx,i,me,nnbrs,myrinfo->nnbrs);
        printf("[%"PRIDX"] R-Neighbors = {",myid);
        for (i=0;i<nparts;++i) {
          if (htable[i] > 0) {
            printf("[p%"PRIDX",w%"PRIDX"]",i,htable[i]);
          }
        }
        printf("}\n");
        printf("[%"PRIDX"] E-Neighbors = {",myid);
        for (i=0;i<myrinfo->nnbrs;++i) {
          printf("[p%"PRIDX",w%"PRIDX"]",mynbrs[i].pid,mynbrs[i].ed);
        }
        printf("}\n");
        return 0;
      } else {
        for (j=0;j<nnbrs;++j) {
          other = mynbrs[j].pid;
          if (htable[other] == -1) {
            printf("[%"PRIDX"] Vertex %"PRIDX"(%"PRIDX")[p%"PRIDX"] thinks its connected to %"PRIDX"/%"PRIDX"",myid,
                gvtx,i,me,other,htable[other]);
            printf(" : {");
            for (k=0;k<nparts;++k) {
              if (htable[k] >= 0) {
                printf("[p%"PRIDX":w%"PRIDX"]",k,htable[k]);
              }
            }
            printf("}\n");
            return 0;
          } else if (htable[other] != mynbrs[j].ed) {
            printf("[%"PRIDX"] Mismatch of neighbor p%"PRIDX" weight of vertex %"PRIDX"(%"PRIDX") of p%"PRIDX""
                "(expected %"PRIDX", but got %"PRIDX")\n",myid,other,gvtx,i,me,htable[other],
                mynbrs[j].ed);
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

idx_t checkWhere(dctrl_t * dctrl, const dgraph_t * graph, 
    const idx_t * const * where) 
{
  idx_t myid,i,mynvtxs;
  const idx_t * mywhere;
  ctrl_t * ctrl = dctrl->ctrl;
  for (myid=0;myid<graph->ndist;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    mywhere = where[myid];
    for (i=0;i<mynvtxs;++i) {
      if (mywhere[i] < 0 || mywhere[i] >= ctrl->nparts) {
        printf("Vertex %"PRIDX" has thinks its in partition %"PRIDX"/%"PRIDX"\n",
            LVTX_2_GVTX(i,myid,graph->dshift),mywhere[i],
          ctrl->nparts);
        return 0;
      }
    }
  }
  return 1;
}

idx_t checkGraph(const dgraph_t * const graph)
{
  const idx_t nvtxs = graph->nvtxs;
  const idx_t nedges = graph->nedges;
  const idx_t rnthreads= graph->ndist;

  idx_t * xadj, * adjncy, * adjwgt;
  idx_t * xudj, * udjncy, * udjwgt;

  idx_t mynvtxs,myid,i,v,j,u,k,nbrid,l,m;

  for (myid=0;myid<graph->ndist;++myid) {
    mynvtxs = graph->mynvtxs[myid];
    xadj = graph->xadj[myid];
    adjncy = graph->adjncy[myid];
    adjwgt = graph->adjwgt[myid];
    for (i=0;i<mynvtxs;++i) {
      v = LVTX_2_GVTX(i,myid,graph->dshift);
      for (j=xadj[i];j<xadj[i+1];++j) {
        u = adjncy[j];
        k = GVTX_2_LVTX(u,graph->dshift);
        nbrid = GVTX_2_THRID(u,graph->dmask);
        ASSERT(LVTX_2_GVTX(k,nbrid,graph->dshift) == u);
        xudj = graph->xadj[nbrid];
        udjncy = graph->adjncy[nbrid];
        udjwgt = graph->adjwgt[nbrid];
        /* find the reverse edge */
        for (l=xudj[k];l<xudj[k+1];++l) {
          if (udjncy[l] == v) {
            if (udjwgt[l] != adjwgt[j]) {
              printf("Adjwgt of edge {%"PRIDX",%"PRIDX"} is uneven\n",v,u);
              return 0;
            } else {
              break;
            }
          }
        }
        if (l == xudj[k+1]) {
          printf("Edge {%"PRIDX",%"PRIDX"} is only in one direction\n",v,u);
          return 0;
        }
      }
    }
  }
  return 1;
}

idx_t checkBND(const dgraph_t * graph) 
{
  idx_t mynvtxs,myid,i,j,k,nbrid,lvtx,gvtx,tid,ted,me,other;
  idx_t * xadj, * adjncy, * adjwgt;

  const idx_t ndist = graph->ndist;
  
  const idx_t * const * const gwhere = (const idx_t **)graph->where;
  

  for (myid=0;myid<ndist;++myid) {
    xadj = graph->xadj[myid];
    adjncy = graph->adjncy[myid];
    adjwgt = graph->adjwgt[myid];
    mynvtxs = graph->mynvtxs[myid];
    for (i=0;i<mynvtxs;++i) {
      tid = ted =0;
      gvtx = LVTX_2_GVTX(i,myid,graph->dshift);
      me = gwhere[myid][i];
      for (j=xadj[i];j<xadj[i+1];++j) {
        lvtx = GVTX_2_LVTX(adjncy[j],graph->dshift);
        nbrid = GVTX_2_THRID(adjncy[j],graph->dmask);
        other = gwhere[nbrid][lvtx];
        if (me != other) {
          ted += adjwgt[j];
        } else {
          tid += adjwgt[j];
        }
      }
      if (ted >= tid) {
        if (graph->bndptr[myid][i] == -1) {
          printf("[%"PRIDX"] vertex %"PRIDX"(%"PRIDX") should be on the border [%"PRIDX":%"PRIDX"]\n",myid,
              gvtx,i,tid,ted);
          return 0;
        }
      } else if (graph->bndptr[myid][i] != -1) {
        printf("[%"PRIDX"] vertex %"PRIDX"(%"PRIDX") should not be on the border [%"PRIDX":%"PRIDX"]\n",myid,
            gvtx,i, tid,ted);
        return 0;
      }
    }
  }

  return 1;
}

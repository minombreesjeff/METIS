/**
 * @file analytics.c
 * @brief Functions for measuring and estimating graph characteristics
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-10-09
 */





#ifndef BOWSTRING_ANALYZE_C
#define BOWSTRING_ANALYZE_C




#include "analyze.h"
#include "graph.h"




/******************************************************************************
* DOMLIB IMPORTS **************************************************************
******************************************************************************/


#define DLPQ_KEY_T vtx_t
#define DLPQ_VAL_T vtx_t
#define DLPQ_PREFIX vv
#define DLPQ_STATIC
#include "dlpq_headers.h"
#undef DLPQ_STATIC
#undef DLPQ_PREFIX
#undef DLPQ_KEY_T
#undef DLPQ_VAL_T




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static void __countingsort_v(
    vtx_t const * const keys,
    vtx_t * const out, 
    vtx_t const min, 
    vtx_t const max, 
    size_t const n)
{
  size_t i;
  size_t const size = (size_t)(max - min)+2;
  size_t * counts = size_calloc(size);
  /* avoid having to do offsets in each iteration */
  for (i=0;i<n;++i) {
    ++counts[(ssize_t)keys[i]];
  }
  size_prefixsum_exc(counts,size);
  for (i=0;i<n;++i) {
    out[counts[(ssize_t)keys[i]]++] = (vtx_t)i;
  }
  dl_free(counts);
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void degree_distribution(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t ** const r_degrees, 
    vtx_t * const r_maxdegree)
{
  vtx_t i, maxdeg, deg;
  vtx_t * degrees;

  maxdeg = 0;

  /* first pass to find maximum degree */
  for (i=0;i<nvtxs;++i) {
    deg = xadj[i+1] - xadj[i];
    if (deg > maxdeg) {
      maxdeg = deg;
    }
  }

  degrees = vtx_calloc(maxdeg+1);

  for (i=0;i<nvtxs;++i) {
    deg = xadj[i+1] - xadj[i];
    ++degrees[deg];
  }

  *r_degrees = degrees;
  *r_maxdegree = maxdeg;
}


void nhop_degree_distribution(
    vtx_t const nvtxs, 
    adj_t const * xadj, 
    vtx_t const * const adjncy, 
    vtx_t const nhops, 
    vtx_t ** const r_degrees, 
    vtx_t * const r_maxdegree)
{
  vtx_t i, k, hop, maxdeg, ndeg;
  adj_t j;
  vtx_t * curdeg, * nextdeg, * degrees;

  if (nhops == 1) {
    /* maybe thi should be reverse */
    degree_distribution(nvtxs,xadj,r_degrees,r_maxdegree);
  } else {
    curdeg = vtx_alloc(nvtxs);
    nextdeg = vtx_alloc(nvtxs);

    /* set initial degrees */
    for (i=0;i<nvtxs;++i) {
      curdeg[i] = xadj[i+1] - xadj[i];
    }

    /* perform passes */
    maxdeg = 0;
    for (hop=1;hop<nhops;++hop) {
      for (i=0;i<nvtxs;++i) {
        ndeg = 0;
        for (j=xadj[i];j<xadj[i+1];++j) {
          k = adjncy[j];
          ndeg += curdeg[k]; 
        }
        nextdeg[i] = ndeg;
        if (ndeg > maxdeg) {
          maxdeg = ndeg;
        }
      }
      dl_swap(curdeg,nextdeg);
    }

    /* drop nextdeg */
    dl_free(nextdeg);

    /* allocate our histogram */
    degrees = vtx_calloc(maxdeg+1);

    /* histogram up our nhop degrees */
    for (i=0;i<nvtxs;++i) {
      ++degrees[curdeg[i]];
    }

    dl_free(curdeg);

    *r_degrees = degrees;
    *r_maxdegree = maxdeg;
  }
}


size_t count_triangles(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy,
    adj_t * radj)
{
  int free_radj;
  vtx_t i, k, kk, maxdeg, p;
  adj_t j,jj,r,ridx;
  size_t ntriangles;
  vtx_t * mark, * perm, * deg, * madjncy;

  adj_t const nedges = xadj[nvtxs];

  ntriangles = 0;

  if (radj) {
    free_radj = 0;
  } else {
    radj = adj_alloc(nedges);
    build_adjncy_index(nvtxs,xadj,adjncy,radj);
    free_radj = 1;
  }

  deg = vtx_alloc(nvtxs);
  mark = vtx_init_alloc(0,nvtxs);
  perm = vtx_init_alloc(0,nvtxs);
  madjncy = vtx_duplicate(adjncy,nedges);

  maxdeg = 0;
  for (i=0;i<nvtxs;++i) {
    deg[i] = xadj[i+1]-xadj[i];
    if (deg[i] > maxdeg) {
      maxdeg = deg[i];
    }
  }

  __countingsort_v(deg,perm,0,maxdeg,nvtxs);

  for (p=0;p<nvtxs;++p) {
    i = perm[p];
    /* mark all vertices adjacent to i */
    for (j=xadj[i];j<xadj[i]+deg[i];++j) {
      k = madjncy[j];
      mark[k] = 1;
    }
    for (j=xadj[i];j<xadj[i]+deg[i];++j) {
      k = madjncy[j];
      if (mark[k]) {
        for (jj=xadj[k];jj<xadj[k]+deg[k];++jj) {
          kk= madjncy[jj];
          if (mark[kk]) {
            ++ntriangles;
          }
        }
        mark[k] = 0;
      }
    }
    /* remove all edges pointing to i */
    for (j=xadj[i];j<xadj[i]+deg[i];++j) {
      k = madjncy[j];
      r = radj[j];
      ridx = xadj[k]+deg[k]-1;
      /* remove the edge */
      madjncy[r] = madjncy[ridx];
      /* update the radj vector */
      radj[r] = radj[ridx];
      /* update the remote radj vector */
      radj[radj[r]] = r;
      --deg[k];
    }
  }

  dl_free(madjncy);
  dl_free(mark);
  dl_free(deg);

  if (free_radj) {
    dl_free(radj);
  }

  return ntriangles;
}


void atomic_cycle_distribution(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    adj_t * radj, 
    size_t ** const r_cycles, 
    vtx_t * const r_maxcycle)
{
  int free_radj;
  vtx_t i, k, v, m, sq, nq, si, ni, curlen, maxlen; 
  adj_t j, l;
  vtx_t * q, * qi, * len;
  int * vvisited, * evisited, * ivisited;
  size_t * cycles;

  adj_t const nedges = xadj[nvtxs];

  q = vtx_alloc(nvtxs);
  qi = vtx_alloc(nvtxs);
  vvisited = int_calloc(nvtxs);
  evisited = int_calloc(nedges);
  ivisited = int_calloc(nvtxs);
  cycles = size_alloc(dl_min(nvtxs,nedges));
  len = vtx_alloc(nvtxs);

  if (radj) {
    free_radj = 0;
  } else {
    radj = adj_alloc(nedges);
    build_adjncy_index(nvtxs,xadj,adjncy,radj);
    free_radj = 1;
  }

  maxlen = 0;
  cycles[0] = 0;

  /* seed the queue */
  sq = nq = 0;
  v = vtx_rand(0,nvtxs);
  q[nq++] = v;
  vvisited[v] = 1;

  /* algorithm from Gashler and Martinez 2012 */
  while (sq < nq) {
    v = q[sq++];
    for (j=xadj[v];j<xadj[v+1];++j) {
      k = adjncy[j];
      if (vvisited[k]) {
        len[k] = 1;
        si = ni = 0;
        qi[ni++] = k;
        ivisited[k] = 1;
        while (si < ni) {
          i = qi[si++];
          for (l=xadj[i];l<xadj[i+1];++l) {
            m = adjncy[l];
            if (!ivisited[m] && evisited[l]) {
              len[m] = len[i]+1;
              qi[ni++] = m;
              ivisited[m] = 1;
              if (m == v) {
                curlen = len[m];
                /* zero out new cycle lengths */
                while (curlen > maxlen) {
                  cycles[++maxlen] = 0;
                }
                ++cycles[curlen];
                /* I might need to break here */
                si = ni;
                break;
              }
            }
          }
        }
        /* clear ivisited */
        if (ni < nvtxs/64) {
          for (i=0;i<ni;++i) {
            ivisited[qi[i]] = 0;
          }
        } else {
          int_set(ivisited,0,nvtxs);
        }
      } else {
        q[nq++] = k;
        vvisited[k] = 1;
      }
      evisited[j] = 1;
      evisited[radj[j]] = 1;
    }
  }

  /* hack to ignore length 2 cycles */
  cycles[2] = 0;
  
  if (r_maxcycle) {
    *r_maxcycle = maxlen;
  }

  if (r_cycles) {
    *r_cycles = cycles;
  }

  if (free_radj) {
    dl_free(radj);
  }
}


void star_distribution(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    vtx_t ** const r_stars, 
    vtx_t * const r_maxstar)
{
  vtx_t i, k, maxstar;
  vtx_t * nleafs, * stars;

  nleafs = vtx_calloc(nvtxs);
  
  /* count leaf vertices */
  maxstar = 0;
  for (i=0;i<nvtxs;++i) {
    if (xadj[i+1] - xadj[i] == 1) {
      k = adjncy[xadj[i]];
      ++nleafs[k];
      if (nleafs[k] > maxstar) {
        maxstar = nleafs[k];
      }
    }
  }

  *r_maxstar = maxstar;

  if (r_stars != NULL) {
    stars = vtx_calloc(maxstar+1);
    for (i=0;i<nvtxs;++i) {
      k = nleafs[i];
      if (k > 1) {
        ++stars[k];
      }
    }
    *r_stars = stars;
  }

  dl_free(nleafs);
}


void calc_domaindegree(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    vlbl_t const nparts, 
    vlbl_t const * const where, 
    wgt_t * const dd)
{
  vtx_t i, k;
  adj_t j;
  vlbl_t me;

  wgt_set(dd,0,nparts);

  for (i=0;i<nvtxs;++i) {
    me = where[i];
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (me != where[k]) {
        if (adjwgt) {
          dd[me] += adjwgt[j];
        } else {
          dd[me] += 1.0;
        }
      }
    }
  }
}


void calc_domainconn(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    vlbl_t const nparts, 
    vlbl_t const * const where, 
    int * const dc)
{
  vtx_t i, k;
  adj_t j;
  vlbl_t me;
  int * dd;

  dd = int_init_alloc(0,nparts*nparts);

  for (i=0;i<nvtxs;++i) {
    me = where[i];
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (me != where[k]) {
        dd[(me*nparts)+where[k]] = 1;
      }
    }
  }

  int_set(dc,0,nparts);

  for (me=0;me<nparts;++me) {
    dc[me] = int_sum(dd+(nparts*me),nparts);
  }

  dl_free(dd);
}


void calc_domaincomvol(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    vlbl_t const nparts, 
    vlbl_t const * const where, 
    wgt_t * const dcvo,
    wgt_t * const dcvi)
{
  vtx_t i;
  adj_t j;
  vlbl_t me, k;
  vtx_t * marker;

  marker = vtx_init_alloc(-1,nparts);
  wgt_set(dcvi,0,nparts);
  wgt_set(dcvo,0,nparts);

  for (i=0;i<nvtxs;++i) {
    me = where[i];
    marker[me] = i;
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = where[adjncy[j]];
      if (marker[k] != i) {
        marker[k] = i;
        if (vwgt) {
          dcvi[k] += vwgt[i];
          dcvo[me] += vwgt[i];
        } else { 
          ++dcvi[k];
          ++dcvo[me];
        }
      }
    }
  }

  dl_free(marker);
}


wgt_t calc_edgecut(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    vlbl_t const * const where)
{
  vtx_t i,k;
  adj_t j;
  wgt_t cut, partcut;

  cut = 0.0;

  for (i=0;i<nvtxs;++i) {
    partcut = 0.0;
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (where[i] != where[k]) {
        if (adjwgt) {
          partcut += adjwgt[j];
        } else {
          partcut += 1.0;
        }
      }
    }
    cut += (partcut/2.0);
  }

  return cut;
}


wgt_t calc_communicationvolume(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  vtx_t i;
  adj_t j;
  vlbl_t k;
  wgt_t comv;
  vtx_t * marker;

  marker = vtx_init_alloc(NULL_VTX,nparts);

  comv = 0;
  for (i=0;i<nvtxs;++i) {
    marker[where[i]] = i;
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = where[adjncy[j]];
      if (marker[k] != i) {
        marker[k] = i;
        if (vwgt) {
          comv += vwgt[i];
        } else {
          ++comv;
        }
      }
    }
  }

  dl_free(marker);

  return comv;
}


wgt_t calc_max_domaincomvol(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const vwgt, 
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  wgt_t max;
  vlbl_t i;
  wgt_t * dcvi = NULL, * dcvo = NULL;

  dcvi = wgt_alloc(nparts);
  dcvo = wgt_alloc(nparts);

  calc_domaincomvol(nvtxs,xadj,adjncy,vwgt,nparts,where,dcvi,dcvo);

  max = 0;
  for (i=0;i<nparts;++i) {
    if (max < dcvi[i]) {
      max = dcvi[i];
    }
  }
  for (i=0;i<nparts;++i) {
    if (max < dcvo[i]) {
      max = dcvo[i];
    }
  }

  dl_free(dcvo);
  dl_free(dcvi);

  return max;
}


wgt_t calc_max_domaindegree(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  wgt_t max;
  wgt_t * dd = NULL;

  dd = wgt_alloc(nparts);

  calc_domaindegree(nvtxs,xadj,adjncy,adjwgt,nparts,where,dd);

  max = wgt_max_value(dd,nparts);

  dl_free(dd);

  return max;
}


double calc_vertex_balance(
    vtx_t const nvtxs,
    adj_t const * const xadj,
    vtx_t const * const adjncy,
    wgt_t const * const vwgt,
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  vtx_t i;
  vlbl_t p;
  wgt_t max, total;
  wgt_t * pwgts;

  pwgts = wgt_calloc(nparts);

  if (vwgt) {
    total = 0;
    for (i=0;i<nvtxs;++i) {
      p = where[i];
      total += vwgt[i];
      pwgts[p] += vwgt[i];
    }
  } else {
    total = (wgt_t)nvtxs;
    for (i=0;i<nvtxs;++i) {
      p = where[i];
      pwgts[p] += 1.0;
    }
  }

  max = 0.0;
  for (p=0;p<nparts;++p) {
    if (max < pwgts[p]) {
      max = pwgts[p];
    }
  }

  dl_free(pwgts);

  return (max*nparts)/(double)total;
}


double calc_edge_balance(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  vtx_t i;
  adj_t j;
  vlbl_t p;
  wgt_t max, total, w;
  wgt_t * pwgts;

  pwgts = wgt_calloc(nparts);

  if (adjwgt) {
    for (i=0;i<nvtxs;++i) {
      p = where[i];
      w = 0;
      for (j=xadj[i];j<xadj[i+1];++j) {
        w += adjwgt[j];
      }
      pwgts[p] += w;
    }
  } else {
    for (i=0;i<nvtxs;++i) {
      p = where[i];
      w = (wgt_t)(xadj[i+1]-xadj[i]);
      pwgts[p] += w;
    }
  }

  max = 0.0;
  total = 0.0;
  for (p=0;p<nparts;++p) {
    w = pwgts[p];
    total += w;
    if (max < w) {
      max = w;
    }
  }

  dl_free(pwgts);

  return (max*nparts)/(double)total;

}


double calc_degree_balance(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  wgt_t max;
  wgt_t total;

  /* not the most efficient */
  max = calc_max_domaindegree(nvtxs,xadj,adjncy,adjwgt,nparts,where);
  total = calc_edgecut(nvtxs,xadj,adjncy,adjwgt,where);

  return (max*nparts)/(double)total;
}


double calc_comvol_balance(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  wgt_t max;
  wgt_t total;

  /* not the most efficient */
  max = calc_max_domaincomvol(nvtxs,xadj,adjncy,adjwgt,nparts,where);
  total = calc_communicationvolume(nvtxs,xadj,adjncy,adjwgt,nparts,where);

  return (max*nparts)/(double)total;
}


double calc_modularity(
    vtx_t const nvtxs, 
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    wgt_t const * const adjwgt, 
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  double q;
  vtx_t i, k;
  adj_t j;
  wgt_t w, myid, myed;
  vlbl_t me, other;
  wgt_t id[nparts], ed[nparts];
  long double td, m;

  for (other=0;other<nparts;++other) {
    id[other] = 0;
    ed[other] = 0;
  }

  m = 0;
  for (i=0;i<nvtxs;++i) {
    me = where[i];
    myid = 0;
    myed = 0;
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      other = where[k];
      if (adjwgt) {
        w = adjwgt[j];
      } else {
        w = 1.0;
      }
      if (me == other) {
        myid += w;
      } else {
        myed += w;
      }
    }
    id[other] += myid;
    ed[other] += myed;
    m += myid + myed;
  }

  q = 0;
  for (other=0;other<nparts;++other) {
    td = ed[other] + id[other];
    q += (1.0 / m) * (id[other] - ((td*td)/m));
  }
  
  return q;
}


double calc_partition_components(
    vtx_t const nvtxs,
    adj_t const * const xadj, 
    vtx_t const * const adjncy, 
    vlbl_t const nparts,
    vlbl_t const * const where)
{
  vlbl_t p, ncom;
  double avg;
  vtx_t * xnvtxs;
  adj_t ** xxadj;
  vtx_t ** xadjncy;

  xnvtxs = vtx_alloc(nparts);
  xxadj = r_adj_alloc(nparts);
  xadjncy = r_vtx_alloc(nparts);

  induce_subgraphs(nvtxs,xadj,adjncy,NULL,NULL,where,nparts,xnvtxs,xxadj, \
      xadjncy,NULL,NULL,NULL,NULL);

  avg = 0;
  for (p=0;p<nparts;++p) {
    label_components(xnvtxs[p],xxadj[p],xadjncy[p],NULL,&ncom);
    avg += ncom;
    dl_free(xxadj[p]);
    dl_free(xadjncy[p]);
  }

  dl_free(xnvtxs);
  dl_free(xxadj);
  dl_free(xadjncy);

  avg /= nparts;

  return avg;
}



#endif

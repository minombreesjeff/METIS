/**
 * @file sparsen.c
 * @brief Edge removal functions
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013-2014, Dominique LaSalle
 * @version 1
 * @date 2013-09-25
 */




#ifndef SPARSEN_C
#define SPARSEN_C




#include "base.h"
#include "sparsen.h"
#include "graph.h"
#include "tree.h"




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


static const adj_t MIN_PSEUDO_SHUFFLE = 256;
static const adj_t MIN_EDGE_SHUFFLE = 16;
static const adj_t MAX_PATH_SEARCH = 1024;
static const adj_t START_ADJ = (adj_t)-1;
static const adj_t MAX_BSEARCH_FACTOR = 4;
static const size_t NNBRHDS = 2;
static const vtx_t MAX_NBRHD_SIZE = 512;




/******************************************************************************
* PRIVATE DOMLIB IMPORTS ******************************************************
******************************************************************************/


#define DLDPQ_PREFIX va 
#define DLDPQ_KEY_T elbl_t 
#define DLDPQ_VAL_T vtx_t
#define DLDPQ_STATIC
#include "dldpq_headers.h"
#undef DLDPQ_STATIC
#undef DLDPQ_PREFIX
#undef DLDPQ_KEY_T
#undef DLDPQ_VAL_T


#define DLDJSET_PREFIX vtx
#define DLDJSET_TYPE_T vtx_t
#define DLDJSET_STATIC
#include "dldjset_headers.h"
#undef DLDJSET_STATIC
#undef DLDJSET_TYPE_T
#undef DLDJSET_PREFIX




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/

static void __countingsort_v(
    const wgt_t * const keys,
    const adj_t * const vals, 
    adj_t * const out, 
    const wgt_t min,
    const wgt_t max, 
    const size_t n)
{
  size_t i, j;
  const size_t size = (size_t)(max - min)+2;
  size_t * counts = size_calloc(size);
  /* avoid having to do offsets in each iteration */
  for (i=0;i<n;++i) {
    ++counts[(ssize_t)keys[i]];
  }
  size_prefixsum_exc(counts,size);
  for (i=0;i<n;++i) {
    j = (size_t)vals[i];
    out[counts[(ssize_t)keys[j]]++] = (adj_t)j;
  }
  dl_free(counts);
}


static int __reweight_graph(
    const vtx_t nvtxs, 
    const adj_t * const xadj,
    const vtx_t * const adjncy, 
    wgt_t * const adjwgt, 
    const adj_t nremoved, 
    const vtx_t * const aradj, 
    const vtx_t * const bradj, 
    const wgt_t * const rwgt)
{
  double pathlen, len;
  size_t npaths;
  vtx_t i, k, m, sq, nq;
  adj_t e, j, l;
  wgt_t w;
  vtx_t * q, * visited;
  adj_t * tree, * radj;

  /* delete me */
  size_t nsearched = 0;
  size_t esearched = 0;

  /* build reverse adjaceny index */
  radj = adj_alloc(xadj[nvtxs]);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  q = vtx_alloc(nvtxs);
  visited = vtx_calloc(nvtxs);
  tree = adj_alloc(xadj[nvtxs]);

  pathlen = 0;
  npaths = 0;

  for (e=0;e<nremoved;++e) {
    /* extract the edge from the removed list */
    i = aradj[e];
    k = bradj[e];
    if (rwgt) {
      w = rwgt[e];
    } else {
      w = 1;
    }
    /* find the path */
    sq = nq = 0;
    visited[i] = 1;
    q[nq++] = START_ADJ;
    while (sq<nq) {
      DL_ASSERT(sq < nq,"Emptied queue!\n");
      l = q[sq++];
      if (l == START_ADJ) {
        i = aradj[e];
      } else {
        i = adjncy[l];
      }
      for (j=xadj[i];j<xadj[i+1];++j) {
        ++esearched;
        tree[j] = l;
        m = adjncy[j];
        if (m == k) {
          /* add weight to the path */
          len = visited[i];
          while (j != START_ADJ) {
            adjwgt[j] += w;
            adjwgt[radj[j]] += w;
            j = tree[j];
          }
          pathlen += len;
          ++npaths;
          goto FOUND_END;
        }
        if (!visited[m]) {
          q[nq++] = j;
          visited[m] = visited[i]+1;
        }
      }
      if (nq > MAX_PATH_SEARCH) {
        goto FOUND_END;
      }
    }

    FOUND_END:

    /* delete me */
    nsearched += nq;
    
    /* clear visited */
    for (i=1;i<nq;++i) {
      DL_ASSERT(visited[adjncy[q[i]]] != 0,"Queued vertex not visited\n");
      visited[adjncy[q[i]]] = 0;
    }
    visited[aradj[e]] = 0;
  }

  dl_free(q);
  dl_free(visited);
  dl_free(tree);
  dl_free(radj);

  printf("Average path length %lf (%zu paths)\n",pathlen/npaths,npaths);

  /* delete me */
  printf("Searched %zu nodes (%lf per removed edge)\n",nsearched,
      nsearched/(double)nremoved);
  printf("Searched %zu edges (%lf per removed edge)\n",esearched,
      esearched/(double)nremoved);

  return BOWSTRING_SUCCESS;
}


static int __approx_reweight_graph(
    const vtx_t nvtxs, 
    const adj_t * const xadj,
    const vtx_t * const adjncy, 
    wgt_t * const adjwgt, 
    const adj_t nremoved, 
    const vtx_t * const aradj, 
    const vtx_t * const bradj, 
    const wgt_t * const rwgt)
{
  double pathlen;
  size_t npaths;
  vtx_t i, l, k, m, sq, nq, n, nm, cm;
  adj_t e, j, len;
  wgt_t w;
  vtx_t * fnbrhd, * q, * visited;
  adj_t * radj, * tree;
  vtx_t * nbrhd;

  /* delete me */
  size_t nsearched = 0;
  size_t esearched = 0;

  /* build reverse adjaceny index */
  radj = adj_alloc(xadj[nvtxs]);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  nbrhd = vtx_alloc(NNBRHDS*nvtxs);  
  fnbrhd = vtx_alloc(nvtxs);
  for (n=0;n<NNBRHDS;++n) {
    neighborhoodify(nvtxs,xadj,adjncy,MAX_NBRHD_SIZE,fnbrhd);
    for (i=0;i<nvtxs;++i) {
      nbrhd[i*NNBRHDS+n] = fnbrhd[i];
    }
  }
  dl_free(fnbrhd);

  q = vtx_alloc(nvtxs);
  visited = vtx_calloc(nvtxs);
  tree = adj_alloc(xadj[nvtxs]);

  pathlen = 0;
  npaths = 0;

  for (e=0;e<nremoved;++e) {
    /* extract the edge from the removed list */
    i = aradj[e];
    k = bradj[e];
    if (rwgt) {
      w = rwgt[e];
    } else {
      w = 1;
    }
    /* find the path */
    sq = nq = 0;
    visited[i] = 1;
    q[nq++] = START_ADJ;
    nm = 1;
    while (sq<nq) {
      l = q[sq++];
      if (l == START_ADJ) {
        i = aradj[e];
      } else {
        i = adjncy[l];
      }
      for (j=xadj[i];j<xadj[i+1];++j) {
        ++esearched;
        m = adjncy[j];
        if (!visited[m]) {
          cm = 0;
          for (n=0;n<NNBRHDS;++n) {
            if (nbrhd[m*NNBRHDS+n] == nbrhd[k*NNBRHDS+n]) {
              ++cm;
            }
            if (NNBRHDS-n + cm < nm) {
              break;
            }
          }
          if (cm > nm) {
            nm = cm;
          }
          /* only search in other nodes neighborhood */
          if (cm == nm) {
            tree[j] = l;
            if (m == k) {
              /* add weight to the path */
              len = visited[i];
              while (j != START_ADJ) {
                adjwgt[j] += w;
                adjwgt[radj[j]] += w;
                j = tree[j];
              }
              pathlen += len;
              ++npaths;
              goto FOUND_END;
            }
            visited[m] = visited[i]+1;
            q[nq++] = j;
          }
        }
      }
      if (nq > MAX_PATH_SEARCH) {
        goto FOUND_END;
      }
    }

    FOUND_END:
    
    /* delete me */
    nsearched += nq;
    
    /* clear visited */
    for (i=1;i<nq;++i) {
      DL_ASSERT(visited[adjncy[q[i]]] != 0,"Queued vertex not visited\n");
      visited[adjncy[q[i]]] = 0;
    }
    visited[aradj[e]] = 0;
  }

  dl_free(nbrhd);

  dl_free(q);
  dl_free(tree);
  dl_free(visited);
  dl_free(radj);

  printf("Average path length %lf (%zu paths)\n",pathlen/npaths,npaths);

  /* delete me */
  printf("Searched %zu nodes (%lf per removed edge)\n",nsearched,
      nsearched/(double)nremoved);
  printf("Searched %zu edges (%lf per removed edge)\n",esearched,
      esearched/(double)nremoved);


  return BOWSTRING_SUCCESS;
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


adj_t prune_ranked_edges(
    const vtx_t nvtxs, 
    const adj_t * const gxadj,
    const vtx_t * const gadjncy, 
    const wgt_t * const gadjwgt, 
    const elbl_t * const rank, 
    const elbl_t maxrank, 
    const double frac, 
    adj_t ** const r_xadj, 
    vtx_t ** const r_adjncy, 
    wgt_t ** const r_adjwgt,
    const int reweight)
{
  unsigned int seed;
  vtx_t i, k;
  adj_t j, knedges, rnedges, nradj, nr;
  elbl_t r, splitrank;
  adj_t * xadj, * prank;
  vtx_t * adjncy, * u, * v, * x = NULL, * y = NULL;
  wgt_t * adjwgt = NULL, * w = NULL, * z = NULL;

  const int do_adjwgt = gadjwgt && r_adjwgt;
  const adj_t nedges = gxadj[nvtxs];
  /* multiple of 2 rounded up */
  const adj_t nnedges = ((adj_t)(ceil(nedges/2)*frac))*2; 

  if (nedges == nnedges) {
    /* we don't need to remove any edges */
    xadj = adj_duplicate(gxadj,nvtxs+1);
    adjncy = vtx_duplicate(gadjncy,nedges);
    if (do_adjwgt) {
      adjwgt = wgt_duplicate(gadjwgt,nedges);
    }
  } else {
    /* we need to remove edges */
    xadj = adj_calloc(nvtxs+1); 
    adjncy = vtx_alloc(nnedges); 
    if (do_adjwgt) {
      adjwgt = wgt_alloc(nnedges);
    }

    /* sum the ranks */
    prank = adj_calloc(maxrank);
    for (j=0;j<nedges;++j) {
      ++prank[rank[j]];
    }

    /* find out what rank we split on */
    knedges = 0;
    splitrank = 0;
    rnedges = 0;
    for (r=0;r<maxrank;++r) {
      knedges += prank[r];
      if (knedges == nnedges) {
        splitrank = r+1;
        rnedges = 0;
        break;
      } else if (knedges > nnedges) {
        splitrank = r;
        /* number of edges of rank splitrank we need to remove */
        rnedges = prank[r] - (knedges - nnedges);
        break;
      }
    }

    DL_ASSERT(rnedges % 2 == 0,"Encountered an odd number of edges "PF_ADJ_T
        " ("PF_ADJ_T"/"PF_ADJ_T") of rank "PF_ELBL_T" to remove\n",rnedges,
        nnedges,nedges,splitrank);

    /* allocate array for removed edges */
    if (r_adjwgt) {
      nr = 0;
      x = vtx_alloc(nedges - nnedges);
      y = vtx_alloc(nedges - nnedges);
      if (gadjwgt) {
        z = wgt_alloc(nedges - nnedges);
      }
    }

    if (rnedges > 0) {
      /* create of list of my splitrank edges */
      v = vtx_alloc(prank[splitrank]);
      u = vtx_alloc(prank[splitrank]);
      if (do_adjwgt) {
        w = wgt_alloc(prank[splitrank]);
      }
      nradj = 0;
      for (i=0;i<nvtxs;++i) {
        for (j=gxadj[i];j<gxadj[i+1];++j) {
          if (rank[j] == splitrank) {
            k = gadjncy[j];
            if (i < k) {
              v[nradj] = i;
              u[nradj] = k;
              if (do_adjwgt) {
                w[nradj] = gadjwgt[j];
              }
              ++nradj;
            }
          } else if (rank[j] < splitrank) {
            /* build partial xadj */
            ++xadj[i];
          } else if (r_adjwgt) {
            k = gadjncy[j];
            if (i < k) {
              x[nr] = i;
              y[nr] = k;
              if (gadjwgt) {
                z[nr] = gadjwgt[j];
              }
              ++nr;
            }
          }
        }
      }

      DL_ASSERT(nradj*2 == prank[splitrank],"Bad nradj: nradj = "PF_ADJ_T", "
          "prank[splitrank] = "PF_ADJ_T"\n",nradj,prank[splitrank]);

      /* permute that list */
      if (nradj > MIN_PSEUDO_SHUFFLE) {
        seed = *(dl_get_rand());
        vtx_pseudo_shuffle_r(v,nradj/8,nradj,&seed);
        seed = *(dl_get_rand());
        vtx_pseudo_shuffle_r(u,nradj/8,nradj,&seed);
        if (do_adjwgt) {
          seed = *(dl_get_rand());
          wgt_pseudo_shuffle_r(w,nradj/8,nradj,&seed);
        }
      } else {
        seed = *(dl_get_rand());
        vtx_shuffle_r(v,nradj,&seed);
        seed = *(dl_get_rand());
        vtx_shuffle_r(u,nradj,&seed);
        if (do_adjwgt) {
          seed = *(dl_get_rand());
          wgt_shuffle_r(w,nradj,&seed);
        }
      }
      dl_set_rand(seed);

      /* build xadj */
      for (j=0;j<rnedges/2;++j) {
        ++xadj[v[j]];
        ++xadj[u[j]]; 
      }
      adj_prefixsum_exc(xadj,nvtxs);

      /* populate edges low rank edges */
      for (i=0;i<nvtxs;++i) {
        for (j=gxadj[i];j<gxadj[i+1];++j) {
          if (rank[j] < splitrank) {
            adjncy[xadj[i]] = gadjncy[j];
            if (do_adjwgt) {
              adjwgt[xadj[i]] = gadjwgt[j];
            }
            ++xadj[i];
          }
        }
      }
      
      /* populate splitrank edges */
      for (j=0;j<rnedges/2;++j) {
        i = v[j];
        k = u[j];
        adjncy[xadj[i]] = k;
        adjncy[xadj[k]] = i;
        if (do_adjwgt) {
          adjwgt[xadj[i]] = adjwgt[xadj[k]] = w[j];
        }
        ++xadj[i];
        ++xadj[k];
      }
      if (r_adjwgt) {
        for (;j<nradj;++j) {
          x[nr] = v[j];
          y[nr] = u[j];
          if (gadjwgt) {
            z[nr] = w[j];
          }
          ++nr;
        }
      }

      /* free intermediate data */
      if (do_adjwgt) {
        dl_free(w);
      }
      dl_free(v);
      dl_free(u);
    } else {
      /* build xadj */
      for (i=0;i<nvtxs;++i) {
        for (j=gxadj[i];j<gxadj[i+1];++j) {
          if (rank[j] < splitrank) {
            /* build partial xadj */
            ++xadj[i];
          }
        }
      }
      DL_ASSERT(adj_sum(xadj,nvtxs) == nnedges,"Incorrect number of edges "
          "counted, found "PF_ADJ_T" but expecting "PF_ADJ_T"\n",
          adj_sum(xadj,nvtxs),nnedges);
      adj_prefixsum_exc(xadj,nvtxs);

      /* populate edges low rank edges */
      for (i=0;i<nvtxs;++i) {
        for (j=gxadj[i];j<gxadj[i+1];++j) {
          if (rank[j] < splitrank) {
            adjncy[xadj[i]] = gadjncy[j];
            if (do_adjwgt) {
              adjwgt[xadj[i]] = gadjwgt[j];
            }
            ++xadj[i];
          } else if (r_adjwgt) {
            k = gadjncy[j];
            if (gxadj[i+1]-gxadj[i] < gxadj[k+1]-gxadj[k]) {
              x[nr] = i;
              y[nr] = k;
              if (gadjwgt) {
                z[nr] = gadjwgt[j]; 
              }
              ++nr;
            }
          }
        }
      }
    }
    dl_free(prank);

    /* shift the xadj */
    for (i=nvtxs;i>0;--i) {
      xadj[i] = xadj[i-1];
    }
    xadj[0] = 0;

    if (r_adjwgt && reweight != BOWSTRING_REWEIGHT_NONE) {
      /* re-weight the edges */
      if (!adjwgt) {
        /* we'll need edge weights to start with */
        adjwgt = wgt_init_alloc(1,xadj[nvtxs]); 
      }

      /* delete me */
      dl_timer_t tmr;
      dl_init_timer(&tmr);
      dl_start_timer(&tmr);

      switch (reweight) {
        case BOWSTRING_REWEIGHT_EXACT:
          __reweight_graph(nvtxs,xadj,adjncy,adjwgt,nr,x,y,w);
          break;
        case BOWSTRING_REWEIGHT_APPROX:
          __approx_reweight_graph(nvtxs,xadj,adjncy,adjwgt,nr,x,y,w);
          break;
        default:
          dl_error("Unknown reweighting type '%d'\n",reweight);
      }

      /* delete me */
      dl_stop_timer(&tmr);
      printf("Rewighting took %lf\n",dl_poll_timer(&tmr));

      dl_free(x);
      dl_free(y);
    }
  }

  DL_ASSERT(xadj[nvtxs] == nnedges, "Bad number of edges after prunning, "
      "expected "PF_ADJ_T" but left with "PF_ADJ_T"\n",nnedges,xadj[nvtxs]);

  /* output assignment */
  if (r_xadj) {
    *r_xadj = xadj;
  }
  if (r_adjncy) {
    *r_adjncy = adjncy;
  }
  if (r_adjwgt) {
    *r_adjwgt = adjwgt;
  }
  
  return nedges-nnedges; 
}


elbl_t build_nirank(
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt, 
    elbl_t * const rank)
{
  vtx_t i, k, nscanned;
  adj_t e, j, l, deg;
  adj_t * radj;
  elbl_t r, maxrank;
  int * es;
  adj_t * perm;
  va_dldpq_t * q;

  const adj_t nedges = xadj[nvtxs];

  radj = adj_alloc(nedges);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  /* allocate necessary structures */
  es = int_calloc(nedges);
  q = va_dldpq_create(0,nedges,0,nvtxs);
  va_dldpq_fill_min_perm(q);
  perm = adj_alloc(nvtxs);

  /* build the spanning forests */
  nscanned = 0;
  maxrank = 0;
  while (nscanned < nvtxs) {
    i = va_dldpq_remove_max(q);
    deg = xadj[i+1]-xadj[i];
    adj_incset(perm,xadj[i],1,deg);
    if (deg > MIN_EDGE_SHUFFLE) {
      adj_pseudo_shuffle(perm,deg/8,deg);
    } else {
      adj_shuffle(perm,deg);
    }
    for (e=0;e<deg;++e) {
      j = perm[e];
      if (!es[j]) {
        k = adjncy[j];
        l = radj[j];
        DL_ASSERT(es[radj[j]] == 0, "Unsymmetrically unscanned edge "PF_ADJ_T
            " with scanned other direction at "PF_ADJ_T"\n",j,radj[j]);
        r = rank[j] = rank[l] = va_dldpq_peek(k,q);
        es[j] = es[l] = 1;
        va_dldpq_inc(k,q);
        if (r+1 > maxrank) {
          maxrank = r+1;
        }
      }
    }
    ++nscanned;
  }

  /* free memory */
  dl_free(perm);
  dl_free(es);
  dl_free(radj);
  va_dldpq_free(q);
  
  return maxrank+1;
}


elbl_t build_mstrank(
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt, 
    elbl_t * const rank)
{
  vtx_t i, k, t, st, et;
  adj_t j, l, lnedges, niters;
  adj_t * radj, * perm;
  vtx_t * alist = NULL, * blist = NULL, * vrank = NULL;
  adj_t * aadj = NULL, * badj = NULL;
  elbl_t maxrank;
  vtx_djset_t * trees;

  const adj_t nedges = xadj[nvtxs];

  radj = adj_alloc(nedges);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  /* allocate necessary structures */
  perm = adj_alloc(nedges/2);

  /* create a random permutation of the edges */
  adj_incset(perm,0,1,nedges/2);
  if (nedges/2 < MIN_PSEUDO_SHUFFLE) {
    adj_shuffle(perm,nedges/2);
  } else {
    adj_pseudo_shuffle(perm,nedges/16,nedges/2);
  }

  /* set up adjancy arrays */
  lnedges = 0;
  alist = vtx_alloc(nedges/2);
  aadj = adj_alloc(nedges/2);
  blist = vtx_alloc(nedges/2);
  badj = adj_alloc(nedges/2);
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (i < k) {
        l = perm[lnedges++];
        alist[l] = i;
        blist[l] = k;
        /* add the forward edge */
        aadj[l] = j;
        badj[l] = radj[j];
      }
    }
  }
  dl_free(radj);
  dl_free(perm);

  /* I think I can do better than this for an upper bound */
  vrank = vtx_calloc(nvtxs);

  /* build the spanning forests */
  maxrank = 0;
  niters = 0;
  trees = vtx_djset_create(0,nedges);
  while (lnedges > 0) {
    l = --lnedges;
    i = alist[l];
    k = blist[l];
    st = 0;
    et = dl_min(vrank[i],vrank[k])+1;
    while (1) { 
      ++niters;
      t = (st + et) / 2;
      if (vtx_djset_find(xadj[i]+t,trees) != vtx_djset_find(xadj[k]+t,trees)) {
        if (t == 0 || vtx_djset_find(xadj[i]+t-1,trees) == 
            vtx_djset_find(xadj[k]+t-1,trees)) {
          break;
        }
        et = t;
      } else {
        if (t == et-1) { 
          t = et;
          break;
        }
        st = t;
      }
    }
    dl_storemax(vrank[i],t);
    dl_storemax(vrank[k],t);
    dl_storemax(maxrank,(elbl_t)(t+1));
    vtx_djset_join(xadj[i]+t,xadj[k]+t,trees);
    rank[aadj[l]] = t;
    rank[badj[l]] = t;
  }

  dl_free(vrank);
  dl_free(alist);
  dl_free(blist);
  dl_free(aadj);
  dl_free(badj);

  /* free my trees */
  vtx_djset_free(trees);

  return maxrank;
}


elbl_t build_astrank(
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt, 
    elbl_t * const rank)
{
  vtx_t i, k, t, st, et;
  adj_t j, l, lnedges, niters, maxniters;
  adj_t * radj, * perm;
  vtx_t * alist = NULL, * blist = NULL, * vrank = NULL;
  adj_t * aadj = NULL, * badj = NULL;
  elbl_t maxrank;
  vtx_djset_t * trees;

  dl_timer_t permtmr, msttmr;

  dl_init_timer(&permtmr);
  dl_init_timer(&msttmr);

  const adj_t nedges = xadj[nvtxs];

  dl_start_timer(&permtmr);

  radj = adj_alloc(nedges);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  /* allocate necessary structures */
  perm = adj_alloc(nedges/2);

  /* create a random permutation of the edges */
  adj_incset(perm,0,1,nedges/2);
  if (nedges/2 < MIN_PSEUDO_SHUFFLE) {
    adj_shuffle(perm,nedges/2);
  } else {
    adj_pseudo_shuffle(perm,nedges/16,nedges/2);
  }

  /* set up adjancy arrays */
  lnedges = 0;
  alist = vtx_alloc(nedges/2);
  aadj = adj_alloc(nedges/2);
  blist = vtx_alloc(nedges/2);
  badj = adj_alloc(nedges/2);
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (i < k) {
        l = perm[lnedges++];
        alist[l] = i;
        blist[l] = k;
        /* add the forward edge */
        aadj[l] = j;
        badj[l] = radj[j];
      }
    }
  }
  dl_free(radj);
  dl_free(perm);

  /* delete me */
  dl_stop_timer(&permtmr);
  printf("Made permutation in %lf s\n",dl_poll_timer(&permtmr));

  dl_start_timer(&msttmr);

  /* I think I can do better than this for an upper bound */
  vrank = vtx_calloc(nvtxs);

  /* build the spanning forests */
  maxrank = 0;
  maxniters = 0;
  niters = 0;
  trees = vtx_djset_create(0,nedges);
  while (lnedges > 0) {
    maxniters += MAX_BSEARCH_FACTOR;
    l = --lnedges;
    i = alist[l];
    k = blist[l];
    st = 0;
    et = dl_min(vrank[i],vrank[k])+1;
    while (1) { 
      if (++niters > maxniters) {
        t = et;
        break;
      } else {
        t = (st + et) / 2;
      }
      if (vtx_djset_find(xadj[i]+t,trees) != vtx_djset_find(xadj[k]+t,trees)) {
        if (t == 0 || vtx_djset_find(xadj[i]+t-1,trees) == 
            vtx_djset_find(xadj[k]+t-1,trees)) {
          break;
        }
        et = t;
      } else {
        if (t == et-1) { 
          t = et;
          break;
        }
        st = t;
      }
    }
    dl_storemax(vrank[i],t);
    dl_storemax(vrank[k],t);
    dl_storemax(maxrank,(elbl_t)(t+1));
    vtx_djset_join(xadj[i]+t,xadj[k]+t,trees);
    rank[aadj[l]] = t;
    rank[badj[l]] = t;
  }

  dl_free(vrank);
  dl_free(alist);
  dl_free(blist);
  dl_free(aadj);
  dl_free(badj);

  /* delete me */
  dl_stop_timer(&msttmr);
  printf("Built "PF_ELBL_T" msts in %lf s\n",maxrank,dl_poll_timer(&msttmr));

  /* free my trees */
  vtx_djset_free(trees);

  return maxrank;
}


elbl_t build_lstrank(
    const vtx_t nvtxs, 
    const adj_t * const xadj, 
    const vtx_t * const adjncy, 
    const wgt_t * const adjwgt, 
    elbl_t * const rank)
{
  vtx_t i, k, t, st, et;
  wgt_t maxwgt;
  adj_t j, l, lnedges, niters, maxniters;
  adj_t * radj, * perm, * rperm;
  vtx_t * alist = NULL, * blist = NULL, * vrank = NULL;
  adj_t * aadj = NULL, * badj = NULL;
  wgt_t * swgt;
  elbl_t maxrank;
  vtx_djset_t * trees;

  dl_timer_t permtmr, msttmr;

  if (!adjwgt) {
    return build_astrank(nvtxs,xadj,adjncy,adjwgt,rank);
  }

  dl_init_timer(&permtmr);
  dl_init_timer(&msttmr);

  const adj_t nedges = xadj[nvtxs];

  dl_start_timer(&permtmr);

  radj = adj_alloc(nedges);
  build_adjncy_index(nvtxs,xadj,adjncy,radj);

  /* allocate necessary structures */
  perm = adj_alloc(nedges/2);
  rperm = adj_alloc(nedges/2);
  swgt = wgt_alloc(nedges/2);

  /* create a random permutation of the edges */
  adj_incset(rperm,0,1,nedges/2);
  if (nedges/2 < MIN_PSEUDO_SHUFFLE) {
    adj_shuffle(rperm,nedges/2);
  } else {
    adj_pseudo_shuffle(rperm,nedges/16,nedges/2);
  }

  l=0;
  maxwgt = 0;
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (i < k) {
        swgt[l++] = adjwgt[j];
        if (adjwgt[j] > maxwgt) {
          maxwgt = adjwgt[j];
        }
      }
    }
  }

  DL_ASSERT(l == nedges/2,"Didn't scan the correct number of edges\n");

  __countingsort_v(swgt,rperm,perm,0,maxwgt,nedges/2);

  dl_free(swgt);
  l = 0;
  for (j=0;j<nedges/2;++j) {
    rperm[perm[j]] = l++;
  }
  dl_free(perm);
  perm = rperm;

  /* set up adjancy arrays */
  lnedges = 0;

  alist = vtx_alloc(nedges/2);
  aadj = adj_alloc(nedges/2);
  blist = vtx_alloc(nedges/2);
  badj = adj_alloc(nedges/2);
  for (i=0;i<nvtxs;++i) {
    for (j=xadj[i];j<xadj[i+1];++j) {
      k = adjncy[j];
      if (i < k) {
        l = perm[lnedges++];
        alist[l] = i;
        blist[l] = k;
        /* add the forward edge */
        aadj[l] = j;
        badj[l] = radj[j];
      }
    }
  }
  dl_free(radj);
  dl_free(perm);

  /* delete me */
  dl_stop_timer(&permtmr);
  printf("Made permutation in %lf s\n",dl_poll_timer(&permtmr));

  dl_start_timer(&msttmr);

  /* I think I can do better than this for an upper bound */
  vrank = vtx_calloc(nvtxs);

  /* build the spanning forests */
  maxrank = 0;
  maxniters = 0;
  niters = 0;
  trees = vtx_djset_create(0,nedges);
  while (lnedges > 0) {
    maxniters += MAX_BSEARCH_FACTOR;
    l = --lnedges;
    
    /* delete me */
    if (adjwgt[aadj[l]] < adjwgt[aadj[0]]) {
      printf("Edges aren't sorted right!\n");
    }
    i = alist[l];
    k = blist[l];
    st = 0;
    et = dl_min(vrank[i],vrank[k])+1;
    while (1) { 
      if (++niters > maxniters) {
        t = et;
        break;
      } else {
        t = (st + et) / 2;
      }
      if (vtx_djset_find(xadj[i]+t,trees) != vtx_djset_find(xadj[k]+t,trees)) {
        if (t == 0 || vtx_djset_find(xadj[i]+t-1,trees) == 
            vtx_djset_find(xadj[k]+t-1,trees)) {
          break;
        }
        et = t;
      } else {
        if (t == et-1) { 
          t = et;
          break;
        }
        st = t;
      }
    }
    dl_storemax(vrank[i],t);
    dl_storemax(vrank[k],t);
    dl_storemax(maxrank,(elbl_t)(t+1));
    vtx_djset_join(xadj[i]+t,xadj[k]+t,trees);
    rank[aadj[l]] = t;
    rank[badj[l]] = t;
  }

  dl_free(vrank);
  dl_free(alist);
  dl_free(blist);
  dl_free(aadj);
  dl_free(badj);

  /* delete me */
  dl_stop_timer(&msttmr);
  printf("Built "PF_ELBL_T" MSTs in %lf s\n",maxrank,dl_poll_timer(&msttmr));

  /* free my trees */
  vtx_djset_free(trees);

  return maxrank;
}




#undef vtx_djset_t
#undef vtx_djset_create 
#undef vtx_djset_find
#undef vtx_djset_join
#undef vtx_djset_free 
#undef r_vtx_djset_calloc




#endif

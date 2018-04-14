/*
 * metis.h 
 *
 * This file contains function prototypes and constant definitions 
 * for METIS
 *
 * Started 8/9/02
 * George
 *
 */

#ifndef METIS_H
#define METIS_H 1


/****************************************************************************
* A set of defines that can be modified by the user
*****************************************************************************/
/*--------------------------------------------------------------------------
 Specifies the width of the elementary data type that will hold information
 about vertices and their adjacency lists.

 Possible values:
   32 : Use 32 bit signed integers
   64 : Use 64 bit signed integers

 A width of 64 should be specified if the number of vertices or the total
 number of edges in the graph exceed the limits of a 32 bit signed integer
 i.e., 2^31-1.
 Proper use of 64 bit integers requires that the c99 standard datatypes
 int32_t and int64_t are supported by the compiler.
 GCC does provides these definitions in stdint.h, but it may require some
 modifications on other architectures.
--------------------------------------------------------------------------*/
#define IDXTYPEWIDTH 32


/*--------------------------------------------------------------------------
 Specifies if the __thread storage directive is available by the compiler
 to indicate thread local storage. This storage directive is available in
 most systems using gcc compiler but it may not be available in other
 systems.

 Possible values:
  0 : Not available and do not use thread local storage
  1 : It is available and the __thread modifier will be used
--------------------------------------------------------------------------*/
#define HAVE_THREADLOCALSTORAGE 0


/****************************************************************************
* In principle, nothing needs to be changed beyond this point, unless the
* int32_t and int64_t cannot be found in the normal places.
*****************************************************************************/


/* Uniform definitions for various compilers */
#if defined(_MSC_VER)
  #define COMPILER_MSC
#endif
#if defined(__ICC)
  #define COMPILER_ICC
#endif
#if defined(__GNUC__)
  #define COMPILER_GCC
#endif



#if defined(COMPILER_GCC)
  #include <stdint.h>
#endif

#if defined(COMPILER_MSC)
  #include <ctrdefs.h>
  #define __thread __declspec( thread )

  typedef __int32                 int32_t;
  typedef __int64                 int64_t;
  typedef unsigned __int32        uint32_t;
  typedef unsigned __int64        uint64_t;
#else
  #include <sys/types.h>
#endif



/*------------------------------------------------------------------------
* Undefine the following #define in order to use short int as the idxtype 
*-------------------------------------------------------------------------*/
#if IDXTYPEWIDTH == 32
  typedef int32_t idxtype;
#elif IDXTYPEWIDTH == 64
  typedef int64_t idxtype;
#else
  #error "Incorrect user-supplied value fo IDXTYPEWIDTH"
#endif





/*------------------------------------------------------------------------
* Function prototypes 
*-------------------------------------------------------------------------*/

#if !defined(__cdecl)
  #define __cdecl
#endif



#ifdef __cplusplus
extern "C" {
#endif

void __cdecl METIS_EstimateMemory(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, 
                   int *optype, int *nbytes);

void __cdecl METIS_PartGraphKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, 
                   int *edgecut, idxtype *part); 

void __cdecl METIS_WPartGraphKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, 
                   int *options, int *edgecut, idxtype *part); 

void __cdecl METIS_PartGraphVKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *vsize, int *wgtflag, int *numflag, int *nparts, int *options, 
                   int *volume, idxtype *part);

void __cdecl METIS_WPartGraphVKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *vsize, int *wgtflag, int *numflag, int *nparts, float *tpwgts, 
                   int *options, int *volume, idxtype *part);

int  __cdecl METIS_MeshToDualCount(int *ne, int *nn, idxtype *elmnts, idxtype *elms, int *etype,
                   int *numflag);

void __cdecl METIS_MeshToDual(int *ne, int *nn, idxtype *elmnts, idxtype *elms, int *etype, 
                   int *numflag,  idxtype *dxadj, idxtype *dadjncy);

int  __cdecl METIS_MixedMeshToDualCount(int *ne, int *nn, idxtype *elmnts, idxtype * elms, 
                   int *etype, int *numflag, idxtype *conmat, int custom);

void __cdecl METIS_MixedMeshToDual(int *ne, int *nn, idxtype *elmnts, idxtype *elms, 
                   int *etype, int *numflag,idxtype *dxadj, idxtype *dadjncy,idxtype *conmat,
                   int custom);

void __cdecl METIS_MeshToNodal(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, 
                   idxtype *dxadj, idxtype *dadjncy);

void __cdecl METIS_MixedMeshToNodal(int *ne, int *nn, idxtype *elmnts, idxtype *etype,
                   int *numflag, idxtype *dxadj, idxtype *dadjncy);

void __cdecl METIS_PartMeshNodal(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag,
                   int *nparts, int *edgecut, idxtype *epart, idxtype *npart);

void __cdecl METIS_PartMixedMeshNodal(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag,
                   int *nparts, int *edgecut, idxtype *epart, idxtype *npart);

void __cdecl METIS_PartMeshDual(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, 
                   int *nparts, int *edgecut, idxtype *epart, idxtype *npart, int wgtflag, 
                   idxtype * vwgt);

void __cdecl METIS_PartMixedMeshDual(int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag,  
                   int *nparts, int *edgecut, idxtype *epart, idxtype *npart, idxtype *conmat, 
                   int custom, int wgtflag, idxtype *vwgt);

void __cdecl METIS_mCPartGraphKway(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, 
                   idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, 
                   float *rubvec, int *options, int *edgecut, idxtype *part);

void __cdecl METIS_mCPartGraphRecursive(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy,
                   idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                   int *options, int *edgecut, idxtype *part);

void __cdecl METIS_mCHPartGraphRecursive(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy,
                   idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                   float *ubvec, int *options, int *edgecut, idxtype *part);

void __cdecl METIS_mCPartGraphRecursiveInternal(int *nvtxs, int *ncon, idxtype *xadj, 
                   idxtype *adjncy, float *nvwgt, idxtype *adjwgt, int *nparts, int *options, 
                   int *edgecut, idxtype *part);

void __cdecl METIS_mCHPartGraphRecursiveInternal(int *nvtxs, int *ncon, idxtype *xadj, 
                   idxtype *adjncy, float *nvwgt, idxtype *adjwgt, int *nparts, float *ubvec, 
                   int *options, int *edgecut, idxtype *part);

void __cdecl METIS_EdgeND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
                   idxtype *perm, idxtype *iperm);

void __cdecl METIS_NodeND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
                   idxtype *perm, idxtype *iperm);


void __cdecl METIS_NodeWND(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag,
                   int *options, idxtype *perm, idxtype *iperm);

void __cdecl METIS_PartGraphKway2(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, 
                   int *edgecut, idxtype *part);

void __cdecl METIS_WPartGraphKway2(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, 
                   int *options, int *edgecut, idxtype *part);

void __cdecl METIS_NodeNDP(int nvtxs, idxtype *xadj, idxtype *adjncy, int npes, int *options, 
                   idxtype *perm, idxtype *iperm, idxtype *sizes);

void __cdecl METIS_NodeComputeSeparator(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, int *options, int *sepsize, idxtype *part);

void __cdecl METIS_EdgeComputeSeparator(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, int *options, int *sepsize, idxtype *part);

void __cdecl METIS_PartGraphRecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, 
                   int *edgecut, idxtype *part);

void __cdecl METIS_WPartGraphRecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, float *tpwgts, 
                   int *options, int *edgecut, idxtype *part);

void __cdecl METIS_PartFillGraph(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, 
                   int *edgecut, idxtype *part);


#ifdef __cplusplus
}
#endif




/*------------------------------------------------------------------------
* Constant definitions 
*-------------------------------------------------------------------------*/

/* Matching Schemes */
#define MTYPE_RM		1
#define MTYPE_HEM		2
#define MTYPE_SHEM		3
#define MTYPE_SHEMKWAY		4
#define MTYPE_SHEBM_ONENORM	5
#define MTYPE_SHEBM_INFNORM	6
#define MTYPE_SBHEM_ONENORM	7
#define MTYPE_SBHEM_INFNORM	8

/* Initial partitioning schemes for PMETIS and ONMETIS */
#define ITYPE_GGPKL		1
#define ITYPE_GGPKLNODE		2
#define ITYPE_RANDOM		2

/* Refinement schemes for PMETIS */
#define RTYPE_FM		1

/* Initial partitioning schemes for KMETIS */
#define ITYPE_PMETIS		1

/* Refinement schemes for KMETIS */
#define RTYPE_KWAYRANDOM	1
#define RTYPE_KWAYGREEDY	2
#define RTYPE_KWAYRANDOM_MCONN	3

/* Refinement schemes for ONMETIS */
#define RTYPE_SEP2SIDED		1
#define RTYPE_SEP1SIDED		2

/* Initial Partitioning Schemes for McKMETIS */
#define ITYPE_McPMETIS		1   	/* Simple McPMETIS */
#define ITYPE_McHPMETIS		2	/* horizontally relaxed McPMETIS */


/* Debug Levels */
#define DBG_TIME	1		/* Perform timing analysis */
#define DBG_OUTPUT	2
#define DBG_COARSEN   	4		/* Show the coarsening progress */
#define DBG_REFINE	8		/* Show info on communication during folding */
#define DBG_IPART	16		/* Show info on initial partition */
#define DBG_MOVEINFO	32		/* Show info on communication during folding */
#define DBG_KWAYPINFO	64		/* Show info on communication during folding */
#define DBG_SEPINFO	128		/* Show info on communication during folding */


/* Metis's version number */
#define METIS_VER_MAJOR         5
#define METIS_VER_MINOR         0
#define METIS_VER_SUBMINOR      0



#endif  /* METIS_H */

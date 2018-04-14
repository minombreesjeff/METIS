/*
 * Copyright 2003, Regents of the University of Minnesota
 *
 * cepic.c
 *
 * This file contains the driving routine for contact/impact simulations
 * for EPIC meshes
 *
 * Started 4/12/03
 * George
 *
 * $Id: cepic.c,v 1.15 2003/05/03 16:10:48 karypis Exp $
 *
 */

#include <metisbin.h>

#define Flip_int32(type)  (((type >>24) & 0x000000ff) | \
                           ((type >> 8) & 0x0000ff00) | \
                           ((type << 8) & 0x00ff0000) | \
                           ((type <<24) & 0xff000000) )

#define Flip_int64(type)  (((type >>56) & 0x00000000000000ff) | \
                           ((type >>40) & 0x000000000000ff00) | \
                           ((type >>24) & 0x0000000000ff0000) | \
                           ((type >>8)  & 0x00000000ff000000) | \
                           ((type <<8)  & 0x000000ff00000000) | \
                           ((type <<24) & 0x0000ff0000000000) | \
                           ((type <<40) & 0x00ff000000000000) | \
                           ((type <<56) & 0xff00000000000000))



/*************************************************************************
* Let the game begin
**************************************************************************/
main(idxtype argc, char *argv[])
{
  idxtype i, j, k, istep, options[10], nn, ne, fstep, lstep, nparts, nboxes, u[3], dim, nsplit, flags=0, NSKIP=1;
  char filename[256];
  idxtype *mien, *mrng, *part, *sflag;
  double *mxyz, *bxyz;
  idxtype *xadj, *adjncy, *cntptr, *cntind;
  idxtype numflag = 0, wgtflag = 0, edgecut, etype=2;
  void *cinfo=NULL;
  FILE *fpin;
  long long int *ltmp;

  if (argc <= 6) {
    fprintf(stderr, "Usage: %s <nn> <ne> <fstep> <lstep> <nparts> [flags] [NSKIP]\n", argv[0]);
    exit(0);
  }

  nn     = atoi(argv[1]);
  ne     = atoi(argv[2]);
  fstep  = atoi(argv[3]);
  lstep  = atoi(argv[4]);
  nparts = atoi(argv[5]);
  if (argc > 6)
    flags = atoi(argv[6]);
  if (argc > 7)
    NSKIP = atoi(argv[7]);

  printf("\n\n------------------------------------------------------------------------------------------\n");
  printf("Reading nn: %d, ne: %d,  fstep: %d, lstep: %d, nparts: %d\n", nn, ne, fstep, lstep, nparts);

  mien = imalloc(4*ne, "main: mien");
  mxyz = dmalloc(3*nn, "main: mxyz");
  mrng = imalloc(4*ne, "main: mrng");
  bxyz = dmalloc(6*ne*4, "main: bxyz");

  part  = imalloc(nn, "main: part");
  sflag = imalloc(nn, "main: sflag");

  xadj   = imalloc(nn+1, "main: xadj");
  adjncy = imalloc(50*nn, "main: adjncy");


  /*========================================================================
   * Read the initial mesh and setup the graph and contact information
   *========================================================================*/
  sprintf(filename, "mien.%04d", fstep);
  fpin = GKfopen(filename, "rb", "main: mien");
  fread(mien, sizeof(int), 4*ne, fpin);
  for (i=0; i<4*ne; i++)
    mien[i] = Flip_int32(mien[i]);
  GKfclose(fpin);

  /*========================================================================
   * Create the nodal graph
   *========================================================================*/
  numflag = mien[idxamin(4*ne, mien)];
  METIS_MeshToNodal(&ne, &nn, mien, &etype, &numflag, xadj, adjncy);



  /*========================================================================
   * Get into the loop in which you go over the different configurations
   *========================================================================*/
  for (k=0, istep=fstep; istep<=lstep; istep++, k++) {
    sprintf(filename, "mxyz.%04d", istep);
    printf("Reading %s...............................................................\n", filename);
    fpin = GKfopen(filename, "rb", "main: mxyz");
    fread(mxyz, sizeof(double), 3*nn, fpin);
    for (i=0; i<3*nn; i++) {
      ltmp = (long long int *)(mxyz+i);
      *ltmp = Flip_int64(*ltmp);
    }
    GKfclose(fpin);

    sprintf(filename, "mrng.%04d", istep);
    fpin = GKfopen(filename, "rb", "main: mrng");
    fread(mrng, sizeof(int), 4*ne, fpin);
    for (i=0; i<4*ne; i++)
      mrng[i] = Flip_int32(mrng[i]);
    GKfclose(fpin);

    /* Determine which nodes are in the surface */
    iset(nn, 0, sflag);
    for (i=0; i<ne; i++) {
      if (mrng[4*i+0] > 0) { /* 1, 2, 3 */
        sflag[mien[4*i+0]-1] = 1;
        sflag[mien[4*i+1]-1] = 1;
        sflag[mien[4*i+2]-1] = 1;
      }
      if (mrng[4*i+1] > 0) { /* 1, 2, 4 */
        sflag[mien[4*i+0]-1] = 1;
        sflag[mien[4*i+1]-1] = 1;
        sflag[mien[4*i+3]-1] = 1;
      }
      if (mrng[4*i+2] > 0) { /* 2, 3, 4 */
        sflag[mien[4*i+1]-1] = 1;
        sflag[mien[4*i+2]-1] = 1;
        sflag[mien[4*i+3]-1] = 1;
      }
      if (mrng[4*i+3] > 0) { /* 1, 3, 4 */
        sflag[mien[4*i+0]-1] = 1;
        sflag[mien[4*i+2]-1] = 1;
        sflag[mien[4*i+3]-1] = 1;
      }
    }

    printf("Contact Nodes: %d of %d\n", isum(nn, sflag), nn);

    /* Compute/Update the partitioning */
    if (k%NSKIP == 0) {
      if (cinfo != NULL)
        METIS_FreeContactInfo(cinfo);

      options[0] = 0;
      cinfo = METIS_PartGraphForContact(&nn, xadj, adjncy, mxyz, sflag, &numflag, &nparts, 
                    options, &edgecut, part);
      for (i=0; i<nn; i++)
        part[i]--;
    }

    switch (flags) {
      case 1:
        if (cinfo != NULL)
          METIS_FreeContactInfo(cinfo);

        cinfo = METIS_SetupContact(&nn, mxyz, sflag, &nparts, part);
        break;
      case 2:
        if (cinfo != NULL)
          METIS_FreeContactInfo(cinfo);

        cinfo = METIS_SetupContact0(&nn, mxyz, sflag, &nparts, part);
        break;
      default:
        METIS_UpdateContactInfo(cinfo, &nn, mxyz, sflag);
    }


    /* Determine the bounding boxes of the surface elements */
    for (nsplit=0, nboxes=0, i=0; i<ne; i++) {
      if (mrng[4*i+0] > 0) { /* 1, 2, 3 */
        u[0] = mien[4*i+0]-1;
        u[1] = mien[4*i+1]-1;
        u[2] = mien[4*i+2]-1;
        bxyz[6*nboxes+0] = bxyz[6*nboxes+3] = mxyz[3*u[0]+0];
        bxyz[6*nboxes+1] = bxyz[6*nboxes+4] = mxyz[3*u[0]+1];
        bxyz[6*nboxes+2] = bxyz[6*nboxes+5] = mxyz[3*u[0]+2];
        for (j=1; j<3; j++) {
          for (dim=0; dim<3; dim++) {
            bxyz[6*nboxes+dim] = (bxyz[6*nboxes+dim] > mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+dim]);
            bxyz[6*nboxes+3+dim] = (bxyz[6*nboxes+3+dim] < mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+3+dim]);
          }
        }
        nboxes++;
        if (part[u[0]] != part[u[1]] || part[u[0]] != part[u[2]])
          nsplit++;
      }
      if (mrng[4*i+1] > 0) { /* 1, 2, 4 */
        u[0] = mien[4*i+0]-1;
        u[1] = mien[4*i+1]-1;
        u[2] = mien[4*i+3]-1;
        bxyz[6*nboxes+0] = bxyz[6*nboxes+3] = mxyz[3*u[0]+0];
        bxyz[6*nboxes+1] = bxyz[6*nboxes+4] = mxyz[3*u[0]+1];
        bxyz[6*nboxes+2] = bxyz[6*nboxes+5] = mxyz[3*u[0]+2];
        for (j=1; j<3; j++) {
          for (dim=0; dim<3; dim++) {
            bxyz[6*nboxes+dim] = (bxyz[6*nboxes+dim] > mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+dim]);
            bxyz[6*nboxes+3+dim] = (bxyz[6*nboxes+3+dim] < mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+3+dim]);
          }
        }
        nboxes++;
        if (part[u[0]] != part[u[1]] || part[u[0]] != part[u[2]])
          nsplit++;
      }
      if (mrng[4*i+2] > 0) { /* 2, 3, 4 */
        u[0] = mien[4*i+1]-1;
        u[1] = mien[4*i+2]-1;
        u[2] = mien[4*i+3]-1;
        bxyz[6*nboxes+0] = bxyz[6*nboxes+3] = mxyz[3*u[0]+0];
        bxyz[6*nboxes+1] = bxyz[6*nboxes+4] = mxyz[3*u[0]+1];
        bxyz[6*nboxes+2] = bxyz[6*nboxes+5] = mxyz[3*u[0]+2];
        for (j=1; j<3; j++) {
          for (dim=0; dim<3; dim++) {
            bxyz[6*nboxes+dim] = (bxyz[6*nboxes+dim] > mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+dim]);
            bxyz[6*nboxes+3+dim] = (bxyz[6*nboxes+3+dim] < mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+3+dim]);
          }
        }
        nboxes++;
        if (part[u[0]] != part[u[1]] || part[u[0]] != part[u[2]])
          nsplit++;
      }
      if (mrng[4*i+3] > 0) { /* 1, 3, 4 */
        u[0] = mien[4*i+0]-1;
        u[1] = mien[4*i+2]-1;
        u[2] = mien[4*i+3]-1;
        bxyz[6*nboxes+0] = bxyz[6*nboxes+3] = mxyz[3*u[0]+0];
        bxyz[6*nboxes+1] = bxyz[6*nboxes+4] = mxyz[3*u[0]+1];
        bxyz[6*nboxes+2] = bxyz[6*nboxes+5] = mxyz[3*u[0]+2];
        for (j=1; j<3; j++) {
          for (dim=0; dim<3; dim++) {
            bxyz[6*nboxes+dim] = (bxyz[6*nboxes+dim] > mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+dim]);
            bxyz[6*nboxes+3+dim] = (bxyz[6*nboxes+3+dim] < mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+3+dim]);
          }
        }
        nboxes++;
        if (part[u[0]] != part[u[1]] || part[u[0]] != part[u[2]])
          nsplit++;
      }
    }


    METIS_FindContacts(cinfo, &nboxes, bxyz, &nparts, &cntptr, &cntind);

    printf("Contacting Elements: %d Indices: %d Nsplit: %d\n", nboxes, cntptr[nboxes]-nboxes, nsplit);

    GKfree((void *)&cntptr, &cntind, LTERM);
  }

}  


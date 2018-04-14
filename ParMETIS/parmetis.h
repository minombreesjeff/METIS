/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * parmetis.h
 *
 * This file contains the prototypes for ParMetis
 *
 * $Id: parmetis.h,v 1.2 1997/07/18 00:32:11 karypis Exp $
 */

void PARKMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARRMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARUAMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARDAMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PAROMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARGKMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
void PARGRMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
void PARGMETIS(idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);

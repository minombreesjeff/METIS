/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * timing.c
 *
 * This file contains routines that deal with timing Metis
 *
 * Started 7/24/97
 * George
 *
 * $Id: timing.c 5993 2009-01-07 02:09:57Z karypis $
 *
 */

#include "metislib.h"


/*************************************************************************
* This function clears the timers
**************************************************************************/
void InitTimers(ctrl_t *ctrl)
{
  gk_clearcputimer(ctrl->TotalTmr);
  gk_clearcputimer(ctrl->InitPartTmr);
  gk_clearcputimer(ctrl->MatchTmr);
  gk_clearcputimer(ctrl->ContractTmr);
  gk_clearcputimer(ctrl->CoarsenTmr);
  gk_clearcputimer(ctrl->UncoarsenTmr);
  gk_clearcputimer(ctrl->RefTmr);
  gk_clearcputimer(ctrl->ProjectTmr);
  gk_clearcputimer(ctrl->SplitTmr);
  gk_clearcputimer(ctrl->SepTmr);
}



/*************************************************************************
* This function prints the various timers
**************************************************************************/
void PrintTimers(ctrl_t *ctrl)
{
  printf("\nTiming Information -------------------------------------------------");
  printf("\n Multilevel: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->TotalTmr));
  printf("\n     Coarsening: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->CoarsenTmr));
  printf("\n            Matching: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->MatchTmr));
  printf("\n            Contract: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->ContractTmr));
  printf("\n     Initial Partition: \t %7.3"PRREAL"", gk_getcputimer(ctrl->InitPartTmr));
  printf("\n   Construct Separator: \t %7.3"PRREAL"", gk_getcputimer(ctrl->SepTmr));
  printf("\n     Uncoarsening: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->UncoarsenTmr));
  printf("\n          Refinement: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->RefTmr));
  printf("\n          Projection: \t\t\t %7.3"PRREAL"", gk_getcputimer(ctrl->ProjectTmr));
  printf("\n     Splitting: \t\t %7.3"PRREAL"", gk_getcputimer(ctrl->SplitTmr));
  printf("\n********************************************************************\n");
}




/*
 * SUMMARY:      InterceptionStorage.c - Calculate interception storage
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate interception storage
 * DESCRIP-END.
 * FUNCTIONS:    InterceptionStorage()
 * COMMENTS:
 * $Id$     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  InterceptionStorage()
*****************************************************************************/
void InterceptionStorage(int NMax, int NAct, float *MaxInt, float *Fract,
			 float *Int, float *Precip, float *KE, float *Height, 
			 unsigned char Understory, float Dt)
{
  float Available;		/* Available storage */
  float Intercepted;		/* Amount of water intercepted during this 
				   timestep */
  int i;			/* counter */
  float OriginalPrecip;

  OriginalPrecip = *Precip;

  
  /* The precipitation is multiplied by the fractional coverage, since if the 
     vegetation covers only 10% of the grid cell, only 10% can be intercepted as a 
     maximum */
  for (i = 0; i < NAct; i++) {
    Available = MaxInt[i] - Int[i];
    if (Available > *Precip * Fract[i])
      Intercepted = (*Precip) * Fract[i];
    else
      Intercepted = Available;
    *Precip -= Intercepted;
    Int[i] += Intercepted;
  }

  /* Find kinetic energy of rainfall for use by the sediment model. */
  if(Understory) {
    
    /* Since the understory is assumed to cover the entire grid cell, all 
       KE is associated with leaf drainage, eq. 23, Morgan et al. (1998) */
    
    if(Height[1] > .14)
      *KE = *Precip * 1000. * (15.8*sqrt(Height[1]) - 5.87); 
    else
      *KE = 0.0;
  }
  else {
    /* If no understory, part of the rainfall reaches the ground as direct throughfall. */
    *KE = 1000.*(OriginalPrecip*(8.95+8.44*log10(OriginalPrecip*1000./Dt))*(1-Fract[0]) +
      *Precip * (15.8*sqrt(Height[0])-5.87));
  }

  /* WORK IN PROGESS */
  /* It should be checked whether the following statement can cause a "loss"
     of water.  When there is interception storage at timestep t, and snow
     will cover this vegetation layer at t = T+1, the amount of water in 
     storage will be lost */
/*   for (i = NAct; i < NMax; i++) */
/*     Int[i] = 0; */
}

/*
 * SUMMARY:      CalcTransmissivity.c - Calculate saturated conductivity
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculates the transmissivity through the saturated part of
 *               the soil profile
 * DESCRIP-END.
 * FUNCTIONS:    CalcTransmissivity()
 * COMMENTS:
 * $Id$     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "functions.h"

/*****************************************************************************
  Function name: CalcTransmissivity()

  Purpose      : Calculates the transmissivity through the saturated part of
                 the soil profile
                 
  Required     : 
    float SoilDepth  - Total soil depth in m
    float WaterTable - Depth of the water table below the soil surface in m
    float LateralKs  - Lateral hydraulic conductivity in m/s
    float KsExponent - Exponent that describes exponential decay of LateralKs
                       with depth below the soil surface

  Returns      : Transmissivity in m2/s

  Modifies     : NA

  Comments     :
    Source:
    Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
      hydrology-vegetation model for complex terrain, Water Resour. Res.,
      30(6), 1665-1679, 1994.

    Based on:
    Beven, K. J., On subsurface stormflow:  An analysis of response times,
      Hydrolog. Sci. J., 4, 505-521, 1982.

    The hydraulic conductivity is assumed exponentially with depth, based on
    material in Beven [1982].
*****************************************************************************/
float CalcTransmissivity(float SoilDepth, float WaterTable, float LateralKs,
			 float KsExponent)
{
  float Transmissivity;		/* Transmissivity (m^2/s) */

  if (fequal(KsExponent, 0.0))
    Transmissivity = LateralKs * (SoilDepth - WaterTable);
  else
    Transmissivity = (LateralKs / KsExponent) *
      (exp(-KsExponent * WaterTable) - exp(-KsExponent * SoilDepth));

  return Transmissivity;
}

/*
 * SUMMARY:      CalcSafetyFactor.c - Calculate the factor of safety
 * USAGE:        Part of MWM
 *
 * AUTHOR:       Laura Bowling and Colleen Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * DESCRIPTION:  Calculate the aerodynamic resistances
 * DESCRIP-END.
 * FUNCTIONS:    CalcSafetyFactor()
 * COMMENTS:
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "DHSVMerror.h"
#include "settings.h"
#include "constants.h"
#include "data.h"

float FindValue(STATSTABLE Stats);

/*****************************************************************************
  Function name: CalcSafetyFactor()

  Purpose      : Calculate the factor of safety for mass wasting failure.
 
                 
  Required     :


  Returns      : float, values between 0 and 1 for failure, 
                               > 1 for stable
                               -0.1 for unconditionally unstable
                               and -999 for not in basin/invalid slope.

  Modifies     : none
   
  Comments     :
*****************************************************************************/

float CalcSafetyFactor(float Slope, int Soil, float SoilDepth, int Veg, SEDTABLE *SedType, VEGTABLE *VType, float M, SOILTABLE *SType)
{
  double root_cohes_kg, soil_cohes_kg, angle_int_frict_rad, slope_angle_rad, fc_soil_density;
  float RootCohesion, FrictionAngle, SoilCohesion, Surcharge;
  double term1;
  float loading;
  float safetyfactor;

  if (Slope >= 0. && Slope <= 90.) { 

    if(SoilDepth<=0.0) SoilDepth=0.001;

    M = M/SoilDepth;
    if(M>=1.0) M=0.99;

    /* Get stochastic parameter values. */
    /* Need to check for valid soil and vegetation types ! */
    RootCohesion = FindValue(VType[Veg - 1].RootCoh);
    FrictionAngle = FindValue(SedType[Soil - 1].Friction);
    SoilCohesion = FindValue(SedType[Soil - 1].Cohesion);
    Surcharge = FindValue(VType[Veg - 1].VegSurcharge);
   
    /* converting root cohesion from kPa to kg/m2 and angles from degrees to radians */
    root_cohes_kg = (RootCohesion * 1000.) / G;
    soil_cohes_kg = (SoilCohesion * 1000.) / G;
    angle_int_frict_rad = RADPDEG * FrictionAngle;
    slope_angle_rad = RADPDEG * Slope;
    
    fc_soil_density = SType[Soil - 1].Dens[0] + 
      (SType[Soil - 1].FCap[0] * WATER_DENSITY);
    
    loading = (Surcharge / (WATER_DENSITY * SoilDepth)) + 
      ((M * SedType[Soil - 1].SatDensity) / WATER_DENSITY) +
      ((1 - M) * (fc_soil_density / WATER_DENSITY));
    
	
    /* Check to see if slope is unconditionally unstable.*/
    term1 = ((SoilCohesion + root_cohes_kg) /
	     (Surcharge + SoilDepth * fc_soil_density) *
	     cos(slope_angle_rad)*cos(slope_angle_rad)) + (tan(angle_int_frict_rad));
		
    if(term1 <= tan(slope_angle_rad)) {
      /* Slope is unconditionally unstable for these values. */
      safetyfactor = -.1;
    }
    else {
      safetyfactor = (((2 * (SoilCohesion + root_cohes_kg)) / (WATER_DENSITY * SoilDepth * (sin(2 * slope_angle_rad)))) + ((loading - M) * ((tan(angle_int_frict_rad)) / (tan(slope_angle_rad))))) / loading; 
    }

  }    /* End of factor of safety calculation loop */
  else if (Slope <   0.0)
    safetyfactor = -999.;
  else  
    safetyfactor = -999.;
  
  return safetyfactor;
}

/*
 * SUMMARY:      MassEnergyBalance.c - Calculate mass and energy balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate mass and energy balance at each pixel
 * DESCRIP-END.
 * FUNCTIONS:    MassEnergyBalance()
 * COMMENTS:
 * $Id$     
 */

/* #define NO_ET */
/* #define NO_SNOW */

#ifdef SNOW_ONLY
#define NO_ET
#define NO_SOIL
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "massenergy.h"
#include "snow.h"
#include "constants.h"
#include "soilmoisture.h"

/*****************************************************************************
  MassEnergyBalance()
*****************************************************************************/
void MassEnergyBalance(int y, int x, float SineSolarAltitude, float DX, 
		       float DY, int Dt, int HeatFluxOption, 
		       int CanopyRadAttOption, int MaxVegLayers, 
		       PIXMET *LocalMet, ROADSTRUCT *LocalNetwork, 
		       PRECIPPIX *LocalPrecip, VEGTABLE *VType, 
		       VEGPIX *LocalVeg, SOILTABLE *SType,
		       SOILPIX *LocalSoil, SNOWPIX *LocalSnow,
		       EVAPPIX *LocalEvap, PIXRAD *TotalRad)
{
  PIXRAD LocalRad;		/* Radiation balance components (W/m^2) */
  float SurfaceWater;		/* Pixel average depth of water before
				   infiltration is calculated (m) */
  float Infiltration;		/* Infiltration into the top soil layer (m)
				 */
  float LowerRa;		/* Aerodynamic resistance for lower layer 
				   (s/m) */
  float LowerWind;		/* Wind for lower layer (m/s) */
  float MaxInfiltration;	/* Maximum infiltration into the top
				   soil layer (m) */
  float MaxRoadbedInfiltration;	/* Maximum infiltration through the road bed
				   soil layer (m) */
  float MeltEnergy;		/* Energy used to melt snow and  change of
				   cold content of snow pack */
  float MoistureFlux;		/* Amount of water transported from the pixel 
				   to the atmosphere (m/timestep) */
  float NetRadiation;		/* Net radiation for a layer (W/m2) */
  float Reference;		/* Reference height for sensible heat 
				   calculation (m) */
  float RoadbedInfiltration;	/* Infiltration through the road bed (m) */
  float Roughness;		/* Roughness length (m) */
  float Rp;			/* radiation flux in visible part of the 
				   spectrum (W/m^2) */
  float UpperRa;		/* Aerodynamic resistance for upper layer 
				   (s/m) */
  float UpperWind;		/* Wind for upper layer (m/s) */
  float SnowLongIn;		/* Incoming longwave radiation at snow surface 
				   (W/m2) */
  float SnowNetShort;		/* Net amount of short wave radiation at the 
				   snow surface (W/m2) */
  float SnowRa;			/* Aerodynamic resistance for snow */
  float SnowWind;		/* Wind 2 m above snow */
  float Tsurf;			/* Surface temperature used in
				   LongwaveBalance() (C) */
  int NVegLActual;		/* Number of vegetation layers above snow */

  /* Calculate the number of vegetation layers above the snow */

  NVegLActual = VType->NVegLayers;
  if (LocalSnow->HasSnow == TRUE && VType->UnderStory == TRUE)
    --NVegLActual;

  /* initialize the total amount of evapotranspiration, and MeltEnergy */

  LocalEvap->ETot = 0.0;
  MeltEnergy = 0.0;
  MoistureFlux = 0.0;

  /* calculate the radiation balance for the ground/snow surface and the
     vegetation layers above that surface */

  RadiationBalance(HeatFluxOption, CanopyRadAttOption, SineSolarAltitude, 
		   LocalMet->Sin, LocalMet->SinBeam, LocalMet->SinDiffuse, 
		   LocalMet->Lin, LocalMet->Tair, LocalVeg->Tcanopy, 
		   LocalSoil->TSurf, SType->Albedo, VType, LocalSnow, 
		   &LocalRad);

  /* calculate the actual aerodynamic resistances and wind speeds */

  UpperWind = VType->U[0] * LocalMet->Wind;
  UpperRa = VType->Ra[0] / LocalMet->Wind;
  if (VType->OverStory == TRUE) {
    LowerWind = VType->U[1] * LocalMet->Wind;
    LowerRa = VType->Ra[1] / LocalMet->Wind;
  }
  else {
    LowerWind = UpperWind;
    LowerRa = UpperRa;
  }

  /* calculate the amount of interception storage, and the amount of 
     throughfall.  Of course this only needs to be done if there is
     vegetation present. */

#ifndef NO_SNOW

  if (VType->OverStory == TRUE &&
      (LocalPrecip->IntSnow[0] || LocalPrecip->SnowFall > 0.0)) {
    SnowInterception(y, x, Dt, VType->Fract[0], VType->LAI[0],
		     VType->MaxInt[0], VType->MaxSnowInt, VType->MDRatio,
		     VType->SnowIntEff, UpperRa, LocalMet->AirDens,
		     LocalMet->Eact, LocalMet->Lv, &LocalRad, LocalMet->Press,
		     LocalMet->Tair, LocalMet->Vpd, UpperWind,
		     &(LocalPrecip->RainFall), &(LocalPrecip->SnowFall),
		     &(LocalPrecip->IntRain[0]), &(LocalPrecip->IntSnow[0]),
		     &(LocalPrecip->TempIntStorage),
		     &(LocalSnow->CanopyVaporMassFlux), &(LocalVeg->Tcanopy),
		     &MeltEnergy);
    MoistureFlux -= LocalSnow->CanopyVaporMassFlux;

    /* Because we now have a new estimate of the canopy temperature we can
       recalculate the longwave balance */
    if (LocalSnow->HasSnow == TRUE)
      Tsurf = LocalSnow->TSurf;
    else if (HeatFluxOption == TRUE)
      Tsurf = LocalSoil->TSurf;
    else
      Tsurf = LocalMet->Tair;
    LongwaveBalance(VType->OverStory, VType->Fract[0], LocalMet->Lin,
		    LocalVeg->Tcanopy, Tsurf, &LocalRad);
  }
  else if (VType->NVegLayers > 0) {
    LocalVeg->Tcanopy = LocalMet->Tair;
    LocalSnow->CanopyVaporMassFlux = 0.0;
    LocalPrecip->TempIntStorage = 0.0;
    InterceptionStorage(VType->NVegLayers, NVegLActual, VType->MaxInt,
			VType->Fract, LocalPrecip->IntRain,
			&(LocalPrecip->RainFall));
  }

  /* if snow is present, simulate the snow pack dynamics */

  if (LocalSnow->HasSnow || LocalPrecip->SnowFall > 0.0) {

    if (VType->OverStory == TRUE) {
      SnowLongIn = LocalRad.LongIn[1];
      SnowNetShort = LocalRad.NetShort[1];
    }
    else {
      SnowLongIn = LocalRad.LongIn[0];
      SnowNetShort = LocalRad.NetShort[0];
    }

    SnowWind = VType->USnow * LocalMet->Wind;
    SnowRa = VType->RaSnow / LocalMet->Wind;

    LocalSnow->Outflow =
      SnowMelt(y, x, Dt, 2. + Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
	       LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
	       LocalMet->Press, LocalPrecip->RainFall, LocalPrecip->SnowFall,
	       LocalMet->Tair, LocalMet->Vpd, SnowWind,
	       &(LocalSnow->PackWater), &(LocalSnow->SurfWater),
	       &(LocalSnow->Swq), &(LocalSnow->VaporMassFlux),
	       &(LocalSnow->TPack), &(LocalSnow->TSurf), &MeltEnergy);

    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */

    LocalPrecip->RainFall = 0.0;
    MoistureFlux -= LocalSnow->VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
       can recalculate the longwave balance */

    Tsurf = LocalSnow->TSurf;
    LongwaveBalance(VType->OverStory, VType->Fract[0], LocalMet->Lin,
		    LocalVeg->Tcanopy, Tsurf, &LocalRad);
  }
  else {
    LocalSnow->Outflow = 0.0;
    LocalSnow->VaporMassFlux = 0.0;
  }

  /* Determine whether a snow pack is still present, or whether everything
     has melted */

  if (LocalSnow->Swq > 0.0)
    LocalSnow->HasSnow = TRUE;
  else
    LocalSnow->HasSnow = FALSE;

  /*do the glacier add */
  if (LocalSnow->Swq < 1.0 && VType->Index == GLACIER) {
    printf("resetting glacier swe of %f to 5.0 meters\n", LocalSnow->Swq);
    LocalSnow->Glacier += (5.0 - LocalSnow->Swq);
    LocalSnow->Swq = 5.0;
    LocalSnow->TPack = 0.0;
    LocalSnow->TSurf = 0.0;
  }
#endif

#ifndef NO_ET
  /* calculate the amount of evapotranspiration from each vegetation layer 
     above the ground/soil surface.  Also calculate the total amount of 
     evapotranspiration from the vegetation */

  if (VType->OverStory == TRUE) {
    Rp = VISFRACT * LocalRad.NetShort[0];
    NetRadiation =
      LocalRad.NetShort[0] +
      LocalRad.LongIn[0] - 2 * VType->Fract[0] * LocalRad.LongOut[0];
    EvapoTranspiration(0, Dt, LocalMet, NetRadiation, Rp, VType, SType,
		       MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[0]),
		       LocalEvap, LocalNetwork->Adjust, UpperRa);
    MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];

    if (LocalSnow->HasSnow != TRUE && VType->UnderStory == TRUE) {
      Rp = VISFRACT * LocalRad.NetShort[1];
      NetRadiation =
	LocalRad.NetShort[1] +
	LocalRad.LongIn[1] - VType->Fract[1] * LocalRad.LongOut[1];
      EvapoTranspiration(1, Dt, LocalMet, NetRadiation, Rp, VType, SType,
			 MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[1]),
			 LocalEvap, LocalNetwork->Adjust, LowerRa);
      MoistureFlux += LocalEvap->EAct[1] + LocalEvap->EInt[1];
    }
    else if (VType->UnderStory == TRUE) {
      LocalEvap->EAct[1] = 0.;
      LocalEvap->EInt[1] = 0.;
    }
  }				/* end if(VType->OverStory == TRUE) */
  else if (LocalSnow->HasSnow != TRUE && VType->UnderStory == TRUE) {
    Rp = VISFRACT * LocalRad.NetShort[0];
    NetRadiation =
      LocalRad.NetShort[0] +
      LocalRad.LongIn[0] - VType->Fract[0] * LocalRad.LongOut[0];
    EvapoTranspiration(0, Dt, LocalMet, NetRadiation, Rp, VType, SType,
		       MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[0]),
		       LocalEvap, LocalNetwork->Adjust, LowerRa);
    MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];
  }
  else if (VType->UnderStory == TRUE) {
    LocalEvap->EAct[0] = 0.;
    LocalEvap->EInt[0] = 0.;
  }

  /* Calculate soil evaporation from the upper soil layer if no snow is 
     present and there is no understory */

  if (LocalSnow->HasSnow != TRUE && VType->UnderStory != TRUE) {

    if (VType->OverStory == TRUE)
      NetRadiation =
	LocalRad.NetShort[1] + LocalRad.LongIn[1] - LocalRad.LongOut[1];
    else
      NetRadiation =
	LocalRad.NetShort[0] + LocalRad.LongIn[0] - LocalRad.LongOut[0];

    LocalEvap->EvapSoil =
      SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
		      LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
		      NetRadiation, LowerRa, MoistureFlux, SType->Porosity[0],
		      SType->Ks[0], SType->Press[0], SType->PoreDist[0],
		      VType->RootDepth[0], &(LocalSoil->Moist[0]),
		      LocalNetwork->Adjust[0]);
  }
  else
    LocalEvap->EvapSoil = 0.0;

  MoistureFlux += LocalEvap->EvapSoil;
  LocalEvap->ETot += LocalEvap->EvapSoil;

#endif

  /* add the water that was not intercepted to the upper soil layer */

#ifndef NO_SOIL

  LocalSoil->SurfaceWater = 0.0;
  SurfaceWater = LocalPrecip->RainFall + LocalSoil->Runoff + LocalSnow->Outflow;

  MaxInfiltration = (1 - VType->ImpervFrac) * LocalNetwork->PercArea[0] * 
    SType->MaxInfiltrationRate * Dt; 

  Infiltration = (1 - VType->ImpervFrac) * LocalNetwork->PercArea[0] * 
    SurfaceWater; 

  if (Infiltration > MaxInfiltration) {
    Infiltration = MaxInfiltration;
  }

  MaxRoadbedInfiltration = (1 - LocalNetwork->PercArea[0]) * 
    LocalNetwork->MaxInfiltrationRate * Dt; 

  RoadbedInfiltration = (1 - LocalNetwork->PercArea[0]) * 
    SurfaceWater; 

  if (RoadbedInfiltration > MaxRoadbedInfiltration) {
    RoadbedInfiltration = MaxRoadbedInfiltration;
  }

  LocalSoil->Runoff = SurfaceWater - Infiltration - RoadbedInfiltration;

  LocalSoil->SurfaceWater = LocalSoil->Runoff;

  /* Calculate unsaturated soil water movement, and adjust soil water 
     table depth */

  UnsaturatedFlow(Dt, DX, DY, Infiltration, RoadbedInfiltration,
		  LocalSoil->SatFlow, SType->NLayers,
		  LocalSoil->Depth, LocalNetwork->Area, VType->RootDepth,
		  SType->Ks, SType->PoreDist, SType->Porosity, SType->FCap,
		  LocalSoil->Perc, LocalNetwork->PercArea,
		  LocalNetwork->Adjust, LocalNetwork->CutBankZone,
		  LocalNetwork->BankHeight, &(LocalSoil->TableDepth),
		  &(LocalSoil->Runoff), LocalSoil->Moist);

  if (HeatFluxOption == TRUE) {
    if (LocalSnow->HasSnow == TRUE) {
      Reference = 2. + Z0_SNOW;
      Roughness = Z0_SNOW;
    }
    else {
      Reference = 2. + Z0_GROUND;
      Roughness = Z0_GROUND;
    }

    SensibleHeatFlux(y, x, Dt, LowerRa, Reference, 0.0f, Roughness,
		     LocalMet, LocalRad.PixelNetShort, LocalRad.PixelLongIn,
		     MoistureFlux, SType->NLayers, VType->RootDepth,
		     SType, MeltEnergy, LocalSoil);
    Tsurf = LocalSoil->TSurf;
    LongwaveBalance(VType->OverStory, VType->Fract[0], LocalMet->Lin,
		    LocalVeg->Tcanopy, Tsurf, &LocalRad);
  }
  else
    NoSensibleHeatFlux(Dt, LocalMet, MoistureFlux, LocalSoil);

#endif

  /* add the components of the radiation balance for the current pixel to 
     the total */
  AggregateRadiation(MaxVegLayers, VType->NVegLayers, &LocalRad, TotalRad);
}

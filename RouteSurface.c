/*
 * SUMMARY:      RouteSurface.c - Route surface flow
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Route surface flow
 * DESCRIP-END.
 * FUNCTIONS:    RouteSurface()
 * COMMENTS:
 * $Id$     
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"


/*****************************************************************************
  RouteSurface()

  If the watertable calculated in WaterTableDepth() was negative, then water is
  ponding on the surface.  At the moment no ponding of water is allowed in
  DHSVM, and the "excess" water is routed to the outlet one pixel per time step
  However, if the pixel contains an impervious fraction, then the surface water
  is immediately routed to the nearest downslope pixel that contains a channel.
  The net effect is that all pixels that have an impervious area are directly
  connected (over the coarse of a single time step) to the channel network, this
  assumption is likely to be true for small urban basins, and perhaps even for
  large rural basins with some urban development

*****************************************************************************/
void RouteSurface(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
		  SOILPIX ** SoilMap, OPTIONSTRUCT *Options,
		  UNITHYDR ** UnitHydrograph,
		  UNITHYDRINFO * HydrographInfo, float *Hydrograph,
		  DUMPSTRUCT *Dump, VEGPIX ** VegMap, VEGTABLE * VType, 
		  SOILTABLE *SType, CHANNEL *ChannelData, SEDPIX **SedMap,
		  PRECIPPIX **PrecipMap, SEDTABLE *SedType)
{
  const char *Routine = "RouteSurface";
  int Lag;			/* Lag time for hydrograph */
  int Step;
  float StreamFlow;
  int TravelTime;
  int WaveLength;
  TIMESTRUCT NextTime;
  TIMESTRUCT VariableTime;

  int i;
  int j;
  int x;
  int y;
  int n;
  float **Runon;

  /* Kinematic wave routing: */
  double slope, alpha, beta;
  double outflow;
  int k;
  int particle_size ;
  float VariableDT;
  float total_in, total_out;
  float **SedIn, SedOut;
  float DR, DS;
  float h, term1, term2, term3;
  float streampower, effectivepower;
  float soliddischarge, TC, floweff;
  float settling;

  if (Options->Sediment) {
    if ((SedIn = (float **) calloc(Map->NY, sizeof(float *))) == NULL) {
      ReportError((char *) Routine, 1);
    }

    for (y = 0; y < Map->NY; y++) {
      if ((SedIn[y] = (float *) calloc(Map->NX, sizeof(float))) == NULL) {
	ReportError((char *) Routine, 1);
      }
    }
  }

  if (Options->HasNetwork) {
    if ((Runon = (float **) calloc(Map->NY, sizeof(float *))) == NULL) {
      ReportError((char *) Routine, 1);
    }

    for (y = 0; y < Map->NY; y++) {
      if ((Runon[y] = (float *) calloc(Map->NX, sizeof(float))) == NULL) {
	ReportError((char *) Routine, 1);
      }
    }
  
    // IExcess in meters/time step
    if(!Options->Routing) {
 
      if(Options->Sediment) {
	fprintf(stderr, "The sediment model must be run with kinematic wave routing.\n");
	exit(0);
      }

      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    SoilMap[y][x].Runoff = SoilMap[y][x].IExcess;
	    SoilMap[y][x].IExcess = 0.0;
	  }
	}
      }
  
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    if (VType[VegMap[y][x].Veg - 1].ImpervFrac > 0.0) {
	      SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].IExcess +=
		SoilMap[y][x].Runoff;
	    }
	    else {
	      for (n = 0; n < NDIRS; n++) {
		int xn = x + xneighbor[n];
		int yn = y + yneighbor[n];
		if (valid_cell(Map, xn, yn)) {
		  SoilMap[yn][xn].IExcess +=
		    SoilMap[y][x].Runoff * ((float) TopoMap[y][x].Dir[n] /
						  (float) TopoMap[y][x].TotalDir);
		}
	      }
	    }
	  }
	}
      }
    }

    /* Begin code for kinematic wave routing. */    
    else { 

      NextTime = *Time;
      VariableTime = *Time;

      /* Holds the value of the next DHSVM time step. */ 
      IncreaseTime(&NextTime);

 
      /* Use the Courant condition to find the maximum stable time step (in seconds). */
      /* Must be an even increment of Dt. */

      VariableDT = FindDT(SoilMap, Map, Time, TopoMap, SedType);
      
      //      fprintf(stderr, "VariableDT = %f\n", VariableDT);

      for (k = 0; k < Map->NumCells; k++) {
	y = Map->OrderedCells[k].y;
	x = Map->OrderedCells[k].x;
	  SoilMap[y][x].Runoff = 0.0;
	  if(Options->Sediment)
	    SedMap[y][x].TotalSediment = 0.0;
      }

      /* Must loop through surface routing multiple times within one DHSVM model time step. */

      while (Before(&(VariableTime.Current), &(NextTime.Current))) {
    
	/* Loop thru all of the cells in descending order of elevation */
	for (k = 0; k < Map->NumCells; k++) {
	  y = Map->OrderedCells[k].y;
	  x = Map->OrderedCells[k].x;
	
	  outflow = SoilMap[y][x].startRunoff;

	  slope = TopoMap[y][x].Slope; /*large grid cell slope.*/
          if (slope == 0) slope=0.0001;
	  else if (slope < 0) {
	    printf("negative slope in RouteSurface.c\n");
	    exit(0);
	  }

	  beta = 3./5.;
	  alpha = pow(SedType[SoilMap[y][x].Soil-1].Manning*pow(Map->DX,2/3)/sqrt(slope),beta);
	  
	  /* Calculate discharge from the grid cell using an explicit finite difference
	     solution of the linear kinematic wave. */

	  if (channel_grid_has_channel(ChannelData->stream_map, x, y) || channel_grid_has_channel(ChannelData->road_map, x, y)) 
	    {
	      outflow = 0.0;
	    }
	  else 
	    {
	      if(Runon[y][x] > 0.0001 || outflow > 0.0001) {

		outflow = ((VariableDT/Map->DX)*Runon[y][x] + 
		       alpha*beta*outflow * pow((outflow+Runon[y][x])/2.0, beta-1.) +
		       SoilMap[y][x].IExcess*Map->DX*VariableDT/Time->Dt) / ((VariableDT/Map->DX) +
									     alpha*beta*pow((outflow+
											     Runon[y][x])/2.0, beta-1.));
	      }
	      else if(SoilMap[y][x].IExcess > 0.0)
		outflow = SoilMap[y][x].IExcess*Map->DX*Map->DY/Time->Dt;
	      else
		outflow = 0.0;
	    }
	  
	  if(outflow < 0.0) 
	    outflow = 0.0;
	  
	  // Update surface water storage.  Make sure calculated outflow doesn't exceed available water.
	  // Otherwise, update surface water storage.
	  
	  if(outflow > (SoilMap[y][x].IExcess*(Map->DX*Map->DY)/VariableDT + Runon[y][x])) {
	    outflow = SoilMap[y][x].IExcess*(Map->DX*Map->DY)/VariableDT + (Runon[y][x]);
	    SoilMap[y][x].IExcess = 0.0;
	  }
	  else
	    SoilMap[y][x].IExcess += (Runon[y][x] - outflow)*VariableDT/(Map->DX*Map->DY);
	  

	  /*************************************************************/
	  /* PERFORM HILLSLOPE SEDIMENT ROUTING.                       */
	  /*************************************************************/
	  
	  if(Options->Sediment) {
	  
	    /* First find potential erosion due to rainfall Morgan et al. (1998). */   
	    /* Kinetic energy of the precip is determined in MassEnergyBalance.c */

	    h = SoilMap[y][x].IExcess; 
	    DR = (SedType[SoilMap[y][x].Soil-1].KIndex/PARTDENSITY) * 
	      PrecipMap[y][x].KineticEnergy*exp(-1.*SEDEXPONENT*h);

	    /* Find transport capacity of the flow out of this grid cell. */
	    
		// Olivier - 2003/07/08 : We need to add a condition on h 
		// if there is no surface water, there should be no transport capacity
		// By the way it prevents a division by 0 in the effectivepower calculation
		if (h!=0)
		{
			streampower= WATER_DENSITY*G*1000.*(outflow/Map->DX)*slope; /* g/s3 */

			effectivepower = sqrt(streampower)/pow(h*100.,2./3.); /* g^1.5 * s^-4.5 * cm^-2/3 */

			soliddischarge = (0.000001977) * pow(effectivepower, 1.044) * pow(SedType[SoilMap[y][x].Soil-1].d50, 0.478);  /* g/cm*s */

			TC = soliddischarge / (10.* (outflow/Map->DX) * PARTDENSITY);
		}
		else TC =0 ;
	    /* Find erosion due to overland flow after Morgan et al. (1998). */
	    
	    floweff = 0.79*exp(-0.85*SedType[SoilMap[y][x].Soil-1].Cohesion.mean);

	    DS = SedType[SoilMap[y][x].Soil-1].d50/1000.;
	    settling = (8.0*VISCOSITY/DS) * (sqrt(1.+ (PARTDENSITY - 1000.)*G*DS*DS*DS/(72.*VISCOSITY*VISCOSITY)) - 
					     1.0)/1000.;
	    /* Calculate sediment mass balance. */
 
	    term1 = (TIMEWEIGHT/Map->DX);
	    term2 = alpha/(2.*VariableDT);
	    term3 = (1.-TIMEWEIGHT/Map->DX);
	
		// Olivier - 2003/07/08 : we need to add a condition on outflow
		// To avoid a division by zero, we should add an if statement,
		// and it makes sense. If there is no outflow, SedOut should be null too.
		if (outflow > 0.0) {

			SedOut = (SedIn[y][x]*(term1*Runon[y][x]-term2*pow(Runon[y][x], beta)) +
			SedMap[y][x].OldSedOut*(term2*pow(SoilMap[y][x].startRunoff, beta) -
						term3*SoilMap[y][x].startRunoff) +
			SedMap[y][x].OldSedIn*(term2*pow(SoilMap[y][x].startRunon, beta) + 
						term3*SoilMap[y][x].startRunon) +
			DR + floweff*Map->DY*settling*TC)/(term2*pow(outflow, beta) + term1*outflow*floweff*Map->DY*settling);

			if(SedOut >= TC)
			SedOut = (SedIn[y][x]*(term1*Runon[y][x]-term2*pow(Runon[y][x], beta)) +
			SedMap[y][x].OldSedOut*(term2*pow(SoilMap[y][x].startRunoff, beta) -
						term3*SoilMap[y][x].startRunoff) +
			SedMap[y][x].OldSedIn*(term2*pow(SoilMap[y][x].startRunon, beta) + 
						term3*SoilMap[y][x].startRunon) +
			DR + Map->DY*settling*TC)/(term2*pow(outflow, beta) + term1*outflow*Map->DY*settling);

	    }
		else SedOut = 0.0 ;
	    

	    SedMap[y][x].OldSedOut = SedOut;
	    SedMap[y][x].OldSedIn = SedIn[y][x];
	    SedMap[y][x].TotalSediment += SedOut;
	    SedMap[y][x].erosion = (SedIn[y][x]*Runon[y][x] - SedOut*outflow)*Time->Dt/(Map->DX*Map->DY);;

		// This is a first attempt to determine overland inflow added to the channel
		// It has to be tested and checked.

		// if grid cell has a channel, sediment is intercepted by the channel
		if (channel_grid_has_channel(ChannelData->stream_map,x,y))
		{

			for (particle_size=0; particle_size<NSEDSIZES; particle_size++)
			{
				// For now, we suppose that the mass of sediments is uniformely 
				// distributed between all NSEDSIZES sizes of particles
				ChannelData->stream_map[x][y]->channel->sediment.overlandinflow[particle_size] += SedMap[y][x].TotalSediment * (1/NSEDSIZES) ;
			}
			SedMap[y][x].OldSedOut = 0.0;
			SedMap[y][x].TotalSediment = 0.0;
		}



	  }

	  /* Save sub-timestep runoff for q(i)(t-1) and q(i-1)(t-1) of next time step. */
	  SoilMap[y][x].startRunoff = outflow;
	  SoilMap[y][x].startRunon = Runon[y][x];
	  

	  /* Calculate total runoff in m/dhsvm timestep. */
	  SoilMap[y][x].Runoff += outflow*VariableDT/(Map->DX*Map->DY); 

	
	  /* Redistribute surface water to downslope pixels. */
	  if(outflow > 0) {  
	    for (n = 0; n < NDIRS; n++) {
	      int xn = x + xneighbor[n];
	      int yn = y + yneighbor[n];

	      /* If a channel cell runoff does not go to downslope pixels. */
	      if (valid_cell(Map, xn, yn)) {
		Runon[yn][xn] +=
		  outflow * ((float) TopoMap[y][x].Dir[n] /
			     (float) TopoMap[y][x].TotalDir);

		if(Options->Sediment)
		  SedIn[yn][xn] +=
		    SedOut * ((float) TopoMap[y][x].Dir[n] /
			     (float) TopoMap[y][x].TotalDir);
	      }
	    } /* end loop thru possible flow directions */
	  }
	  
	  /* Initialize runon for next timestep. */
	  Runon[y][x] = 0.0;

	  total_in += SoilMap[y][x].IExcess;
	  total_out += SoilMap[y][x].Runoff;
	  
	} /* end loop thru ordered basin cells */
	
        /* Increases time by VariableDT. */
	IncreaseVariableTime(&VariableTime, VariableDT, &NextTime);
	
	/*************************************************************/

      } /* End of internal time step loop. */

    }/* End of code added for kinematic wave routing. */

    if(Options->Sediment) {
      for (y = 0; y < Map->NY; y++) {
	free(SedIn[y]);
      }
      free(SedIn);
    }

    for (y = 0; y < Map->NY; y++) {
      free(Runon[y]);
    }
    free(Runon);
  }

  /* MAKE SURE THIS WORKS WITH A TIMESTEP IN SECONDS */
  else {			/* No network, so use unit hydrograph 
				   method */
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  TravelTime = (int) TopoMap[y][x].Travel;
	  if (TravelTime != 0) {
	    WaveLength = HydrographInfo->WaveLength[TravelTime - 1];
	    for (Step = 0; Step < WaveLength; Step++) {
	      Lag = UnitHydrograph[TravelTime - 1][Step].TimeStep;
	      Hydrograph[Lag] += SoilMap[y][x].Runoff *
		UnitHydrograph[TravelTime - 1][Step].Fraction;
	    }
	    SoilMap[y][x].Runoff = 0.0;
	  }
	}
      }
    }

    StreamFlow = 0.0;
    for (i = 0; i < Time->Dt; i++)
      StreamFlow += (Hydrograph[i] * Map->DX * Map->DY) / Time->Dt;

    /* Advance Hydrograph */
    for (i = 0; i < Time->Dt; i++) {
      for (j = 0; j < HydrographInfo->TotalWaveLength - 1; j++) {
	Hydrograph[j] = Hydrograph[j + 1];
      }
    }

    /* Set the last elements of the hydrograph to zero */
    for (i = 0; i < Time->Dt; i++)
      Hydrograph[HydrographInfo->TotalWaveLength - (i + 1)] = 0.0;

    PrintDate(&(Time->Current), Dump->Stream.FilePtr);
    fprintf(Dump->Stream.FilePtr, " %g\n", StreamFlow);
  }
}


/*****************************************************************************
  FindDT()
  Find the variable time step that will satisfy the courant condition for stability 
  in overland flow routing.
*****************************************************************************/

float FindDT(SOILPIX **SoilMap, MAPSIZE *Map, TIMESTRUCT *Time, 
	     TOPOPIX **TopoMap, SEDTABLE *SedType)
{
  int x, y;
  float slope, beta, alpha;
  float Ck;
  float DT, minDT;
  float numinc;
  float maxRunoff;

  maxRunoff = -99.;
  minDT = 36000.;

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {

      slope = TopoMap[y][x].Slope;
      beta = 3./5.;
      alpha = pow(SedType[SoilMap[y][x].Soil-1].Manning*pow(Map->DX,2/3)/sqrt(slope),beta);

      /* Calculate flow velocity from discharge  using manning's equation. */
      Ck = 1./(alpha*beta*pow(SoilMap[y][x].Runoff, beta -1.));

      if(SoilMap[y][x].Runoff > maxRunoff) {
	maxRunoff = SoilMap[y][x].Runoff;
      }

      if(Map->DX/Ck < minDT)
	minDT = Map->DX/Ck;

    }
  }

  /* Find the time step that divides evenly into Time->DT */
  
  numinc = (float) ceil((double)Time->Dt/minDT);
  DT = Time->Dt/numinc;
	
  if(DT > Time->Dt)
    DT = (float) Time->Dt;

  return DT;
}


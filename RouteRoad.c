/*
 * SUMMARY:      RouteRoad.c - Route surface flow and erosion for Forest Roads
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Colleen O. Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Nov-03
 * DESCRIPTION:  Route surface flow and erosion for Forest Roads
 * DESCRIP-END.
 * FUNCTIONS:    RouteRoad()
 * COMMENTS:
 * $Id$     
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "channel_grid.h"

/*****************************************************************************
  RouteRoad()
*****************************************************************************/
void RouteRoad(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap, 
	       SOILPIX ** SoilMap, ROADSTRUCT ** Network, SOILTABLE * SType, 
	       CHANNEL * ChannelData) 
{
  const char *Routine = "RouteRoad";
  int i,x,y;                     /* Counters */
  float dx, dy;                  /* Road grid cell dimensions (m)*/
  float cells;                   /* Number of road grid cells in a basin grid 
				    cell */
  float roadwidth, roadlength;   /* Road dimensions in current basin grid 
				    cell (m)*/
  float roadwater;               /* Depth (m) of water on the road surface */
  float slope;                   /* Slope of road surface the water travels 
				    over (m/m)*/    
  float *Runon;                  /* Water on moving downslope along the road (m) */
  float beta, alpha;             /* For Mannings. alpha is channel parameter 
				    including wetted perimeter, n, and slope. */
  float outflow;                 /* Flow out of a road grid cell in the 
				    subtimestep (m3/s) */
    
  TIMESTRUCT NextTime;
  TIMESTRUCT VariableTime;
  float VariableDT;              /* Maximum stable time step (s) */  

  if ((Runon = (float *) calloc(CELLFACTOR, sizeof(float))) == NULL)
      ReportError((char *) Routine, 1);
  
  NextTime = *Time;
  VariableTime = *Time;
  
  /* Holds the value of the next DHSVM time step. */ 
  IncreaseTime(&NextTime);
  
  /* Since Network.IExcess stays in the cell it is generated in,
     route each basin grid cell, with a road, separately */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
	  
	  /* Discretizing road into a grid for finite difference solution.
	     This assume the cells are oriented with the direction of flow. */
	  
	  dx = Network[y][x].FlowLength/(float) CELLFACTOR;
	  dy = dx; /* road grid cells are square */
	  roadlength = channel_grid_cell_length(ChannelData->road_map, x, y);
	  roadwidth = channel_grid_cell_width(ChannelData->road_map, x, y);
	  cells = roadlength*roadwidth/(dx*dy);

	  slope = Network[y][x].FlowSlope;
	  if (slope == 0) slope=0.0001;
	  else if (slope < 0) {
	    printf("RouteRoad.c: negative slope\n");
	    exit(0);
	  }
	  
	  beta = 3./5.;
	  alpha = pow(SType[SoilMap[y][x].Soil-1].Manning*pow(dx,2./3.)/sqrt(slope),beta);

	  /* Evenly distribute water over road surface. */
	  roadwater = (Network[y][x].IExcess * Map->DX * Map->DY)/
	    (roadlength*roadwidth);
	  for (i = 0; i < CELLFACTOR; i++)
	      Network[y][x].h[i] += roadwater;
	  
	  /* Use the Courant condition to find the maximum stable time step 
	     (in seconds). Must be an even increment of Dt. */
	  VariableDT = FindDTRoad(Network, Time, y, x, dx, beta, alpha);  
	  
	  /* Must loop through road segment routing multiple times within 
	     one DHSVM model time step. */
	  while (Before(&(VariableTime.Current), &(NextTime.Current))) {
	    
	    /* Loop through road grid cells starting at crown or road edge*/
	    for (i = 0; i < CELLFACTOR; i++){
	      
	      outflow = Network[y][x].startRunoff[i]; 
	      
	      /* Calculate discharge from the road segment using an explicit
		 finite difference solution of the linear kinematic wave. */
	      if(Runon[i] > 0.0001 || outflow > 0.0001) {
		outflow = ((VariableDT/dx)*Runon[i] + alpha*beta*outflow * 
			   pow((outflow+Runon[i])/2.,beta-1.) +
			   Network[y][x].h[i]*dx*VariableDT/Time->Dt)/
		  ((VariableDT/dx) + alpha*beta*pow((outflow+Runon[i])/2., beta-1.));
		
	      }
	      else if(Network[y][x].h[i] > 0.0)
		outflow = Network[y][x].h[i]*dy*dx/Time->Dt;
	      else
		outflow = 0.0;
	      
	      if(outflow < 0.0) outflow = 0.0;
	      
	      /*Update surface water storage.  Make sure calculated outflow 
		doesn't exceed available water.  Otherwise, update surface 
		water storage. */
	      if(outflow > (Network[y][x].h[i]*dy*dx)/Time->Dt + Runon[i]) 
		outflow = (Network[y][x].h[i]*dy*dx)/Time->Dt + Runon[i];
	      
	      Network[y][x].h[i] += (Runon[i]-outflow)*VariableDT/(dy*dx);
	      
	      /* Save sub-timestep runoff for q(i)(t-1) and q(i-1)(t-1) of next time step. */
	      Network[y][x].startRunoff[i] = outflow;
	      Network[y][x].startRunon[i] = Runon[i];
	      
	      /* Redistribute surface water to downslope pixel. */
	      if(outflow > 0.){
		if(i < (CELLFACTOR-1))
		  Runon[i+1] += outflow;
		/* If last pixel send outflow off road */
		else {
		  if(ChannelData->road_map[x][y]->channel->class->crown == CHAN_OUTSLOPED)
		    SoilMap[y][x].IExcess += ((outflow*VariableDT)/(Map->DX * Map->DY))*
		      (cells/(float)CELLFACTOR);
		  /* If the road is crowned, then the same amount of outflow goes to 
		     the ditch and off the road edge into the same pixel. This is similar to 
		     culvert flow. 0.5 accounts for 1/2 the cells on are one side of the
		     crown. */
		  else if(ChannelData->road_map[x][y]->channel->class->crown == CHAN_CROWNED){
		    channel_grid_inc_inflow(ChannelData->road_map, x, y, 
					    outflow*VariableDT*0.5*(cells/(float)CELLFACTOR));
		    SoilMap[y][x].RoadInt += ((outflow*VariableDT)/(Map->DX * Map->DY))*0.5*
		      (cells/(float)CELLFACTOR);
		    SoilMap[y][x].IExcess += ((outflow*VariableDT)/(Map->DX * Map->DY))*0.5*
		      (cells/(float)CELLFACTOR);
		  }
		  else { /* INSLOPED and all goes to ditch */
		    channel_grid_inc_inflow(ChannelData->road_map, x, y,
					    outflow*VariableDT*(cells/(float)CELLFACTOR));
		    SoilMap[y][x].RoadInt += ((outflow*VariableDT)/(Map->DX * Map->DY))*
		      (cells/(float)CELLFACTOR);
		  }
		}
	      }
	      /* Initialize for next timestep. */
	      Runon[i] = 0.0;
	    }
	    /* Increases time by VariableDT. */
	    IncreaseVariableTime(&VariableTime, VariableDT, &NextTime);
	    
	  } /* End of internal time step loop. */
	  /* Initialize for next DHSVM time step */
	  Network[y][x].IExcess = 0.0;
	}
      }
    }
  }/* End loop through basin grid cells */
  free(Runon);
  
}

/*****************************************************************************
  FindDTRoad()
  Find the variable time step that will satisfy the courant condition for stability 
  in overland flow routing.
*****************************************************************************/

float FindDTRoad(ROADSTRUCT **Network, TIMESTRUCT *Time, int y, int x, 
	     float dx, float beta, float alpha)
{
  int i;           /* counter */
  float Ck;        /* Flow velocity based on discharge using Manning's equation. */
  float DT, minDT; /* seconds */
  float numinc;
  float runoff, maxRunoff;
  
  maxRunoff = -99.;
  minDT = Time->Dt; 
  
  for (i = 0; i < CELLFACTOR; i++){
    
    runoff = Network[y][x].startRunoff[i];
    if (runoff == 0) runoff = 0.000000001;
    
    Ck = 1./(alpha*beta*pow(runoff, beta -1.)); 
    
    if(runoff > maxRunoff)
      maxRunoff = runoff;
    
    if(dx/Ck < minDT)
      minDT = dx/Ck;
    
  }
  
  /* Find the time step that divides evenly into Time->DT */
  numinc = (float) ceil((double)Time->Dt/minDT);
  DT = Time->Dt/numinc;
  
  if(DT > Time->Dt)
    DT = (float) Time->Dt;
  
  return DT;
}


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
		  DUMPSTRUCT *Dump, VEGPIX ** VegMap, VEGTABLE * VType)
{
  const char *Routine = "RouteSurface";
  int Lag;			/* Lag time for hydrograph */
  int Step;
  float StreamFlow;
  int TravelTime;
  int WaveLength;

  int i;
  int j;
  int x;
  int y;
  int n;
  float **surface;

  int dumpflag; 
  int count;
  void *Array;
  char buffer[32];
  char FileName[100];
  FILE *FilePtr;

  count = 0;

  if (Options->HasNetwork) {
    if ((surface = (float **) malloc(Map->NY * sizeof(float *))) == NULL) {
      ReportError((char *) Routine, 1);
    }

    for (y = 0; y < Map->NY; y++) {
      if ((surface[y] = (float *) malloc(Map->NX * sizeof(float))) == NULL) {
	ReportError((char *) Routine, 1);
      }
      else {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    surface[y][x] = SoilMap[y][x].Runoff;
	    SoilMap[y][x].Runoff = 0.0;
	    if(surface[y][x] > .01)
	      count += 1;
	  }
	}
      }
    }
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  if (VType[VegMap[y][x].Veg - 1].ImpervFrac > 0.0) {
	    SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].Runoff +=
	      surface[y][x];
	  }
	  else {
	    for (n = 0; n < NDIRS; n++) {
	      int xn = x + xneighbor[n];
	      int yn = y + yneighbor[n];
	      if (valid_cell(Map, xn, yn)) {
		SoilMap[yn][xn].Runoff +=
		  surface[y][x] * ((float) TopoMap[y][x].Dir[n] /
				   (float) TopoMap[y][x].TotalDir);
	      }
 	    }
	  }
	}
      }
    }

    /*************************************************************/
    /* Hack code added to dump surface runoff maps. */
    if (Options->Sediment) {
      dumpflag = 0;
    
      if(count > 97) dumpflag = 1;

      if(dumpflag == 1) {

	if (!(Array = calloc(Map->NY * Map->NX, sizeof(float))))
	  ReportError((char *) Routine, 1);
	for (y = 0; y < Map->NY; y++)
	  for (x = 0; x < Map->NX; x++)
	    ((float *)Array)[y*Map->NX + x] = surface[y][x];
	SPrintDate(&(Time->Current), buffer);
	
	sprintf(FileName, "%s/Map.%s.Runoff.bin",Dump->Path, buffer);
	
	if (!(FilePtr = fopen(FileName, "wb")))
	  ReportError(FileName, 3);

	fwrite(Array, sizeof(float), Map->NY*Map->NX, FilePtr);
	fclose(FilePtr);
	free(Array);
      }
    }

     /*************************************************************/
    /* End added code. */


    for (y = 0; y < Map->NY; y++) {
      free(surface[y]);
    }
    free(surface);
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

/*
 * SUMMARY:      MassBalance.c - calculate basin-wide mass balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Battelle - Pacific Northwest National Laboratory
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Oct-96
 * DESCRIPTION:  Calculate water and sediment mass balance errors
 *               
 * DESCRIP-END.
 * FUNCTIONS:    FinalMassBalance()
 * COMMENTS:
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  Aggregate()
  
  Calculate the average values for the different fluxes and state variables
  over the basin.  Only the runoff is calculated as a total volume instead
  of an average.  In the current implementation the local radiation
  elements are not stored for the entire area.  Therefore these components
  are aggregated in AggregateRadiation() inside MassEnergyBalance().

  The aggregated values are set to zero in the function RestAggregate,
  which is executed at the beginning of each time step.
*****************************************************************************/
void FinalMassBalance(FILES * Out, AGGREGATED * Total, WATERBALANCE * Mass,
		      OPTIONSTRUCT * Options, float roadarearatio)
{
  float NewWaterStorage;	/* water storage at the end of the time step */
  float Output;			/* total water flux leaving the basin;  */
  float MassError;		/* mass balance error m  */
  float Input;
  float MWMMassError;
  float SedInput, SedOutput, SedMassError;

  NewWaterStorage = Total->Soil.IExcess + Total->Road.IExcess + 
    Total->CanopyWater + Total->SoilWater +
    Total->Snow.Swq + Total->Soil.SatFlow;

  Output = Mass->CumChannelInt + ( Mass->CumRoadInt  -
    Mass->CumCulvertReturnFlow )+ Mass->CumET - Mass->CumSnowVaporFlux;

  Input = Mass->CumPrecipIn;
  MassError = (NewWaterStorage - Mass->StartWaterStorage) +
   Output - Input;

/*   fprintf(Out->FilePtr, "\n Final Mass Balance\n"); */
/*   fprintf(Out->FilePtr, " Precip (mm): %f\n", Mass->CumPrecipIn * 1000.);  */
/*   fprintf(Out->FilePtr, " ET (mm): %f\n", Mass->CumET * 1000.);  */
/*   fprintf(Out->FilePtr, " SnowVaporFlux (mm): %f\n", Mass->CumSnowVaporFlux *
     1000.);  */
/*   fprintf(Out->FilePtr, " Runoff (mm): %f\n", Total->Runoff * 1000.);  */
/*   fprintf(Out->FilePtr, " RoadInt (mm): %f\n", Mass->CumRoadInt * 1000.);  */
/*   fprintf(Out->FilePtr, " ChannelInt (mm): %f\n", Mass->CumChannelInt *
     1000.);  */
/*   fprintf(Out->FilePtr, " Mass Error (mm): %f\n", MassError * 1000.);  */

  fprintf(stderr, "\nFinal Mass Balance\n");
  fprintf(stderr, " \nInflow, %.2f:\n", Input*1000.);
  fprintf(stderr, " Precip (mm): %f\n", Mass->CumPrecipIn * 1000.);
  fprintf(stderr, " \nOutflow %.2f:\n", Output*1000.);
  fprintf(stderr, " ET (mm): %f\n", Mass->CumET * 1000.);
  fprintf(stderr, " SnowVaporFlux (mm): %f\n", -1.*Mass->CumSnowVaporFlux * 1000.);
  fprintf(stderr, " CulvertToChannel (mm): %f\n", Mass->CumCulvertToChannel *
	  1000.);
  fprintf(stderr, " ChannelInt (mm): %f\n", Mass->CumChannelInt * 1000.);
  fprintf(stderr, " \nStorage:\n");
  fprintf(stderr, " Initial Storage (mm): %f\n", Mass->StartWaterStorage * 1000.);
  fprintf(stderr, " End of Run Storage (mm): %f\n", NewWaterStorage * 1000.);
  fprintf(stderr, " \tFinal SWQ (mm): %f\n", Total->Snow.Swq * 1000.);
  fprintf(stderr, " \tFinal Soil Moisture (mm): %f\n", (Total->SoilWater + Total->Soil.SatFlow) * 1000.);
  fprintf(stderr, " \tFinal Surface (mm): %f\n", (Total->Soil.IExcess + Total->Road.IExcess + 
						  Total->CanopyWater) * 1000.);
  fprintf(stderr, " \tFinal Road Surface (mm): %f\n", Total->Road.IExcess *1000.);


  fprintf(stderr, " \nOther:\n");
  fprintf(stderr, " RoadInt (mm): %f\n", Mass->CumRoadInt * 1000.);
  fprintf(stderr, " CulvertReturnFlow (mm): %f\n", Mass->CumCulvertReturnFlow *
	  1000.);
  fprintf(stderr, " Mass Error (mm): %f\n", MassError * 1000.);
  fprintf(stderr, " Mass added to glacier (mm) %f\n",
	  Total->Snow.Glacier * 1000.);

  if(Options->Sediment){
    
    fprintf(stderr, "\nFinal Sediment Mass Balance\n");
    
    if (Options->MassWaste){
      
      MWMMassError = Mass->CumMassWasting - Mass->CumSedimentToChannel - 
	Mass->CumMassDeposition;
      
      fprintf(stderr, " \nTotal Mass Wasting\n");
      fprintf(stderr, " MassWasted (m3): %.2e\n", Mass->CumMassWasting);
      fprintf(stderr, " SedimentToChannel (m3): %.2e\n", 
	      Mass->CumSedimentToChannel);
      fprintf(stderr, " MassDepostion (m3): %.2e\n", Mass->CumMassDeposition);
      fprintf(stderr, " Mass Error (m3): %e\n", MWMMassError);
    }
    
    if (Options->InitSedFlag){
      fprintf(stderr, " \nAverage Surface Erosion\n");
      fprintf(stderr, " Surface Erosion (mm): %.2e\n", 
	      Mass->CumSedimentErosion);
      fprintf(stderr, " Surface Erosion (kg/hectare): %.2e\n", 
	      Mass->CumSedimentErosion * PARTDENSITY * (float)MMTOM * 10000.);
    }
    
    if (Options->RoadRouting){
      fprintf(stderr, " \nBasin Average Road Surface Erosion\n");
      fprintf(stderr, " Road Surface Erosion (mm): %.2e\n", 
	      Mass->CumRoadErosion * 1000.);
      fprintf(stderr, " Road Surface Erosion (kg/hectare): %.2e\n", 
	      Mass->CumRoadErosion * PARTDENSITY * 10000.);
      fprintf(stderr, "Road Sediment to Hillslope (mm): %.2e\n", 
	      Mass->CumRoadSedHill * 1000.);

      /* roadarearatio is used to convert basin average road erosion to 
	 road erosion averaged over the road surface area */

      fprintf(stderr, " \nAverage Road Surface Erosion\n");
      fprintf(stderr, " Road Surface Erosion (mm): %.2e\n", 
	      Mass->CumRoadErosion/roadarearatio * 1000.);
      fprintf(stderr, " Road Surface Erosion (kg/hectare): %.2e\n", 
	      Mass->CumRoadErosion/roadarearatio * PARTDENSITY * 10000.);
      fprintf(stderr, "Road Sediment to Hillslope (mm): %.2e\n", 
	      Mass->CumRoadSedHill/roadarearatio * 1000.);
    }

    SedInput = Mass->CumDebrisInflow + 
      (Mass->CumSedOverlandInflow - Mass->CumCulvertSedToChannel) + 
      Mass->CumSedOverroadInflow;

    /* NOTE: CulvertSedToChannel + CulvertReturnSedFlow = CulvertSedFlow */
    SedOutput = Mass->CumSedimentOutflow - 
      (Mass->CumCulvertSedToChannel + Mass->CumCulvertReturnSedFlow) ; 
    
    SedMassError = (Total->ChannelSedimentStorage + 
		    Total->ChannelSuspendedSediment - 
		    Mass->StartChannelSedimentStorage) + 
      SedOutput - SedInput;  
    
    fprintf(stderr, " \nChannel Erosion");
    fprintf(stderr, " \nInflow %.2e:\n", SedInput); 
    fprintf(stderr, " DebrisInflow (kg): %e\n", Mass->CumDebrisInflow); 
    fprintf(stderr, " OverlandInflow (kg): %e\n", 
	    Mass->CumSedOverlandInflow);
    fprintf(stderr, " OverroadInflow (kg): %e\n", 
	     Mass->CumSedOverroadInflow);
    
    if (Options->ChannelRouting){
      fprintf(stderr, " \nOutflow %.2e:\n",  SedOutput);
      fprintf(stderr, " SedimentOutflow (kg): %e\n", Mass->CumSedimentOutflow);
      fprintf(stderr, " CulvertReturnSedFlow (kg): %e\n", 
	      Mass->CumCulvertReturnSedFlow); 
      fprintf(stderr, " CulvertSedToChannel (kg): %e\n", 
	      Mass->CumCulvertSedToChannel);    
      
      fprintf(stderr, " \nStorage:\n");
      fprintf(stderr, " Initial Storage (kg): %e\n", 
	      Mass->StartChannelSedimentStorage); 
      fprintf(stderr, " End of Run Storage (kg): %e\n", 
	      Total->ChannelSedimentStorage + Total->ChannelSuspendedSediment); 
      fprintf(stderr, " \tFinal Bed Storage (kg): %e\n", 
	      Total->ChannelSedimentStorage); 
      fprintf(stderr, " \tFinal Suspended Sediment (kg): %e\n", 
	      Total->ChannelSuspendedSediment); 
      fprintf(stderr, " \nMass Error (kg): %e\n", SedMassError); 
    }
  }
}


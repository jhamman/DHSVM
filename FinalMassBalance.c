/*
 * SUMMARY:      MassBalance.c - calculate basin-wide mass balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Battelle - Pacific Northwest National Laboratory
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Oct-96
 * DESCRIPTION:  Calculate water mass balance errors
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
void FinalMassBalance(FILES * Out, AGGREGATED * Total, WATERBALANCE * Mass)
{
  float NewWaterStorage;	/* water storage at the end of the time step */
  float Output;			/* total water flux leaving the basin;  */
  float MassError;		/* mass balance error m  */
  float Input;

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

  fprintf(stderr, " \nOther:\n");
  fprintf(stderr, " RoadInt (mm): %f\n", Mass->CumRoadInt * 1000.);
  fprintf(stderr, " CulvertReturnFlow (mm): %f\n", Mass->CumCulvertReturnFlow *
	  1000.);
  fprintf(stderr, " Mass Error (mm): %f\n", MassError * 1000.);
  fprintf(stderr, " Mass added to glacier (mm) %f\n",
	  Total->Snow.Glacier * 1000.);
}

/*
 * SUMMARY:      RouteChannelSediment
 * USAGE:        
 * * AUTHOR:       Edwin P. Maurer
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Sep-02
 * Last Change:  Thu Jun 19 09:27:02 2003 by Ed Maurer <edm@u.washington.edu>
 * DESCRIPTION:  
 * DESCRIP-END.   
 * FUNCTIONS:    main()
 * COMMENTS:
 */

/******************************************************************************/
/*			     INCLUDES                                  */
/******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "data.h"
#include "functions.h"
#include "DHSVMChannel.h"
#include "DHSVMerror.h"

/*****************************************************************************
  Alloc_Chan_Sed_Mem
*****************************************************************************/
/* void Alloc_Chan_Sed_Mem(float ** DummyVar) */
/* { */
/*    if (!(*DummyVar = (float *) calloc(NSEDSIZES, sizeof(float)))) */
/*     ReportError(" Alloc_Chan_Sed_Mem", 1); */
/*  } */
/*****************************************************************************
  InitChannelSediment)

  Assign initial colluvium mass to each unique channel ID (amount
  of storage, kg)
*****************************************************************************/
int InitChannelSediment(Channel * Head, AGGREGATED * Total)
{
  if (Head != NULL){ 
    Channel *Current = NULL;
    int i;
    float InitialDepth = 0.010; /* initial depth of sediment in the channel, m */
    float bulkporosity, initvol;
    
    bulkporosity = 0.245+0.14*pow(DEBRISd50,-0.21); /* Komura, 1961 relation */
   
    /* Assign the storages to the correct IDs */
    Current = Head;
    while (Current) {
      
      initvol = Current->length * InitialDepth * Current->class->width;
      for(i=0;i<NSEDSIZES;i++) {
	Current->sediment.mass[i] = 
	  initvol*(1.-bulkporosity)*((float) PARTDENSITY)*(1./((float) NSEDSIZES));
	Current->sediment.debrisinflow[i]=0.0; 
	Current->sediment.overlandinflow[i]=0.0;
	Current->sediment.overroadinflow[i]=0.0;
/* 	Current->sediment.inflow[i]=0.0; */
	Current->sediment.inflowrate[i]=0.0;
	Current->sediment.last_inflowrate[i]=0.0; 
	Current->sediment.outflow[i]=0.0;
	Current->sediment.last_outflow[i]=0.0; 
	Current->sediment.outflowrate[i]=0.0;
	Current->sediment.last_outflowrate[i]=0.0;
	Total->ChannelSedimentStorage += Current->sediment.mass[i];
      }
      Current = Current->next;
    }
  }
  return (0);
}
/*****************************************************************************
  InitChannelSedInflow
  
  Assign initial colluvium mass to each unique channel ID (amount
  of storage, kg)
*****************************************************************************/
int InitChannelSedInflow(Channel * Head)
{
  if (Head != NULL){ 
    Channel *Current = NULL;
    int i;
    
    Current = Head;
    while (Current) {
      for(i=0;i<NSEDSIZES;i++) {
	Current->sediment.inflow[i] = 0.0;
      }
      Current->sediment.outflowconc = 0.0;
      Current->sediment.totalmass = 0.;
      Current = Current->next;
    }
  }
  return (0);
}
/*****************************************************************************
  SaveChannelSedInflow
  
  For FinalMassBalance output
*****************************************************************************/
int SaveChannelSedInflow(Channel * Head, AGGREGATED * Total)
{
  if (Head != NULL){ 
    Channel *Current = NULL;
    int i;
    
    Current = Head;
    while (Current) {
      for(i=0;i<NSEDSIZES;i++) {
	Total->DebrisInflow += Current->sediment.debrisinflow[i];
	Current->sediment.debrisinflow[i] = 0.;
	Total->SedimentOverlandInflow += Current->sediment.overlandinflow[i];
	Current->sediment.overlandinflow[i] = 0.;
	Total->SedimentOverroadInflow += Current->sediment.overroadinflow[i];
	Current->sediment.overroadinflow[i] = 0.;
      }
      Current = Current->next;
    }
  }
  return (0);
}


/*****************************************************************************
  RouteChannelSediment()

  Read in DHSVM sediment mass and inflows for each channel segment, and 
  route sediment downstream. Sorts by particle size, transports finer material
  first, as done by Williams (1980).

*****************************************************************************/
void RouteChannelSediment(Channel * Head, TIMESTRUCT Time, 
			  DUMPSTRUCT *Dump, AGGREGATED * Total,
			  float *SedDiams)
{
  Channel *Current = NULL;
  float DS,DT_sed,numinc;
  float flowdepth,Qavg,V,dIdt,dOdt,dMdt,Vsed,Vshear,Vshearcrit;
  float minDT_sed,TotalCapacityUp,TotalCapacityDown;
  float lateral_sed_inflow_rate;
  float TotalCapacity, CapacityUsed;
  float Qup,Qdown;
  float phi=0.55, theta=0.55,term3,term4; /*space and time weighting factors*/
  int i,tstep;
  int order;
  int order_count;
 
  /* the next 5 lines are from channel_route_network - used to order streams */
  for (order = 1;; order += 1) {
    order_count = 0;
    Current = Head;
    
    while (Current != NULL) {
      if (Current->order == order) {
/* 	Current->sediment.outflowconc = 0.0; */
	CapacityUsed = 0.0;
	
/* 	if (Current->outlet != NULL) */
/* 	  Current->outlet->sediment.totalmass = 0.; */
	
	/* rate of inflow and outflow change over model time step*/
	dIdt = (Current->inflow - Current->last_inflow)/(float) Time.Dt;
	dOdt = (Current->outflow - Current->last_outflow)/(float) Time.Dt;
	
	/****************************************/
	/* Estimate sub-time step for the reach */
	/****************************************/
	minDT_sed = 3600.;
	/* Estimate flow velocity from discharge using manning's equation. */
	Qavg = (Current->inflow+Current->outflow)/(2.0*(float) Time.Dt);

	/* If there is no flow (true for roads), move on to the next segment */
	if(Qavg > 0){
	  if(Current->slope>0.0) {
	    flowdepth = pow(Qavg*Current->class->friction/(Current->class->width*sqrt(Current->slope)),0.6);
	    V = Qavg/(flowdepth*Current->class->width);
	  }
	  else V=0.01;
	  if(Current->length/V < minDT_sed) minDT_sed = 1.0*Current->length/V;
	  numinc = (float) ceil((double)Time.Dt/minDT_sed);
	  if(numinc<1) numinc=1;
	  DT_sed = (float) Time.Dt/numinc;
	  

	  /* Initialize sediment.outflow for this segment 
	     and calculate inflow from upstream reach */
	  
	  for(i=0;i<NSEDSIZES;i++) {
	    Current->sediment.outflow[i]=0.0;
	    Current->sediment.inflowrate[i] = Current->sediment.inflow[i]/(float) Time.Dt;
	  }
	  
	  /****************************************/
	  /* Loop for each sub-timestep           */
	  /****************************************/
	  for(tstep=0;tstep<numinc;tstep++) {
	    
	    CapacityUsed=0.0;
	    
	    Qup = Current->last_inflow + dIdt*tstep*DT_sed;
	    Qdown = Current->last_outflow + dOdt*tstep*DT_sed;
	    
	    /****************************************/
	    /* Loop for each particle size          */
	    /****************************************/
	    /*DO NOT USE BAGNOLD's EQ. FOR D<0.015 mm - this is wash load anyway*/
	    for(i=0;i<NSEDSIZES;i++) {
	      DS = SedDiams[i]*((float) MMTOM); /* convert from mm to m */
	      dMdt=0;
	      
	      /* lateral inflow for the reach per second kg/s */
	      lateral_sed_inflow_rate = (Current->sediment.debrisinflow[i] + 
					 Current->sediment.overlandinflow[i] +
					 Current->sediment.overroadinflow[i])/(float) Time.Dt;
	      
	      /****************************************/
	      /* Find rate of bed change and new mass */
	      /****************************************/
	      
	      /* Set theta to 1.0 to prevent instabilities during mass wasting inflow  */ 
	      if(Current->sediment.debrisinflow[i]>0)
		theta=1.0;
	      
	      /*  Set theta to 1.0 to prevent instabilities during large differences 
		  between current and previous steps */
	      if(Current->sediment.inflowrate[i]>0 || Current->sediment.last_inflowrate[i]>0 ){
		if(abs(1-Current->sediment.last_inflowrate[i]/Current->sediment.inflowrate[i])>0.75 || abs(1-Current->sediment.inflowrate[i]/Current->sediment.last_inflowrate[i])>0.75 || abs(1-Current->sediment.outflowrate[i]/Current->sediment.inflowrate[i])>0.7  )
		  theta=1.0;
		else theta=0.55;} /* this should be .55 */
	      else theta=1.0;
	      /* TotalCapacity is in kg/s */
	      if(SedDiams[i] < 0.062) { /* per Wicks and Bathurst, wash load */
		TotalCapacity = 
		  Current->sediment.inflowrate[i]+Current->sediment.mass[i]/DT_sed;
	      }
	      else {
		TotalCapacityUp = CalcBagnold(DS,&Time,Qup,Current->class->width,
					      Current->class->friction,Current->slope);
		TotalCapacityDown = CalcBagnold(DS,&Time,Qdown,Current->class->width,
						Current->class->friction,Current->slope);
		TotalCapacity=phi*TotalCapacityDown + (1.0-phi)*TotalCapacityUp;
		TotalCapacity -= CapacityUsed; /* Avoid mult use of streampower */
	      }
	      
	      if(TotalCapacity<=0) TotalCapacity=0.0;
	      
	      if(TotalCapacity*DT_sed > Current->sediment.mass[i]) {	   
		dMdt= -Current->sediment.mass[i]/DT_sed;
		Current->sediment.mass[i] = 0.;	     
	      }
	      
	      else {
		dMdt =-TotalCapacity;
		Current->sediment.mass[i] -=  TotalCapacity*DT_sed;	     
	      }
	      
	      /****************************************/
	      /* Calculate reach sed outflow rate     */
	      /****************************************/
	      /* limit it to the total available sediment transport capacity */
	      term3 = (1.-theta) * 
		(Current->sediment.last_outflowrate[i]-Current->sediment.last_inflowrate[i]);
	      term4 = theta * Current->sediment.inflowrate[i];
	      
	      Current->sediment.outflowrate[i] = 
		(1./theta)*(lateral_sed_inflow_rate-dMdt-term3+term4);
	      if(Current->sediment.outflowrate[i]<0.0){
		Current->sediment.outflowrate[i]=0.0;}
	      
	      if(Current->sediment.outflowrate[i]>=TotalCapacity) {
		Current->sediment.mass[i] += 
		  (Current->sediment.outflowrate[i]-TotalCapacity)*DT_sed;
		dMdt += Current->sediment.outflowrate[i]-TotalCapacity;
		Current->sediment.outflowrate[i]=TotalCapacity;
	      }
	      
	      
	      /****************************************/
	      /* Assign new values to next step old   */
	      /****************************************/
	      Current->sediment.last_outflowrate[i]=Current->sediment.outflowrate[i];
	      Current->sediment.last_inflowrate[i]=Current->sediment.inflowrate[i];
	      
	      /****************************************/
	      /* Accumulate reach sed outflow mass    */
	      /****************************************/
	      Current->sediment.outflow[i] += Current->sediment.outflowrate[i]*DT_sed;
	      
	      CapacityUsed += Current->sediment.outflowrate[i];
	      
	      Current->sediment.totalmass += Current->sediment.mass[i];
	      
	      /* Calculate outflow concentration in mg/l 
		 Check for division by zero (concentration should be 0 if 
		 there is no outflow) */
	      if (Current->outflow >= 0.0) {
		if(Current->slope>0.0) 
		  flowdepth = pow(Current->outflow*Current->class->friction/(Current->class->width*sqrt(Current->slope)),0.6);
		else flowdepth = 0.005;
		
		V = Current->outflow/(flowdepth*Current->class->width);
		
		if (Current->slope >= 0.0)
		  Vshear = sqrt(G*flowdepth*Current->slope); /* approximate */
		else Vshear = 0.0; /* could calculate with friction slope */ 
		
		if(SedDiams[i] < 15) Vshearcrit = -0.0003*SedDiams[i]*SedDiams[i]+0.0109*SedDiams[i]+0.0106;
		else Vshearcrit = 0.0019*SedDiams[i]+0.0926; /* both in m/s */
	    
	      /* 		if (Vshearcrit > 0.95*Vshear) */
	      /*	  Vsed = 8.5*Vshear*sqrt(0.05); arbitrary - substituted Vshearcrit = 0.95*Vshear in equation below */
	      /*	else sed velocity per Wicks and Bathurst, fine particles travel at flow velocity */
	      /* 	  Vsed = 8.5*Vshear*sqrt(1-(Vshearcrit/Vshear)); */
	    
	      /* 	if(Vsed>V || SedDiams[i]<=0.062) Vsed=V; */
	    
	      /* 	if(Vsed<(0.2*V)) Vsed=0.2*V; */
		Vsed = V;	
		Current->sediment.outflowconc += 
		  (1000.0*Current->sediment.outflow[i]/(Current->outflow))*(V/Vsed);
	      }	  
	    } /* close loop for each sediment size */	  	  
	  } /* end of sub-time step loop */
	  
	  for(i=0;i<NSEDSIZES;i++) {
	    
	    /* pass the sediment mass outflow to the next downstream reach */
	    if(Current->outlet != NULL){
	      Current->outlet->sediment.inflow[i] += Current->sediment.last_outflow[i];
	      Current->sediment.last_outflow[i] =  Current->sediment.outflow[i];
	      /* Needed for last time step to balance mass */
	      Total->ChannelSuspendedSediment += Current->sediment.outflow[i];
	    }
	    /* If no stream segment outlet, there is a road sink or a basin outlet.
	       Track this for the sediment mass balance. */
	    else{
	      Total->SedimentOutflow += Current->sediment.outflow[i];
	    }
	    Total->ChannelSedimentStorage += Current->sediment.mass[i];	  
	  }
	} /* end 	if(Qavg > 0){ */
	else {/* if Qvag < 0 */
	  for(i=0;i<NSEDSIZES;i++) {
	    Current->sediment.mass[i] += Current->sediment.debrisinflow[i] + 
	      Current->sediment.overlandinflow[i] + Current->sediment.overroadinflow[i];
	    Total->ChannelSedimentStorage += Current->sediment.mass[i];	  
	  }
	}
	
	/* the next 7 lines are from channel_route_network -- closes the loop above */
	order_count += 1;
      } /* close if statement checking for stream order */
      Current = Current->next;  
    } /* close while statement checking that CURRENT != NULL */
    if (order_count == 0)
      break;
  } /* close loop for the stream order */
}


/*****************************************************************************
  RouteCulvertSediment()

*****************************************************************************/
void RouteCulvertSediment(CHANNEL * ChannelData, MAPSIZE * Map,
			  TOPOPIX ** TopoMap, SEDPIX ** SedMap, 
			  AGGREGATED * Total, float *SedDiams)
{
  int x, y;
  float CulvertSedFlow;   /* culvert flow of sediment, kg */
  int i;
  float temp = 0.;
  Total->CulvertReturnSedFlow = 0.0;

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
 
	for(i=0; i<NSEDSIZES; i++) {
	    
	  CulvertSedFlow = ChannelCulvertSedFlow(y, x, ChannelData, i);
	  CulvertSedFlow /= Map->DX * Map->DY;

	  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
	    /* Percent delivery to streams is conservative and based on particle size */
	    if (SedDiams[i] <= 0.063){
	      ChannelData->stream_map[x][y]->channel->sediment.overlandinflow[i] += CulvertSedFlow;
	      Total->CulvertSedToChannel += CulvertSedFlow;
	      CulvertSedFlow = 0.;
	    }
	    if ((SedDiams[i] > 0.063) && (SedDiams[i] <= 0.5)){
	      ChannelData->stream_map[x][y]->channel->sediment.overlandinflow[i] += 0.3*CulvertSedFlow;
	      Total->CulvertSedToChannel += 0.3*CulvertSedFlow;
	      Total->CulvertReturnSedFlow += 0.7*CulvertSedFlow;
	      CulvertSedFlow = 0.;
	    }
	    if ((SedDiams[i] > 0.5) && (SedDiams[i] <= 2.)){
	      ChannelData->stream_map[x][y]->channel->sediment.overlandinflow[i] += 0.1*CulvertSedFlow;
	      Total->CulvertSedToChannel += 0.1*CulvertSedFlow;
	      Total->CulvertReturnSedFlow += 0.9*CulvertSedFlow;
	      CulvertSedFlow = 0.;
	    }
	    Total->CulvertReturnSedFlow += CulvertSedFlow;
	  }
	  else {
	    Total->CulvertReturnSedFlow += CulvertSedFlow;
	  }
	}
      }
    }
  }
}


/*****************************************************************************
   ChannelCulvertSedFlow ()   
   computes sediment outflow (kg) of channel/road network to a grid cell, if it
   contains a sink (sink check is in channel_grid_sed_outflow)

*****************************************************************************/
double ChannelCulvertSedFlow(int y, int x, CHANNEL * ChannelData, int i)
{
  if (channel_grid_has_channel(ChannelData->road_map, x, y)){
    return channel_grid_sed_outflow(ChannelData->road_map, x, y, i);
  }
  else {
    return 0;
  }
}

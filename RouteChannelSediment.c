/*
 * SUMMARY:      RouteChannelSediment
 * USAGE:        
 *
 * AUTHOR:       Edwin P. Maurer
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Sep-02
 * Last Change:  Thu Jun 19 09:27:02 2003 by Ed Maurer <edm@u.washington.edu>
 * DESCRIPTION:  
 * DESCRIP-END.   
 cd
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
void InitChannelSediment(Channel * Head, AGGREGATED * Total)
{
  Channel *Current = NULL;
  int i;
  float InitialDepth = 0.010; /* initial depth of sediment in the channel, m */
  float bulkporosity, initvol;

  printf("Initializing channel sediment\n\n");

  bulkporosity = 0.245+0.14*pow(DEBRISd50,-0.21); /* Komura, 1961 relation */

  /* Assign the storages to the correct IDs */
  Current = Head;
  while (Current) {
    
    initvol = Current->length * InitialDepth * Current->class->width;
    for(i=0;i<NSEDSIZES;i++) {
      Current->sediment.debrisinflow[i]=0.0; 
      Current->sediment.overlandinflow[i]=0.0;
      Current->sediment.inflow[i]=0.0;
      Current->sediment.inflowrate[i]=0.0;
      Current->sediment.last_inflowrate[i]=0.0; 
      Current->sediment.outflow[i]=0.0;
      Current->sediment.outflowrate[i]=0.0;
      Current->sediment.last_outflowrate[i]=0.0; 
      Current->sediment.mass[i] = 
	initvol*(1-bulkporosity)*((float) PARTDENSITY)*(1/((float) NSEDSIZES));
      Total->ChannelSedimentStorage += Current->sediment.mass[i];
    }
    Current = Current->next;
  }
}
/*****************************************************************************
  InitChannelSedInflow
  
  Assign initial colluvium mass to each unique channel ID (amount
  of storage, kg)
*****************************************************************************/
void InitChannelSedInflow(Channel * Head)
{
  Channel *Current = NULL;
  int i;

  Current = Head;
  while (Current) {
    for(i=0;i<NSEDSIZES;i++) {
      Current->sediment.inflow[i] = 0.0;
    }
    Current = Current->next;
  }
}

/*****************************************************************************
  RouteChannelSediment()

  Read in DHSVM sediment mass and inflows for each channel segment, and 
  route sediment downstream. Sorts by particle size, transports finer material
  first, as done by Williams (1980).

*****************************************************************************/
void RouteChannelSediment(Channel * Head, Channel *RoadHead, TIMESTRUCT Time, 
			  DUMPSTRUCT *Dump, AGGREGATED * Total)
{
  Channel *Current = NULL;
  float DS,DT_sed,numinc;
  float flowdepth,Qavg,V,dIdt,dOdt,dMdt,Vsed,Vshear,Vshearcrit;
  float minDT_sed,TotalCapacityUp,TotalCapacityDown;
  float lateral_sed_inflow_rate;
  float TotalCapacity, CapacityUsed;
  float SedDiams[NSEDSIZES];
  float Qup,Qdown;
  float phi=0.55, theta=0.55,term3,term4; /*space and time weighting factors*/
  int i,tstep;
  int order;
  int order_count;

  /* For each of the sediment diameters, calculate the mass balance */
  DistributeSedimentDiams(SedDiams);  /* find diameter for each portion */
 
  /* the next 5 lines are from channel_route_network - used to order streams */
  for (order = 1;; order += 1) {
    order_count = 0;
    Current = Head;
    
    while (Current != NULL) {
      if (Current->order == order) {
	
	Current->sediment.outflowconc = 0.0;
	CapacityUsed = 0.0;
	
	if (Current->outlet != NULL)
	  Current->outlet->sediment.totalmass = 0.;
	
	/* rate of inflow and outflow change over model time step*/
	dIdt = (Current->inflow - Current->last_inflow)/(float) Time.Dt;
	dOdt = (Current->outflow - Current->last_outflow)/(float) Time.Dt;
	
	/****************************************/
	/* Estimate sub-time step for the reach */
	/****************************************/
	minDT_sed = 3600.;
	/* Estimate flow velocity from discharge using manning's equation. */
	Qavg = (Current->inflow+Current->outflow)/(2.0*(float) Time.Dt);
	if(Current->slope>0.0) {
	  flowdepth = pow(Qavg*Current->class->friction/(Current->class->width*sqrt(Current->slope)),0.6);
	  V = Qavg/(flowdepth*Current->class->width);
	}
	else V=0.01;
	if(Current->length/V < minDT_sed) minDT_sed = Current->length/V;
	numinc = (float) ceil((double)Time.Dt/minDT_sed);
	if(numinc<1) numinc=1;
	DT_sed = (float) Time.Dt/numinc;

	/****************************************/
	/* Loop for each particle size          */
	/****************************************/
	/*DO NOT USE BAGNOLD's EQ. FOR D<0.015 mm - this is wash load anyway*/
	for(i=0;i<NSEDSIZES;i++) {
	  Current->sediment.outflow[i]=0.0;
	  DS = SedDiams[i]*((float) MMTOM); /* convert from mm to m */

	  /****************************************/
    	  /* Calculate segment sed inflows        */
	  /****************************************/
	  /*exclude the slug of sediment from the first time step */
	  if(IsEqualTime(&(Time.Current), &(Time.Start))) 
	    lateral_sed_inflow_rate = 0.0;
     	  /* lateral inflow for the reach per second kg/s */
	  else lateral_sed_inflow_rate = (Current->sediment.debrisinflow[i] + Current->sediment.overlandinflow[i])/(float) Time.Dt;

	  /* inflow from upstream reach */
	  Current->sediment.inflowrate[i] = Current->sediment.inflow[i]/(float) Time.Dt;

	  /****************************************/
	  /* Loop for each sub-timestep           */
	  /****************************************/
	  for(tstep=0;tstep<numinc;tstep++) {
	    
	    Qup = Current->last_inflow + dIdt*tstep*DT_sed;
	    Qdown = Current->last_outflow + dOdt*tstep*DT_sed;
	    
	    /****************************************/
	    /* Find rate of bed change and new mass */
	    /****************************************/
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
	    if(TotalCapacity*DT_sed >= Current->sediment.mass[i]) {
	      dMdt = Current->sediment.mass[i]/DT_sed;
	      Current->sediment.mass[i] = 0.;
	    }
	    else {
	      dMdt = TotalCapacity;
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
	    if(Current->sediment.outflowrate[i]<0.0) Current->sediment.outflowrate[i]=0.0;
	    if(Current->sediment.outflowrate[i]>TotalCapacity) {
	      Current->sediment.mass[i] += 
		(Current->sediment.outflowrate[i]-TotalCapacity)*DT_sed;
	      dMdt -= Current->sediment.outflowrate[i]-TotalCapacity;
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
	    
	    CapacityUsed += dMdt;
	    
	  } /* end of sub-time step loop */
	  
	  Current->sediment.totalmass += Current->sediment.mass[i];
	  
	  /* calculate outflow concentration in mg/l */
	  
	  // Olivier - 2003/07/08 : We should add a condition on current->outflow
	  // to avoid division by zero in the calculation of Current->sediment.outflowconc
	  // It makes sense as the concentration should be 0 if there is no outflow
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
		Vsed=V;	
		Current->sediment.outflowconc += 
		(1000.0*Current->sediment.outflow[i]/(Current->outflow))*(V/Vsed);
	  }
	  
	  /* pass the sediment mass outflow to the next downstream reach */
	  if(Current->outlet != NULL){
	    Current->outlet->sediment.inflow[i] += Current->sediment.outflow[i];
	    /* Needed for last time step to balance mass */
	    Total->ChannelSuspendedSediment += Current->sediment.outflow[i];
	  }
	  /* If no stream segment outlet, there is a road sink or a basin outlet.
	     Track this for the sediment mass balance. */
	  else{
	    Total->SedimentOutflow += Current->sediment.outflow[i];
	  }
	  /* Mass Balance Variables */
	  Total->DebrisInflow += Current->sediment.debrisinflow[i];
	  Current->sediment.debrisinflow[i] = 0.;
	  
	  Total->SedimentOverlandInflow += Current->sediment.overlandinflow[i];
	  Current->sediment.overlandinflow[i] = 0.;
	  
	  Total->ChannelSedimentStorage += Current->sediment.mass[i];
	} /* close loop for each sediment size */
	/* the next 7 lines are from channel_route_network -- closes the loop above */
	order_count += 1;
//fprintf(stdout,"order %d order_count %d id %d outflowconc %f\n",order,order_count,Current->id,Current->sediment.outflowconc);
      } /* close if statement checking for stream order */
      Current = Current->next;  
    } /* close while statement checking that CURRENT != NULL */
    if (order_count == 0)
      break;
  } /* close loop for the stream order */
}

/*
 * SUMMARY:      MainMWM.c - Mass Wasting Module
 * USAGE:        MWM
 *
 * AUTHOR:       Colleen O. Doten/Laura Bowling
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Sep-02
 * Last Change: Thu Jun 19 09:27:02 2003 by Ed Maurer <edm@u.washington.edu>
 * DESCRIPTION:  Main routine to drive MWM - the Mass Wasting Module for DHSVM 
 * DESCRIP-END.cd
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
#include "settings.h"
#include "Calendar.h"
#include "getinit.h"
#include "DHSVMerror.h"
#include "data.h"
#include "fileio.h"
#include "constants.h"
#include "DHSVMChannel.h"

#define BUFSIZE      255
#define empty(s) !(s)

float ElevationSlope(MAPSIZE *Map, FINEPIX ***FineMap, int y, int x, int *nexty, int *nextx, int prevy, int prevx, float *Aspect);
float Bagnold(float DS,TIMESTRUCT * Time, float, float, float, float);
float CalcSafetyFactor(float Slope, int Soil, float SoilDepth, int Veg, SEDTABLE *SedType, VEGTABLE *VType, float M, SOILTABLE *SType);
void enqueue(node **head, node **tail, int y, int x);
void dequeue(node **head, node **tail, int *y, int *x);
void RouteDebrisFlow(float *SedimentToChannel, int prevy, int prevx, float SlopeAspect, CHANNEL *ChannelData, MAPSIZE *Map); 
void InitChannelSediment(Channel * Head);
void InitChannelSedInflow(Channel * Head);
void DistributeSedimentDiams(float SedDiam[NSEDSIZES]);
void OutputChannelSediment(Channel * Head, TIMESTRUCT Time, DUMPSTRUCT *Dump);
void RouteChannelSediment(Channel * Head, Channel *RoadHead, TIMESTRUCT Time, DUMPSTRUCT *Dump);
void Alloc_Chan_Sed_Mem(float ** DummyVar);

/******************************************************************************/
/*			       MAIN                                    */
/******************************************************************************/
void MainMWM(SEDPIX **SedMap, FINEPIX *** FineMap, VEGTABLE *VType, SEDTABLE *SedType,
	    CHANNEL *ChannelData, char *DumpPath, SOILPIX **SoilMap, TIMESTRUCT *Time,
	    MAPSIZE *Map, TOPOPIX **TopoMap, SOILTABLE *SType, VEGPIX **VegMap,
	     int MaxStreamID) 
{
  int x,y,i,j,ii,jj,iter;  /* Counters. */
  int coursei, coursej;
  int nextx, nexty;
  int prevx, prevy;
  int numfailures;
  char buffer[32];
  char sumoutfile[100], outfile[100];  /* Character arrays to hold file name. */ 
  int **failure;
  float factor_safety;
  float LocalSlope;
  FILE *fs, *fo;          /* File pointers. */
  int numpixels;
  int cells, count;
  int massitertemp;       /* if massiter is 0, sets the counter to 1 here */
  float TotalVolume;
  node *head, *tail;
  float SlopeAspect, SedimentToChannel;
  float *SegmentSediment;  /* The cumulative sediment content over all stochastic
			      iterations for each channel segment. */
  float *InitialSegmentSediment; /* Placeholder of segment sediment load at 
				    beginning of time step. */
  float **SedThickness;    /* Cumulative sediment depth over all stochastic
			      iterations for each pixel.  */
  float **InitialSediment; /* Place holder of pixel sediment load at beginning of
			      time step. */
  float *SedDiams, *SedDist;
  head = NULL;
  tail = NULL;

  /*****************************************************************************
   Allocate memory 
  ****************************************************************************/

  if (!(failure = (int **)calloc(Map->NYfine, sizeof(int *))))
    ReportError("MainSEM", 1);
  for(i=0; i<Map->NYfine; i++) {
      if (!(failure[i] = (int *)calloc(Map->NXfine, sizeof(int))))
	ReportError("MainSEM", 1);
  }

  if (!(SedThickness = (float **)calloc(Map->NYfine, sizeof(float *))))
    ReportError("MainSEM", 1);
  for(i=0; i<Map->NYfine; i++) {
      if (!(SedThickness[i] = (float *)calloc(Map->NXfine, sizeof(float))))
	ReportError("MainSEM", 1);
  }

  if (!(InitialSediment = (float **)calloc(Map->NYfine, sizeof(float *))))
    ReportError("MainSEM", 1);
  for(i=0; i<Map->NYfine; i++) {
      if (!(InitialSediment[i] = (float *)calloc(Map->NXfine, sizeof(float))))
	ReportError("MainSEM", 1);
  }

  if (!(SegmentSediment = (float *)calloc(MaxStreamID, sizeof(float ))))
    ReportError("MainSEM", 1);

  if (!(SedDiams = (float *)calloc(NSEDSIZES, sizeof(float))))
    ReportError("MainSEM", 1);

  if (!(SedDist = (float *)calloc(NSEDSIZES, sizeof(float))))
    ReportError("MainSEM", 1);

  if (!(InitialSegmentSediment = (float *)calloc(NSEDSIZES, sizeof(float ))))
    ReportError("MainSEM", 1);

  /*****************************************************************************
  Perform Calculations 
  *****************************************************************************/
    
  /* Redistribute soil moisture from course grid to fine grid. 
     Currently assuming equal distribution; this will be replaced by Colleen's routine. */

  for (i = 0; i < Map->NY; i++) {
    for (j  = 0; j < Map->NX; j++) {

      if (INBASIN(TopoMap[i][j].Mask)) {
	for(ii=0; ii< Map->DY/Map->DMASS; ii++) {
	  for(jj=0; jj< Map->DX/Map->DMASS; jj++) {
	    y = (int) i*Map->DY/Map->DMASS + ii;
	    x = (int) j*Map->DX/Map->DMASS + jj;

	    (*FineMap)[y][x].SatThickness = SoilMap[i][j].Depth - SoilMap[i][j].TableDepth;
	  }
	}
      }
    }
  }

  /* Initialize arrays. */
  for (i = 0; i < Map->NYfine; i++) {
    for (j  = 0; j < Map->NXfine; j++) {
    
      InitialSediment[i][j] = (*FineMap)[i][j].sediment;
      (*FineMap)[i][j].probability = 0.;
    }
  }

  update_sediment_array(ChannelData->streams, InitialSegmentSediment);
	    
  /************************************************************************/
  /* Begin iteration for multiple ensembles. */
  /************************************************************************/

  /* sloppy fix for when MASSITER=0 -- this only needs to be checked once */
  if(MASSITER==0) massitertemp=1;
  else massitertemp=MASSITER;

  for(iter=0; iter < massitertemp; iter++) {

    fprintf(stderr,"iter=%d\n",iter);

    /************************************************************************/
    /* Begin factor of safety code. */
    /************************************************************************/
    for(i=0; i<Map->NY; i++) {
      for(j=0; j<Map->NX; j++) {

	/* Check to make sure region is in the basin. */
	if (INBASIN(TopoMap[i][j].Mask)) {		
		

	  /* Step over each fine resolution cell within the model grid cell. */
	  for(ii=0; ii<Map->DY/Map->DMASS; ii++) {
	    for(jj=0; jj<Map->DX/Map->DMASS; jj++) {
	      y = i*Map->DY/Map->DMASS + ii; /* Fine resolution counters. */
	      x = j*Map->DX/Map->DMASS + jj;
	      prevy = y;
	      prevx = x;
	      coursei = i;
	      coursej = j;

	      numpixels = 0;
	      SedimentToChannel = 0.0;
	      /* First check for original failure. */
	      if((*FineMap)[y][x].SatThickness/SoilMap[i][j].Depth > MTHRESH && failure[y][x] == 0) {

		LocalSlope = ElevationSlope(Map, FineMap, y, x, &nexty, &nextx, y, x, &SlopeAspect);

		if(LocalSlope >= 10. && LocalSlope <= 90.) {
		  factor_safety = CalcSafetyFactor(LocalSlope, SoilMap[i][j].Soil, 
						   (*FineMap)[y][x].sediment, 
						   VegMap[i][j].Veg, SedType, VType, 
						   (*FineMap)[y][x].SatThickness, SType);

		  /* check if fine pixel fails */
		  if (factor_safety < 1 && factor_safety > 0) {
		    numpixels = 1;
		    failure[y][x] = 1;	      
			  
		    if(!channel_grid_has_channel(ChannelData->stream_map, j, i)) {

		      /* Update sediment depth. All sediment leaves failed fine pixel */
		      SedThickness[nexty][nextx] += (*FineMap)[y][x].sediment;
		      (*FineMap)[y][x].sediment = 0.0;
		      //		      fprintf(stderr, "Original failure did not intersect a channel.\n");
		    }
		    else {
		      //	      fprintf(stderr, "Original failure intersected a channel:i=%d, j=%d, y=%d, x=%d\n", i, j, y, x);
		    }
		  }
		  
		}		      
	      }

	      /* Follow failures down slope; skipped if no original failure. */
	      while(failure[y][x] == 1 && numpixels >= 1 && !channel_grid_has_channel(ChannelData->stream_map, coursej, coursei)) {
		      
		/* Update counters. */
		prevy = y;
		prevx = x;
		y = nexty;
		x = nextx;
		coursei = floor(y*Map->DMASS/Map->DY);
		coursej = floor(x*Map->DMASS/Map->DX);

		LocalSlope = ElevationSlope(Map, FineMap, y, x, &nexty, &nextx, prevy, prevx, &SlopeAspect);

		/*  Check that not a sink and that we have not encountered 
		    a stream segment. */
		if(LocalSlope >= 0) {
	
		  factor_safety = CalcSafetyFactor(LocalSlope, SoilMap[coursei][coursej].Soil, 
						   (*FineMap)[y][x].sediment, 
						   VegMap[coursei][coursej].Veg, SedType, VType, 
						   (*FineMap)[y][x].SatThickness, SType);

			
		  if (factor_safety < 1 && factor_safety > 0) {
		    numpixels += 1;
		    failure[y][x] = 1;
	
		    if(!channel_grid_has_channel(ChannelData->stream_map, coursej, coursei)) {

		      /* Update sediment depth. All sedminet leaves failed fine pixel */
		      SedThickness[nexty][nextx] += (*FineMap)[y][x].sediment;
		      (*FineMap)[y][x].sediment = 0.0;
		      //      fprintf(stderr, "Secondary failure did not encounter a channel!\n");
		    }
		    else {
		      //      fprintf(stderr, "Secondary failure encountered a channel!\n");
		    }
		  }
		}
	      }  /* End of while loop. */
	      
	      /* Failure has stopped, now calculate runout distance and 
		 redistribute sediment. */
	      if(numpixels >= 1) {

		TotalVolume = (*FineMap)[y][x].sediment;
		(*FineMap)[y][x].sediment = 0.;
		cells = 1;

		/* queue begins with initial unfailed pixel. */
		enqueue(&head, &tail, y, x); 

		while(LocalSlope > 4. && !channel_grid_has_channel(ChannelData->stream_map, coursej, coursei)) {
		  /* Redistribution stops if last pixel was a channel. */
		  /* Update counters. */
		  prevy = y;
		  prevx = x;
		  y = nexty;
		  x = nextx;
		  coursei = floor(y*Map->DMASS/Map->DY);
		  coursej = floor(x*Map->DMASS/Map->DX);

		  LocalSlope = ElevationSlope(Map, FineMap, y, x, &nexty, &nextx, prevy, prevx, &SlopeAspect);
			
		  enqueue(&head, &tail, y, x);
		  cells++;
		}
		prevy = y;
		prevx = x;

		for(count=0; count < cells; count++) {
		  dequeue(&head, &tail, &y, &x);
		   coursei = floor(y*Map->DMASS/Map->DY);
		  coursej = floor(x*Map->DMASS/Map->DX);

		  if(channel_grid_has_channel(ChannelData->stream_map, coursej, coursei) && count ==0) {
		    /* TotalVolume at this point is a depth in m over one fine
		       map grid cell - convert to m3 */
		    SedimentToChannel = TotalVolume*(Map->DMASS*Map->DMASS)/(float) cells;
		  }
		  else {
		    /* Redistribute sediment equally among all hillslope cells. */
		    (*FineMap)[y][x].sediment += TotalVolume/cells;
		  }
		}
		      
		if(SedimentToChannel > 0.0) {
		  if(SlopeAspect < 0.) {
		    fprintf(stderr, "No Valid aspect.\n");
		    exit(0);
		  }
		  else
		    RouteDebrisFlow(&SedimentToChannel, coursei, coursej, SlopeAspect, ChannelData, Map); 
		}
		
	      }
	    }  /* End of jj loop. */
	  }
	}
	
	
      }       
    }    /* End of course resolution loop. */

    /* Record failures and Reset failure map for new iteration. */
    for (i = 0; i < Map->NYfine; i++) {
      for (j  = 0; j < Map->NXfine; j++) {

	y = (int) floor(i/(Map->DY/Map->DMASS));
	x = (int) floor(j/(Map->DX/Map->DMASS));

	if (INBASIN(TopoMap[y][x].Mask)) {
	  (*FineMap)[i][j].probability += (float) failure[i][j];
	  
	  /* Record cumulative sediment volume. */
	  SedThickness[i][j] += (*FineMap)[i][j].sediment;
	  
	  /* Reset sediment thickness for each iteration; otherwise there is 
	     a decreasing probability of failure for subsequent iterations. */
	  /* If not in stochastic mode, then allow a history of past failures
	     and do not reset sediment depth */
	  if(massitertemp>1) {
	    (*FineMap)[i][j].sediment = InitialSediment[i][j];
	    failure[i][j] = 0;
	  }
	}
	else {
	  (*FineMap)[i][j].probability = -999.;
	  (*FineMap)[i][j].DeltaDepth = -999.;
	}
      }
    }

    /* Record cumulative stream sediment volumes. */
    initialize_sediment_array(ChannelData->streams, SegmentSediment);

    /* Reset channel sediment volume for each iteration. */
    update_sediment_array(ChannelData->streams, InitialSegmentSediment);
   
    
  }    /* End iteration loop */

  /*************************************************************************/
  /* Create output files...currently hard-coded, should be moved to dump   */
  /* functions for user specification. Creates the following files in the  */
  /* specified output directory:                                           */ 
  /* failure_summary.txt - no. of failed pixels for each date that the mwm */
  /*                       algorithm is run.  Failed pixel = probability   */
  /*                       of failure > 0.5                                */ 
  /* "date"_failure.txt - map of probability of failure for each pixel     */
  /* "date"_Deltasoildepth.txt - map of cumulative change in soil depth    */
  /*                             since beginning of model run.             */
  /*************************************************************************/

  numfailures = 0;
  for (i = 0; i < Map->NYfine; i++) {
    for (j  = 0; j < Map->NXfine; j++) {

      	y = (int) floor(i/(Map->DY/Map->DMASS));
	x = (int) floor(j/(Map->DX/Map->DMASS));

	if (INBASIN(TopoMap[y][x].Mask)) {

	  (*FineMap)[i][j].probability /= (float)massitertemp;
	  (*FineMap)[i][j].sediment = SedThickness[i][j]/(float)massitertemp;

	  if((*FineMap)[i][j].probability > .5)
	    numfailures +=1;
	  
	  (*FineMap)[i][j].DeltaDepth = (*FineMap)[i][j].sediment - SoilMap[y][x].Depth;
	}
	else {
	  (*FineMap)[i][j].probability = OUTSIDEBASIN;
	  (*FineMap)[i][j].DeltaDepth = OUTSIDEBASIN;
	  (*FineMap)[i][j].sediment = OUTSIDEBASIN;
	}
    }
  }

  /*average sediment delivery to each stream segment***************/
  for(i=0; i<MaxStreamID; i++) {
    SegmentSediment[i] /= (float)massitertemp;
    if(SegmentSediment[i] < 0.0) SegmentSediment[i]=0.0;
    /*    fprintf(stderr, "segment %d: segmentsediment = %f\n", i, 
	  SegmentSediment[i]-InitialSegmentSediment[i]); */
  }
  update_sediment_array(ChannelData->streams, SegmentSediment);
  /* Take new sediment inflow and distribute it by representative diameters*/
  /* and convert to mass */
  sed_vol_to_distrib_mass(ChannelData->streams, SegmentSediment);
  /* back to map analysis******************************************/

  sprintf(sumoutfile, "%sfailure_summary.txt", DumpPath);
  
  if((fs=fopen(sumoutfile,"a")) == NULL)
    {
      fprintf(stderr,"Cannot open factor of safety summary output file.\n");
      exit(0);
    }

  SPrintDate(&(Time->Current), buffer);
  fprintf(fs, "%-20s %7d\n", buffer, numfailures); 
  fprintf(stderr, "%d pixels failed\n", numfailures);
  fclose(fs);

  /* If any pixels failed, output map of failure probabilities & deltasoildepth. 
     Note: no failures does not mean that the probability of failure is zero. */ 
     
  if(numfailures > 0) {
    sprintf(outfile, "%s%s_failure.txt", DumpPath, buffer);

    if((fo=fopen(outfile,"w")) == NULL)
      {
	fprintf(stderr,"Cannot open factor of safety output file.\n");
	exit(0);
      }
	
    /* Printing header to output file */
    fprintf(fo,"ncols %11d\n",Map->NXfine);
    fprintf(fo,"nrows %11d\n",Map->NYfine);
    fprintf(fo,"xllcorner %.1f\n",Map->Xorig);
    fprintf(fo,"yllcorner %.1f\n",Map->Yorig - Map->NY*Map->DY);
    fprintf(fo,"cellsize %.0f\n",Map->DMASS);
    fprintf(fo,"NODATA_value %.2f\n",-999.);

    for (i = 0; i < Map->NYfine; i++) {
      for (j  = 0; j < Map->NXfine; j++) {
      
	y = (int) floor(i/(Map->DY/Map->DMASS));
	x = (int) floor(j/(Map->DX/Map->DMASS));

	/* Check to make sure region is in the basin. */
	if (INBASIN(TopoMap[y][x].Mask)) 		
	  fprintf(fo, "%.2f ", (*FineMap)[i][j].probability);
	else
	  fprintf(fo, "-999. ");
  
      }
      fprintf(fo, "\n");
    }
	
  sprintf(sumoutfile, "%s%s_Deltasoildepth.txt", DumpPath, buffer);
  if((fs=fopen(sumoutfile,"w")) == NULL)
    {
      fprintf(stderr,"Cannot open soil depth output file.\n");
      exit(0);
    }

  /* Printing header to output file */
  fprintf(fs,"ncols %11d\n",Map->NXfine);
  fprintf(fs,"nrows %11d\n",Map->NYfine);
  fprintf(fs,"xllcorner %.1f\n",Map->Xorig);
  fprintf(fs,"yllcorner %.1f\n",Map->Yorig - Map->NY*Map->DY);
  fprintf(fs,"cellsize %.0f\n",Map->DMASS);
  fprintf(fs,"NODATA_value %.2f\n",-999.);

  for (i = 0; i < Map->NYfine; i++) {
    for (j  = 0; j < Map->NXfine; j++) {
      
      y = (int) floor(i/(Map->DY/Map->DMASS));
      x = (int) floor(j/(Map->DX/Map->DMASS));

      /* Check to make sure region is in the basin. */
      if (INBASIN(TopoMap[y][x].Mask)) 		
	fprintf(fs, "%.2f ", (*FineMap)[i][j].DeltaDepth);
      else
	fprintf(fs, "-999. ");
    }
    fprintf(fs, "\n");
  }
  /* Close files. */
  fclose(fo);
  fclose(fs);
  }

  for(i=0; i<Map->NYfine; i++) { 
    free(failure[i]);
    free(SedThickness[i]);
    free(InitialSediment[i]);
  }
  free(failure);
  free(SedThickness);
  free(InitialSediment);
  free(SegmentSediment);
  }

/*****************************************************************************
  End of Main
*****************************************************************************/

  void enqueue(node **head, node **tail, int y, int x)
{
  node *new;
  new = (node *) malloc(sizeof(node));
  new->x = x;
  new->y = y;
  new->next = NULL;

  if(empty(*head)) {
    *head = new;
    *tail = new;
  }
  else {
    //    if(!empty(*tail)) {
    //   fprintf(stderr,"New node is not at end of queue\n");
    //  exit(0);
    //   }
    (*tail)->next = new;
    *tail = new;
  }
}

void dequeue(node **head, node **tail, int *y, int *x)
{
  //  if(!empty(*head))
  //    {
  //      fprintf(stderr,"Node is not at head of queue\n");
  //      exit(0);
  //      }
  *y = (*head)->y;
  *x = (*head)->x;
  *head = (*head)->next;
  if(head == NULL) tail = NULL;
}

/*****************************************************************************
  Alloc_Chan_Sed_Mem
*****************************************************************************/
void Alloc_Chan_Sed_Mem(float ** DummyVar)
{
  printf("in sub\n");
  if (!(*DummyVar = (float *) calloc(NSEDSIZES, sizeof(float))))
    ReportError("sed_chan_alloc", 1);
  printf("did it\n");
}
/*****************************************************************************
  InitChannelSediment)

  Assign initial colluvium mass to each unique channel ID (amount
  of storage, kg)
*****************************************************************************/
void InitChannelSediment(Channel * Head)
{
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
  DistributeSedimentDiams()

  For new lateral sediment inflow, Find the particle diameters for each portion
  Assumes a lognormal distribution
*****************************************************************************/
void DistributeSedimentDiams(float SedDiam[NSEDSIZES])
{
  /* to be consistent with FindValue, use the Tukey (1960) approx to normal
     CDF -- y is probability, function returns Z(y) */
#ifndef NORMALDIST
#define NORMALDIST(mean, stdev, y) (4.91 * stdev * (pow(y,.14) - pow(( 1 - y ),.14)) + mean )
#endif

  int i;
  float mn,std,z;
  float pctfiner;

  /* slope of the lognormal curve, log(d) on y-axis, Z on x-axis, is stdev */
  mn = log10(DEBRISd50);
  std = log10(DEBRISd90)-log10(DEBRISd50)/(NORMALDIST(0,1,0.9)-NORMALDIST(0,1,0.5));

  pctfiner = 1.0/(2.0*NSEDSIZES); /* midpoint of finest interval */

  for(i=0;i<NSEDSIZES;i++) {
    z = NORMALDIST(mn,std,pctfiner);
    SedDiam[i] = pow(10,mn+std*z);
    pctfiner += 1.0/NSEDSIZES;
  }
}
/*****************************************************************************
  OutputChannelSediment()

  Output colluvium volumes (m3) for each unique channel ID for each time step.
*****************************************************************************/
void OutputChannelSediment(Channel * Head, TIMESTRUCT Time, DUMPSTRUCT *Dump)
{
  Channel *Current = NULL;
  char buffer[20], FileName[100];
  int i;
  FILE *fo;

  SPrintDate(&(Time.Current), buffer);
  sprintf(FileName, "%s/Channel.sediment.%s.asc", Dump->Path, buffer);

  if((fo=fopen(FileName,"w")) == NULL)
    ReportError("OutputChannelSediment", 3);
  
  /* Assign the storages to the correct IDs */
  Current = Head;
  while (Current) {
    fprintf(fo, "%d\t", Current->id);
    for(i=0;i<NSEDSIZES;i++) fprintf(fo, "%.3f\t", Current->sediment.mass[i]);
    fprintf(fo, "\n");
    Current = Current->next;
  }
  fclose(fo);

}

/*****************************************************************************
  RouteChannelSediment()

  Read in DHSVM sediment mass and inflows for each channel segment, and 
  route sediment downstream. Sorts by particle size, transports finer material
  first, as done by Williams (1980).

*****************************************************************************/
void RouteChannelSediment(Channel * Head, Channel *RoadHead, TIMESTRUCT Time, DUMPSTRUCT *Dump)
{
  Channel *Current = NULL;
  float DS,DT_sed,Q,numinc;
  float flowdepth,Qavg,V,dIdt,dOdt,dMdt,Vsed,Vshear,Vshearcrit;
  float minDT_sed,TotalCapacityUp,TotalCapacityDown;
  float lateral_sed_inflow_rate;
  float TotalCapacity, CapacityUsed;
  float SedDiam[NSEDSIZES];
  float Qup,Qdown;
  float phi=0.55, theta=0.55,term3,term4; /*space and time weighting factors*/
  int i,tstep;
  int order;
  int order_count;

  /* For each of the sediment diameters, calculate the mass balance */
  DistributeSedimentDiams(SedDiam); /* find diameter for each portion */

  /* the next 5 lines are from channel_route_network - used to order streams */
  for (order = 1;; order += 1) {
    order_count = 0;
    Current = Head;

    while (Current != NULL) {
      if (Current->order == order) {

	Current->sediment.outflowconc=0.0;
	CapacityUsed = 0.0;
	Current->outlet->sediment.totalmass=0;
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
	  V = Q/(flowdepth*Current->class->width);
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
	  DS = SedDiam[i]*((float) MMTOM); /* convert from mm to m */

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
	    if(SedDiam[i] < 0.062) { /* per Wicks and Bathurst, wash load */
	      TotalCapacity = 
		Current->sediment.inflowrate[i]+Current->sediment.mass[i]/DT_sed;
	    }
	    else {
	      TotalCapacityUp = Bagnold(DS,&Time,Qup,Current->class->width,
					Current->class->friction,Current->slope);
	      TotalCapacityDown = Bagnold(DS,&Time,Qdown,Current->class->width,
					  Current->class->friction,Current->slope);
	      TotalCapacity=phi*TotalCapacityDown + (1.0-phi)*TotalCapacityUp;
	      TotalCapacity -= CapacityUsed; /* Avoid mult use of streampower */
	    }
	    if(TotalCapacity<=0) TotalCapacity=0.0;
	    if(TotalCapacity*DT_sed >= Current->sediment.mass[i]) {
	      dMdt = Current->sediment.mass[i]/DT_sed;
	      Current->sediment.mass[i] = 0;
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
	  if (Current->outflow==0.0)
			Current->sediment.outflowconc = 0.0;
	  else {
		if(Current->slope>0.0) 
			flowdepth = pow(Current->outflow*Current->class->friction/(Current->class->width*sqrt(Current->slope)),0.6);
		else flowdepth = 0.005;
		V = Current->outflow/(flowdepth*Current->class->width);
		Vshear = sqrt(G*flowdepth*Current->slope); /* approximate */
		if(SedDiam[i] < 15) Vshearcrit = -0.0003*SedDiam[i]*SedDiam[i]+0.0109*SedDiam[i]+0.0106;
		else Vshearcrit = 0.0019*SedDiam[i]+0.0926; /* both in m/s */
		if (Vshearcrit > Vshear) Vshearcrit = 0.95*Vshear; /*arbitrary */
		/* sed velocity per Wicks and Bathurst, fine particles travel at flow velocity */
		Vsed = 8.5*Vshear*sqrt(1-(Vshearcrit/Vshear));
		if(Vsed>V || SedDiam[i]<=0.062) Vsed=V;
		if(Vsed<(0.2*V)) Vsed=0.2*V;
		Current->sediment.outflowconc += 
			(1000.0*Current->sediment.outflow[i]/(Current->outflow))*(V/Vsed);
	  }
	  
	  /* pass the sediment mass outflow to the next downstream reach */
	  if(Current->outlet != NULL)
	    Current->outlet->sediment.inflow[i] += Current->sediment.outflow[i];
	  
	} /* close loop for each sediment size */
	
      } /* close if statement checking for stream order */
      
	/* the next 7 lines are from channel_route_network -- closes the loop above */
      order_count += 1;
    } /* close while statement checking that CURRENT != NULL */
    Current = Current->next;

    if (order_count == 0)
      break;
  } /* close loop for the stream order */
}


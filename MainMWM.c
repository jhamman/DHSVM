/*
 * SUMMARY:      MainMWM.c - Mass Wasting Module
 * USAGE:        MWM
 *
 * AUTHOR:       Colleen O. Doten/Laura Bowling
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Sep-02
 * Last Change:  
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

float CalcSafetyFactor(float Slope, int Soil, float SoilDepth, int Veg, SEDTABLE *SedType, VEGTABLE *VType, float M, SOILTABLE *SType);
void enqueue(node **head, node **tail, int y, int x);
void dequeue(node **head, node **tail, int *y, int *x);
void RouteDebrisFlow(float *SedimentToChannel, int prevy, int prevx, float SlopeAspect, 
		     CHANNEL *ChannelData, MAPSIZE *Map); 
void InitChannelSediment(Channel * Head);
void OutputChannelSediment(Channel * Head, TIMESTRUCT Time, DUMPSTRUCT *Dump);
void RouteChannelSediment(Channel * Head, Channel *RoadHead, TIMESTRUCT Time, DUMPSTRUCT *Dump);

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

  if (!(SegmentSediment = (float *)calloc(MaxStreamID, sizeof(int))))
    ReportError("MainSEM", 1);

   if (!(InitialSegmentSediment = (float *)calloc(MaxStreamID, sizeof(int))))
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

  initialize_sediment_array(ChannelData->streams, InitialSegmentSediment);
	    
  /************************************************************************/
  /* Begin iteration for multiple ensembles. */
  /************************************************************************/
  for(iter=0; iter < MASSITER; iter++) {

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

		  if (factor_safety < 1 && factor_safety > 0) {
		    numpixels = 1;
		    failure[y][x] = 1;	      
			  
		    if(!channel_grid_has_channel(ChannelData->stream_map, j, i)) {

		      /* Update sediment depth. */
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

		      /* Update sediment depth. */
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
		    SedimentToChannel = TotalVolume/cells;
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
	(*FineMap)[i][j].sediment = InitialSediment[i][j];
	failure[i][j] = 0;
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

	  (*FineMap)[i][j].probability /= (float)MASSITER;
	  (*FineMap)[i][j].sediment = SedThickness[i][j]/(float)MASSITER;

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

  for(i=0; i<MaxStreamID; i++) {
    SegmentSediment[i] /= (float)MASSITER;
    //    fprintf(stderr, "segment %d: segmentsediment = %f\n", i, SegmentSediment[i]-InitialSegmentSediment[i]);
  }
  update_sediment_array(ChannelData->streams, SegmentSediment);

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
  InitChannelSediment)

  Assign initial colluvium volumes to each unique channel ID (amount
  of storage (m3))
*****************************************************************************/
void InitChannelSediment(Channel * Head)
{
  Channel *Current = NULL;
  float InitialDepth = 0.25;

  /* Assign the storages to the correct IDs */
  Current = Head;
  while (Current) {
    Current->sediment = Current->length * InitialDepth * Current->class->width;
    Current = Current->next;
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
  FILE *fo;

  SPrintDate(&(Time.Current), buffer);
  sprintf(FileName, "%s/Channel.sediment.%s.asc", Dump->Path, buffer);

  if((fo=fopen(FileName,"w")) == NULL)
    ReportError("OutputChannelSediment", 3);
  
  /* Assign the storages to the correct IDs */
  Current = Head;
  while (Current) {
    fprintf(fo, "%d\t%f\n", Current->id, Current->sediment);
    Current = Current->next;
  }
  fclose(fo);

}

/*****************************************************************************
  RouteChannelSediment()

  Read in DHSVM discharge volumes for each channel segment, and route sediment 
  downstream.
.
*****************************************************************************/
void RouteChannelSediment(Channel * Head, Channel *RoadHead, TIMESTRUCT Time, DUMPSTRUCT *Dump)
{
  Channel *Current = NULL;
  char buffer[20], FileName[100];
  FILE *fo;

  SPrintDate(&(Time.Current), buffer);
  sprintf(FileName, "%s/Channel.sediment.%s.asc", Dump->Path, buffer);

  if((fo=fopen(FileName,"w")) == NULL)
    ReportError("OutputChannelSediment", 3);
  
  /* Assign the storages to the correct IDs */
  Current = Head;
  while (Current) {
    fprintf(fo, "%d\t%f\n", Current->id, Current->sediment);
    Current = Current->next;
  }
  fclose(fo);

} 

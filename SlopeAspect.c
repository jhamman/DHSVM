/*
 * SUMMARY:      SlopeAspect.c - Calculate slope and aspect of each pixel
 * USAGE:        Part of DHSVM/MWM
 *
 * AUTHOR:       William A Perkins
 * ORG:          Battelle Memorial Institute Pacific Northwest Laboratory
 * E-MAIL:       perk@clio.muse.pnl.gov
 * ORIG-DATE:    21-May-96
 * DESCRIPTION:  This module contains two routines to compute "slope" and
 *               "aspect"  direction of slope): one which uses only terrain
 *               elevations and another which uses water table elevations.
 * DESCRIP-END.
 * FUNCTIONS:    valid_cell()
 *               valid_cell_fine()
 *               slope_aspect()
 *               flow_fractions()
 *               ElevationSlopeAspect()
 *               HeadSlopeAspect()
 *               ElevationSlope()
 *               ElevationSlopeAspectfine()
 * COMMENTS:
 * $Id$     
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"

/* These indices are so neighbors can be looked up quickly */

int xneighbor[NDIRS] = {
#if NDIRS == 4
  0, 1, 0, -1
#elif NDIRS == 8
  -1, 0, 1, 1, 1, 0, -1, -1
#endif
};
int yneighbor[NDIRS] = {
#if NDIRS == 4
  -1, 0, 1, 0
#elif NDIRS == 8
  1, 1, 1, 0, -1, -1, -1, 0
#endif
};

float temp_aspect[NDIRSfine] = {
#if  NDIRSfine == 4
  180., 270., 0., 90.
#elif NDIRSfine == 8
 135., 180., 225., 270., 315., 0., 45., 90.
#endif
};

int xneighborfine[NDIRSfine] = {
#if NDIRSfine == 4
  0, 1, 0, -1
#elif NDIRSfine == 8
  -1, 0, 1, 1, 1, 0, -1, -1
#endif
};
int yneighborfine[NDIRSfine] = {
#if NDIRSfine == 4
  -1, 0, 1, 0
#elif NDIRSfine == 8
  1, 1, 1, 0, -1, -1, -1, 0
#endif
};

/* -------------------------------------------------------------
   valid_cell
   Checks to see if grid indices, x and y, are within the grid 
   defined by the specified Map
   ------------------------------------------------------------- */
int valid_cell(MAPSIZE * Map, int x, int y)
{
  return (x >= 0 && y >= 0 && x < Map->NX && y < Map->NY);
}

/******************************************************************************/
/*   valid_cell_fine                                                          */
/*   Checks to see if grid indices, x and y, are within the grid              */
/*   defined by the specified Map                                             */ 
/******************************************************************************/

int valid_cell_fine(MAPSIZE *Map, int x, int y) 
{
  return (x >= 0 && y >= 0 && x < Map->NXfine && y < Map->NYfine);
}


/* -------------------------------------------------------------
   slope_aspect
   Calculation of slope and aspect given elevations of cell and neighbors
   ------------------------------------------------------------- */
static void slope_aspect(float dx, float dy, float celev, float
			 nelev[NDIRS], float *slope, float *aspect)
{
  int n;
  float dzdx, dzdy;
  switch (NDIRS) {
  case 8:
    /* for eight neighbors, this is
       exactly the same algorithm that
       Arc/Info uses */

    for (n = 0; n < NDIRS; n++) {
      if (nelev[n] == OUTSIDEBASIN) {
	nelev[n] = celev;
      }
    }

    dzdx = ((nelev[7] + 2 * nelev[6] + nelev[5]) -
	    (nelev[1] + 2 * nelev[2] + nelev[3])) / (8 * dx);
    dzdy = ((nelev[7] + 2 * nelev[0] + nelev[1]) -
	    (nelev[5] + 2 * nelev[4] + nelev[3])) / (8 * dy);
    break;

  case 4:
    if (nelev[1] == (float) OUTSIDEBASIN && nelev[3] == (float) OUTSIDEBASIN) {
      dzdx = 0.0;
    }
    else if (nelev[1] == (float) OUTSIDEBASIN) {
      dzdx = (nelev[3] - celev) / dx;
    }
    else if (nelev[3] == (float) OUTSIDEBASIN) {
      dzdx = (celev - nelev[1]) / dx;
    }
    else {
      dzdx = (nelev[3] - nelev[1]) / (2 * dx);
    }

    if (nelev[0] == (float) OUTSIDEBASIN && nelev[2] == (float) OUTSIDEBASIN) {
      dzdy = 0.0;
    }
    else if (nelev[2] == (float) OUTSIDEBASIN) {
      dzdy = (celev - nelev[0]) / dy;
    }
    else if (nelev[0] == (float) OUTSIDEBASIN) {
      dzdy = (nelev[2] - celev) / dy;
    }
    else {
      dzdy = (nelev[2] - nelev[0]) / (2 * dy);
    }
    break;
  default:
    ReportError("slope_aspect",65);
    assert(0);			/* nothing else works */
    break;
  }

  *slope = sqrt(dzdx * dzdx + dzdy * dzdy);
  if (fequal(dzdx, 0.0) && fequal(dzdy, 0.0)) {
    *aspect = 0.0;
  }
  else {
    *aspect = atan2(dzdx, dzdy);
  }

/*  if (*slope>2.0) printf("slope >2 = %g\n",*slope); */
  return;
}

/* -------------------------------------------------------------
   flow_fractions
   Computes subsurface flow fractions given the slope and aspect 
------------------------------------------------------------- */
static void flow_fractions(float dx, float dy, float slope, float aspect,
			   float nelev[NDIRS], float *grad,
			   unsigned char dir[NDIRS], unsigned int *total_dir)
{
  float cosine = cos(aspect);
  float sine = sin(aspect);
  float total_width, effective_width;
  int n;

  switch (NDIRS) {
  case 4:

    /* fudge any cells which flow outside
       the basin by just pointing the
       aspect in the opposite direction */

    if ((cosine > 0 && nelev[0] == (float) OUTSIDEBASIN) ||
	(cosine < 0 && nelev[2] == (float) OUTSIDEBASIN))
      cosine = -cosine;
    if ((sine > 0 && nelev[1] == (float) OUTSIDEBASIN) ||
	(sine < 0 && nelev[3] == (float) OUTSIDEBASIN))
      sine = -sine;

    /* compute flow widths */

    total_width = fabs(cosine) * dx + fabs(sine) * dy;
    *grad = slope * total_width;
    *total_dir = 0;
    for (n = 0; n < NDIRS; n++) {
      switch (n) {
      case 0:
	effective_width = (cosine > 0 ? cosine * dx : 0.0);
	break;
      case 2:
	effective_width = (cosine < 0 ? -cosine * dx : 0.0);
	break;
      case 1:
	effective_width = (sine > 0 ? sine * dx : 0.0);
	break;
      case 3:
	effective_width = (sine < 0 ? -sine * dx : 0.0);
	break;
      default:
	ReportError("flow_fractions",65);
	assert(0);		/* How can this happen? */
      }
      dir[n] = (int) ((effective_width / total_width) * 255.0 + 0.5);
      *total_dir += dir[n];
    }
    break;
  case 8:
    ReportError("flow_fractions",65);
    assert(0);			/* can't handle this */
    break;
  default:
    ReportError("flow_fractions",65);
    assert(0);			/* other cases don't work either */
  }
  return;
}

/* -------------------------------------------------------------
   ElevationSlopeAspect
   ------------------------------------------------------------- */
void ElevationSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap)
{
  const char *Routine = "ElevationSlopeAspect";
  int x;
  int y;
  int n;
  int k;
  float neighbor_elev[NDIRS];

  /* fill neighbor array */
  
  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (TopoMap[y][x].Mask) {
	
	/* Count the number of cells in the basin.  Need this to allocate memory for
	   the new, smaller Elev[] and Coords[][].  */
	if (INBASIN(TopoMap[y][x].Mask)) 
	  Map->NumCells++;
	
	for (n = 0; n < NDIRS; n++) {
	  int xn = x + xneighbor[n];
	  int yn = y + yneighbor[n];
	  
	  if (valid_cell(Map, xn, yn)) {
	    neighbor_elev[n] = ((TopoMap[yn][xn].Mask) ? TopoMap[yn][xn].Dem : (float) OUTSIDEBASIN);
	  }
	  else {
	    neighbor_elev[n] = OUTSIDEBASIN;
	  }
	}
	
	slope_aspect(Map->DX, Map->DY, TopoMap[y][x].Dem, neighbor_elev,
		     &(TopoMap[y][x].Slope), &(TopoMap[y][x].Aspect));
	
	/* fill Dirs in TopoMap too */
	
	flow_fractions(Map->DX, Map->DY, TopoMap[y][x].Slope,
		       TopoMap[y][x].Aspect,
		       neighbor_elev, &(TopoMap[y][x].FlowGrad),
		       TopoMap[y][x].Dir, &(TopoMap[y][x].TotalDir));
	
      } // end of if statement (TopoMap[y][x].Mask)
      
    }
  } // end of for loops
  
  /* Create a structure to hold elevations of only those cells
     within the basin and the y,x of those cells.*/
 
  if (!(Map->OrderedCells = (ITEM *) calloc(Map->NumCells, sizeof(ITEM))))
    ReportError((char *) Routine, 1);
  
  k = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      /* Save the elevation, y, and x in the ITEM structure. */
      if (INBASIN(TopoMap[y][x].Mask)) {
        Map->OrderedCells[k].Rank = TopoMap[y][x].Dem;
        Map->OrderedCells[k].y = y;
        Map->OrderedCells[k].x = x;
        k++;
      }
    }
  }

  /* Sort Elev in descending order-- Elev.x and Elev.y hold indices. */
  
  quick(Map->OrderedCells, Map->NumCells);

  /* End of modifications to create ordered cell coordinates.  SRW 10/02, LCB 03/03 */

  return;
}

/* -------------------------------------------------------------
   QuickSort
   ------------------------------------------------------------- */

/**********************************************************************
        this subroutine starts the quick sort
**********************************************************************/

void quick(ITEM *OrderedCells, int count)
{
  qs(OrderedCells,0,count-1);
}

void qs(ITEM *item, int left, int right)
/**********************************************************************
        this is the quick sort subroutine - it returns the values in
        an array from high to low.
**********************************************************************/
{
  register int i,j;
  ITEM x,y;

  i=left;
  j=right;
  x=item[(left+right)/2];

  do {
    while(item[i].Rank<x.Rank && i<right) i++;
    while(x.Rank<item[j].Rank && j>left) j--;

    if (i<=j) {
      y=item[i];
      item[i]=item[j];
      item[j]=y;
      i++;
      j--;
    }
  } while (i<=j);

  if(left<j) qs(item,left,j);
  if(i<right) qs(item,i,right);

}

/* -------------------------------------------------------------
   HeadSlopeAspect
   This computes slope and aspect using the water table elevation. 
   ------------------------------------------------------------- */
void HeadSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SOILPIX ** SoilMap)
{
  int x;
  int y;
  int n;
  float neighbor_elev[NDIRS];

  /* let's assume for now that WaterLevel is the SOILPIX map is
     computed elsewhere */

  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (TopoMap[y][x].Mask) {

	float slope, aspect;

	for (n = 0; n < NDIRS; n++) {
	  int xn = x + xneighbor[n];
	  int yn = y + yneighbor[n];

	  if (valid_cell(Map, xn, yn)) {
	    neighbor_elev[n] =
	      ((TopoMap[yn][xn].Mask) ? SoilMap[yn][xn].WaterLevel :
	       OUTSIDEBASIN);
	  }
	  else {
	    neighbor_elev[n] = OUTSIDEBASIN;
	  }
	}

	slope_aspect(Map->DX, Map->DY, SoilMap[y][x].WaterLevel, neighbor_elev,
		     &slope, &aspect);
	flow_fractions(Map->DX, Map->DY, slope, aspect, neighbor_elev,
		       &(TopoMap[y][x].FlowGrad), TopoMap[y][x].Dir,
		       &(TopoMap[y][x].TotalDir));
      }
    }
  }
  return;
}

/******************************************************************************/
/*			     ElevationSlope                            */
/* Part of MWM, should probably be merged w/ ElevationSlopeAspect function.   */
/******************************************************************************/

float ElevationSlope(MAPSIZE *Map, FINEPIX ***FineMap, int y, int x, int *nexty, 
		     int *nextx, int prevy, int prevx, float *Aspect) 
{
  int n;
  float neighbor_elev[NDIRSfine];
  float Slope;
  float temp_slope[NDIRSfine];
  double length_diagonal;
  float dx, dy, celev;
  float elevsink=0.0;

  /* fill neighbor array */
  
  for (n = 0; n < NDIRSfine; n++) {
    int xn = x + xneighborfine[n];
    int yn = y + yneighborfine[n];
          
    if (valid_cell_fine(Map, xn, yn)) 
      neighbor_elev[n] = (*FineMap)[yn][xn].bedrock + (*FineMap)[yn][xn].sediment;
    else 
      neighbor_elev[n] = OUTSIDEBASIN;
  }
        
  dx = Map->DMASS;
  dy = Map->DMASS;
  celev = (*FineMap)[y][x].bedrock + (*FineMap)[y][x].sediment;

  length_diagonal = sqrt((pow(dx, 2)) + (pow(dy, 2))); 

  switch (NDIRSfine){
  case 8:     
    for (n = 0; n < NDIRSfine; n++) {
      if (neighbor_elev[n] == OUTSIDEBASIN) 
	neighbor_elev[n] = celev;
      if(n==0 || n==2 || n==4 || n==6)
	temp_slope[n] = (atan((celev - neighbor_elev[n]) / length_diagonal))
	  * DEGPRAD;
      else if(n==1 || n==5)
	temp_slope[n] = (atan((celev - neighbor_elev[n]) / dy)) * DEGPRAD;
      else
	temp_slope[n] = (atan((celev - neighbor_elev[n]) / dx)) * DEGPRAD;
    }
    break;
  
  case 4:
    ReportError("ElevationSlope",65);
    assert(0); 
    
    /* Removes sinks since preprocessing does not remove sinks in 4 directions. 
       This is WORK IN PROGRESS, thus the assert(0) prior to it. */
       if(celev < neighbor_elev[0] && celev < neighbor_elev[1] && celev < neighbor_elev[2] && celev < neighbor_elev[3]){
      elevsink=neighbor_elev[0];
      for (n = 1; n < NDIRSfine; n++) {
	if (neighbor_elev[n] < elevsink && neighbor_elev[n] != OUTSIDEBASIN)
	  elevsink= neighbor_elev[n];
      }
      celev=elevsink;
    } 
  
    for (n = 0; n < NDIRSfine; n++) {
      if (neighbor_elev[n] == OUTSIDEBASIN) 
	neighbor_elev[n] = celev;
      if(n==1 || n==3)
	temp_slope[n] = (atan((celev - neighbor_elev[n]) / dy)) * DEGPRAD;
      else
	temp_slope[n] = (atan((celev - neighbor_elev[n]) / dx)) * DEGPRAD;
    }
    break;
  default:
    ReportError("ElevationSlope",65);
    assert(0);	/* nothing else works */
    break;
  }
  
  Slope = -999.;
  *Aspect = -99.;

  for (n = 0; n < NDIRSfine; n++){
    if (temp_slope[n] >= 0.) {
      if(temp_slope[n] > Slope) {
	if((y + yneighborfine[n]) == prevy && (x + xneighborfine[n]) == prevx) {
	 /*  fprintf(stderr, "%d %d %d Not allowed to back track!\n", n, y, x); */
	}
	else {
	  Slope = temp_slope[n];
	  *Aspect = temp_aspect[n] * PI / 180.0;
	  *nexty = y + yneighborfine[n];
	  *nextx = x + xneighborfine[n];
	}
      }

    }
  }
  
  if(Slope == -999.) {
    fprintf(stderr, "Sink encountered in cell y= %d x= %d, all routes from here go up!\n", y,x);
  }

  return Slope;
}

/* -------------------------------------------------------------
   ElevationSlopeAspectfine
   ------------------------------------------------------------- */
void ElevationSlopeAspectfine(MAPSIZE * Map, FINEPIX ** FineMap)
{
  const char *Routine = "ElevationSlopeAspectfine";
  int x;
  int y;
 /*  int n; */
  int k;
 /*  float neighbor_elev[NDIRSfine]; */

  /*fill neighbor array */

  for (x = 0; x < Map->NXfine; x++) { 
    for (y = 0; y < Map->NYfine; y++) {
      if (FineMap[y][x].Mask) { 
	
	/* Count the number of cells in the basin.  Need this to allocate memory for
	   the new, smaller Elev[] and Coords[][].  */
	if (INBASIN(FineMap[y][x].Mask)) 
	  Map->NumCellsfine++;
	    	
/* 	for (n = 0; n < NDIRSfine; n++) { */
/* 	  int xn = x + xneighborfine[n]; */
/* 	  int yn = y + yneighborfine[n]; */
	  
/* 	  if (valid_cell_fine(Map, xn, yn)) { */
/* 	    neighbor_elev[n] = ((FineMap[yn][xn].Mask) ? FineMap[yn][xn].Dem : (float) OUTSIDEBASIN); */
/* 	  } */
/* 	  else { */
/* 	    neighbor_elev[n] = OUTSIDEBASIN; */
/* 	  } */
/* 	} */
	
/* 	slope_aspect(Map->DMASS, Map->DMASS, FineMap[y][x].Dem, neighbor_elev, */
/* 		     &(FineMap[y][x].Slope), &(FineMap[y][x].Aspect)); */
	
      }/*  end of (FineMap[y][x].Mask */
    }
  }  /* end of for loops */
  
  /* Create a structure to hold elevations of only those cells
     within the basin and the y,x of those cells.*/
  
  if (!(Map->OrderedCellsfine = (ITEM *) calloc(Map->NumCellsfine, sizeof(ITEM))))
    ReportError((char *) Routine, 1);
  
  k = 0;
  for (y = 0; y < Map->NYfine; y++) {
    for (x = 0; x < Map->NXfine; x++) {
      /* Save the elevation, y, and x in the ITEM structure. */
      if (INBASIN(FineMap[y][x].Mask)) {
	Map->OrderedCellsfine[k].Rank = FineMap[y][x].Dem;
	Map->OrderedCellsfine[k].y = y;
	Map->OrderedCellsfine[k].x = x;
	k++;
      }
    }
  }
 
  /* Sort Elev in descending order-- Elev.x and Elev.y hold indices. */
  
  quick(Map->OrderedCellsfine, Map->NumCellsfine);

  return;
}

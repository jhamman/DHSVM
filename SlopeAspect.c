/*
 * SUMMARY:      SlopeAspect.c - Calculate slope and aspect of each pixel
 * USAGE:        Part of DHSVM/MWM
 *
 * AUTHOR:       William A Perkins
 * ORG:          Battelle Memorial Institute Pacific Northwest Laboratory
 * E-MAIL:       perk@clio.muse.pnl.gov
 * ORIG-DATE:    21-May-96
 * DESCRIPTION:  This module contains two routines to compute "slope" and
 *               "aspect"  (direction of slope): one which uses only terrain
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

int xdirection[NDIRS] = {
  0, 1, 0, -1
};
int ydirection[NDIRS] = {
  -1, 0, 1, 0
};

float temp_aspect[NNEIGHBORS] = {
  225., 180., 135., 90., 45., 0., 315., 270.
};

/* NNEIGHBORS used for redistribution of subsurface flow in topoindex
 * and for slope/aspect calculations, must equal 8. */
int xneighbor[NNEIGHBORS] = {
-1, 0, 1, 1, 1, 0, -1, -1
};

int yneighbor[NNEIGHBORS] = {
 1, 1, 1, 0, -1, -1, -1, 0
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
			 nelev[NNEIGHBORS], float *slope, float *aspect)
{
  int n;
  float dzdx, dzdy;
 
  for (n = 0; n < NNEIGHBORS; n++) {
    if (nelev[n] == OUTSIDEBASIN) {
      nelev[n] = celev;
    }
  }

  dzdx = ((nelev[0] + 2 * nelev[7] + nelev[6]) -
	  (nelev[2] + 2 * nelev[3] + nelev[4])) / (8 * dx);
  dzdy = ((nelev[0] + 2 * nelev[1] + nelev[2]) -
	  (nelev[4] + 2 * nelev[5] + nelev[6])) / (8 * dy);
 
  *slope = sqrt(dzdx * dzdx + dzdy * dzdy);
  
  /* Aspect is in radians cw from north, in the range -pi to pi. */
  if (fequal(dzdx, 0.0) && fequal(dzdy, 0.0)) {
    *aspect = 0.0;
  }
  else {
    *aspect = atan2(dzdx, dzdy);
  }

  //fprintf(stdout, "slope = %f\, aspect=%f\n",*slope, *aspect*180/3.14159);
  return;
}

/* -------------------------------------------------------------
   flow_fractions
   Computes subsurface flow fractions given the slope and aspect 
------------------------------------------------------------- */
static void flow_fractions(float dx, float dy, float slope, float aspect,
			   float nelev[NNEIGHBORS], float *grad,
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

    if ((cosine > 0 && nelev[5] == (float) OUTSIDEBASIN) ||
	(cosine < 0 && nelev[1] == (float) OUTSIDEBASIN))
      cosine = -cosine;
    if ((sine > 0 && nelev[3] == (float) OUTSIDEBASIN) ||
	(sine < 0 && nelev[7] == (float) OUTSIDEBASIN))
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
	effective_width = (sine > 0 ? sine * dy : 0.0);
	break;
      case 3:
	effective_width = (sine < 0 ? -sine * dy : 0.0);
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
  float neighbor_elev[NNEIGHBORS];
  int tempdir[NDIRS];
  int steepestdirection;
  unsigned int sum;
  float min;
  int xn, yn;

  /* fill neighbor array */
  
  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {

	/* Count the number of cells in the basin.  
	   Need this to allocate memory for
	   the new, smaller Elev[] and Coords[][].  */
	Map->NumCells++;

	for (n = 0; n < NNEIGHBORS; n++) {
	  xn = x + xneighbor[n];
	  yn = y + yneighbor[n];
	  
	  
	  if (valid_cell(Map, xn, yn)) {
	    neighbor_elev[n] = ((TopoMap[yn][xn].Mask) ? TopoMap[yn][xn].Dem : (float) OUTSIDEBASIN);
	  }
	  else {
	    neighbor_elev[n] = (float) OUTSIDEBASIN;
	  }
	}
	
	slope_aspect(Map->DX, Map->DY, TopoMap[y][x].Dem, neighbor_elev,
		     &(TopoMap[y][x].Slope), &(TopoMap[y][x].Aspect));
	
	/* fill Dirs in TopoMap too */
	
	flow_fractions(Map->DX, Map->DY, TopoMap[y][x].Slope,
		       TopoMap[y][x].Aspect,
		       neighbor_elev, &(TopoMap[y][x].FlowGrad),
		       TopoMap[y][x].Dir, &(TopoMap[y][x].TotalDir));

	/* Check that upslope neighbors && outsidebasin Dir = 0. */
	
	for (n = 0; n < NDIRS; n++)
	  tempdir[n] = 0; /* initialize */
	   
	sum = 0;
	for (n = 0; n < NDIRS; n++) {
	  xn = x + xdirection[n];
	  yn = y + ydirection[n];
	  
	  if ((TopoMap[y][x].Dir[n] > 0) && valid_cell(Map, xn, yn)) {
	    
	    if (!INBASIN(TopoMap[yn][xn].Mask))  { 
	      /* Can never have flow in this direction.*/
	      (TopoMap[y][x].TotalDir) -= TopoMap[y][x].Dir[n];
	      TopoMap[y][x].Dir[n] = 0; 
	    }
	    else if(TopoMap[yn][xn].Dem >= TopoMap[y][x].Dem) {
	      /* Will put flow in this direction if no other choice to
		 preserve water balance. */
	      (TopoMap[y][x].TotalDir) -= TopoMap[y][x].Dir[n];
	      tempdir[n]=TopoMap[y][x].Dir[n];
	      TopoMap[y][x].Dir[n] = 0; 
	    }
	  }
	  else if((TopoMap[y][x].Dir[n] > 0) && valid_cell(Map, xn, yn)) {
	    /* Can never have flow in this direction, should not be 
	       possible to get here. */
	    (TopoMap[y][x].TotalDir) -= TopoMap[y][x].Dir[n];
	      TopoMap[y][x].Dir[n] = 0; 
	    }
	  sum+=TopoMap[y][x].Dir[n];
	}
	
	/* If there is a sink, check again to see if there 
	   is a direction of steepest descent. Does not account 
	   for ties.*/
	steepestdirection = -99;
	if(sum==0) {
	 
	  min = DHSVM_HUGE;
	       
	  for (n = 0; n < NDIRS; n++) {
	    xn = x + xdirection[n];
	    yn = y + ydirection[n];
	  
	    if (valid_cell(Map, xn, yn)) {
	      if (INBASIN(TopoMap[yn][xn].Mask)) {
		if(TopoMap[yn][xn].Dem < min)
		{ 
		  min = TopoMap[yn][xn].Dem;
		  steepestdirection = n;
		}}
	    }
	  }
	  
	  if(min < TopoMap[y][x].Dem) {
	    TopoMap[y][x].Dir[steepestdirection] = (int)(255.0 + 0.5);
	    TopoMap[y][x].TotalDir = (int)(255.0 + 0.5);
	  }
	  else {
	    /*  Last resort: set the Dir of the cell to the cell that is
		closest in elevation. This should only happen for the 
		basin outlet, unless the Dem wasn't filled. */
	  
	    TopoMap[y][x].Dir[steepestdirection] = (int)(255.0 + 0.5);
	    TopoMap[y][x].TotalDir = (int)(255.0 + 0.5);
	    
	    xn = x + xdirection[steepestdirection];
	    yn = y + ydirection[steepestdirection];
	  }
	}
	
      } // end if (INBASIN(TopoMap[y][x].Mask)) {
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
void HeadSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SOILPIX ** SoilMap,
		     float **FlowGrad, unsigned char ***Dir, unsigned int **TotalDir)
{
  int x;
  int y;
  int n;
  float neighbor_elev[NNEIGHBORS];

  /* let's assume for now that WaterLevel is the SOILPIX map is
     computed elsewhere */

  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (INBASIN(TopoMap[y][x].Mask)) {

	float slope, aspect;

	for (n = 0; n < NNEIGHBORS; n++) {
	  int xn = x + xneighbor[n];
	  int yn = y + yneighbor[n];

	  if (valid_cell(Map, xn, yn)) {
	    neighbor_elev[n] =
	      ((TopoMap[yn][xn].Mask) ? SoilMap[yn][xn].WaterLevel :
	       (float) OUTSIDEBASIN);
	  }
	  else {
	    neighbor_elev[n] = (float) OUTSIDEBASIN;
	  }
	}

	slope_aspect(Map->DX, Map->DY, SoilMap[y][x].WaterLevel, neighbor_elev,
		     &slope, &aspect);
	flow_fractions(Map->DX, Map->DY, slope, aspect, neighbor_elev,
		       &(FlowGrad[y][x]), Dir[y][x], &(TotalDir[y][x])); 
	/* these were changed from TopoMap */
      }
    }
  }
  return;
}

/******************************************************************************/
/*			     ElevationSlope                            */
/* Part of MWM, should probably be merged w/ ElevationSlopeAspect function.   */
/******************************************************************************/

float ElevationSlope(MAPSIZE *Map, TOPOPIX **TopoMap, FINEPIX ***FineMap, int y, int x, int *nexty, 
		     int *nextx, int prevy, int prevx, float *Aspect) 
{
  int n, direction;
  float soil_elev[NNEIGHBORS];
  float bedrock_elev[NNEIGHBORS];
  float Slope;
  float temp_slope[NNEIGHBORS];
  double length_diagonal;
  float dx, dy, celev;
  int coarsej, coarsei;

  /* fill neighbor array */
  
  for (n = 0; n < NNEIGHBORS; n++) {

    int xn = x + xneighbor[n];
    int yn = y + yneighbor[n];
   
    // Initialize soil_elev and bedrock_elev
    soil_elev[n] = (float) OUTSIDEBASIN;
    bedrock_elev[n] = (float) OUTSIDEBASIN;

    // Check whether yn, xn are within FineMap array bounds
    if (valid_cell_fine(Map,xn,yn)){

      coarsej = floor(yn*Map->DMASS/Map->DY);
      coarsei = floor(xn*Map->DMASS/Map->DX);

      // Check whether FineMap element has been allocated for this cell
      // (equivalent to checking whether parent coarse grid cell is within coarse mask)
      if (INBASIN(TopoMap[coarsej][coarsei].Mask)) { 

	bedrock_elev[n] = (((*FineMap[yn][xn]).Mask) ? (*FineMap[yn][xn]).bedrock : (float) OUTSIDEBASIN);
	soil_elev[n] = (((*FineMap[yn][xn]).Mask) ? (*FineMap[yn][xn]).bedrock+(*FineMap[yn][xn]).sediment : (float) OUTSIDEBASIN);
	
      }
    }
    
  }       
  /*  Find bedrock slope in all directions. Negative slope = ascent, positive slope = descent.  */     
  dx = Map->DMASS;
  dy = Map->DMASS;
  celev = (*FineMap[y][x]).bedrock;


  length_diagonal = sqrt((pow(dx, 2)) + (pow(dy, 2))); 

  for (n = 0; n < NNEIGHBORS; n++) {
    if (bedrock_elev[n] == OUTSIDEBASIN) 
      bedrock_elev[n] = DHSVM_HUGE;
    
    if(n==0 || n==2 || n==4 || n==6)
      temp_slope[n] = (atan((celev - bedrock_elev[n]) / length_diagonal))
	* DEGPRAD;
    else if(n==1 || n==5)
      temp_slope[n] = (atan((celev - bedrock_elev[n]) / dy)) * DEGPRAD;
    else
      temp_slope[n] = (atan((celev - bedrock_elev[n]) / dx)) * DEGPRAD;
  }
    
 /* Find largest (positive) slope, this is the direction of failure along bedrock plain.  
     Backtracking isn't a problem if using the bedrock, but sinks may exist. */ 
   
  Slope = -999.;
  *Aspect = -99.;

  for (n = 0; n < NNEIGHBORS; n++){
    if(temp_slope[n] > Slope) {
      Slope = temp_slope[n];
      *Aspect = temp_aspect[n] * PI / 180.0;
      direction = n;
      *nexty = y + yneighbor[n];
      *nextx = x + xneighbor[n];
    }
  }

  /* If no positive slope found, a bedrock sink was encountered.  Assuming the sink should be filled to 
     the lowest "pour elevation", aspect should have already been assigned correctly. */

  /* Find dynamic slope in direction of steepest descent. */
  
  celev = (*FineMap[y][x]).bedrock + (*FineMap[y][x]).sediment;
  if(direction==0 || direction==2 || direction==4 || direction==6)
    Slope = (atan((celev - soil_elev[direction]) / length_diagonal))
      * DEGPRAD;
  else if(direction==1 || direction==5)
    Slope = (atan((celev - soil_elev[direction]) / dy)) * DEGPRAD;
  else
    Slope = (atan((celev - soil_elev[direction]) / dx)) * DEGPRAD;

  /* It is possible that a "soil" sink could be encountered at this point.  
     This is not really an error, and is checked for in MainMWM. */
  // if(Slope < 0.0 ) {
  // fprintf(stderr, "Sink encountered in cell y= %d x= %d, all routes from here go up!\n", y,x);
  // }

  if(Slope == -999. || *Aspect == -99.) {
    fprintf(stderr, "Aspect not assigned, this shouldn't have happened.\n");
    exit(0);
  }

  return Slope;
}

/* -------------------------------------------------------------
   ElevationSlopeAspectfine
   ------------------------------------------------------------- */
void ElevationSlopeAspectfine(MAPSIZE * Map, FINEPIX ***FineMap, TOPOPIX **TopoMap)
{
  const char *Routine = "ElevationSlopeAspectfine";
  int x;
  int y;
 /*  int n; */
  int k;
  int coarsei, coarsej;
 /*  float neighbor_elev[NNEIGHBORS]; */

  /* Create a structure to hold elevations of all cells within the coarse
     resolution mask and the y,x of those cells.*/
  
  if (!(Map->OrderedCellsfine = (ITEM *) calloc(Map->NumCellsfine, sizeof(ITEM))))
    ReportError((char *) Routine, 1);
  
  k = 0;
  for (y = 0; y < Map->NYfine; y++) {
    for (x = 0; x < Map->NXfine; x++) {
      
      coarsei = floor(y*Map->DMASS/Map->DY);
      coarsej = floor(x*Map->DMASS/Map->DX);
      
      /* Save the elevation, y, and x in the ITEM structure. */
      if (INBASIN(TopoMap[coarsei][coarsej].Mask)) {
	Map->OrderedCellsfine[k].Rank = (*FineMap[y][x]).Dem;
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

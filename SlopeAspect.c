/*
 * SUMMARY:      SlopeAspect.c - Calculate slope and aspect of each pixel
 * USAGE:        Part of DHSVM
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
 *               slope_aspect()
 *               flow_fractions()
 *               ElevationSlopeAspect()
 *               HeadSlopeAspect()
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

/* -------------------------------------------------------------
   valid_cell
   Checks to see if grid indices, x and y, are within the grid 
   defined by the specified Map
   ------------------------------------------------------------- */
int valid_cell(MAPSIZE * Map, int x, int y)
{
  return (x >= 0 && y >= 0 && x < Map->NX && y < Map->NY);
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
	assert(0);		/* How can this happen? */
      }
      dir[n] = (int) ((effective_width / total_width) * 255.0 + 0.5);
      *total_dir += dir[n];
    }
    break;
  case 8:
    assert(0);			/* can't handle this */
    break;
  default:
    assert(0);			/* other cases don't work either */
  }
  return;
}

/* -------------------------------------------------------------
   ElevationSlopeAspect
   ------------------------------------------------------------- */
void ElevationSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap)
{
  int x;
  int y;
  int n;
  float neighbor_elev[NDIRS];

  /* fill neighbor array */

  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (TopoMap[y][x].Mask) {

	for (n = 0; n < NDIRS; n++) {
	  int xn = x + xneighbor[n];
	  int yn = y + yneighbor[n];

	  if (valid_cell(Map, xn, yn)) {
	    neighbor_elev[n] =
	      ((TopoMap[yn][xn].Mask) ? TopoMap[yn][xn].
	       Dem : (float) OUTSIDEBASIN);
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

      }
    }
  }
  return;
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

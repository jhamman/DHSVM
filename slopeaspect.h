/*
 * SUMMARY:      slopeaspect.h - header file for SlopeAspect.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for SlopeAspect.c
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id$     
 */

#ifndef SLOPEASPECT_H
#define SLOPEASPECT_H

#include "settings.h"
#include "data.h"

/* -------------------------------------------------------------
   available variables
   ------------------------------------------------------------- */
extern int xneighbor[NDIRS];
extern int yneighbor[NDIRS];

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
void ElevationSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap);
void HeadSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SOILPIX ** SoilMap);
int valid_cell(MAPSIZE * Map, int x, int y);

#endif

/*
 * SUMMARY:      InitFineMaps() - Initialize fine resolution coverages
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Laura C. Bowling/Colleen O. Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       colleen@hydro.washington.edu
 * ORIG-DATE:    Oct-03
 * Last Change:  Mon Oct 27 14:40:44 2003 by Colleen O. Doten <colleen@hydro.washington.edu>
 * DESCRIPTION:  Initialize terrain coverages for fine resolution map
 * DESCRIP-END.
 * FUNCTIONS:    InitFineMaps()
 * COMMENTS:     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "constants.h"
#include "getinit.h"
#include "varid.h"
#include "sizeofnt.h"
#include "slopeaspect.h"

void CalcTopoIndex (MAPSIZE *Map, FINEPIX **FineMap);

/*****************************************************************************
  InitFineMaps()
*****************************************************************************/
void InitFineMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, 
		     LAYER *Soil, TOPOPIX ***TopoMap, SOILPIX ***SoilMap, 
		     FINEPIX ***FineMap)
{
 
  const char *Routine = "InitFineMaps";
  char VarName[BUFSIZE+1];	/* Variable name */
  int i, k, x, y;		/* Counters */
  int ii, jj, xx, yy;            /* Counters */
  int NumberType;		/* Number type of data set */
  float *Elev;                   /* Surface elevation */
   STRINIENTRY StrEnv[] = {
    {"FINEDEM", "DEM FILE"        , ""  , ""},
    {NULL       , NULL            , ""  , NULL}
  };
  
  printf("Initializing mass wasting resolution maps\n");

  /* Process the [FINEDEM] section in the input file */

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
                  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }
  
  /* Read the elevation dataset */ 
  
  GetVarName(001, 0, VarName);
  GetVarNumberType(001, &NumberType);
  if (!(Elev = (float *) calloc(Map->NXfine * Map->NYfine, 
				SizeOfNumberType(NumberType)))) 
    ReportError((char *) Routine, 1);
  Read2DMatrix(StrEnv[demfile].VarStr, Elev, NumberType, Map->NYfine, Map->NXfine, 0,
	       VarName);
  
  /* Assign the attributes to the correct map pixel */
  if (!(*FineMap = (FINEPIX **) calloc(Map->NYfine, sizeof(FINEPIX *))))
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NYfine; y++) {
    if (!((*FineMap)[y] = (FINEPIX *) calloc(Map->NXfine, sizeof(FINEPIX))))
      ReportError((char *) Routine, 1);
  }
  
  for (y = 0, i = 0; y < Map->NYfine; y++) 
    for (x = 0; x < Map->NXfine; x++, i++) 
      (*FineMap)[y][x].Dem  = Elev[i]; 
  free(Elev);
  
  /* Create fine resolution mask, sediment and bedrock maps. */
  
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      for (ii=0; ii< Map->DY/Map->DMASS; ii++) {
	for (jj=0; jj< Map->DX/Map->DMASS; jj++) {
	  yy = (int) y*Map->DY/Map->DMASS + ii;
	  xx = (int) x*Map->DX/Map->DMASS + jj;
	  (*FineMap)[yy][xx].Mask = (*TopoMap)[y][x].Mask;
	  (*FineMap)[yy][xx].bedrock = (*FineMap)[yy][xx].Dem - (*SoilMap)[y][x].Depth;
	  (*FineMap)[yy][xx].sediment = (*SoilMap)[y][x].Depth;
	  
	}
      }
    }
  }
  
  Map->NumFineIn = (Map->DX/Map->DMASS) * (Map->DY/Map->DMASS);
  
  /* Calculate slope, aspect, magnitude of subsurface flow gradient, and 
     fraction of flow flowing in each direction based on the land surface 
     slope. */
  ElevationSlopeAspectfine(Map, *FineMap); 
  
  printf("Basin has %d active pixels in the mass wasting resolution map\n",
	 Map->NumCellsfine);
  
  /* Calculate the topographic index */
  CalcTopoIndex(Map, *FineMap);
  
  for (y = 0; y < Map->NY; y++) {
    for (x  = 0; x < Map->NX; x++) {
      if (INBASIN((*TopoMap)[y][x].Mask)) {
	if (!((*TopoMap)[y][x].OrderedTopoIndex = (ITEM *) calloc(Map->NumFineIn, sizeof(ITEM))))
	  ReportError((char *) Routine, 1);
      }
    }
  }
  
  for (y = 0; y < Map->NY; y++) {
    for (x  = 0; x < Map->NX; x++) {
      if (INBASIN((*TopoMap)[y][x].Mask)) {
	k = 0;
	for(ii=0; ii< Map->DY/Map->DMASS; ii++) {
	  for(jj=0; jj< Map->DX/Map->DMASS; jj++) {
	    yy = (int) y*Map->DY/Map->DMASS + ii;
	    xx = (int) x*Map->DX/Map->DMASS + jj;
	    (*TopoMap)[y][x].OrderedTopoIndex[k].Rank = (*FineMap)[yy][xx].TopoIndex;
	    (*TopoMap)[y][x].OrderedTopoIndex[k].y = yy;
	    (*TopoMap)[y][x].OrderedTopoIndex[k].x = xx;
	    k++;
	  }
	}
       	quick((*TopoMap)[y][x] .OrderedTopoIndex, Map->NumFineIn);
      }
    }
  }
}
  

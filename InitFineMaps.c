/*
 * SUMMARY:      InitFineMaps() - Initialize fine resolution coverages
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * Last Change: Thu Feb 11 23:45:44 1999 by Bart Nijssen <nijssen@u.washington.edu>
 * DESCRIPTION:  Initialize terrain coverages
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


/*****************************************************************************
  InitFineMaps()
*****************************************************************************/
void InitFineMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, 
		     LAYER *Soil, TOPOPIX ***TopoMap, SOILPIX ***SoilMap, 
		     FINEPIX ***FineMap)
{
 
  const char *Routine = "InitFineMaps";
  char VarName[BUFSIZE+1];	/* Variable name */
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  int ii, jj, xx, yy;
  int NumberType;		/* Number type of data set */
  float *Elev;                  /* Surface elevation */
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

  /* Create fine resolution sediment and bedrock maps. */

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      for(ii=0; ii< Map->DY/Map->DMASS; ii++) {
	for(jj=0; jj< Map->DX/Map->DMASS; jj++) {
	  yy = (int) y*Map->DY/Map->DMASS + ii;
	  xx = (int) x*Map->DX/Map->DMASS+jj;
	  (*FineMap)[yy][xx].bedrock = (*FineMap)[yy][xx].Dem - (*SoilMap)[y][x].Depth;
	  (*FineMap)[yy][xx].sediment = (*SoilMap)[y][x].Depth;
	}}
    }
  }
}


/*
 * SUMMARY:      InitParameters.c - Initialize constants for MWM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96 
 * DESCRIPTION:  Initialize constants for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    InitParameters()
 * COMMENTS:
 * $Id$     
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "Calendar.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "getinit.h"
#include "constants.h"
#include "rad.h"

/*****************************************************************************
  Function name: InitParameterss()

  Purpose      : Initialize constants and settings for DHSVM run
                 Processes the following sections in InFile:
                 [PARAMETERS]

  Required     :
    LISTPTR Input          - Linked list with input strings
    OPTIONSTRUCT *Options   - Structure with different program options
    MAPSIZE *Map            - Coverage and resolution of model area

  Returns      : void

  Modifies     : (see list of required above)

  Comments     :
*****************************************************************************/
void InitParameters(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map)
{
  int i;			/* counter */

  STRINIENTRY StrEnv[] = {
    {"PARAMETERS", "MASS WASTING SPACING", "", ""},
    {"PARAMETERS", "MAXIMUM ITERATIONS", "", ""},
    {"PARAMETERS", "CHANNEL PARENT D50", "", ""},
    {"PARAMETERS", "CHANNEL PARENT D90", "", ""},
    {"PARAMETERS", "DEBRIS FLOW D50", "", ""},
    {"PARAMETERS", "DEBRIS FLOW D90", "", ""},
    {NULL, NULL, "", NULL}
  };

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default, 
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
  }

  /**************** Determine model options ****************/

  if (!CopyFloat(&(Map->DMASS), StrEnv[mass_spacing].VarStr, 1))
    ReportError(StrEnv[mass_spacing].KeyName, 51);
  
  Map->NYfine = Map->NY * (Map->DY/Map->DMASS);
  Map->NXfine = Map->NX * (Map->DY/Map->DMASS);
  Map->NumCellsfine = 0;

  if (!CopyFloat(&MASSITER, StrEnv[max_iterations].VarStr, 1))
  ReportError(StrEnv[max_iterations].KeyName, 51);

  if (!CopyFloat(&CHANNELd50, StrEnv[channeld50].VarStr, 1))
    ReportError(StrEnv[channeld50].KeyName, 51);

  if (!CopyFloat(&CHANNELd90, StrEnv[channeld90].VarStr, 1))
    ReportError(StrEnv[channeld90].KeyName, 51);

  if (!CopyFloat(&DEBRISd50, StrEnv[debrisd50].VarStr, 1))
    ReportError(StrEnv[debrisd50].KeyName, 51);
  
  if (!CopyFloat(&DEBRISd90, StrEnv[debrisd90].VarStr, 1))
    ReportError(StrEnv[debrisd90].KeyName, 51);
}

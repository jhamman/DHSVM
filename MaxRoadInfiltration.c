/*
 * SUMMARY:      MaxRoadInfiltration.c - Calculate area averaged road infiltration 
 *               rate
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * E-MAIL:       nijssen@u.arizona.edu
 * ORIG-DATE:    
 * DESCRIPTION:  This subroutine calculates an area averaged maximum
 *               infiltration rate for the roads in a grid cell.
 * DESCRIP-END.
 * FUNCTIONS:    MaxRoadInfiltration()
 * COMMENTS:
 * $Id$     
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMChannel.h"
#include "functions.h"

/*****************************************************************************
  Function name: MaxRoadInfiltration()

*****************************************************************************/
float MaxRoadInfiltration(ChannelMapPtr **map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  float area;
  float infiltration = 0.0;
  float tot_area = 0.;

  while (cell != NULL) {
    area = cell->length * cell->cut_width;
    infiltration += area * cell->channel->class->infiltration;
    tot_area += area;
    cell = cell->next;
  }

  infiltration /= tot_area;

  return infiltration;
}

/* -------------------------------------------------------------
   file: DHSVMChannel.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created August 30, 1996 by  William A Perkins
   $Id$
   ------------------------------------------------------------- */

#ifndef _DHSVMChannel_h_
#define _DHSVMChannel_h_

#include "settings.h"		/* for data.h */
#include "data.h"
#include "getinit.h"
#include "channel.h"
#include "channel_grid.h"

/* -------------------------------------------------------------
   struct CHANNEL
   ------------------------------------------------------------- */
typedef struct {
  ChannelClass *stream_class;
  ChannelClass *road_class;
  Channel *streams;
  Channel *roads;
  ChannelMapPtr **stream_map;
  ChannelMapPtr **road_map;
  FILE *streamout;
  FILE *roadout;
  FILE *streamflowout;
  FILE *roadflowout;
} CHANNEL;

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
void InitChannel(LISTPTR Input, MAPSIZE *Map, int deltat, CHANNEL *channel,
		 SOILPIX **SoilMap, int *MaxStreamID, int *MaxRoadID);
void InitChannelDump(CHANNEL *channel, char *DumpPath);
double ChannelCulvertFlow(int y, int x, CHANNEL *ChannelData);
void RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map,
		  TOPOPIX **TopoMap, SOILPIX **SoilMapi, AGGREGATED *Total);
void ChannelCut(int y, int x, CHANNEL *ChannelData, ROADSTRUCT *Network);
uchar ChannelFraction(TOPOPIX *topo, ChannelMapRec *rds);

#endif

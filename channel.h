/* -------------------------------------------------------------
   file: channel.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created October 24, 1995 by  William A Perkins
   $Id$
   ------------------------------------------------------------- */

#ifndef _channel_h_
#define _channel_h_

typedef unsigned short int SegmentID, ClassID;

/* -------------------------------------------------------------
   struct ChannelClass
   ------------------------------------------------------------- */
typedef enum {
  CHAN_OUTSLOPED, CHAN_CROWNED, CHAN_INSLOPED
} ChannelCrownType;

typedef struct _channel_class_rec_ {
  ClassID id;			/* unique identifier */

  float width;			/* ``channel'' width */
  float bank_height;		/* bank height for streams (or cut height for roads) */
  float friction;		/* Manning's n */
  float infiltration;		/* infiltration through ditch surface.  Note,
				   this may not be what you think it is, so be
				   sure to read the documentation before you use
				   it.  It is ONLY used for road networks and if
				   the option ROAD INFILTRATION is set to TRUE. */
  ChannelCrownType crown;

  struct _channel_class_rec_ *next;

} ChannelClass;

/* -------------------------------------------------------------
   struct Channel
   This is the basic unit of channel information.
   ------------------------------------------------------------- */
struct _channel_rec_ {
  SegmentID id;

  unsigned order;		/* determines computation order */
  char *record_name;		/* The name this segment is to have in
				   the output, if output is recorded */
  char record;			/* TRUE if outflow values are to be
				   saved by channel_save_outflow */

  float length;			/* Parameters */
  float slope;
  float K;
  float X;

  ChannelClass *class;		/* ChannelClass identifier */

  /* necessary routing terms */

  float lateral_inflow;		/* cubic meters */
  float last_inflow;		/* cubic meters */
  float last_outflow;		/* cubic meters */
  float last_storage;		/* cubic meters */
  float inflow;			/* cubic meters */
  float outflow;		/* cubic meters */
  float storage;		/* cubic meters */
  float last_lateral_inflow;	/* cubic meters */

  struct _channel_rec_ *outlet;	/* NULL if does not drain to another segment */

  struct _channel_rec_ *next;
};
typedef struct _channel_rec_ Channel, *ChannelPtr;

/* -------------------------------------------------------------
   externally available routines
   ------------------------------------------------------------- */

				/* ChannelClass */

ChannelClass *channel_read_classes(const char *file);
void channel_free_classes(ChannelClass * head);

				/* Channel */

Channel *channel_read_network(const char *file, ChannelClass * class_list);
void channel_routing_parameters(Channel * net, int deltat);
Channel *channel_find_segment(Channel * net, SegmentID id);
int channel_step_initialize_network(Channel * net);
int channel_incr_lat_inflow(Channel * segment, float linflow);
int channel_route_network(Channel * net, int deltat);
int channel_save_outflow(double time, Channel * net, FILE * file, FILE * file2);
int channel_save_outflow_text(char *tstring, Channel * net, FILE * out,
			      FILE * out2, int flag);
void channel_free_network(Channel * net);

				/* Module */

void channel_init(void);
void channel_done(void);

#endif

# 
# SUMMARY:      makefile for DHSVM
# USAGE:        make DHSVM

#	$Id$	

OBJS = AdjustStorage.o Aggregate.o AggregateRadiation.o		     \
CalcAerodynamic.o CalcAvailableWater.o CalcDistance.o		     \
CalcEffectiveKh.o CalcKhDry.o CalcSnowAlbedo.o CalcSolar.o	     \
CalcTotalWater.o CalcTransmissivity.o CalcWeights.o Calendar.o	     \
CanopyResistance.o ChannelState.o CheckOut.o CutBankGeometry.o	     \
DHSVMChannel.o Desorption.o Draw.o EvalExponentIntegral.o	     \
EvapoTranspiration.o ExecDump.o FileIOBin.o FileIONetCDF.o Files.o   \
FinalMassBalance.o GetInit.o GetMetData.o InArea.o InitAggregated.o  \
InitArray.o InitConstants.o InitDump.o InitFileIO.o		     \
InitInterpolationWeights.o InitMetMaps.o InitMetSources.o	     \
InitModelState.o InitNetwork.o InitNewMonth.o InitSnowMap.o	     \
InitTables.o InitTerrainMaps.o InitUnitHydrograph.o InitXGraphics.o  \
InterceptionStorage.o IsStationLocation.o LapseT.o LookupTable.o     \
MainDHSVM.o MakeLocalMetData.o MassBalance.o MassEnergyBalance.o     \
MassRelease.o MaxRoadInfiltration.o NoEvap.o RadiationBalance.o	     \
ReadMetRecord.o ReadRadarMap.o ReportError.o ResetAggregate.o	     \
RootBrent.o Round.o RouteSubSurface.o RouteSurface.o		     \
SatVaporPressure.o SensibleHeatFlux.o SeparateRadiation.o SizeOfNT.o \
SlopeAspect.o SnowInterception.o SnowMelt.o SnowPackEnergyBalance.o  \
SoilEvaporation.o StabilityCorrection.o StoreModelState.o	     \
SurfaceEnergyBalance.o UnsaturatedFlow.o VarID.o WaterTableDepth.o   \
channel.o channel_grid.o equal.o errorhandler.o globals.o tableio.o

SRCS = $(OBJS:%.o=%.c)

HDRS = Calendar.h DHSVMChannel.h DHSVMerror.h brent.h channel.h	    \
channel_grid.h constants.h data.h errorhandler.h fifoNetCDF.h	    \
fifobin.h fileio.h functions.h getinit.h lookuptable.h massenergy.h \
rad.h settings.h sizeofnt.h slopeaspect.h snow.h soilmoisture.h	    \
tableio.h varid.h

OTHER = makefile tableio.lex

RCSDIR= RCS/
GET= co
REL=

 
DEFS =  -DHAVE_X11 -DHAVE_NETCDF
#possible DEFS -DHAVE_NETCDF -DHAVE_X11 -DSHOW_MET_ONLY -DSNOW_ONLY
CFLAGS =  -O3 -I/usr/include/sys -I/usr/X11R6/include $(DEFS)
CC = cc
FLEX = /usr/bin/flex
LIBS = -lm -L/usr/X11R6/lib -lX11 -L/sw/lib -lnetcdf

DHSVM: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) -o DHSVM $(LIBS)

clean::
	rm -f DHSVM

library: libBinIO.a

BINIOOBJ = \
FileIOBin.o Files.o InitArray.o SizeOfNT.o Calendar.o \
ReportError.o

BINIOLIBOBJ = $(BINIOOBJ:%.o=libBinIO.a(%.o))

libBinIO.a: $(BINIOLIBOBJ)
	-ranlib $@

clean::
	rm -f libBinIO.a

# -------------------------------------------------------------
# rules for individual objects (created with make depend)
# -------------------------------------------------------------
AdjustStorage.o : settings.h soilmoisture.h 
Aggregate.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h 
AggregateRadiation.o : settings.h data.h Calendar.h massenergy.h 
CalcAerodynamic.o : DHSVMerror.h settings.h constants.h functions.h data.h \
  Calendar.h DHSVMChannel.h getinit.h channel.h channel_grid.h 
CalcAvailableWater.o : settings.h soilmoisture.h 
CalcDistance.o : settings.h data.h Calendar.h functions.h DHSVMChannel.h \
  getinit.h channel.h channel_grid.h 
CalcEffectiveKh.o : settings.h constants.h DHSVMerror.h functions.h data.h \
  Calendar.h DHSVMChannel.h getinit.h channel.h channel_grid.h 
CalcKhDry.o : settings.h functions.h data.h Calendar.h DHSVMChannel.h \
  getinit.h channel.h channel_grid.h 
CalcSnowAlbedo.o : settings.h constants.h data.h Calendar.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h 
CalcSolar.o : constants.h settings.h Calendar.h functions.h data.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h rad.h 
CalcTotalWater.o : settings.h soilmoisture.h 
CalcTransmissivity.o : settings.h functions.h data.h Calendar.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h 
CalcWeights.o : constants.h settings.h data.h Calendar.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h 
Calendar.o : settings.h functions.h data.h Calendar.h DHSVMChannel.h \
  getinit.h channel.h channel_grid.h DHSVMerror.h 
CanopyResistance.o : settings.h massenergy.h data.h Calendar.h constants.h 
ChannelState.o : settings.h data.h Calendar.h DHSVMerror.h fileio.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h sizeofnt.h 
CheckOut.o : DHSVMerror.h settings.h data.h Calendar.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h 
CutBankGeometry.o : settings.h soilmoisture.h 
DHSVMChannel.o : constants.h getinit.h DHSVMChannel.h settings.h data.h \
  Calendar.h channel.h channel_grid.h DHSVMerror.h functions.h \
  errorhandler.h fileio.h 
Desorption.o : settings.h massenergy.h data.h Calendar.h constants.h 
Draw.o : settings.h data.h Calendar.h functions.h DHSVMChannel.h getinit.h \
  channel.h channel_grid.h snow.h 
EvalExponentIntegral.o : settings.h data.h Calendar.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h 
EvapoTranspiration.o : settings.h data.h Calendar.h DHSVMerror.h \
  massenergy.h constants.h 
ExecDump.o : settings.h data.h Calendar.h fileio.h sizeofnt.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h 
FileIOBin.o : fifobin.h fileio.h sizeofnt.h settings.h DHSVMerror.h 
FileIONetCDF.o : 
Files.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h fileio.h 
FinalMassBalance.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h 
GetInit.o : DHSVMerror.h fileio.h getinit.h 
GetMetData.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h rad.h 
InArea.o : constants.h settings.h data.h Calendar.h 
InitAggregated.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h 
InitArray.o : functions.h data.h settings.h Calendar.h DHSVMChannel.h \
  getinit.h channel.h channel_grid.h 
InitConstants.o : settings.h data.h Calendar.h fileio.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h rad.h 
InitDump.o : settings.h data.h Calendar.h DHSVMerror.h fileio.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h sizeofnt.h \
  varid.h 
InitFileIO.o : fileio.h fifobin.h fifoNetCDF.h DHSVMerror.h 
InitInterpolationWeights.o : settings.h data.h Calendar.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h 
InitMetMaps.o : settings.h constants.h data.h Calendar.h DHSVMerror.h \
  fileio.h functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  rad.h sizeofnt.h 
InitMetSources.o : settings.h data.h Calendar.h fileio.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h rad.h 
InitModelState.o : settings.h data.h Calendar.h DHSVMerror.h fileio.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h sizeofnt.h soilmoisture.h varid.h 
InitNetwork.o : constants.h settings.h data.h Calendar.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  soilmoisture.h 
InitNewMonth.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h fifobin.h \
  fileio.h rad.h slopeaspect.h sizeofnt.h 
InitSnowMap.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h 
InitTables.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h fileio.h 
InitTerrainMaps.o : settings.h data.h Calendar.h DHSVMerror.h fileio.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h sizeofnt.h slopeaspect.h varid.h 
InitUnitHydrograph.o : settings.h constants.h data.h Calendar.h \
  DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
  channel_grid.h fileio.h sizeofnt.h 
InitXGraphics.o : settings.h data.h Calendar.h DHSVMerror.h 
InterceptionStorage.o : settings.h data.h Calendar.h DHSVMerror.h \
  massenergy.h constants.h 
IsStationLocation.o : settings.h data.h Calendar.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h 
LapseT.o : settings.h data.h Calendar.h functions.h DHSVMChannel.h \
  getinit.h channel.h channel_grid.h 
LookupTable.o : lookuptable.h DHSVMerror.h 
MainDHSVM.o : settings.h constants.h data.h Calendar.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h fileio.h 
MakeLocalMetData.o : settings.h data.h Calendar.h snow.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h rad.h 
MassBalance.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h 
MassEnergyBalance.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h massenergy.h snow.h \
  constants.h soilmoisture.h 
MassRelease.o : constants.h settings.h massenergy.h data.h Calendar.h \
  snow.h 
MaxRoadInfiltration.o : settings.h data.h Calendar.h DHSVMChannel.h \
  getinit.h channel.h channel_grid.h functions.h 
NoEvap.o : settings.h data.h Calendar.h massenergy.h 
RadiationBalance.o : settings.h data.h Calendar.h DHSVMerror.h massenergy.h \
  constants.h 
ReadMetRecord.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h 
ReadRadarMap.o : settings.h data.h Calendar.h DHSVMerror.h fileio.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h sizeofnt.h 
ReportError.o : settings.h data.h Calendar.h DHSVMerror.h 
ResetAggregate.o : settings.h data.h Calendar.h functions.h DHSVMChannel.h \
  getinit.h channel.h channel_grid.h 
RootBrent.o : settings.h brent.h massenergy.h data.h Calendar.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h DHSVMerror.h 
Round.o : functions.h data.h settings.h Calendar.h DHSVMChannel.h getinit.h \
  channel.h channel_grid.h DHSVMerror.h 
RouteSubSurface.o : settings.h data.h Calendar.h DHSVMerror.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h \
  soilmoisture.h slopeaspect.h 
RouteSurface.o : settings.h data.h Calendar.h slopeaspect.h DHSVMerror.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h 
SatVaporPressure.o : lookuptable.h 
SensibleHeatFlux.o : settings.h data.h Calendar.h DHSVMerror.h massenergy.h \
  constants.h brent.h functions.h DHSVMChannel.h getinit.h channel.h \
  channel_grid.h 
SeparateRadiation.o : settings.h rad.h 
SizeOfNT.o : DHSVMerror.h sizeofnt.h 
SlopeAspect.o : constants.h settings.h data.h Calendar.h functions.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h slopeaspect.h 
SnowInterception.o : brent.h constants.h settings.h massenergy.h data.h \
  Calendar.h snow.h functions.h DHSVMChannel.h getinit.h channel.h \
  channel_grid.h 
SnowMelt.o : brent.h constants.h settings.h massenergy.h data.h Calendar.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h snow.h 
SnowPackEnergyBalance.o : settings.h constants.h massenergy.h data.h \
  Calendar.h snow.h functions.h DHSVMChannel.h getinit.h channel.h \
  channel_grid.h 
SoilEvaporation.o : settings.h DHSVMerror.h massenergy.h data.h Calendar.h \
  constants.h 
StabilityCorrection.o : settings.h massenergy.h data.h Calendar.h \
  constants.h 
StoreModelState.o : settings.h data.h Calendar.h DHSVMerror.h fileio.h \
  functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
  constants.h sizeofnt.h varid.h 
SurfaceEnergyBalance.o : settings.h massenergy.h data.h Calendar.h \
  constants.h 
UnsaturatedFlow.o : constants.h settings.h functions.h data.h Calendar.h \
  DHSVMChannel.h getinit.h channel.h channel_grid.h soilmoisture.h 
VarID.o : settings.h data.h Calendar.h DHSVMerror.h sizeofnt.h varid.h 
WaterTableDepth.o : settings.h soilmoisture.h 
channel.o : errorhandler.h channel.h tableio.h settings.h 
channel_grid.o : channel_grid.h channel.h settings.h data.h Calendar.h \
  tableio.h errorhandler.h 
equal.o : functions.h data.h settings.h Calendar.h DHSVMChannel.h getinit.h \
  channel.h channel_grid.h 
errorhandler.o : errorhandler.h 
globals.o : 
tableio.o : tableio.h errorhandler.h settings.h 

tableio.c: tableio.lex
	$(FLEX) -Ptable_yy -o$@ $<

#clean::
#	rm -f tableio.c

# -------------------------------------------------------------
# sources
# -------------------------------------------------------------
sources: $(SRCS) $(HDRS) $(OTHER)

clean::
	rm -f $(OBJS)
	rm -f *~

# -------------------------------------------------------------
# tags 
# so we can find our way around
# -------------------------------------------------------------
tags: TAGS
TAGS: $(SRCS) $(HDRS)
	etags $(SRCS) $(HDRS)

clean::
	\rm -f TAGS


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
depend: .depend
.depend: $(SRCS)
	$(CC) $(CFLAGS) -MM $(SRCS) > $@

clean:: 
	rm -f .depend


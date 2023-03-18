#name:		  	mh370inmarsatprocessing.py
#created:		March 2022
#by:			paul.kennedy@guardiangeomatics.com
#description:   Read the Inamrsat flight log file and convert into spatial data.

import sys
import os
import math
# import numpy as np
import os.path
import time
from datetime import datetime
from datetime import timedelta
# import savgol
from argparse import ArgumentParser
import matplotlib.pyplot as plt
# from scipy import stats
# from scipy import signal
# from scipy.ndimage.filters import uniform_filter1d
import timeseries
import fileutils
import geodetic
import simplekml

###############################################################################
def main():
	parser = ArgumentParser(description='Read an Inmarsat log file and convert to spatial data.')
	parser.add_argument('-i', dest='inputfile', action='store', default="inmarsatlog.csv", help='Input filename to process.')

	args = parser.parse_args()

	process(args)

###############################################################################
def process(args):

	#setup the basics....
	flightpath = []

	#07/03/2014 16:42:00 UTC departs KL
	date_object = datetime.strptime("07/03/2014 16:42:00", '%d/%m/%Y %H:%M:%S')
	timestamp = to_timestamp(date_object)
	flightpath.append(Ccoordinate(timestamp, 101.698056, 2.743333, 21.0, "Depart KL Airport"))

	date_object = datetime.strptime("07/03/2014 17:07:00", '%d/%m/%Y %H:%M:%S')
	timestamp = to_timestamp(date_object)
	flightpath.append(Ccoordinate(timestamp, 101.698056, 2.743333, 21.0, "Final ACARS Transmission"))

	# read the inmarsat log file....
	print ("File to process: %d" % (len(args.inputfile)))

	
	inmarsat = Cinmarsat(args.inputfile)

	engine = CMH370Engine()

	kml = simplekml.Kml()
	kml.AltitudeMode = 'absolute'

	kml.newpoint(name= engine.groundstation.name, coords=[(engine.groundstation.longitude, engine.groundstation.latitude, engine.groundstation.elevation)])  # lon, lat, optional height

	ls = kml.newlinestring(name="Pathway", description="linePerthtoIOR",
						coords=[(engine.groundstation.longitude, engine.groundstation.latitude, engine.groundstation.elevation),
						 (engine.IORnominalposition.longitude, engine.IORnominalposition.latitude, engine.IORnominalposition.elevation)])
	ls.altitudemode = simplekml.AltitudeMode.absolute
	kml.newpoint(name= engine.IORnominalposition.name, altitudemode='absolute', coords=[(engine.IORnominalposition.longitude, engine.IORnominalposition.latitude, engine.IORnominalposition.elevation)])  # lon, lat, optional height
	kml.newpoint(name= engine.KLIAposition.name,  coords=[(engine.KLIAposition.longitude, engine.KLIAposition.latitude, engine.KLIAposition.elevation)])  # lon, lat, optional height

	ls = kml.newlinestring(name="Pathway", description="lineKLIAtoIOR",
						coords=[(engine.KLIAposition.longitude, engine.KLIAposition.latitude, engine.KLIAposition.elevation),
						 (engine.IORnominalposition.longitude, engine.IORnominalposition.latitude, engine.IORnominalposition.elevation)])
	ls.altitudemode = simplekml.AltitudeMode.absolute

	kml.save("pk1.kml")

	return
###############################################################################
def	exportxyz(records, inputfilename):

	inputfolder = os.getcwd()
	name = os.path.basename(inputfilename)
	parts = os.path.splitext(name)
	outname = os.path.join(inputfolder, parts[0] + "_I" + ".ixyz")
	outname = fileutils.createOutputFileName(outname, ext="")
	with open(outname, 'w') as fileptr:
		for record in records:
			msg = "%.10f %.10f %.10f %.10f\n" % (record[0], record[1], record[2], record[3])
			fileptr.write(msg)
	
###############################################################################
def to_timestamp(dateObject):
	return (dateObject - datetime(1970, 1, 1)).total_seconds()

###############################################################################
def to_DateTime(recordDate, recordTime):
	'''return a python date object from a split date and time record. works with kongsberg date and time structures'''
	date_object = datetime.strptime(str(recordDate), '%Y%m%d') + timedelta(0,recordTime)
	return date_object

###############################################################################
def from_timestamp(unixtime):
	return datetime.utcfromtimestamp(unixtime)

# ###############################################################################
# def distancecrosscourse(pointx, pointy, startx, starty, endx, endy):
# 	p1=np.array([startx, starty])
# 	p2=np.array([endx, endy])
# 	p3=np.array([pointx, pointy])
# 	d=np.cross(p2-p1,p3-p1)/np.linalg.norm(p2-p1)
# 	return d

###############################################################################
class Ccoordinate:
	#class to hold a timestamp coordinate.
###############################################################################
	def __init__(self, timestamp, longitude, latitude, elevation=0, name=""):
		self.timestamp	= timestamp
		self.longitude	= longitude
		self.latitude	= latitude
		self.elevation 	= elevation
		self.name 		= name

###############################################################################
class Cinmarsatlogrecord:

###############################################################################
	def __init__(self, timestamp, channelname, ocean, groundstation, channelid, channeltype, description, bfo, bto):

		self.timestamp 		= timestamp
		self.channelname 	= channelname
		self.ocean 			= ocean
		self.groundstation 	= groundstation
		self.channelid 		= channelid
		self.channeltype 	= channeltype
		self.description 	= description
		self.bfo 			= bfo
		self.bto 			= bto

		self.speedlightmetrespersec		= 299792458 

###############################################################################
class Cinmarsat:
	'''# a class to read and parse the inmarsat log file'''
	# Time,Channel,NameOcean,Region,GESID(octal),ChannelUnitID,ChannelType,SUType,BurstFrequencyOffset(Hz)BFO,BurstTimingOffset(microseconds)BTO
	# 7/03/2014,16:00:13.406,IOR-R1200-0-36D3,IOR,305,8,R-ChannelRX,0x15-Log-on/Log-offAcknowledge,103,14820
	# 7/03/2014,16:00:13.906,IOR-P10500-0-3859,IOR,305,10,R-ChannelTX,0x15-Log-on/Log-offAcknowledge

###############################################################################
	def __init__(self, filename=""):
		self.aircraftlogrecords = []

		counter = 0
		with open(filename, 'r') as fileptr:
			for line in fileptr:
				counter = counter + 1
				if len(line) < 40:
					continue
				if line.startswith("#"):
					continue
				words = line.strip().split(",")
				# print (words)
				if len(words) < 9:
					continue
				recdate  = datetime.strptime(words[0] + " " + words[1], '%d/%m/%Y %H:%M:%S.%f') 
				timestamp = to_timestamp(recdate)
				rec = Cinmarsatlogrecord(recdate, words[2], words[3], words[4], words[5], words[6], words[7], words[8], words[9])
				self.aircraftlogrecords.append(rec)


###############################################################################
class CMH370Engine:

###############################################################################
	def __init__(self):
		#http://wikimapia.org/1647465/Inmarsat-Land-Earth-Station-Perth
		# WGS84 
		# 31°48'16"S   
		# 115°52'57"E
		# 22.16 this is the altitude of the ESA ground station, very close by (within 100m)
		#GES (Perth) 
		#115.88250000, -31.8044444, 22.16

		# https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm
		# X : -2368393m
		# Y : 4881307m
		# Z : -3342034m

		# x,y,z = GeographictoECEF (115.88250000, -31.8044444, 22.16)
		# longitude, latitude, elevation = ECEFtoGeographic (x,y,z)

		perthlon = 115.887162
		perthlat = -31.805231

		self.groundstation = Ccoordinate(0, perthlon, perthlat, 22.16)
		self.IORnominalposition = Ccoordinate(0, 64.500, 0.00, 35793600)

		#AES (KLIA) ECEF −1293000,6238300,303500
		# KLIA Lat:2.7 Long:101.7
		long, lat, elev = geodetic.ECEFtoGeographic(-1293000,6238300,303500)
		recdate  = datetime.strptime("07/03/2014 16:00:00", '%d/%m/%Y %H:%M:%S') 
		timestamp = to_timestamp(recdate)
		self.KLIAposition = Ccoordinate(timestamp, long, lat, elev)

		###########################################################################################################
		# read the inmarsat satellite postion file.  it drifts over time
		# now load the satellite positions.  the satellit is not geostationary so we need to take this into account
		# ECEF x,y,z,RangeToGES(m),RangeToKLIA(m)
		# 7/03/2014,00:00:00,18118900,38081800,706700,  39222700,37296000
		self.satellitepositions = []
		filename = "IORSatellitePositions.csv"
		with open(filename, 'r') as fileptr:
			for line in fileptr:
				if len(line) < 40:
					continue
				if line.startswith("#"):
					continue
				words = line.strip().split(",")
				recdate  = datetime.strptime(words[0] + " " + words[1], '%d/%m/%Y %H:%M:%S') 
				timestamp = to_timestamp(recdate)
				long, lat, elev = geodetic.ECEFtoGeographic(float(words[2]), float(words[3]), float(words[4]))
				rec = Ccoordinate(recdate, long, lat, elev)
				self.satellitepositions.append(rec)

###############################################################################
if __name__ == "__main__":
		main()
###############################################################################

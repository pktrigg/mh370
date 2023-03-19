#name:		  	mh370inmarsatprocessing.py
#created:		March 2022
#by:			paul.kennedy@guardiangeomatics.com
#description:   Read the Inamrsat flight log file and convert into spatial data.

import sys
import os
import math
import numpy as np
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
	inmarsatlog = Cinmarsatlogreader(args.inputfile)

	#create an instance of the engine.  this holds all data.
	engine = CMH370Engine()

	for rec in inmarsatlog.aircraftlogrecords:
		r = engine.calculateATO_BTORangemetres(rec.timestamp, rec.bto)

	#lets use the KLIA log record to demonstrate the maths for spherical intersection works...
	calcintersectionofspheres(engine.IORPosition, IORRange)




	#write out the results...
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
# https://archive.lib.msu.edu/crcmath/math/math/s/s563.htm
def calcintersectionofspheres(IORPosition, IORRange):

	# d = ecef height of satellite
	d = IORPosition.z

	# r = range aircraft to satellite

	# R = radius of earth at location directly below the IOR satellite
	x,y,z = geodetic.GeographictoECEF(IORPosition.longitude, IORPosition.latitude, 0.0)
	R = z
	rangeFromECEFnadir = (1/(2*d) * math.sqrt((-d+r-R) * (-d-r+R) * (-d+r+R) * (d+r+R)))

###############################################################################
def  position(a,r,u,v):
	return a + r*np.array([np.cos(u)*np.sin(v),np.sin(u)*np.sin(v),np.cos(v)])

###############################################################################
def f(args):
	# u1,v1,u2,v2,u3,v3,_,_,_ = args
	u1,v1,u2,v2 = args
	pos1 = position(a1,r1,u1,v1)
	pos2 = position(a2,r2,u2,v2)
	# pos3 = position(a3,r3,u3,v3)
	return np.array([pos1 - pos2, pos1 - pos3, pos2 - pos3]).flatten()
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
def getrangeSatelliteToPerth(timestamp, satellitepositions, earthposition):

	# given a timestamp compute the satellite position at that time
	ts 		= []
	long 	= []
	lat 	= []
	elev 	= []
	for rec in satellitepositions:
		ts.append(rec.timestamp)
		long.append(rec.longitude)
		lat.append(rec.latitude)
		elev.append(rec.elevation)

	tslongitude = timeseries.cTimeSeries(ts, long)
	tslatitude = timeseries.cTimeSeries(ts, lat)
	tselev = timeseries.cTimeSeries(ts, elev)

	latitudeattimestamp = (tslongitude.getValueAt(timestamp))
	longitudeattimestamp = (tslatitude.getValueAt(timestamp))
	elevationattimestamp = (tselev.getValueAt(timestamp))

	# we now know where the satellite is, so we can now compute the range to the earth satation
	range1, brg1, brg2 = geodetic.calculateRangeBearingFromGeographicals(longitudeattimestamp, latitudeattimestamp, earthposition.longitude, earthposition.latitude)
	x1,y1,z1 = geodetic.GeographictoECEF(earthposition.longitude, earthposition.latitude, earthposition.elevation)
	x2,y2,z2 = geodetic.GeographictoECEF(longitudeattimestamp, latitudeattimestamp, elevationattimestamp)
	range = geodetic.calculateECEFRange(x1, y1, z1, x2, y2, z2)

	return range

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
		self.x,self.y,self.z	= geodetic.GeographictoECEF(longitude, latitude, elevation)

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

###############################################################################
class Cinmarsatlogreader:
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
				bto = float(words[9]) / 1000000 # microseconds to seconds
				rec = Cinmarsatlogrecord(timestamp, words[2], words[3], words[4], words[5], words[6], words[7], words[8], bto)
				self.aircraftlogrecords.append(rec)

		# the first record represents the aircraft on the ground in KLIA...
		# self.kliarange
###############################################################################
class CMH370Engine:

###############################################################################
	def __init__(self):
		self.speedlightmetrespersec		= 299792458 

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
				rec = Ccoordinate(timestamp, long, lat, elev)
				self.satellitepositions.append(rec)


	###############################################################################
	def calculateATO_BTORangemetres(self, timestamp, BTO):

		# we need to convert the BTO to a range in metres...
		# file:///C:/Users/pkdesktop/Guardian%20Geomatics/Technical%20-%20Documents/MH370/the-search-for-mh370.pdf
		# Seventeen measurements were taken during this 30 minute period, and these were
		# processed to estimate the fixed timing bias. The mean bias of −495,679 μs was then
		# used to predict the path length from the measured data

		btobias = 495679 / 1000000 # microseconds to seconds conversion pkpkpk
		bto = float(BTO)
		rangeSatelliteToPerth = getrangeSatelliteToPerth(timestamp, self.satellitepositions, self.groundstation)

		ATOrangeSatelliteToAircraft = ((self.speedlightmetrespersec * (bto + btobias))/ 2) - rangeSatelliteToPerth
		return ATOrangeSatelliteToAircraft

###############################################################################
if __name__ == "__main__":
		main()
###############################################################################

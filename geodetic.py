#!/usr/bin/python
#
# ---------------------------------------------------------------------
# |																	 |
# |	geodetic.cc -  a collection of geodetic functions			   |
# |	Paul Kennedy May 2016										   |
# |	Jim Leven  - Dec 99											 |
# |																	 |
# | originally from:													|
# | http://wegener.mechanik.tu-darmstadt.de/GMT-Help/Archiv/att-8710/Geodetic_py |
# |ftp://pdsimage2.wr.usgs.gov/pub/pigpen/Python/Geodetic_py.py		|
# |																	 |
# ---------------------------------------------------------------------
#
#
# ------------------------------------------------------------------------------
# | Algrothims from Geocentric Datum of Australia Technical Manual			|
# | 												|
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html					|
# | 												|
# | This page last updated 11 May 1999 								|
# | 												|
# | Computations on the Ellipsoid									|
# | 												|
# | There are a number of formulae that are available		   		|
# | to calculate accurate geodetic positions, 							|
# | azimuths and distances on the ellipsoid.							|
# | 												|
# | Vincenty's formulae (Vincenty, 1975) may be used 						|
# | for lines ranging from a few cm to nearly 20,000 km, 					|
# | with millimetre accuracy. 									|
# | The formulae have been extensively tested 								|
# | for the Australian region, by comparison with results	   		|
# | from other formulae (Rainsford, 1955 & Sodano, 1965). 					|
# |												|
# | * Inverse problem: azimuth and distance from known 					|
# |			latitudes and longitudes 					|
# | * Direct problem: Latitude and longitude from known 					|
# |			position, azimuth and distance. 				|
# | * Sample data 										|
# | * Excel spreadsheet 											|
# | 												|
# | Vincenty's Inverse formulae										|
# | Given: latitude and longitude of two points				 		|
# |			(latitude1, longitude1 and latitude2, longitude2), 	|
# | Calculate: the ellipsoidal distance (s) and 						|
# | forward and reverse azimuths between the points (alpha1Tp2, alpha21).	|
# |											|
# ------------------------------------------------------------------------------

import math
import numpy as np
import sys
import os.path
import pyproj
from pyproj import Transformer

###############################################################################
def GeographictoECEF (longitude, latitude, elevation):
	ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
	lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
	x,y,z = pyproj.transform(lla, ecef, longitude, latitude, elevation, radians=False)
	return x,y,z

###############################################################################
def ECEFtoGeographic (x,y,z):
	ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
	lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
	lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)
	return lon, lat, alt

###############################################################################
def toCartesian (longitude, latitude, elevation):
	from pyproj import Transformer
	trans_GPS_to_XYZ = Transformer.from_crs(4979, 4978, always_xy=False)
	x,y,z = trans_GPS_to_XYZ.transform(longitude, latitude, elevation)
	return x,y,z

###############################################################################
def fromCartesian (x,y,z):
	from pyproj import Transformer
	trans_XYZ_to_GPS = Transformer.from_crs(4978,4979, always_xy=False)
	longitude, latitude, elevation = trans_XYZ_to_GPS.transform(x,y,z)
	return longitude, latitude, elevation

###############################################################################
def medfilt (x, k):
	"""Apply a length-k median filter to a 1D array x.
	Boundaries are extended by repeating endpoints.
	"""
	assert k % 2 == 1, "Median filter length must be odd."
	assert x.ndim == 1, "Input must be one-dimensional."
	k2 = (k - 1) // 2
	y = np.zeros ((len (x), k), dtype=x.dtype)
	y[:,k2] = x
	for i in range (k2):
		j = k2 - i
		y[j:,i] = x[:-j]
		y[:j,i] = x[0]
		y[:-j,-(i+1)] = x[j:]
		y[-j:,-(i+1)] = x[-1]
	return np.median (y, axis=1)

###############################################################################
# from: http://mathforum.org/library/drmath/view/62034.html
def calculateRangeBearingFromGridPosition(easting1, northing1, easting2, northing2):
	"""given 2 east, north, pairs, compute the range and bearing"""

	dx = easting2-easting1
	dy = northing2-northing1

	bearing = 90 - (180/math.pi)*math.atan2(northing2-northing1, easting2-easting1)
	return (math.sqrt((dx*dx)+(dy*dy)), bearing)

###############################################################################
# taken frm http://gis.stackexchange.com/questions/76077/how-to-create-points-based-on-the-distance-and-bearing-from-a-survey-point
def calculateGridPositionFromRangeBearing(easting, northing, distance, bearing):
	"""given an east, north, range and bearing, compute a new coordinate on the grid"""
	point =   (easting, northing)
	angle =   90 - bearing
	bearing = math.radians(bearing)
	angle =   math.radians(angle)

	# polar coordinates
	dist_x = distance * math.cos(angle)
	dist_y = distance * math.sin(angle)

	xfinal = point[0] + dist_x
	yfinal = point[1] + dist_y

	# direction cosines
	cosa = math.cos(angle)
	cosb = math.cos(bearing)
	xfinal = point[0] + (distance * cosa)
	yfinal = point[1] + (distance * cosb)

	return [xfinal, yfinal]

def calculateRangeBearingFromGeographicals(longitude1, latitude1,  longitude2,  latitude2 ) :
		"""
		Returns s, the distance between two geographic points on the ellipsoid
		and alpha1, alpha2, the forward and reverse azimuths between these points.
		lats, longs and azimuths are in decimal degrees, distance in metres

		Returns ( s, alpha1Tp2,  alpha21 ) as a tuple
		"""
		f = 1.0 / 298.257223563		# WGS84
		a = 6378137.0 			# metres

		if (abs( latitude2 - latitude1 ) < 1e-8) and ( abs( longitude2 - longitude1) < 1e-8 ) :
				return 0.0, 0.0, 0.0

		piD4   = math.atan( 1.0 )
		two_pi = piD4 * 8.0

		latitude1	= latitude1 * piD4 / 45.0
		longitude1 = longitude1 * piD4 / 45.0		# unfortunately lambda is a key word!
		latitude2	= latitude2 * piD4 / 45.0
		longitude2 = longitude2 * piD4 / 45.0

		b = a * (1.0 - f)

		TanU1 = (1-f) * math.tan( latitude1 )
		TanU2 = (1-f) * math.tan( latitude2 )

		U1 = math.atan(TanU1)
		U2 = math.atan(TanU2)

		lembda = longitude2 - longitude1
		last_lembda = -4000000.0		# an impossibe value
		omega = lembda

		# Iterate the following equations,
		#  until there is no significant change in lembda

		while ( last_lembda < -3000000.0 or lembda != 0 and abs( (last_lembda - lembda)/lembda) > 1.0e-9 ) :

				sqr_sin_sigma = pow( math.cos(U2) * math.sin(lembda), 2) + \
						pow( (math.cos(U1) * math.sin(U2) - \
						math.sin(U1) *  math.cos(U2) * math.cos(lembda) ), 2 )

				Sin_sigma = math.sqrt( sqr_sin_sigma )

				Cos_sigma = math.sin(U1) * math.sin(U2) + math.cos(U1) * math.cos(U2) * math.cos(lembda)

				sigma = math.atan2( Sin_sigma, Cos_sigma )

				Sin_alpha = math.cos(U1) * math.cos(U2) * math.sin(lembda) / math.sin(sigma)
				alpha = math.asin( Sin_alpha )

				Cos2sigma_m = math.cos(sigma) - (2 * math.sin(U1) * math.sin(U2) / pow(math.cos(alpha), 2) )

				C = (f/16) * pow(math.cos(alpha), 2) * (4 + f * (4 - 3 * pow(math.cos(alpha), 2)))

				last_lembda = lembda

				lembda = omega + (1-C) * f * math.sin(alpha) * (sigma + C * math.sin(sigma) * \
						(Cos2sigma_m + C * math.cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2) )))

		u2 = pow(math.cos(alpha),2) * (a*a-b*b) / (b*b)

		A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))

		B = (u2/1024) * (256 + u2 * (-128+ u2 * (74 - 47 * u2)))

		delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4) * \
				(Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2) ) - \
				(B/6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) * \
				(-3 + 4 * pow(Cos2sigma_m,2 ) )))

		s = b * A * (sigma - delta_sigma)

		alpha1Tp2 = math.atan2( (math.cos(U2) * math.sin(lembda)), \
				(math.cos(U1) * math.sin(U2) - math.sin(U1) * math.cos(U2) * math.cos(lembda)))

		alpha21 = math.atan2( (math.cos(U1) * math.sin(lembda)), \
				(-math.sin(U1) * math.cos(U2) + math.cos(U1) * math.sin(U2) * math.cos(lembda)))

		if ( alpha1Tp2 < 0.0 ) :
				alpha1Tp2 =  alpha1Tp2 + two_pi
		if ( alpha1Tp2 > two_pi ) :
				alpha1Tp2 = alpha1Tp2 - two_pi

		alpha21 = alpha21 + two_pi / 2.0
		if ( alpha21 < 0.0 ) :
				alpha21 = alpha21 + two_pi
		if ( alpha21 > two_pi ) :
				alpha21 = alpha21 - two_pi

		alpha1Tp2	= alpha1Tp2	* 45.0 / piD4
		alpha21	= alpha21	* 45.0 / piD4
		return s, alpha1Tp2,  alpha21

   # END of Vincenty's Inverse formulae


#-------------------------------------------------------------------------------
# Vincenty's Direct formulae							|
# Given: latitude and longitude of a point (latitude1, longitude1) and 			|
# the geodetic azimuth (alpha1Tp2) 						|
# and ellipsoidal distance in metres (s) to a second point,			|
# 										|
# Calculate: the latitude and longitude of the second point (latitude2, longitude2) 	|
# and the reverse azimuth (alpha21).						|
# 										|
#-------------------------------------------------------------------------------

def calculateGeographicalPositionFromRangeBearing(latitude1, longitude1, alpha1To2, s ):
		"""
		Returns the lat and long of projected point and reverse azimuth
		given a reference point and a distance and azimuth to project.
		lats, longs and azimuths are passed in decimal degrees

		Returns ( latitude2,  longitude2,  alpha2To1 ) as a tuple

		"""
		f = 1.0 / 298.257223563		# WGS84
		a = 6378137.0 			# metres

		piD4 = math.atan( 1.0 )
		two_pi = piD4 * 8.0

		latitude1	= latitude1	* piD4 / 45.0
		longitude1 = longitude1 * piD4 / 45.0
		alpha1To2 = alpha1To2 * piD4 / 45.0
		if ( alpha1To2 < 0.0 ) :
				alpha1To2 = alpha1To2 + two_pi
		if ( alpha1To2 > two_pi ) :
				alpha1To2 = alpha1To2 - two_pi

		b = a * (1.0 - f)

		TanU1 = (1-f) * math.tan(latitude1)
		U1 = math.atan( TanU1 )
		sigma1 = math.atan2( TanU1, math.cos(alpha1To2) )
		Sinalpha = math.cos(U1) * math.sin(alpha1To2)
		cosalpha_sq = 1.0 - Sinalpha * Sinalpha

		u2 = cosalpha_sq * (a * a - b * b ) / (b * b)
		A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * \
				(320 - 175 * u2) ) )
		B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2) ) )

		# Starting with the approximation
		sigma = (s / (b * A))

		last_sigma = 2.0 * sigma + 2.0	# something impossible

		# Iterate the following three equations
		#  until there is no significant change in sigma

		# two_sigma_m , delta_sigma
		while ( abs( (last_sigma - sigma) / sigma) > 1.0e-9 ) :
				two_sigma_m = 2 * sigma1 + sigma

				delta_sigma = B * math.sin(sigma) * ( math.cos(two_sigma_m) \
						+ (B/4) * (math.cos(sigma) * \
						(-1 + 2 * math.pow( math.cos(two_sigma_m), 2 ) -  \
						(B/6) * math.cos(two_sigma_m) * \
						(-3 + 4 * math.pow(math.sin(sigma), 2 )) *  \
						(-3 + 4 * math.pow( math.cos (two_sigma_m), 2 ))))) \

				last_sigma = sigma
				sigma = (s / (b * A)) + delta_sigma

		latitude2 = math.atan2 ( (math.sin(U1) * math.cos(sigma) + math.cos(U1) * math.sin(sigma) * math.cos(alpha1To2) ), \
				((1-f) * math.sqrt( math.pow(Sinalpha, 2) +  \
				pow(math.sin(U1) * math.sin(sigma) - math.cos(U1) * math.cos(sigma) * math.cos(alpha1To2), 2))))

		lembda = math.atan2( (math.sin(sigma) * math.sin(alpha1To2 )), (math.cos(U1) * math.cos(sigma) -  \
				math.sin(U1) *  math.sin(sigma) * math.cos(alpha1To2)))

		C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq ))

		omega = lembda - (1-C) * f * Sinalpha *  \
				(sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) + \
				C * math.cos(sigma) * (-1 + 2 * math.pow(math.cos(two_sigma_m),2) )))

		longitude2 = longitude1 + omega

		alpha21 = math.atan2 ( Sinalpha, (-math.sin(U1) * math.sin(sigma) +  \
				math.cos(U1) * math.cos(sigma) * math.cos(alpha1To2)))

		alpha21 = alpha21 + two_pi / 2.0
		if ( alpha21 < 0.0 ) :
				alpha21 = alpha21 + two_pi
		if ( alpha21 > two_pi ) :
				alpha21 = alpha21 - two_pi

		latitude2	   = latitude2	   * 45.0 / piD4
		longitude2	= longitude2	* 45.0 / piD4
		alpha21	= alpha21	* 45.0 / piD4

		return latitude2,  longitude2,  alpha21
# END of Vincenty's Direct formulae

#*******************************************************************

def est_dist(  latitude1,  longitude1,  latitude2,  longitude2 ) :
		"""

		Returns an estimate of the distance between two geographic points
		This is a quick and dirty vinc_dist
		which will generally estimate the distance to within 1%
		Returns distance in metres

		"""
		f = 1.0 / 298.257223563		# WGS84
		a = 6378137.0 			# metres

		piD4   = 0.785398163397

		latitude1	= latitude1 * piD4 / 45.0
		longitude1 = longitude1 * piD4 / 45.0
		latitude2	= latitude2 * piD4 / 45.0
		longitude2 = longitude2 * piD4 / 45.0

		c = math.cos((latitude2+latitude1)/2.0)

		return math.sqrt( pow(math.fabs(latitude2-latitude1), 2) + \
				pow(math.fabs(longitude2-longitude1)*c, 2) ) * a * ( 1.0 - f + f * c )
   	# END of rough estimate of the distance.

def getPRJFromEPSG(EPSGCode):
	'''read through the SRID.csv file from Pyproj to find the correct PRJ string for a given EPSG code.  This is used to write out a sensible PRJ file alongside a shape file. '''
	localpath = os.path.dirname(os.path.realpath(__file__))
	sys.path.append(localpath)
	filename = os.path.join(localpath, "srid.csv")
	# filename = 'srid.csv'
	if os.path.isfile(filename):
		datafile = open(filename)
		for line in datafile:
			if line.startswith(str(EPSGCode)):
				return line.split(";")[1]
				break
	return ""
def loadProj(EPSGCode):
	'''load a pyproj object using the supplied code'''
	# wgs84=pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
	#note: anaconda conda install has a bug when installing.  It is stupid and forgets to install the proj data folder.
	#to fix this, you need to
	#copy the data folder from c:\ggtools\python\pyproj to # C:\ProgramData\Anaconda3\Lib\site-packages\pyproj
	#rename the datadir.py to datadir.bak and then copy the datadir.py from the c:\ggtools\python\pyproj into the folder
	try:
		projection = pyproj.Proj("+init=EPSG:" + str(EPSGCode))
	except:
		return None
	return projection

def	writePRJ(filename, EPSGCode):
	'''try and find a matching PRJ string from the Proj CSV file.  If we find one, write it as a PRJ file so the shape file opens nicely in GIS'''
	prjstring = getPRJFromEPSG(EPSGCode)
	prj = open(filename, 'w')
	if len(prjstring) > 0:
		prj.write(prjstring) # python will convert \n to os.linesep
	else:
		prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]') # python will convert \n to os.linesep
	prj.close() # you can omit in most cases as the destructor will call it


#--------------------------------------------------------------------------
# Notes:
#
# * "The inverse formulae may give no solution over a line
# 	between two nearly antipodal points. This will occur when
# 	lembda ... is greater than pi in absolute value". (Vincenty, 1975)
#
# * In Vincenty (1975) L is used for the difference in longitude,
# 	however for consistency with other formulae in this Manual,
# 	omega is used here.
#
# * Variables specific to Vincenty's formulae are shown below,
# 	others common throughout the manual are shown in the Glossary.
#
#
# alpha = Azimuth of the geodesic at the equator
# U = Reduced latitude
# lembda = Difference in longitude on an auxiliary sphere (longitude1 & longitude2
# 		are the geodetic longitudes of points 1 & 2)
# sigma = Angular distance on a sphere, from point 1 to point 2
# sigma1 = Angular distance on a sphere, from the equator to point 1
# sigma2 = Angular distance on a sphere, from the equator to point 2
# sigma_m = Angular distance on a sphere, from the equator to the
# 		midpoint of the line from point 1 to point 2
# u, A, B, C = Internal variables
#
#
# Sample Data
#
# Flinders Peak
# -37 57'03.72030"
# 144 25'29.52440"
# Buninyong
# -37 39'10.15610"
# 143 55'35.38390"
# Ellipsoidal Distance
# 54,972.271 m
#
# Forward Azimuth
# 306 52'05.37"
#
# Reverse Azimuth
# 127 10'25.07"
#
#
#*******************************************************************

# Test driver

if __name__ == "__main__" :

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

	x,y,z = GeographictoECEF (115.88250000, -31.8044444, 22.16)
	longitude, latitude, elevation = ECEFtoGeographic (x,y,z)

	# pkpk 18/03/2023 the ECEF conversions work fine in both directions

	# x,y,z = toCartesian (64.500, 0.00, 357881226)
	# longitude, latitude, elevation = fromCartesian (x,y,z)

	f = 1.0 / 298.257223563		# WGS84
	a = 6378137.0 			# metres

	print  ("\n Ellipsoidal major axis =  %12.3f metres\n" % ( a ))
	print  ("\n Inverse flattening	 =  %15.9f\n" % ( 1.0/f ))

	print ("\n Test Flinders Peak to Buninyon")
	print ("\n ****************************** \n")
	latitude1 = -(( 3.7203 / 60. + 57) / 60. + 37 )
	longitude1 = ( 29.5244 / 60. + 25) / 60. + 144
	print ("Flinders Peak = %12.6f, %13.6f \n" % ( latitude1, longitude1 ))
	deg = int(latitude1)
	min = int(abs( ( latitude1 - deg) * 60.0 ))
	sec = abs(latitude1 * 3600 - deg * 3600) - min * 60
	print (" Flinders Peak =   %3i\xF8%3i\' %6.3f\",  " % ( deg, min, sec ),)
	deg = int(longitude1)
	min = int(abs( ( longitude1 - deg) * 60.0 ))
	sec = abs(longitude1 * 3600 - deg * 3600) - min * 60
	print (" %3i\xF8%3i\' %6.3f\" \n" % ( deg, min, sec ))

	latitude2 = -(( 10.1561 / 60. + 39) / 60. + 37 )
	longitude2 = ( 35.3839 / 60. + 55) / 60. + 143
	print ("\n Buninyon	  = %12.6f, %13.6f \n" % ( latitude2, longitude2 ))

	deg = int(latitude2)
	min = int(abs( ( latitude2 - deg) * 60.0 ))
	sec = abs(latitude2 * 3600 - deg * 3600) - min * 60
	print (" Buninyon	  =   %3i\xF8%3i\' %6.3f\",  " % ( deg, min, sec ),)
	deg = int(longitude2)
	min = int(abs( ( longitude2 - deg) * 60.0 ))
	sec = abs(longitude2 * 3600 - deg * 3600) - min * 60
	print (" %3i\xF8%3i\' %6.3f\" \n" % ( deg, min, sec ))

	dist   = est_dist(latitude1, longitude1, latitude2,  longitude2 )

	print ("\n Ellipsoidal Distance = %15.3f metres\n			should be		 54972.271 m\n" % ( dist ))
	# print ("\n Forward and back azimuths = %15.6f, %15.6f \n" % ( alpha1Tp2, alpha21 ))
	# deg = int(alpha1Tp2)
	# min = int( abs(( alpha1Tp2 - deg) * 60.0 ) )
	# sec = abs(alpha1Tp2 * 3600 - deg * 3600) - min * 60
	# print (" Forward azimuth = %3i\xF8%3i\' %6.3f\"\n" % ( deg, min, sec ))
	# deg = int(alpha21)
	# min = int(abs( ( alpha21 - deg) * 60.0 ))
	# sec = abs(alpha21 * 3600 - deg * 3600) - min * 60
	# print (" Reverse azimuth = %3i\xF8%3i\' %6.3f\"\n" % ( deg, min, sec ))


	# Test the direct function */
	latitude1 = -(( 3.7203 / 60. + 57) / 60. + 37 )
	longitude1 = ( 29.5244 / 60. + 25) / 60. + 144
	dist = 54972.271
	alpha1Tp2 = ( 5.37 / 60. + 52) / 60. + 306
	latitude2 = longitude2 = 0.0
	alpha21 = 0.0

	latitude2, longitude2, alpha21 = calculateGeographicalPositionFromRangeBearing (latitude1, longitude1, alpha1Tp2, dist )

	print ("\n Projected point =%11.6f, %13.6f \n" % ( latitude2, longitude2 ))
	deg = int(latitude2)
	min = int(abs( ( latitude2 - deg) * 60.0 ))
	sec = abs( latitude2 * 3600 - deg * 3600) - min * 60
	print (" Projected Point = %3i\xF8%3i\' %6.3f\", " % ( deg, min, sec ),)
	deg = int(longitude2)
	min = int(abs( ( longitude2 - deg) * 60.0 ))
	sec = abs(longitude2 * 3600 - deg * 3600) - min * 60
	print ("  %3i\xF8%3i\' %6.3f\"\n" % ( deg, min, sec ))
	print (" Should be Buninyon \n" )
	print ("\n Reverse azimuth = %10.6f \n" % ( alpha21 ))
	deg = int(alpha21)
	min = int(abs( ( alpha21 - deg) * 60.0 ))
	sec = abs(alpha21 * 3600 - deg * 3600) - min * 60
	print (" Reverse azimuth = %3i\xF8%3i\' %6.3f\"\n\n" % ( deg, min, sec ))


	# test the EPSG transformer
	from pyproj import Transformer
	from pyproj import CRS

	dlon = 33.0
	dlat = 40.0

	projWGS84=CRS.from_epsg(4326)
	projUTM36N=CRS.from_epsg(32636)
	proj22780="+proj=sterea +lat_0=34.2 +lon_0=39.15 +k=0.9995341 +x_0=0 +y_0=0 +a=6378249.2 +b=6356515 +towgs84=-190.421,8.532,238.69,0,0,0,0 +units=m +no_defs"
	#sourced from https://epsg.io/22780 
	#alternatively use:
	#proj22780=CRS.from_epsg(22780)

	print ("WGS84",dlon,dlat)

	transformer = Transformer.from_crs(projWGS84, projUTM36N)
	x,y=transformer.transform(dlon,dlat)
	print ("UTM36N",x,y)

	transformer = Transformer.from_crs(projUTM36N, proj22780)
	x,y=transformer.transform(x,y)
	print ("EPSG22780",x,y)
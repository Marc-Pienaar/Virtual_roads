#!/usr/bin/python

from osgeo import gdal, ogr, osr,gdalconst,gdal_array
import time
import numpy as np
import os
from os import path
import math 
from sklearn.metrics.pairwise import euclidean_distances
import geopandas as gpd
import skimage
from skimage.graph import MCP_Geometric
import sys



gdal.UseExceptions()
gdal.SetConfigOption("GDAL_DATA", "gdal-data/")

''' Small function to derive slope from a DEM using gdal''' 

def calculate_slope(DEM,slope_path):
	start_time = time.time()
	slope=None
	if os.path.exists(slope_path):
		ds = gdal.Open(slope_path)
		slope = np.array(ds.GetRasterBand(1).ReadAsArray())
		ds = None
	else:
		gdal.DEMProcessing(slope_path, DEM, 'slope',scale = 111120)
		ds = gdal.Open(slope_path)
		slope = np.array(ds.GetRasterBand(1).ReadAsArray())
		ds = None
	print("--- %s seconds for slope creation and array extraction---" % (time.time() - start_time))#~ 3 seconds ---
	return slope

def assign_threshold(cost_surface,slope,slope_threshold,na):
	"""
	@ Author Marc Pienaar
	Replaces all values >= to the slope_threshold with NA values
	----------
	cost_surface: Numpy array
		Array of cost surface values extracted from a tif

	slope: Numpy array
		Array of slope values ither extracted from a tif or derived using the calculate_slope function

	slope_threshold: Real
		The threshold in degrees (e.g 15) to use as a cutoff value
	
	na:  Real
		The na value to assign all values in the cost surface that are >= to the slope threshold 

	Returns
	-------
	An array of values where number >= to the slope threshold have been replaced by a NA value
	"""
	start_time = time.time()
	cost_surface[cost_surface<0]=na
	slope[slope<0]=0
	slope[slope>=slope_threshold]=na
	cost_surface[slope==na]=na#add slope NAs to cost surface
	print("--- %s seconds for adding slope values to the cost surface---" % (time.time() - start_time))
	return slope
	
def set_obstacles(obstacle_path,cost_surface,ds,na):
	"""
	@ Author Marc Pienaar
	Rasterises obstacles in the obstacle path folder and adds them as NA values to the cost surface.
	----------
	obstacle_path: String
		Directory path for obstacle shapefiles

	cost_surface: Numpy array
		Array of cost surface values extracted from a tif

	ds: GDAL raster datasource object 
		The original Digital Elevation Model (DEM) to extract parameters from
	
	na:  Real
		The na value to assign all values in the cost surface that are >= to the slope threshold 

	Returns
	-------
	An array of values where number obstacles have been rasterised and corresponding values in the cost surface replaced by a NA value
	"""
	# get some paramteres like extent, projection, and resolution
	GeoTransform = ds.GetGeoTransform()
	Projection = ds.GetProjection()
	x_min, xres, xskew, y_max, yskew, yres = GeoTransform
	x_max = x_min + (ds.RasterXSize * xres)
	y_min = y_max + (ds.RasterYSize * yres)
	x_res = ds.RasterXSize
	y_res = ds.RasterYSize
	pixel_width = xres
	obstacles = []
	for root, dirs, files in os.walk(obstacle_path):
		for file in files:
			if file.endswith(".shp"):
				obstacles.append(os.path.join(root, file))				
	#add some obstacles to the cost surface
	accum=0
	for i in obstacles:
		obs_out="outputs/obs_"+str(accum)+".tif"
		accum= accum+1
		shpDriver = ogr.GetDriverByName("ESRI Shapefile")
		source_ds = ogr.Open(i)  # open the original shp
		mb_l = source_ds.GetLayer()
		target_ds = gdal.GetDriverByName('GTiff').Create(obs_out, x_res, y_res, 1, gdal.GDT_Float32)
		target_ds.SetGeoTransform(GeoTransform)
		target_ds.SetProjection(Projection)
		bandlist = target_ds.GetRasterBand(1)
		bandlist.SetNoDataValue(na)
		gdal.RasterizeLayer(target_ds, [1], mb_l, burn_values=[0])
		obss=np.array(target_ds.GetRasterBand(1).ReadAsArray())
		cost_surface[obss==0]=na
		source_ds = None
		target_ds = None
		try:#remove the obstacle files
			if os.path.isfile(obs_out):
				os.remove(obs_out)
		except:
			pass
	return cost_surface


def BSU_template(inraster,outraster,res,res2):
	"""
	@ Author Marc Pienaar
	creates a nested version of a blank input grid 
	----------
	inraster: String
		Directory path for a raster

	outraster: String
		Directory path for an output raster

	res: Real 
		a resolution in meters to set the bsu templat to
	
	res2: Real
		a resolution in meters to set the bsu templat to

	Returns
	-------
	A nested version of an equal area intput grid
	"""
	start_time = time.time()
	ds = gdal.Open(inraster)
	GeoTransform = ds.GetGeoTransform()
	Projection = ds.GetProjection()
	x_min, xres, xskew, y_max, yskew, yres = GeoTransform
	x_res = ds.RasterXSize
	y_res = ds.RasterYSize	
	x_max = x_min + (ds.RasterXSize * xres)
	y_min = y_max + (ds.RasterYSize * yres)
	gtLst = list(GeoTransform)
	gtLst[1] = res2
	gtLst[5] = -res2
	out_ds = gdal.Translate(outraster, ds,format='GTiff', width=math.ceil(x_res/res), height=math.ceil(y_res/res)+1)#e.g. 100m resolution re sampling
	x = np.array([[0,1],[1,0]])
	data=np.tile(x, (int(math.ceil(y_res/res)/2),int(math.ceil(x_res/res)/2)))
	out_ds=None
	del out_ds	
	#write binary values so that it is easier to see the grid
	out_ds=gdal.Open(outraster,gdalconst.GF_Write)
	out_ds.GetRasterBand(1).WriteArray( data )
	out_ds.SetGeoTransform(gtLst)
	out_ds=None
	del out_ds	
	ds=None
	print("--- %s seconds to create template---" % (time.time() - start_time))#~ 60 seconds ---	

def rasterise_shp(inputraster,shapein, rasterout,na):
	"""
	@ Author Marc Pienaar
	Rasterises a shapefile into a raster grid
	----------
	inputraster: String
		Directory path for a raster template to use

	shapein: String
		Directory path for a shapefile to rasterise 

	rasterout: String
		Directory path for the output rasterised shape
	
	na:  Real
		The na value to assign all values in the cost surface that are >= to the slope threshold 

	Returns
	-------
	A rastersied shapefile saved as rasterout
	"""	
	start_time = time.time()
	ds = gdal.Open(inputraster)
	### get extent
	GeoTransform = ds.GetGeoTransform()
	Projection = ds.GetProjection()
	x_min, xres, xskew, y_max, yskew, yres = GeoTransform
	x_max = x_min + (ds.RasterXSize * xres)
	y_min = y_max + (ds.RasterYSize * yres)
	x_res = ds.RasterXSize
	y_res = ds.RasterYSize
	pixel_width = xres
	targetprj = osr.SpatialReference(wkt=Projection)
	####add the roads (create a seondary cost surface)
	shpDriver = ogr.GetDriverByName("ESRI Shapefile")
	source_ds = ogr.Open(shapein)  # open the original shp
	mb_l = source_ds.GetLayer()
	target_ds = gdal.GetDriverByName('GTiff').Create(rasterout, x_res, y_res, 1, gdal.GDT_Float32)	
	target_ds.SetGeoTransform(GeoTransform)
	target_ds.SetProjection(Projection)
	bandlist = target_ds.GetRasterBand(1)
	bandlist.SetNoDataValue(na)
	gdal.RasterizeLayer(target_ds, [1], mb_l, burn_values=[0])
	source_ds = None
	target_ds = None
	ds=None
	print("rasterisation took --- %s seconds ---" % (time.time() - start_time))
	
def raster_to_point(inraster, shpout):
	"""
	@ Author Marc Pienaar
	Converts valid raster cells to a point shapefile
	----------
	inraster: String
		The input raster object to use

	shpout: String
		Directory path for a shapefile output file 	

	Returns
	-------
	A point shapefile from a raster input file 
	"""		
	start_time = time.time()
	array =raster2array(inraster)
	gridlatitude, gridlongitude, gridvalues, gridvaluesX, gridvaluesY, grid_coords = array2shp(array, inraster)
	#create a new shapefile
	# set up the shapefile driver
	driver = ogr.GetDriverByName("ESRI Shapefile")
	## create the data source
	data_source = driver.CreateDataSource(shpout)
	### create the spatial reference, WGS84
	srs = osr.SpatialReference()
	ds = gdal.Open(inraster)
	srs.ImportFromWkt( ds.GetProjectionRef() )
	ds=None
	##srs.ImportFromEPSG(4326)
	## create the layer
	layer = data_source.CreateLayer("points", srs, ogr.wkbPoint)
	layer.CreateField(ogr.FieldDefn("Latitude", ogr.OFTReal))
	layer.CreateField(ogr.FieldDefn("Longitude", ogr.OFTReal))
	for i in range(0,len(gridlatitude)):
		feature = ogr.Feature(layer.GetLayerDefn())
		feature.SetField("Latitude", gridlatitude[i])
		feature.SetField("Longitude", gridlongitude[i])
		wkt = "POINT(%f %f)" %  (float(gridlongitude[i]) , float(gridlatitude[i]))
		# Create the point from the Well Known Txt
		point = ogr.CreateGeometryFromWkt(wkt)
		# Set the feature geometry using the point
		feature.SetGeometry(point)
		# Create the feature in the layer (shapefile)
		layer.CreateFeature(feature)
		# Dereference the feature
		feature = None
	# Save and close the data source
	data_source = None
	print("point layer took --- %s seconds ---" % (time.time() - start_time))


def raster2array(rasterfn):	
	##return an array of values from a raster	
	raster = gdal.Open(rasterfn)
	band = raster.GetRasterBand(1)
	array = band.ReadAsArray()
	raster = None
	return array

def array2shp(array, rasterfn):
	##return coordinate and other values from and array input	
	gridlatitude = []
	gridlongitude = []
	gridvalues = []
	gridvaluesX = []
	gridvaluesY = []
	grid_coords = []
	field_ids = []
	raster = gdal.Open(rasterfn)
	row_count = array.shape[0]
	for ridx, row in enumerate(array):		
		for cidx, value in enumerate(row):
			if value >= 0:
				Xcoord, Ycoord =pixelOffset2coord(raster, cidx, ridx)
				gridlatitude.append(Ycoord)
				gridlongitude.append(Xcoord)
				gridvalues.append(int(value))
				gridvaluesX.append(cidx)
				gridvaluesY.append(ridx)
				grid_coords.append([ridx, cidx])
			else:
				pass
	raster = None
	return gridlatitude, gridlongitude, gridvalues, gridvaluesX, gridvaluesY, grid_coords


def pixelOffset2coord(raster, xOffset, yOffset):
	#Get the center of raster cells
	geotransform = raster.GetGeoTransform()
	originX = geotransform[0]
	originY = geotransform[3]
	pixelWidth = geotransform[1]
	pixelHeight = geotransform[5]
	coordX = (originX + pixelWidth / 2) + pixelWidth * xOffset
	coordY = (originY + pixelHeight / 2) + pixelHeight * yOffset
	return coordX, coordY

def reprojectshape_point(inShp,inshp2,outShp,layername):
	"""
	@ Author Marc Pienaar
	eprojects a shapefile 
	----------
	inShp: String
		The input point shapefile projection to use

	inShp2: String
		The input point shapefile to reproject

	layername: String
		arbritary layer name to assign to the output shapefile

	outShp: String
		Directory path for a shapefile output file 

	Returns
	-------
	Reprojects a shapefile 
	"""		
	start_time = time.time()
	driver = ogr.GetDriverByName('ESRI Shapefile')	
	
	# Remove output shapefile if it already exists
	if os.path.exists(outShp):
		driver.DeleteDataSource(outShp)
		
	# get the input layer
	inDataSet = driver.Open(inshp2)
	inLayer = inDataSet.GetLayer()	
	srs = inLayer.GetSpatialRef()	
	
	inDataSet2 = driver.Open(inShp)
	inLayer2 = inDataSet2.GetLayer()	
	srs2 = inLayer2.GetSpatialRef()
	# create the CoordinateTransformation
	coordTrans = osr.CoordinateTransformation(srs, srs2)
	# create the output layer
	outputShapefile = outShp
	if os.path.exists(outputShapefile):
		driver.DeleteDataSource(outputShapefile)
	outDataSet = driver.CreateDataSource(outputShapefile)
	outLayer = outDataSet.CreateLayer(layername,srs2,geom_type=ogr.wkbPoint)
#	outLayer = outDataSet.CreateLayer(outShp,srs)
	# add fields
	inLayerDefn = inLayer.GetLayerDefn()
	for i in range(0, inLayerDefn.GetFieldCount()):
		fieldDefn = inLayerDefn.GetFieldDefn(i)
		outLayer.CreateField(fieldDefn)
#		print(fieldDefn)
	# get the output layer's feature definition
	outLayerDefn = outLayer.GetLayerDefn()
	# loop through the input features
	inFeature = inLayer.GetNextFeature()
	while inFeature:
		# get the input geometry
		geom = inFeature.GetGeometryRef()
		# reproject the geometry
		geom.Transform(coordTrans)
		# create a new feature
		outFeature = ogr.Feature(outLayerDefn)
		# set the geometry and attribute
#		bufferDistance = 500
#		poly = geom.Buffer(geom)
		outFeature.SetGeometry(geom)
		for i in range(0, outLayerDefn.GetFieldCount()):
			outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
#			print(outLayerDefn.GetFieldDefn(i).GetNameRef())
		# add the feature to the shapefile
		outLayer.CreateFeature(outFeature)
		# dereference the features and get the next input feature
		outFeature = None
		inFeature = inLayer.GetNextFeature()
		
	# Save and close the shapefiles
	inDataSet2 = None
	inDataSet = None
	outDataSet = None
	print("--- %s seconds to reproject---" % (time.time() - start_time))#~ 60 seconds ---
	

def distance_in_meters(x1,y1,x2,y2):
	#small function to convert distance between degree points to meters
	def rad(x):
		return(x*math.pi/180.)	
	# approximate radius of earth in km
	R = 6373.0
	lat1 = rad(y1)
	lon1 = rad(x1)
	lat2 = rad(y2)
	lon2 = rad(x2)
	dlon = lon2 - lon1
	dlat = lat2 - lat1	
	a = math.sin(dlat / 2)**2 +  math.cos(lat1) *  math.cos(lat2) *  math.sin(dlon / 2)**2
	c = 2 *  math.atan2( a ** 0.5,  (1 - a)**0.5)	
	distance = R * c
	return(distance)

def copy_point_shapefile(inShp,outShp):
	"""
	@ Author Marc Pienaar
	makes a copy of a point shapefile
	----------
	inShp: String
		The input point shapefile to use
	
	outShp: String
		Directory path for a shapefile output file copy

	Returns
	-------
	A point shapefile copy 
	"""			
	driver = ogr.GetDriverByName('ESRI Shapefile')
	# Remove output shapefile if it already exists
	if os.path.exists(outShp):
		driver.DeleteDataSource(outShp)
	# get the input layer
	inDataSource = driver.Open(inShp)
	inLayer = inDataSource.GetLayer()	
	srs = inLayer.GetSpatialRef()
	# Create the output Layer
	outShapefile = outShp
	outDriver = ogr.GetDriverByName("ESRI Shapefile")
	# Create the output shapefile
	outDataSource = outDriver.CreateDataSource(outShapefile)
	outLayer = outDataSource.CreateLayer("points", srs, geom_type=ogr.wkbPoint)
	# Add input Layer Fields to the output Layer
	inLayerDefn = inLayer.GetLayerDefn()
	for i in range(0, inLayerDefn.GetFieldCount()):
		fieldDefn = inLayerDefn.GetFieldDefn(i)
		outLayer.CreateField(fieldDefn)
	# Get the output Layer's Feature Definition
	outLayerDefn = outLayer.GetLayerDefn()
	# Add features to the ouput Layer
	for i in range(0, inLayer.GetFeatureCount()):
		# Get the input Feature
		inFeature = inLayer.GetFeature(i)
		# Create output Feature
		outFeature = ogr.Feature(outLayerDefn)
		# Add field values from input Layer
		for i in range(0, outLayerDefn.GetFieldCount()):
			outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
		# Set geometry as centroid
		geom = inFeature.GetGeometryRef()
#		centroid = geom.Centroid()
		outFeature.SetGeometry(geom)
		# Add new feature to output Layer
		outLayer.CreateFeature(outFeature)
	# Close DataSources
	inDataSource.Destroy()
	outDataSource.Destroy()
	
def copy_polygon_shapefile(inShp,outShp):
	"""
	@ Author Marc Pienaar
	makes a copy of a polygon shapefile
	----------
	inShp: String
		The input point shapefile to use
	
	outShp: String
		Directory path for a shapefile output file copy

	Returns
	-------
	A polygon shapefile copy 
	"""		
	driver = ogr.GetDriverByName('ESRI Shapefile')
	# Remove output shapefile if it already exists
	if os.path.exists(outShp):
		driver.DeleteDataSource(outShp)
	# get the input layer
	inDataSource = driver.Open(inShp)
	inLayer = inDataSource.GetLayer()	
	srs = inLayer.GetSpatialRef()	
	# Create the output Layer
	outShapefile = outShp
	outDriver = ogr.GetDriverByName("ESRI Shapefile")	
	# Create the output shapefile
	outDataSource = outDriver.CreateDataSource(outShapefile)
	outLayer = outDataSource.CreateLayer("points", srs, geom_type=ogr.wkbPolygon)	
	# Add input Layer Fields to the output Layer
	inLayerDefn = inLayer.GetLayerDefn()
	for i in range(0, inLayerDefn.GetFieldCount()):
		fieldDefn = inLayerDefn.GetFieldDefn(i)
		outLayer.CreateField(fieldDefn)		
	# Get the output Layer's Feature Definition
	outLayerDefn = outLayer.GetLayerDefn()	
	# Add features to the ouput Layer
	for i in range(0, inLayer.GetFeatureCount()):
		# Get the input Feature
		inFeature = inLayer.GetFeature(i)
		# Create output Feature
		outFeature = ogr.Feature(outLayerDefn)
		# Add field values from input Layer
		for i in range(0, outLayerDefn.GetFieldCount()):
			outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
		# Set geometry as centroid
		geom = inFeature.GetGeometryRef()
#		centroid = geom.Centroid()
		outFeature.SetGeometry(geom)
		# Add new feature to output Layer
		outLayer.CreateFeature(outFeature)		
	# Close DataSources
	inDataSource.Destroy()
	outDataSource.Destroy()
	
def clip_point(inShp,ROI,outShp):
	"""
	@ Author Marc Pienaar
	Clips a point shapefile to a ROI
	----------
	inShp: String
		The input point shapefile to use

	ROI: String
		The input shapefile path to use for clipping
	
	outShp: String
		Directory path for a shapefile output file copy

	Returns
	-------
	A clipped shapefile 
	"""		
	start_time = time.time()
	driver = ogr.GetDriverByName('ESRI Shapefile')
	# Remove output shapefile if it already exists
	if os.path.exists(outShp):
		driver.DeleteDataSource(outShp)
	# get the input layer
	inDataSet = driver.Open(inShp)
	inLayer = inDataSet.GetLayer()	
	srs = inLayer.GetSpatialRef()	
	inDataSet2 = driver.Open(ROI)
	inLayer2 = inDataSet2.GetLayer()	
	srs2 = inLayer2.GetSpatialRef()
	outputShapefile = outShp
	if os.path.exists(outputShapefile):
		driver.DeleteDataSource(outputShapefile)
	outDataSet = driver.CreateDataSource(outputShapefile)
	outLayer = outDataSet.CreateLayer("points",srs,geom_type=ogr.wkbPoint)	
	inLayer.Intersection(inLayer2, outLayer)
	inDataSet2 = None
	inDataSet = None
	outDataSet = None	
	#now update the inlayer and delete the temp outlayer
	copy_point_shapefile(outShp,inShp)
	# Remove output shapefile if it already exists
	if os.path.exists(outShp):
		driver.DeleteDataSource(outShp)
	print("--- %s seconds to reproject---" % (time.time() - start_time))
	
def clip_poly(inShp,ROI,outShp):
	"""
	@ Author Marc Pienaar
	Clips a polygon shapefile to a ROI
	----------
	inShp: String
		The input point shapefile to use

	ROI: String
		The input shapefile path to use for clipping
	
	outShp: String
		Directory path for a shapefile output file copy

	Returns
	-------
	A clipped shapefile 
	"""		
	start_time = time.time()
	driver = ogr.GetDriverByName('ESRI Shapefile')
	# Remove output shapefile if it already exists
	if os.path.exists(outShp):
		driver.DeleteDataSource(outShp)
	# get the input layer
	inDataSet = driver.Open(inShp)
	inLayer = inDataSet.GetLayer()	
	srs = inLayer.GetSpatialRef()	
	
	inDataSet2 = driver.Open(ROI)
	inLayer2 = inDataSet2.GetLayer()	
	srs2 = inLayer2.GetSpatialRef()
	
	outputShapefile = outShp
	if os.path.exists(outputShapefile):
		driver.DeleteDataSource(outputShapefile)
	outDataSet = driver.CreateDataSource(outputShapefile)
	outLayer = outDataSet.CreateLayer("poly",srs,geom_type=ogr.wkbPolygon)	
	inLayer.Intersection(inLayer2, outLayer)
	inDataSet2 = None
	inDataSet = None
	outDataSet = None
	
	#now update the inlayer and delete the temp outlayer
#	copy_polygon_shapefile(outShp,inShp)
#	# Remove output shapefile if it already exists
#	if os.path.exists(outShp):
#		driver.DeleteDataSource(outShp)
		
		
	print("--- %s seconds to reproject---" % (time.time() - start_time))#~ 60 seconds ---
	
def get_current_window(centoids_in,shapeout_temp,numpoints):
	"""
	@ Author Marc Pienaar
	Create a temporary window (clipped shapefile) from a number of points using a Eudlidean disstance function, for use during analysis
	----------
	centoids_in: String
		The input point shapefile to use

	shapeout_temp: String
		The temporary shapefile path to use for a temp window
	
	numpoints: int
		the number of points from the centoids_in layer to use to define the window size

	Returns
	-------
	A temporary clipped shapefile
	"""		
	#get all the points into an array
	x1=[]
	y1=[]
	shpDriver = ogr.GetDriverByName("ESRI Shapefile")
	source_ds = ogr.Open(centoids_in)  # open the original shp
	layer = source_ds.GetLayer()
	srs = layer.GetSpatialRef()
	featureCount = layer.GetFeatureCount()
	for feature in layer:
		point = feature.GetGeometryRef()
		x1.append(point.GetX())
		y1.append(point.GetY())
	source_ds=None
	
	X=np.column_stack((x1,y1)) #create a numpy array
	#get distance to first centoid
	dist=euclidean_distances(X[0].reshape(1, -1),X)
	#sort the array
	dist2 = np.argsort(dist.flatten())
	x2=[]
	y2=[]
	if len(x1)<numpoints:
		for i in range(0,math.ceil(len(x1))):
			x2.append(x1[dist2[i]])
			y2.append(y1[dist2[i]])
	else:
		for i in range(0,math.ceil(numpoints)):
			x2.append(x1[dist2[i]])
			y2.append(y1[dist2[i]])
	#now create a new shapefile of the closest points 
	driver = ogr.GetDriverByName('ESRI Shapefile')
	
	if os.path.exists(shapeout_temp):
		driver.DeleteDataSource(shapeout_temp)
		
	outDataSet = driver.CreateDataSource(shapeout_temp)
	outLayer = outDataSet.CreateLayer("points",srs,geom_type=ogr.wkbPoint)	
	# Add an ID field
	idField = ogr.FieldDefn("id", ogr.OFTInteger)
	outLayer.CreateField(idField)
	for i in range(0,len(x2)):
		# Create the feature and set values
		featureDefn = outLayer.GetLayerDefn()
		feature = ogr.Feature(featureDefn)
		feature.SetField("id", i)
		wkt = "POINT(%f %f)" %  (float(x2[i]) , float(y2[i]))
		# Create the point from the Well Known Txt
		point = ogr.CreateGeometryFromWkt(wkt)
		# Set the feature geometry using the point
		feature.SetGeometry(point)
		outLayer.CreateFeature(feature)
		feature = None	
	#close the resource
	outDataSet=None
	
def extent_to_shape(inShapefile,outShapefile,buffer):
	"""
	@ Author Marc Pienaar
	creates a shapefile from an extent
	----------
	inShapefile: String
		The input shapefile to use for the extent

	outShapefile: String
		The temporary shapefile path to use for a the extent shapefile
	
	buffer: real
		the buffer in degrees to add to the extent bounds

	Returns
	-------
	A temporary polygin shapefile representing an extent
	"""		
	# Get a Layer's Extent
#	inShapefile = "states.shp"
	inDriver = ogr.GetDriverByName("ESRI Shapefile")
	inDataSource = inDriver.Open(inShapefile, 0)
	inLayer = inDataSource.GetLayer()
	extent = inLayer.GetExtent()
	
	srs = inLayer.GetSpatialRef()
#	outDataSource = inDriver.CreateDataSource(outShapefile)
#	outLayer = outDataSource.CreateLayer( out_lyr_name, srs, geom_type=ogr.wkbMultiLineString )
#
	# Create a Polygon from the extent tuple
	ring = ogr.Geometry(ogr.wkbLinearRing)
	ring.AddPoint(extent[0]-buffer,extent[2]-buffer)
	ring.AddPoint(extent[1]+buffer, extent[2]-buffer)
	ring.AddPoint(extent[1]+buffer, extent[3]+buffer)
	ring.AddPoint(extent[0]-buffer, extent[3]+buffer)
	ring.AddPoint(extent[0]-buffer,extent[2]+buffer)
	poly = ogr.Geometry(ogr.wkbPolygon)
	poly.AddGeometry(ring)
	
	# Save extent to a new Shapefile
#	outShapefile = "states_extent.shp"
	outDriver = ogr.GetDriverByName("ESRI Shapefile")
	
	# Remove output shapefile if it already exists
	if os.path.exists(outShapefile):
		outDriver.DeleteDataSource(outShapefile)
		
	# Create the output shapefile
	outDataSource = outDriver.CreateDataSource(outShapefile)
	outLayer = outDataSource.CreateLayer("extent", srs,geom_type=ogr.wkbPolygon)
	
	# Add an ID field
	idField = ogr.FieldDefn("id", ogr.OFTInteger)
	outLayer.CreateField(idField)
	
	# Create the feature and set values
	featureDefn = outLayer.GetLayerDefn()
	feature = ogr.Feature(featureDefn)
	feature.SetGeometry(poly)
	feature.SetField("id", 1)
	outLayer.CreateFeature(feature)
	
	# Close DataSource
	inDataSource.Destroy()
	outDataSource.Destroy()
	
def crop_raster_to_shape(raster_in,shape,raster_out):
	"""
	@ Author Marc Pienaar
	Crops a raster using a input shapefile
	----------
	raster_in: String
		The input rater path to use 

	shape: String
		The shapefile path to use 
	
	raster_out: String
		The raster output path
	Returns
	-------
	A cropped raster 
	"""		
	start_time = time.time()
#	if os.path.exists(raster_out):
#		pass
#	else:
	driverName = "ESRI Shapefile"
	drv = ogr.GetDriverByName(driverName)
	dataSource_in = drv.Open(shape,0)  # 0 means read-only. 1 means writeable.
	layer = dataSource_in.GetLayer(0)
	featureCount = layer.GetFeatureCount()
	targetprj = layer.GetSpatialRef()
	minXl, maxXl, minYl, maxYl = layer.GetExtent()
	dataSource_in.Destroy()  # source
	#readraster
	ds=gdal.Open(raster_in)
	geotransform = ds.GetGeoTransform()
	prj=ds.GetProjection()
	x_res = ds.RasterXSize
	y_res = ds.RasterYSize	
	
	sourceprj=osr.SpatialReference(wkt=prj)
	ds = None
	transform = osr.CoordinateTransformation(targetprj,sourceprj)
	#perform_crop
	out_ds=gdal.Warp(raster_out,raster_in,dstSRS=sourceprj,cutlineDSName=shape,cropToCutline=True)
	out_ds=None
	print("cropping took --- %s seconds ---" % (time.time() - start_time))
	
def point_to_buffer(inShp,outshapetemp,layername,bufferDistance):
	"""
	@ Author Marc Pienaar
	converts a point shapefile into a polygon shapefile using a small buffer
	----------
	inShp: String
		The input shapefile path to use 

	outshapetemp: String
		The shapefile path to write out to 
	
	layername: String
		arbritarty layer anme to assign 

	bufferDistance: Real
		the distance in degrees to add a buffer to a point
	Returns
	-------
	A polygon shapefile from a point shapefile
	"""		
	start_time = time.time()
	driver = ogr.GetDriverByName('ESRI Shapefile')
	# Remove output shapefile if it already exists
	if os.path.exists(outshapetemp):
		driver.DeleteDataSource(outshapetemp)
	# get the input layer
	inDataSet = driver.Open(inShp)
	inLayer = inDataSet.GetLayer()	
	
	trg = osr.SpatialReference()
	trg.ImportFromEPSG(4326)
	srs = inLayer.GetSpatialRef()
	# create the CoordinateTransformation
	coordTrans = osr.CoordinateTransformation(srs, trg)
#	coordTrans.MorphToESRI()
		# create the output layer
	outputShapefile = outshapetemp
	if os.path.exists(outputShapefile):
		driver.DeleteDataSource(outputShapefile)
	outDataSet = driver.CreateDataSource(outputShapefile)
	outLayer = outDataSet.CreateLayer(layername,srs,geom_type=ogr.wkbMultiPolygon)
#	outLayer = outDataSet.CreateLayer(outShp,srs)
	# add fields
	inLayerDefn = inLayer.GetLayerDefn()
	for i in range(0, inLayerDefn.GetFieldCount()):
		fieldDefn = inLayerDefn.GetFieldDefn(i)
		outLayer.CreateField(fieldDefn)
	# get the output layer's feature definition
	outLayerDefn = outLayer.GetLayerDefn()
	# loop through the input features
	inFeature = inLayer.GetNextFeature()
	while inFeature:
		# get the input geometry
		geom = inFeature.GetGeometryRef()
		# reproject the geometry
#		geom.Transform(coordTrans)
		# create a new feature
		outFeature = ogr.Feature(outLayerDefn)
		# set the geometry and attribute
#		bufferDistance = 500
		poly = geom.Buffer(bufferDistance)
		outFeature.SetGeometry(poly)
		for i in range(0, outLayerDefn.GetFieldCount()):
			outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
		# add the feature to the shapefile
		outLayer.CreateFeature(outFeature)
		# dereference the features and get the next input feature
		outFeature = None
		inFeature = inLayer.GetNextFeature()
		
	# Save and close the shapefiles
	inDataSet = None
	outDataSet = Non
	print("--- %s seconds to reproject---" % (time.time() - start_time))#~ 60 seconds ---
	
def point_shapefile_from_array(shapfile_out,x_coords,y_coords):
	"""
	@ Author Marc Pienaar
	converts an array to a point shapefile 
	----------
	shapfile_out: String
		The output shapefile to write to

	x_coords: array
		array of x coordinates
	
	y_coords: array
		array of y coordinates 
	
	Returns
	-------
	A point shapefile from an array of coordinates
	"""		
	# set up the shapefile driver
	driver = ogr.GetDriverByName("ESRI Shapefile")
	# create the data source
	data_source = driver.CreateDataSource(shapfile_out)
	# create the spatial reference, WGS84
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(4326)
	# create the layer
	layer = data_source.CreateLayer("points", srs, ogr.wkbPoint)
	#layer.CreateField(ogr.FieldDefn("Latitude", ogr.OFTReal))
	layer.CreateField(ogr.FieldDefn("Lat", ogr.OFTReal))
	layer.CreateField(ogr.FieldDefn("Lon", ogr.OFTReal))
	for i in range(0,len(x_coords)):
		feature = ogr.Feature(layer.GetLayerDefn())
		feature.SetField("Lat", y_coords[i])
		feature.SetField("Lon", x_coords[i])
		wkt = "POINT(%f %f)" %  (float(x_coords[i]) , float(y_coords[i]))
		# Create the point from the Well Known Txt
		point = ogr.CreateGeometryFromWkt(wkt)
		# Set the feature geometry using the point
		feature.SetGeometry(point)
		# Create the feature in the layer (shapefile)
		layer.CreateFeature(feature)
		# Dereference the feature
		feature = None
	# Save and close the data source
	data_source = None
	
def combined_point_shapefile_from_array(shapefile_out,x_coords1,y_coords1,x_coords2,y_coords2):
	"""
	@ Author Marc Pienaar
	converts an two array's to a point shapefile 
	----------
	shapfile_out: String
		The output shapefile to write to

	x_coords1: array
		array of x coordinates
	
	y_coords1: array
		array of y coordinates 

	x_coords2: array
		array of x coordinates
	
	y_coords2: array
		array of y coordinates 
	
	Returns
	-------
	A point shapefile from an array of coordinates
	"""		
	
	#combine the two to get an extent
	#create a new shapefile
	# set up the shapefile driver
	driver = ogr.GetDriverByName("ESRI Shapefile")
	# create the data source
	data_source = driver.CreateDataSource(shapefile_out)
	# create the spatial reference, WGS84
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(4326)
	# create the layer
	layer = data_source.CreateLayer("points", srs, ogr.wkbPoint)
	layer.CreateField(ogr.FieldDefn("Lat", ogr.OFTReal))
	layer.CreateField(ogr.FieldDefn("Lon", ogr.OFTReal))
	for i in range(0,len(x_coords1)):
		feature = ogr.Feature(layer.GetLayerDefn())
		feature.SetField("Lat", y_coords1[i])
		feature.SetField("Lon", x_coords1[i])
		wkt = "POINT(%f %f)" %  (float(x_coords1[i]) , float(y_coords1[i]))
		# Create the point from the Well Known Txt
		point = ogr.CreateGeometryFromWkt(wkt)
		# Set the feature geometry using the point
		feature.SetGeometry(point)
		# Create the feature in the layer (shapefile)
		layer.CreateFeature(feature)
		# Dereference the feature
		feature = None	
	for i in range(0,len(x_coords2)):
		feature = ogr.Feature(layer.GetLayerDefn())
		feature.SetField("Lat", y_coords2[i])
		feature.SetField("Lon", x_coords2[i])
		wkt = "POINT(%f %f)" %  (float(x_coords2[i]) , float(y_coords2[i]))
		# Create the point from the Well Known Txt
		point = ogr.CreateGeometryFromWkt(wkt)
		# Set the feature geometry using the point
		feature.SetGeometry(point)
		# Create the feature in the layer (shapefile)
		layer.CreateFeature(feature)
		# Dereference the feature
		feature = None	
	# Save and close the data source
	data_source = None
	
def first_run(J,road_centoids_in,centoids_in2,numpoints_roads,temp_target,shapein,shapecombined,shapeout2,rasterout,rasterout2,cost_path,shapeout,shapein2,ROI):
	"""
	@ Author Marc Pienaar
	The first entry point for the virtual roads algorithm
	
	"""		
	start_time = time.time()
	#get all values from road network
	x_roads=[];y_roads=[]
	shpDriver = ogr.GetDriverByName("ESRI Shapefile")
	source_ds = ogr.Open(road_centoids_in)  # open the original shp
	layer = source_ds.GetLayer()
	for feature in layer:
		geom = feature.GetGeometryRef()
		x_roads.append(geom.GetX())
		y_roads.append(geom.GetY())
	source_ds=None
	X_roads=np.column_stack((x_roads,y_roads)) #create a numpy array
	
	#get origin points
	x1=[];y1=[]
	shpDriver = ogr.GetDriverByName("ESRI Shapefile")
	source_ds = ogr.Open(centoids_in2)  # open the original shp
	layer = source_ds.GetLayer()
	featureCount = layer.GetFeatureCount()
	for feature in layer:
		point = feature.GetGeometryRef()
		x1.append(point.GetX())
		y1.append(point.GetY())	
	source_ds=None
	
	X2=np.column_stack((x1,y1)) #create a numpy array
	dist=euclidean_distances(X2,X_roads)
	#
	x3=[];y3=[]
	for i in range(0,len(dist)):
		dist3 = np.argsort(dist[i].flatten())
		for j in range(0,int(numpoints_roads)):
			x3.append(x_roads[dist3[j]])
			y3.append(y_roads[dist3[j]])		
	
	#create temporary destination points
	point_shapefile_from_array(temp_target,x3,y3)
	#create temporary origin points
	point_shapefile_from_array(shapein,x1,y1)
	#create a combined points layer
	combined_point_shapefile_from_array(shapecombined,x1,y1,x3,y3)
	
	extent_to_shape(shapecombined,shapeout,0.01)#1km buffer
	#make sure the extent is within the roi
	#crop the layers to ROI using a temporary file to copy
	#func.clip_poly(shapeout,ROI,temp)
	shapeout22=gpd.read_file(shapeout)
	ROI2=gpd.read_file(ROI)
	clipped = gpd.clip(shapeout22, ROI2)
	clipped.to_file(shapeout)
	clipped.to_file(shapeout+str(J)+".shp")
	##crop the raster
	crop_raster_to_shape(cost_path,shapeout,rasterout)
	crop_raster_to_shape(cost_path,shapeout,rasterout2)
	crop_raster_to_shape(cost_path,shapeout,rasterout+str(J)+".tif")
	#create a buffer layer
	point_to_buffer(shapein2,shapeout2,"temp",0.002)
	
	
def second_run(numpoints_window,cost_path,shapein,temp_target,rasterout,rasterout2,NGI_road_path,outputpath,n):
	"""
	@ Author Marc Pienaar
	The second entry point for the virtual roads algorithm
	
	"""	
	for k in range(0,int(numpoints_window)):
		try:
			cost_path=rasterout
			cost_path2=rasterout2
			#Generate a cost surface
			ds = gdal.Open(cost_path)
			# get extent
			GeoTransform = ds.GetGeoTransform()
			Projection = ds.GetProjection()
			x_min, xres, xskew, y_max, yskew, yres = GeoTransform
			x_max = x_min + (ds.RasterXSize * xres)
			y_min = y_max + (ds.RasterYSize * yres)
			x_res = ds.RasterXSize
			y_res = ds.RasterYSize
			pixel_width = xres
			targetprj = osr.SpatialReference(wkt=Projection)
			meters = degrees_to_meters(pixel_width)
			cost_surface = np.array(ds.GetRasterBand(1).ReadAsArray())
			
			gridlatitude = [];gridlongitude = [];gridvalues = [];gridvaluesX = [];gridvaluesY = [];grid_coords = []
			#roads
			roadlatitude = [];roadlongitude = [];roadvalues = [];roadvaluesX = [];roadvaluesY = [];road_coords = []
			source_centoids_in=shapein
			target_centoids_in=temp_target
			
			shpDriver = ogr.GetDriverByName("ESRI Shapefile")
			source_ds = ogr.Open(source_centoids_in)  # open the original shp
			layer = source_ds.GetLayer()
			for feature in layer:
				geom = feature.GetGeometryRef()
				gridvaluesX.append(geom.GetX())
				gridvaluesY.append(geom.GetY())
			source_ds=None
			
			source_ds = ogr.Open(target_centoids_in)  # open the original shp
			layer = source_ds.GetLayer()
			for feature in layer:
				geom = feature.GetGeometryRef()
				roadvaluesX.append(geom.GetX())
				roadvaluesY.append(geom.GetY())
			source_ds=None
			
			# create temp cost surface
			path_surface2 = np.copy(cost_surface)
			
			#get the roads file and rasterise to get a point layer
			roads=NGI_road_path
			obs_out="temp/obs_roads_test.tif"
			na=-999999
			source_ds = ogr.Open(roads)  # open the original shp
			mb_l = source_ds.GetLayer()
			target_ds = gdal.GetDriverByName('GTiff').Create(obs_out, x_res, y_res, 1, gdal.GDT_Float32)
			target_ds.SetGeoTransform(GeoTransform)
			target_ds.SetProjection(Projection)
			bandlist = target_ds.GetRasterBand(1)
			bandlist.SetNoDataValue(na)
			gdal.RasterizeLayer(target_ds, [1], mb_l, burn_values=[0])
			obss=np.array(target_ds.GetRasterBand(1).ReadAsArray())
			source_ds=None
			target_ds = None
			
			
			source=source_centoids_in
			path=""
			na=-999999
			source_ds = ogr.Open(source)  # open the original shp
			mb_l = source_ds.GetLayer()
			obs_out2="temp/obs_source_test.tif"
			target_ds = gdal.GetDriverByName('GTiff').Create(obs_out2, x_res, y_res, 1, gdal.GDT_Float32)
			target_ds.SetGeoTransform(GeoTransform)
			target_ds.SetProjection(Projection)
			bandlist = target_ds.GetRasterBand(1)
			bandlist.SetNoDataValue(na)
			gdal.RasterizeLayer(target_ds, [1], mb_l, burn_values=[0])
			source_ds=None
			target_ds = None
			row,col,data=coord_to_pixel(ds,gridvaluesY,gridvaluesX,k)	
			costsi = []
			costsa = []
			coordsa = []
			# create temp cost surface
			path_surface2 = np.copy(cost_surface)
			out_ds=gdal.Open(cost_path2,gdalconst.GF_Write)
			out_ds.GetRasterBand(1).WriteArray( path_surface2 )
			out_ds=None
			#obs_out
			array =raster2array(obs_out)
			roadlatitude, roadlongitude, roadvalues, roadvaluesX, \
			roadvaluesY, road_coords =array2shp(array, obs_out)
			maxx = -1000000000000
			mcp = skimage.graph.MCP_Geometric(path_surface2, fully_connected=True)
			cumulative_costs, traceback = mcp.find_costs([[row, col]])  # start points
			cities = np.array(road_coords)  # end points
			ncities = cities.shape[0]
			paths = np.empty(path_surface2.shape)
			paths.fill(-1)
			costs2 = 0
			optimal_route = []
			i_val = 0
			x = sys.maxsize
			raster=gdal.Open(cost_path2)
			for i in range(ncities):
		#				try:
				cost3 = cumulative_costs[cities[i][0], cities[i][1]]
				if cost3 < x:
					x = min(x, cost3)
					i_val = i
		#				except:
		#					pass
			route = mcp.traceback([cities[i_val, :][0], cities[i_val, :][1]])
			optimal_route = route
			a=[]
			for j in range(0,len(route)):
				costsa.append(costs2)
				coordsa.append(pixel_to_coord(raster, route[j][1], route[j][0]))
			raster=None
			maxx=max(maxx,costs2)
			print(maxx/1000)
			num="00000"
			print(k+n)
			numList2 = [outputpath, "/",num,"_", str(int(k+n)), ".shp"]
			numList3 = [outputpath, "/",num,"_", str(int(k+n)), ".prj"]
			filename2 = ''.join(numList2)
			filename3 = ''.join(numList3)
			driver = ogr.GetDriverByName('Esri Shapefile')
			ds2 = driver.CreateDataSource(filename2)
			srs = osr.SpatialReference()
			srs.ImportFromEPSG(4326)
			srs.MorphToESRI()
			file = open(filename3, 'w')
			file.write(srs.ExportToWkt())
			file.close()
			# ds.SetProjection(srs.ExportToWkt())
			layer = ds2.CreateLayer('path', geom_type=ogr.wkbLineString)
			line = ogr.Geometry(ogr.wkbLineString)
			for i in range(0, len(coordsa)):
				line.AddPoint(coordsa[i][0], coordsa[i][1])
			wkt = line.ExportToWkt()
			geom = ogr.CreateGeometryFromWkt(wkt)
			field_testfield = ogr.FieldDefn("dist_m", ogr.OFTReal)
			field_testfield.SetWidth(50)
			layer.CreateField(field_testfield)
			
			feature = ogr.Feature(layer.GetLayerDefn())
			feature.SetField("dist_m", costs2)
			feature.SetGeometry(geom)
			layer.CreateFeature(feature)
			feature = None
			ds2 = None
			ds=None
		except:
			pass
		

def degrees_to_meters(lat):
	#small utility to get degrees to meters estimate
	import math
	radius = 6371.0088
	lat1 = 0. * math.pi / 180.
	lat2 = 0. * math.pi / 180.
	lon1 = lat * math.pi / 180.
	lon2 = 0. * math.pi / 180.
	# dlat = lat2 - lat1
	dlon = lon2 - lon1
	rad = math.acos(
		(math.sin(lat1) * math.sin(lat2))
		+
		(math.cos(lat2) * math.cos(dlon))
	)
	answerlon = (radius * 1000) * rad;
	if lat < 0:
		answerlon = answerlon - answerlon - answerlon;
	return (answerlon)

def coord_to_pixel(raster, gridvaluesY,gridvaluesX,k):
	#small utility to convert a coordinate to a raster pixel coordinate
	cols = raster.RasterXSize
	rows = raster.RasterYSize
	transform = raster.GetGeoTransform()
	xOrigin = transform[0]
	yOrigin = transform[3]
	pixelWidth = transform[1]
	pixelHeight = -transform[5]
	band = raster.GetRasterBand(1)
	data = band.ReadAsArray(0, 0, cols, rows)
#	print(gridvaluesX[k],xOrigin, "mmmmmm")
	col = int((gridvaluesX[k] - xOrigin) / pixelWidth)
	row = int((yOrigin -gridvaluesY[k] ) / pixelHeight)
	return row,col,data[row][col]

def coord_to_pixel2(raster, gridvaluesY,gridvaluesX):
	#small utility to convert a coordinate to a raster pixel coordinate
	cols = raster.RasterXSize
	rows = raster.RasterYSize
	transform = raster.GetGeoTransform()
	xOrigin = transform[0]
	yOrigin = transform[3]
	pixelWidth = transform[1]
	pixelHeight = -transform[5]
	band = raster.GetRasterBand(1)
	data = band.ReadAsArray(0, 0, cols, rows)
	
	col = int((gridvaluesX - xOrigin) / pixelWidth)
	row = int((yOrigin -gridvaluesY ) / pixelHeight)
	return row,col,data[row][col]



def pixel_to_coord(raster, xOffset, yOffset):
	#small utility to convert araster pixel coordinate to a world coordinate. 
	
	geotransform = raster.GetGeoTransform()
	originX = geotransform[0]
	originY = geotransform[3]
	pixelWidth = geotransform[1]
	pixelHeight = geotransform[5]
	coordX = (originX + pixelWidth / 2) + pixelWidth * xOffset
	coordY = (originY + pixelHeight / 2) + pixelHeight * yOffset
	return coordX, coordY
#!/usr/bin/python

import Virtual_roads_functions as func
from osgeo import gdal, ogr, osr,gdalconst,gdal_array
import time
import numpy as np
import os
from os import path
import pathlib
import math 
from sklearn.metrics.pairwise import euclidean_distances
import geopandas as gpd
import skimage
from skimage.graph import MCP_Geometric
import sys

''' use gdal exceptoions '''
gdal.UseExceptions()
gdal.SetConfigOption('GDAL_DATA', 'gdal-data/')
#creat a start time instance to measure model run performace 
start_time = time.time()
#inputs 
dempath = 'inputs/SRTM_clip.tif' #use a merged SRTM DEM covering the ROI to create this
obstacle_path='inputs/obstacles/'#Collection of obstacle shapefiles to burn into the cost surface. e.g. water bodies, restricted areas, etc. 
NGI_road_path='inputs/roads/NGI_2019_road_clip.shp'#the road layer
inputraster='inputs/BSU_100m_clip.tif'#An equal area rasster to use to get hi-res start and and points
ROI='inputs/ROI.shp'#a ROI for this example
#Create output dirctories
pathlib.Path('outputs').mkdir(parents=True, exist_ok=True) 
pathlib.Path('output_paths').mkdir(parents=True, exist_ok=True) 
pathlib.Path('temp').mkdir(parents=True, exist_ok=True) 
#get a reference to the output_paths directory
outputpath='output_paths'
temp='outputs/temp.shp'
#define outputs
cost_path='outputs/cost.tif'
slope_path = 'outputs/slope.tif'
out_raster50m='outputs/RSA_BSU_50m.tif'
road_out50m='outputs/road_50m_target.tif'
target_out_50='outputs/road_50m_target.shp'
source_out_1000='outputs/1000m_source.shp'
target_out_50wgs='outputs/road_50m_target_wgs84.shp'
source_out_1000wgs='outputs/1000m_source_wgs84.shp'
#define parameters for the model
#set wahtever values you desire here	
na=-999999#assign missing values with this number
slope_threshold=15#10 degrees is the legal limit for trucks in new zealand, we'll use 15 degree max slope threshold assuming it would be navagable by a tracktor trailor combination

#begin the code
#open the DEM file 
ds = gdal.Open(dempath)
# get extent and other info
GeoTransform = ds.GetGeoTransform()
Projection = ds.GetProjection()
x_min, xres, xskew, y_max, yskew, yres = GeoTransform
x_max = x_min + (ds.RasterXSize * xres)
y_min = y_max + (ds.RasterYSize * yres)
x_res = ds.RasterXSize
y_res = ds.RasterYSize
pixel_width = xres
##create a cost surface 
out_ds = gdal.Translate(cost_path, ds,format='GTiff', width=x_res, height=y_res)#1
out_ds=None
#extract the raster data into a cost surface array
cost_surface = np.array(ds.GetRasterBand(1).ReadAsArray())
#generate slope array
slope=func.calculate_slope(dempath,slope_path)
#set slope thresholds
cost_surface=func.assign_threshold(cost_surface,slope,slope_threshold,na)
#get list of obstacles and add them to the cost surface
cost_surface = func.set_obstacles(obstacle_path,cost_surface,ds,na)
#write out the cost surface 
out_ds=gdal.Open(cost_path,gdalconst.GF_Write)
out_ds.GetRasterBand(1).WriteArray( cost_surface )
out_ds=None
ds = None
print('--- %s seconds for creating the cost surface---' % (time.time() - start_time))
print()

#create a directory for the virtual road outputs
start_time = time.time()
func.BSU_template(inputraster,out_raster50m,0.5,50)#e.g 50m res for target points

#rasterise the road layer (e.g. 50m resolution)
func.rasterise_shp(out_raster50m,NGI_road_path,road_out50m,na)
#get target points i.e. 50m road points
func.raster_to_point(road_out50m,target_out_50)

##last step is to re-project to the same projection as the cost surface i.e. WGS84, before we can start running the algorithm
func.reprojectshape_point(NGI_road_path,target_out_50,target_out_50wgs,'DN')
func.reprojectshape_point(NGI_road_path,source_out_1000,source_out_1000wgs,'DN')

#crop the layers to ROI using a temporary file to copy
func.clip_point(source_out_1000wgs,ROI,temp)
func.clip_point(target_out_50wgs,ROI,temp)

#first we duplicate the source points - we are going to itterate over these and remove them
source_out_copy='inputs/updated_points.shp'
func.copy_point_shapefile(source_out_1000wgs, source_out_copy)
#Now we get a subset and process those
shapeout_temp='outputs/window.shp'


shapecombined='temp/temp_combined.shp'
shapein='temp/temp_source.shp'
shapein2='temp/temp_source.shp'
temp_target='temp/temp_target.shp'
shapeout='temp/extent.shp'
shapeout2='temp/temp3.shp'
rasterout='temp/cost_clip.tif'
rasterout2='temp/cost_clip2.tif'

numpoints_roads=250#this is the number of points to process at a time to demonstrate the functionality
#numpoints_window=77#this is the total pointss divided by 3 to use in the outer itteration 
numpoints_window2=11#this is the number of points to process at a time in the inner loop to demonstrate the functionality if it were a larger ccost surface for instance


n=0

for i in range(0,24):
	try:
		func.get_current_window(source_out_copy,shapeout_temp,numpoints_window2)
		
		centoids_in1=target_out_50wgs
		centoids_in2=shapeout_temp
		func.first_run(i,centoids_in1,centoids_in2,numpoints_roads,temp_target,shapein,shapecombined,shapeout2,rasterout,rasterout2,cost_path,shapeout,shapein2,ROI)
		
		func.second_run(numpoints_window2,cost_path,shapein,temp_target,rasterout,rasterout2,NGI_road_path,outputpath,n)
		
		start_time = time.time()
		g1 = gpd.GeoDataFrame.from_file(source_out_copy)
		g2 = gpd.GeoDataFrame.from_file('temp/temp3.shp')
		res_difference = gpd.overlay(g1, g2, how='difference')
		res_difference.to_file(driver = 'ESRI Shapefile', filename= source_out_copy)
		print((time.time() - start_time),'---seconds to update source')#~ 60 seconds ---
		
		n=n+numpoints_window2
	except:
		print("reached the end of the run")





print('--- %s seconds total time---' % (time.time() - start_time))

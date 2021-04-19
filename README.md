# Virtual roads example
@ Main authors: Marc Pienaar (marc@saeon.ac.za, marc.pienaar@gmail.com)
Run the virtual_roads_test.py file to see what the model does.
Outputs are generated in an output folder

Quick description and example
-----------------
The virtual roads algorithm was developed as part of the South African Environmental Observation Network BioEnergy-Atlas (SAEON BEA) generalisable feasibility model (see https://github.com/Marc-Pienaar/BEA_Feasibility_model_test for other examples). It is written in python and uses the scikit-image open source image processing library for its minimum cost path routing (https://github.com/scikit-image/scikit-image; Van der Walt, et al., 2014). Each virtual road represents a theoretical set of paths (following a slope < 15 degrees) that could be navigated by an off-road vehicle (tractor trailer combination). These virtual roads have been created at every _n_km equally spaced interval and represent paths from origin points (areas where no roads exist) to the closest road segment according to South Africa’s National Geo-spatial Information (NGI) 2019 road layer (http://www.ngi.gov.za). The cost surface is based on 1-arcsecond (~30 meter) resolution Shuttle Radar Topography Mission (SRTM) version 3 (also known as SRTM plus) digital elevation model (DEM) (https://lpdaac.usgs.gov/news/nasa-shuttle-radar-topography-mission-srtm-version-30-srtm-plus-product-release/). All major hydrological obstacles (Rivers, dams, lakes, etc.) and areas with an average slope < 15 degrees have been removed from the cost surface and treated as obstacles. Due to the massive size of the cost surface when merged (2 711 921 602 pixels), the algorithm uses a variable window approach to avoid edge effects from individual SRTM tiles, by cropping sections from the larger merged cost surface into computationally manageable ‘chunks’. Each window or ‘cost surface crop’ is determined by first selecting an origin point and finding a small subset of its nearest neighbours using a standard Euclidean distance function. A further set of the closest destination points (road points in other words) associated with each of these origin points is also extracted using a Euclidean distance function. These origin and destination points (road points) are then combined to get an extent (with some overlap) which is used to crop the larger cost surface and process the virtual roads in a batch operation. Once a ‘chunk’ or batch of paths have been calculated, they are saved as shapefiles, and the origin points used to create these paths are removed from further analysis. The process is repeated in this manner until all origin points have been depleted. An overview and eample of the algorithm is provided below. 

![Virtual_roads_thread_marc_pienaar](https://user-images.githubusercontent.com/50328370/115218425-d5932c80-a106-11eb-8478-689ca5a7cbb9.png)
**Fgure 1**. overview of the SAEON BEA virtual roads algorithm 

The SAEON BEA virtual roads used a merged SRTM v3 DEM covering South Africa at 1-arcsecond (~30m) resolution (2 711 921 602 pixels) for its cost surface. In this example, a much smaller area is used to illustrate more detailed functionality. The region of Interest (ROI) used in this example is near the town of Upington on the Orange River, South Africa (see figure 2) 

<img width="307" alt="figure2" src="https://user-images.githubusercontent.com/50328370/115217345-b3e57580-a105-11eb-92b6-e35a820ccf3e.png">
**Figure 2**. The geographic area for the worked example on the orange river, near the town of Upington, South Africa [21.3031944449999990,-28.4962499999999999 : 21.4462499999999991,-28.3559722219999983– crs: EPSG:4326 - WGS 84 – Geographic]. The ROI is outlined in black and represents a  ~225km^2 area (~15X15 km). A Bing aerial image is shown as the backdrop. 

The methodology is summarised briefly in section 1 (and illustrated schematically in Figure 1), here it is described according to the following sections

1.	Creation of the cost surface.
2.	Creation of the source and destination points.
3.	Virtual roads variable window and minimum cost path algorithm 
4.	Creation of simplified lookup layers  

**Creation of the cost surface**

Figure 3 illustrates the creation of the cost surface. Here, a cost surface (DEM) is first used to calculate the slope using the ‘gdaldem’ processing algorithm (https://gdal.org/programs/gdaldem.html), then various obstacles are removed from the cost surface prior to calculating the minimum cost path. The process can be summarised by the following two steps:

1.	The slope is calculated from the DEM using the ‘gdaldem’ function and used as the cost surface. It is then converted into an array of values (NumPy array). Any values (the index of the values) corresponding to a slope value of less than 15 degrees are replaced by a NA or missing value.  
2.	An ‘obstacle’ layer – hence, shapefiles that have various obstacle features, in this case water features are rasterised (with the GDAL rasterise function) using the cost surface as the grid template. The new raster is also converted to an array of values (0 for a pixel / feature value, NA or missing value otherwise), and similar to step 1 above any values (the index of the values) corresponding to an obstacle value (i.e. a value of 0) in the cost surface is replaced by a NA or missing value in the cost surface array. 

![figure3](https://user-images.githubusercontent.com/50328370/115219554-f14b0280-a107-11eb-93a1-428ce3ef76fa.png)

**Figure 3**. Creation of the cost surface. A is an image from the SRTM v3 DEM for the ROI. B is the slope derived from A using the GDAL ‘gdaldem’ algorithm. C is a shapefile of water features (purple) with a Bing aerial image as a backdrop. D is the cost surface (slope) in B with the corresponding water features in C and any pixel with a slope value <15 degrees removed – the red lines represent the road layer for this area. 

**Creation of the source and destination points** 

The SANBI and STATS SA 100m Basic Spatial Unit or BSU raster [’+proj=aea +lat_1=-22 +lat_2=-38 +lat_0=-30 +lon_0=25 +x_0=1400000 +y_0=1300000 +datum=WGS84 +units=m +no_defs] is used to derive both the origin and destination points in meter units using various nested aggregations of this grid. In the case of the origin points, a 1km aggregate of the BSU grid is generated, and the centroid of each grid cell is used as an origin points. For the destination points, a vector-based road layer (here the NGI road layer) is rasterised into a higher resolution aggregate of the BSU grid, in this example 50m, and the centroids of all valid pixels extracted as a destination point layer. Both of these point layers are then reprojected to the same projection as the cost surface to use for various operations such as window selection and of course routing. Figure 4 illustrates the creation of the origin and destination points, summarised with the following steps:

1.	Rasterise the destination point layer (here, the NGI road layer) to an aggregate of the BSU grid (here, 50m) using GDAL rasterise function, with valid pixels having a value of 0, NA or missing otherwise. 
2.	Convert the rasterised layer into a point layer using the centroid values of the grid cells with valid pixel (note any resolution aggregate can be used, where a finer resolution would results in greater computational time, but more accuracy, and a coarser resolution in less computation time, but less accuracy – a 100m grid was use for the Virtual roads dataset for South Africa to reduce the computation time 
3.	Reproject the destination point layer (road points) to the same projection as the cost surface layer.
4.	Create an aggregate 1km version of the StatsSA BSU grid and convert all pixels (centroids) to a point layer in a similar manner to the creation of destination points in steps 2 and 3 above and reproject these to the same projection as the cost surface layer.

<img width="451" alt="figure4" src="https://user-images.githubusercontent.com/50328370/115217609-f6a74d80-a105-11eb-8701-926ddc800769.png">
**Figure 4**. Worked example creation of the source and destination points. A is a cropped BSU 50m aggregate grid (with binary values [0,1] to render the grid cells) over the ROI (note the difference in projection). B is the cost surface overlaid with an 1km aggregate of the BSU grid (with binary values [0,1] to render the grid cells). The centroids values of the 1km BSU grid have been converted to origin points (green), and the destination points (red dots) are derived by rasterising the road layer to the 50m grid shown in A and extracted to points in B. 

**Virtual roads variable window and minimum cost path algorithm **

A variable window approach is used to allow reasonable computational chunks, but also to minimise edge effects, that would otherwise exist when using standard SRTM tiles as an input. The window size is determined by a selection of origin points and associated ‘closest’ destination points, determined using a simple nested Euclidean distance function. The extent from the subset of origin and potential destination points is used to determine the window size (along with a small buffer). The number of origin points, associated destination points, relative to each origin point, and buffer size increase the probability of avoiding edge effects, but at the same tome increase the computation time. In other words, the large the number of origin points and associated destination points in the subset, as well as, a larger buffer, will result in a reduced chance of paths with edge effects during the routing operation. Figure 5 illustrates the variable window and minimum cost path algorithm approach for this example, summarised according to the following steps: 

1.	Divide the total number of origin points into reasonable sub-sets. Here there are a total of 220 divided into 20 equal steps (11 origin points at a time). 
2.	The algorithm iterates through each chunk  (i.e. 20 iterations) by selecting an origin point, then finds its closets neighbours (11 in this case – see Figure 5B). For each of these origin points, it finds the closets n number of destination points (250 was arbitrarily chosen). The combined origin and destination points (see Figure 5B) are used to create an extent and crop the cost surface (with a small buffer – 1km in this example). The extent is then cropped to the ROI, and used to crop the cost surface into a workable window size. 
3.	The minimum cost path is calculated in a loop for each of these origin subset points (Figure 5C), saved as a shapefile, then the origin points in the sub-set are converted to polygons (i.e given a small buffer), and a difference operation is performed to remove these values from the original origin points layer. 
4.	The next iteration is performed (as in step 3 above) on a new subset until all origin points have been depleted (Figure 5D)    

<img width="451" alt="figure5" src="https://user-images.githubusercontent.com/50328370/115217720-12125880-a106-11eb-895f-6a1b56dde1b3.png">
**Figure 5**. Variable window and minimum cost path algorithm. A is a copy of Figure 4B and shows the origin points (green dots), destination points (red dots) and a a 1km BSU grid and cost surface underlain. B is the first subset of points (both origin and associated destination [brown dots]), as well as the extend (with a 1km buffer) cropped to fit the ROI (green square). C is a close up of the first set of points show in B, with paths (green dashed lines) and buffered values (used to remove these points in the next iteration for the origin points). D show all the paths (green dashed lines) and the first 5 windows used (shaded blocks) to perform the minimum cost path routing.  

**Refs**

Farr, T.G., Rosen, P.A., Caro, E., Crippen, R., Duren, R., Hensley, S., Kobrick, M., Paller, M., Rodriguez, E., Roth, L., Seal, D., Shaffer, S., Shimada, J., Umland, J., Werner, M., Oskin, M., Burbank, D., and Alsdorf, D.E., 2007, The shuttle radar topography mission: Reviews of Geophysics, v. 45, no. 2, RG2004, at https://doi.org/10.1029/2005RG000183.

Van der Walt, S., Sch"onberger, Johannes L, Nunez-Iglesias, J., Boulogne, Franccois, Warner, J. D., Yager, N., … Yu, T. (2014). scikit-image: image processing in Python. PeerJ, 2, e453.




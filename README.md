# CTI South East Alaska project

In this first python project, the goal was to produce some code able to process the full Topographic Wetness Index (TWI), 
also called Compound Topographic Index (CTI) over the watersheds of south east Alaska.
At the moment, provided with a NED grid and a watershed shapefile, it has able to extract the corresponding DEM from the NED database
and run the full process; including merging, reprojecting in 3338 and calculating the CTI over the watershed. At the moment no tools are
easily available for computing a decent sink filling of the DEM, some groups are working on it as Pygeoprocessing and pydem.

#This code was as well used as training and is, at the moment, not totally functionnal.

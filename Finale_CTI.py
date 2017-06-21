
#This script was written in order to take advantage of  multicores linux machine as a faster alternative to the use of ARCGIS and TAUDEM
#The version of TAUDEM used in this code is : V5.3.5

def mask_raster(tif_file, shps , output_file , crop = False , invert = False) :
	'''Using rasterio to mask the rasters with a shapefile'''

	import rasterio.tools.mask
	import fiona
	print 'Masking in progress'

	if type(shps) != list : shps = [shps]

	features = [feature["geometry"] for i in shps for feature in fiona.open(i)]

	with rasterio.open(tif_file) as src:
	    out_image, out_transform = rasterio.tools.mask.mask( src, features , nodata = nodata_value , crop=crop ,invert = invert)
	    out_meta = src.meta.copy()

	out_meta.update({"driver": "GTiff",
	                 "height": out_image.shape[1],
	                 "width": out_image.shape[2],
	                 "transform": out_transform,
	                 "crs" : crs,
	                 "nodata" : nodata_value
	               })

	with rasterio.open(output_file, "w", **out_meta) as dest:
	    dest.write(out_image)

def slope_reclasser(base_DEM , slope_file , output , origin_value , target_value ):
	''' This function just reclass the slope's 0 data with 0.001 so it doesn't impact the following steps
	'''
	print 'Reclassing slope'
	with rasterio.drivers() :

		#Open both the slope and contributing area maps
		slope = rasterio.open(slope_file)
		sl_arr = slope.read(1)
		sl_arr[sl_arr == origin_value ] = target_value
		sl_arr[sl_arr<= -1] = nodata_value
		meta = rasterio.open( base_DEM ).meta
		
		meta.update( crs=crs, nodata=nodata_value )

		with rasterio.open(output, "w", **meta) as dst :
			dst.write_band(1,sl_arr.astype(rasterio.float32))


def CTI(sca,slp,cti) :
	'''Does the actual CTI calculation : log(sca/slope)'''
	print 'Computing CTI'
	with rasterio.drivers() :
		sca_arr = rasterio.open(sca).read(1)
		slp_arr = rasterio.open(slp).read(1)
		meta = rasterio.open(slp).meta
		sca_arr[sca_arr == nodata_value ] = np.nan
		slp_arr[slp_arr == nodata_value ] = np.nan

		tmp = np.log(np.divide(sca_arr,slp_arr))
		
		tmp[np.isnan(tmp)] = nodata_value
		meta.update( crs=crs, nodata=nodata_value )

		with rasterio.open(cti, "w", **meta) as dst :
			dst.write_band(1,tmp.astype(rasterio.float32))


def processing() :
	print 'Preprocessing the DEM'

	#Compute pit filling algorithm from TAUDEM
	os.system('mpirun -n %d ./pitremove -z %s -fel %s'%(cpu, base_DEM , os.path.join(out,'AKPCTR_DEMfel.tif') ))

	#################################### first part #######################################

	print 'Working on CTI with no masking'
	#first working on a version without any glacier masking

	#Run the flow direction taudem algorithm 
	os.system('mpirun -n %d ./dinfflowdir -fel %s -slp %s -ang %s' %(cpu , os.path.join(out,'AKPCTR_DEMfel.tif') , os.path.join(out,'AKPCTR_Nomask_slp.tif') , os.path.join(out,'AKPCTR_Nomask_ang.tif') ))

	#Reclass the slope raster from 0 to 0.0001
	slope_reclasser(base_DEM,os.path.join(out,'AKPCTR_Nomask_slp.tif'),os.path.join(out,'AKPCTR_Nomask_No0slp.tif'),0,0.001)

	#Run the contributing area taudem algorithm
	os.system('mpirun -n %d ./areadinf -ang %s -sca %s'%(cpu,os.path.join(out,'AKPCTR_Nomask_ang.tif') , os.path.join(out,'AKPCTR_Nomask_sca.tif') ))

	#Run the actual TWI calculation
	CTI(os.path.join(out,'AKPCTR_Nomask_sca.tif'),os.path.join(out,'AKPCTR_Nomask_No0slp.tif'),os.path.join(out,'CTI_Nomask_full.tif'))

	# #Mask the CTI raster by land mask
	mask_raster(os.path.join(out,'CTI_Nomask_full.tif'), glacier,os.path.join(out,'CTI_Nomask_full.tif'), invert=True)

	#Crop the CTI with the AOI shapefile
	mask_raster(os.path.join(out,'CTI_Nomask_full.tif'), buff_raster4k,os.path.join(out,'CTI_Nomask_AOIcropped.tif'),crop=True)


	#################################### Second part #######################################

	print 'Working on CTI without Ice'
	#Working on a CTI version with glacier masking, downslope pixels are masked by TAUDEM during the process

	# Mask out the glacier
	mask_raster(os.path.join(out,'AKPCTR_DEMfel.tif'), glacier ,os.path.join(out,'AKPCTR_NoIce.tif'),invert=True)

	# Run the flow direction taudem algorithm 
	os.system('mpirun -n %d ./dinfflowdir -fel %s -slp %s -ang %s'%(cpu,os.path.join(out,'AKPCTR_NoIce.tif') , os.path.join(out,'AKPCTR_NoIce_slp.tif') , os.path.join(out,'AKPCTR_NoIce_ang.tif') ))

	# Reclass the slope raster from 0 to 0.0001
	slope_reclasser(base_DEM,os.path.join(out,'AKPCTR_NoIce_slp.tif'),os.path.join(out,'AKPCTR_NoIce_No0slp.tif'),0,0.001)

	# Run the contributing area taudem algorithm
	os.system('mpirun -n %d ./areadinf -ang %s -sca %s' %(cpu,os.path.join(out,'AKPCTR_NoIce_ang.tif') , os.path.join(out,'AKPCTR_NoIce_sca.tif') ))

	# Run the actual TWI calculation
	CTI(os.path.join(out,'AKPCTR_NoIce_sca.tif'),os.path.join(out,'AKPCTR_NoIce_No0slp.tif'),os.path.join(out,'CTI_NoIce_full.tif'))

	# Crop the CTI with the AOI shapefile
	mask_raster(os.path.join(out,'CTI_NoIce_full.tif'), buff_raster4k,os.path.join(out,'CTI_NoIce_AOIcropped.tif'),crop=True)


if __name__ == '__main__':
	# some setup
	import os, glob, rasterio
	import numpy as np

	#set some reference files
	buff_raster4k = '/workspace/Shared/Users/jschroder/CTI/data/CTI_ref/AKPCTR_buff.shp'
	base_DEM = '/workspace/Shared/Users/jschroder/CTI/data/CTI_ref/AKPCTR_DEM_ExtraNoDataRemoved/AKPCTR_DEM_FR.tif'
	out = '/atlas_scratch/jschroder/CTI_finale2/'
	if not os.path.exists( out ):
		os.mkdir( out )

	cpu = 96
	nodata_value = -9999.0
	glacier = '/workspace/Shared/Users/jschroder/CTI/data/CTI_ref/RGIv5_AKPCTR.shp'
	with rasterio.open(base_DEM) as src :
		crs = src.meta['crs']
	
	#set path for Taudem EXE
	os.chdir('/home/UA/jschroder/src/TauDEM-Develop1/')

	processing()





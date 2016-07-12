
#This script was written in order to take advantage of 32 cores linux machine as a faster alternative to the use of ARCGIS and TAUDEM
#The version of TAUDEM used in this code is : V5.0

def mask_rasters(base_raster , mask_raster , mask_value , output , new_value=None):
	''' Function built to mask a raster with another raster. The Function takes a mask_raster or a list
	of mask_raster as input. Mask value is the value that will be masked in the base_raster and the new value
	will be the new value taken by the base_raster under the mask_value.
	'''
	def replace_masked_values( base_arr, mask_arr, mask_value , new_value):
		'''
		modify masked values in a raster using other rasters
		'''
		base_arr[ mask_arr == mask_value ] = new_value
		return base_arr
	print 'Masking raster'
	mask_arr_list = [ rasterio.open( path ).read( 1 ) for path in mask_raster ]
	base_arr = rasterio.open( base_raster ).read( 1 )
	meta = rasterio.open( base_raster ).meta

	meta.update( crs={'init':'epsg:26931'}, nodata=nodata_value )

	for mask_arr in mask_arr_list:
		base_arr = replace_masked_values( base_arr, mask_arr, mask_value, new_value )

	base_arr[base_arr < -1000] = nodata_value

	with rasterio.open( output, 'w', **meta ) as out:
		out.write( base_arr, 1 )

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
		
		meta.update( crs={'init':'epsg:26931'}, nodata=nodata_value )

		with rasterio.open(output, "w", **meta) as dst :
			dst.write_band(1,sl_arr.astype(rasterio.float32))


def smoothing_data(raster , origin_value , new_value , method='<='):
	'''Simple reclassing factory used to reclass values based on the argument method,
	if method is <= all pixels <= to origin value will be reclassed as new_value'''

	print 'Smoothing data on : %s' %raster


	with rasterio.drivers() :

		rast_arr = rasterio.open( raster ).read( 1 )
		if method == '=':
			rast_arr[rast_arr == origin_value] = new_value
		elif method  == '<=' : 
			rast_arr[rast_arr <= origin_value] = new_value
		elif method  == '>=' : 
			rast_arr[rast_arr >= origin_value] = new_value

		meta = rasterio.open(raster).meta
		meta.update( nodata = new_value )

		with rasterio.open(os.path.join(out,'Base_DEM.tif'), "w", **meta) as dst :
			dst.write_band(1,rast_arr.astype(rasterio.float32))

def CTI(sca,slp,cti) :
	print 'Computing CTI'
	with rasterio.drivers() :
		sca_arr = rasterio.open(sca).read(1)
		slp_arr = rasterio.open(slp).read(1)
		meta = rasterio.open(slp).meta

		tmp = np.log(sca_arr/slp_arr)
		tmp[tmp < -1] = nodata_value
		tmp[np.isnan(tmp)] = nodata_value
		meta.update( crs={'init':'epsg:26931'}, nodata=nodata_value )

		with rasterio.open(cti, "w", **meta) as dst :
			dst.write_band(1,tmp.astype(rasterio.float32))

def processing() :
	print 'Preprocessing the DEM'
	
	#Smoothing, making sure that we don't have two nodata values as this was an issue
	smoothing_data(base_DEM , -1000, nodata_value , '<=')

	#Compute pit filling algorithm from TAUDEM
	os.system('mpirun -n 32 ./pitremove -z %s -fel %s'%( os.path.join(out,'Base_DEM.tif') , os.path.join(out,'AKPCTR_DEMfel.tif') ))

	#################################### first part #######################################

	print 'Working on CTI without Seam or Ice'

	# Mask the glaciers as proposed by Frances
	mask_rasters(os.path.join(out,'AKPCTR_DEMfel.tif'), [ glacier , seam ],1,os.path.join(out,'AKPCTR_NoIceSeam.tif'),nodata_value)

	#Run the flow direction taudem algorithm 
	os.system('mpirun -n 32 ./dinfflowdir -fel %s -slp %s -ang %s' %(os.path.join(out,'AKPCTR_NoIceSeam.tif') , os.path.join(out,'AKPCTR_NoIceSeam_slp.tif') , os.path.join(out,'AKPCTR_NoIceSeam_ang.tif') ))

	#Reclass the slope raster from 0 to 0.0001
	slope_reclasser(base_DEM,os.path.join(out,'AKPCTR_NoIceSeam_slp.tif'),os.path.join(out,'AKPCTR_NoIceSeam_No0slp.tif'),0,0.0001)

	#Run the contributing area taudem algorithm
	os.system('mpirun -n 32 ./areadinf -ang %s -sca %s'%(os.path.join(out,'AKPCTR_NoIceSeam_ang.tif') , os.path.join(out,'AKPCTR_NoIceSeam_sca.tif') ))

	#Smoothing, making sure that we don't have two nodata values as this was an issue
	smoothing_data(os.path.join(out,'AKPCTR_NoIceSeam_sca.tif'),-1, nodata_value)

	#Run the actual TWI calculation
	CTI(os.path.join(out,'AKPCTR_NoIceSeam_sca.tif'),os.path.join(out,'AKPCTR_NoIceSeam_No0slp.tif'),os.path.join(out,'CTI_NoIceSeam_full.tif'))

	# # #Mask the CTI raster by land mask
	mask_rasters(os.path.join(out,'CTI_NoIceSeam_full.tif'), [buff_raster4k],1,os.path.join(out,'CTI_NoIceSeam_Fullext.tif'),nodata_value)


	#################################### Second part #######################################

	print 'Working on CTI without Ice'

	#Mask out the glacier
	mask_rasters(os.path.join(out,'AKPCTR_DEMfel.tif'), [ glacier ],1,os.path.join(out,'AKPCTR_NoIce.tif'),nodata_value)

	#Run the flow direction taudem algorithm 
	os.system('mpirun -n 32 ./dinfflowdir -fel %s -slp %s -ang %s'%(os.path.join(out,'AKPCTR_NoIce.tif') , os.path.join(out,'AKPCTR_NoIce_slp.tif') , os.path.join(out,'AKPCTR_NoIce_ang.tif') ))

	#Reclass the slope raster from 0 to 0.0001
	slope_reclasser(base_DEM,os.path.join(out,'AKPCTR_NoIce_slp.tif'),os.path.join(out,'AKPCTR_NoIce_No0slp.tif'),0,0.0001)

	#Run the contributing area taudem algorithm
	os.system('mpirun -n 32 ./areadinf -ang %s -sca %s' %(os.path.join(out,'AKPCTR_NoIce_ang.tif') , os.path.join(out,'AKPCTR_NoIce_sca.tif') ))

	#Smoothing, making sure that we don't have two nodata values as this was an issue
	smoothing_data(os.path.join(out,'AKPCTR_NoIce_sca.tif'),-1, nodata_value)

	#Run the actual TWI calculation
	CTI(os.path.join(out,'AKPCTR_NoIce_sca.tif'),os.path.join(out,'AKPCTR_NoIce_No0slp.tif'),os.path.join(out,'CTI_NoIce_full.tif'))

	#Mask the CTI raster by land mask
	mask_rasters(os.path.join(out,'CTI_NoIce_full.tif'), [buff_raster4k],1,os.path.join(out,'CTI_NoIce_Fullext.tif'),nodata_value)


if __name__ == '__main__':
	# some setup
	import os, glob, rasterio
	import numpy as np

	#set some reference files
	#buff_raster4k = '/workspace/Shared/Users/jschroder/TMP/AKPCTR_4km.tif'
	buff_raster4k = '/workspace/Shared/Users/jschroder/CTI/data/CTI_ref/AKPCTR_4km.tif'
	base_DEM = '/workspace/Shared/Users/jschroder/CTI/data/CTI_ref/AKPCTR_DEMcor.tif'
	out = '/workspace/Shared/Users/jschroder/TRY_AKPCTR_WS_Outline_4kfinale_newTAUDEM/'
	if not os.path.exists( out ):
		os.mkdir( out )

	nodata_value = -9999.0
	glacier = '/workspace/Shared/Users/jschroder/CTI/data/CTI_ref/RGI_v5_AKPCTR_NewExtent.tif'
	seam = '/workspace/Shared/Users/jschroder/CTI/data/CTI_ref/Seam3000_NewExtent.tif'

	#set path for Taudem EXE
	os.chdir('/home/UA/jschroder/src/TauDEM-Develop1/')

	processing()





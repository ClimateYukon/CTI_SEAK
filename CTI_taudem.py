def crop( large, small, output_name, output_path= None,compress=False, *args, **kwargs ):
	'''
	crop a larger raster by a smaller overlapping one.
	this function assumes the 2 rasters are in the same CRS
	arguments:
	----------
	large_rst = rasterio.raster object of the larger raster to be cropped
	small_rst = rasterio.raster object of the smaller raster used to crop the larger
	output_path = [optional] string path to directory to output the new GTiff created. 
		if None, it will be the directory of the large_rst.
	compress = boolean. Default = True.  If True lzw-compress. if False do not.
	returns:
	--------
	path to the newly generated cropped raster.  With the side-effect of outputting the 
	raster to disk.
	notes:
	------
	potential gotcha here is that it will only work on a single banded raster.
	'''
	large_rst = rasterio.open( large )
	small_rst = rasterio.open( small )

	window = large_rst.window( *small_rst.bounds )
	crop_affine = large_rst.window_transform( window )
	window_arr = large_rst.read( 1, window=window )
	# make a new meta object to pass to the cropped raster
	height, width = window_arr.shape
	meta = small_rst.meta
	meta.update( affine=crop_affine,
				height=height,
				width=width,
				nodata=large_rst.nodata,
				crs=large_rst.crs,
				count=1,
				dtype=large_rst.dtypes[0] )

	if output_path:
		output_filename = os.path.join( output_path, output_name )
	else:
		output_path = os.path.dirname( large_rst.name )
		output_filename = os.path.join( output_path, output_name )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( window_arr, 1 )
	return output_filename

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

	mask_arr_list = [ rasterio.open( path ).read( 1 ) for path in mask_raster ]
	base_arr = rasterio.open( base_raster ).read( 1 )
	meta = rasterio.open( base_raster ).meta

	if new_value == None :
		new_value = meta['nodata']
	else :
		meta.update(nodata=new_value)
		pass

	meta.update( crs={'init':'epsg:26931'} )

	for mask_arr in mask_arr_list:
		base_arr = replace_masked_values( base_arr, mask_arr, mask_value, new_value )

	print 'Writing masked array'
	print meta['nodata']

	with rasterio.open( output, 'w', **meta ) as out:
		out.write( base_arr, 1 )

def slope_reclasser(base_DEM,slope_file,output,origin_value,target_value):
	''' This function just reclass the slope's 0 data with 0.001 so it doesn't impact the following steps
	'''
	with rasterio.drivers() :

		#Open both the slope and contributing area maps
		slope = rasterio.open(slope_file)
		sl_arr = slope.read(1)
		sl_arr[sl_arr == origin_value ] = target_value
		sl_arr[sl_arr<= -1] = nodata_value
		meta = rasterio.open( base_DEM ).meta
				
		print 'Writing reclassed slope'
		meta.update( crs={'init':'epsg:26931'}, nodata=nodata_value )

		with rasterio.open(output, "w", **meta) as dst :
			dst.write_band(1,sl_arr.astype(rasterio.float32))

def quality_mask(base_DEM, fel , NoIce , output):

	fel_arr = rasterio.open(fel).read(1)
	sca_arr = rasterio.open(NoIce).read(1)
	fel_arr[sca_arr != -1] = nodata_value
	fel_arr[fel_arr <= 0] = nodata_value

	meta = rasterio.open( base_DEM ).meta
	meta.update( crs={'init':'epsg:26931'} )

	with rasterio.open(output, "w", **meta) as dst :
		dst.write_band(1,fel_arr.astype(rasterio.float32))


def nodata(raster):
	with rasterio.drivers() :

		rast = rasterio.open( raster ).read( 1 )
		meta = rasterio.open(raster).meta
		meta.update(nodata=nodata_value)
		with rasterio.open(i, "w", **meta) as dst :
			dst.write_band(1,rast.astype(rasterio.float32))

def CTI(sca,slp,cti) :
	with rasterio.drivers() :
		sca_arr = rasterio.open(sca).read(1)
		slp_arr = rasterio.open(slp).read(1)
		meta = rasterio.open(slp).meta

		tmp = np.log(sca_arr/slp_arr)
		tmp[tmp < -1] = nodata_value
		meta.update( crs={'init':'epsg:26931'}, nodata=nodata_value )

		with rasterio.open(cti, "w", **meta) as dst :
			dst.write_band(1,tmp.astype(rasterio.float32))

def processing() :


	#Crop with a 10k buffer so we dont use Atlas for 24h calculation flat area on the sea
	mask_rasters(base_DEM, [buff_raster10k],1,os.path.join(out,'AKPCTR_DEM4k.tif'))

	#Compute pit filling Taudem algorithm
	os.system('mpirun -n 32 pitremove -z %s -fel %s'%( os.path.join(out,'AKPCTR_DEM4k.tif') , os.path.join(out,'AKPCTR_DEMfel.tif') ))

	# Mask the glaciers as proposed by Frances
	mask_rasters(os.path.join(out,'AKPCTR_DEMfel.tif'), [ glacier , seam ],1,os.path.join(out,'AKPCTR_NoIceSeam.tif'))

	#Run the flow direction taudem algorithm 
	os.system('mpirun -n 32 dinfflowdir -fel %s -slp %s -ang %s' %(os.path.join(out,'AKPCTR_NoIceSeam.tif') , os.path.join(out,'AKPCTR_NoIceSeam_slp.tif') , os.path.join(out,'AKPCTR_NoIceSeam_ang.tif') ))

	#Reclass the slope raster from 0 to 0.0001
	#slope_reclasser(base_DEM,os.path.join(out,'AKPCTR_NoIceSeam_slp.tif'),os.path.join(out,'AKPCTR_NoIceSeam_slp.tif'),-1,nodata_value)
	slope_reclasser(base_DEM,os.path.join(out,'AKPCTR_NoIceSeam_slp.tif'),os.path.join(out,'AKPCTR_NoIceSeam_No0slp.tif'),0,0.0001)

	#Run the contributing area taudem algorithm
	os.system('mpirun -n 32 areadinf -ang %s -sca %s'%(os.path.join(out,'AKPCTR_NoIceSeam_ang.tif') , os.path.join(out,'AKPCTR_NoIceSeam_sca.tif') ))

	#Run the actual TWI calculation
	#os.system('mpirun -n 32 twi -sca %s -slp %s -twi %s'%(os.path.join(out,'AKPCTR_NoIceSeam_sca.tif'), os.path.join(out,'AKPCTR_NoIceSeam_No0slp.tif') , os.path.join(out,'CTI_NoIceSeam_1.tif') ))
	CTI(os.path.join(out,'AKPCTR_NoIceSeam_sca.tif'),os.path.join(out,'AKPCTR_NoIceSeam_No0slp.tif'),os.path.join(out,'CTI_NoIceSeam_full.tif'))

	# #Clip the CTI raster by land mask
	#crop( os.path.join(out,'CTI_NoIceSeam_full.tif'), buff_raster4k, os.path.join(out,'CTI_NoIceSeam_1.tif'))
	mask_rasters(os.path.join(out,'CTI_NoIceSeam_full.tif'), [buff_raster4k],1,os.path.join(out,'CTI_NoIceSeam.tif'),nodata_value)


	#################################### Second part #######################################

	print 'Second part'
	#Run the flow direction taudem algorithm 
	
	mask_rasters(os.path.join(out,'AKPCTR_DEMfel.tif'), [ glacier ],1,os.path.join(out,'AKPCTR_NoIce.tif'))
	os.system('mpirun -n 32 dinfflowdir -fel %s -slp %s -ang %s'%(os.path.join(out,'AKPCTR_NoIce.tif') , os.path.join(out,'AKPCTR_NoIce_slp.tif') , os.path.join(out,'AKPCTR_NoIce_ang.tif') ))

	#Reclass the slope raster from 0 to 0.0001

	slope_reclasser(base_DEM,os.path.join(out,'AKPCTR_NoIce_slp.tif'),os.path.join(out,'AKPCTR_NoIce_No0slp.tif'),0,0.0001)

	#Run the contributing area taudem algorithm
	os.system('mpirun -n 32 areadinf -ang %s -sca %s' %(os.path.join(out,'AKPCTR_NoIce_ang.tif') , os.path.join(out,'AKPCTR_NoIce_sca.tif') ))

	#Run the actual TWI calculation
	#os.system('mpirun -n 32 twi -sca %s -slp %s -twi %s'%(os.path.join(out,'AKPCTR_sca.tif') , os.path.join(out,'AKPCTR_No0slp.tif'), os.path.join(out,'CTI_1.tif') ))
	CTI(os.path.join(out,'AKPCTR_NoIce_sca.tif'),os.path.join(out,'AKPCTR_NoIce_No0slp.tif'),os.path.join(out,'CTI_NoIce_full.tif'))

	#Clip the CTI raster by land mask
	#crop( os.path.join(out,'CTI_NoIce_full.tif'), buff_raster4k, os.path.join(out,'CTI_NoIce_1.tif'))
	mask_rasters(os.path.join(out,'CTI_NoIce_full.tif'), [buff_raster4k],1,os.path.join(out,'CTI_NoIce.tif'),nodata_value)



if __name__ == '__main__':
	# some setup
	import os, glob, rasterio
	import numpy as np

	#set some reference files

	buff_raster10k = '/workspace/Shared/Users/jschroder/CTI_ref/Buffer_raster.tif'
	buff_raster4k = '/workspace/Shared/Users/jschroder/CTI_ref/Buffer_rasterFR.tif'
	#buff_raster4k = '/workspace/Shared/Users/jschroder/CTI_ref/AKPCTR_WS_Outline_4kmBuffcor.tif'
	base_DEM = '/workspace/Shared/Users/jschroder/CTI_ref/AKPCTR_DEMcor.tif'
	out = '/workspace/Shared/Users/jschroder/New_approach_finalebuilt/'
	if not os.path.exists( out ):
		os.mkdir( out )

	nodata_value = -9999
	glacier = '/workspace/Shared/Users/jschroder/CTI_ref/RGI_v5_AKPCTR_NewExtent.tif'
	seam = '/workspace/Shared/Users/jschroder/CTI_ref/Seam3000_NewExtent.tif'

	#set path for Taudem EXE
	os.chdir('/home/UA/jschroder/TauDEM/')

	processing()





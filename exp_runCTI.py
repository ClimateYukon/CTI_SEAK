def slope_reclasser(base_DEM , slope_file , output , origin_value , target_value ):
	''' This function just reclass the slope's 0 data with 0.001 so it doesn't impact the following steps
	'''

	with rasterio.drivers() :

		#Open both the slope and contributing area maps
		slope = rasterio.open(slope_file)
		sl_arr = slope.read(1)
		sl_arr[sl_arr == origin_value ] = target_value
		sl_arr[sl_arr<= -1] = nodata_value
		meta = rasterio.open( base_DEM ).meta
		
		meta.update( crs={'init':'epsg:3338'}, nodata=nodata_value )

		with rasterio.open(output, "w", **meta) as dst :
			dst.write_band(1,sl_arr.astype(rasterio.float32))

def CTI(sca,slp,cti) :

	with rasterio.drivers() :
		sca_arr = rasterio.open(sca).read(1)
		slp_arr = rasterio.open(slp).read(1)
		meta = rasterio.open(slp).meta

		tmp = np.log(sca_arr/slp_arr)
		tmp[tmp < -1] = nodata_value
		tmp[np.isnan(tmp)] = nodata_value
		meta.update( crs={'init':'epsg:3338'}, nodata=nodata_value )

		with rasterio.open(cti, "w", **meta) as dst :
			dst.write_band(1,tmp.astype(rasterio.float32))

def processing(folder) :

	HUC = folder.split('/')[-1]
	files =  glob.glob(os.path.join(folder,'*.img')) 
	tiles = ' '.join(map(str,files))
	merged = os.path.join(folder, 'raw' + HUC + '.tif')
	shp = os.path.join(folder,HUC + '_buff.shp')
	process = folder

	print 'merging'
	#Merge all the tiles together within the HUC8 folder
	merged = os.path.join(folder, 'RAW' + HUC + '.tif')
	print 'merged is %s'%merged
	os.system('gdal_merge.py -co COMPRESS=LZW -of GTiff -o  %s %s' % (merged, tiles))
	#Project the mosaic to 26931, will have to check into the nodata values, weird issues there

	base  = os.path.join(folder, 'base_DEM' + HUC + '.tif')
	print 'wraping'

	os.system('gdalwarp -overwrite -s_srs EPSG:4269 -t_srs EPSG:3338 -r bilinear -wm 999 -multi -dstnodata -9999 -q -cutline %s %s %s'% ( shp,merged , base) ) 

	#Compute pit filling algorithm from TAUDEM


def CTIprocess2(process):
	HUC = process.split('/')[-1]
	base  = os.path.join(process, 'base_DEM' + HUC + '.tif')
	os.system('mpirun -n %d ./pitremove -z %s -fel %s'%(proc , base , os.path.join(process, HUC +'AKPCTR_DEMfel.tif') ))
	
	#Run the flow direction taudem algorithm 
	os.system('mpirun -n %d ./dinfflowdir -fel %s -slp %s -ang %s'%(proc , os.path.join(process,HUC +'AKPCTR_DEMfel.tif') , os.path.join(process, HUC +'AKPCTR_NoIce_slp.tif') , os.path.join(process, HUC +'AKPCTR_NoIce_ang.tif') ))

	#Reclass the slope raster from 0 to 0.0001
	slope_reclasser(base ,os.path.join(process, HUC +'AKPCTR_NoIce_slp.tif'),os.path.join(process, HUC +'AKPCTR_NoIce_No0slp.tif'),0,0.0001)

	#Run the contributing area taudem algorithm
	os.system('mpirun -n %d ./areadinf -ang %s -sca %s' %(proc , os.path.join(process, HUC +'AKPCTR_NoIce_ang.tif') , os.path.join(process, HUC +'AKPCTR_NoIce_sca.tif') ))

	# #Smoothing, making sure that we don't have two nodata values as this was an issue
	# smoothing_data(os.path.join(process, HUC +'AKPCTR_NoIce_sca.tif'),-1, nodata_value)

	#Run the actual TWI calculation
	CTI(os.path.join(process, HUC +'AKPCTR_NoIce_sca.tif'),os.path.join(process, HUC +'AKPCTR_NoIce_No0slp.tif'),os.path.join(process, HUC + 'CTI_NoIce_full.tif'))




if __name__ == '__main__':
	# some setupexit
	import os, glob, rasterio
	import numpy as np
	from pathos import multiprocessing as mp

	base_folder = '/atlas_scratch/jschroder/Processing/'
	folders = [os.path.join(base_folder,i) for i in os.listdir(base_folder) if len(glob.glob(os.path.join(base_folder,str(i),'*CTI*'))) == 0]
	proc = 64
	nodata_value = -9999.0


	#set path for Taudem EXE
	os.chdir('/home/UA/jschroder/src/TauDEM-Develop1/')
	# pool = mp.Pool( 32 )
	# pool.map(processing,folders)
	# pool.close()
	# pool.join()

	for folder in folders:
		CTIprocess2(folder)
		print '%s is done' %(folder)


	


	# process_folder = [os.path.join(base_folder,str(i)) for i in folder_alt if len(glob.glob(os.path.join(base_folder,str(i),'*CTI*.tif'))) == 0 ]

	# #set path for Taudem EXE
	# os.chdir('/home/UA/jschroder/src/TauDEM-Develop1/')
	# # pool = mp.Pool( 32 )
	# # pool.map(processing,process_folder)
	# # pool.close()
	# # pool.join()

	# for folder in process_folder :

	#     CTIprocess(folder)
	#     print '%s is done' %(folder)




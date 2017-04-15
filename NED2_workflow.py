def dl(url):
	outpath = DL_NED2
	wget.download(url,out=outpath)
	file = os.path.join(outpath , os.path.basename(url))
	zfile = zipfile.ZipFile(file)
	zfile.extractall(outpath)
	print '%s done' %( file)

def create_index(in_path,out_path):
	if not os.path.exists(out_path): os.makedirs(out_path)
	ls = [ os.path.join(in_path,i) for i in glob.glob(os.path.join(in_path,'*.img'))    ]
	tiles = ' '.join(map(str,ls))
	outtiles = os.path.join(out_path,'tiles.shp')

	os.system('gdaltindex %s %s'%(outtiles , tiles))
	return outtiles

def dissolve(watersheds, process_folder):
	''' This function takes a watershed shapefile and group its feature by HUC_8 in order to reduce the number of DEM to process
	but still have a manageable sized DEM throughout the process, base_watershed is the watershed source and watersheds is the HUC_8 grouped shapefile'''

	from shapely.geometry import shape, mapping
	from shapely.ops import unary_union
	import fiona ,sys, os
	import itertools
	from osgeo import ogr


	'''from https://gist.github.com/jgomezdans/808194 create buffered watersheds, this part creates buffer around watersheds in order to account for inaccuracy.'''
	water_buff = add_suffix(watersheds,'buff')
	ds=ogr.Open(watersheds)
	drv=ds.GetDriver()
	if os.path.exists(water_buff):
		drv.DeleteDataSource(water_buff)
	drv.CopyDataSource(ds,water_buff)
	ds.Destroy()

	ds=ogr.Open(add_suffix(watersheds,'buff'),1)
	lyr=ds.GetLayer(0)
	for i in range(0,lyr.GetFeatureCount()):
		feat=lyr.GetFeature(i)
		lyr.DeleteFeature(i)
		geom=feat.GetGeometryRef()
		feat.SetGeometry(geom.Buffer(float(0.02)))
		lyr.CreateFeature(feat)
	ds.Destroy()

	#For later it is good to have the unit (HUC_8) in every folder for further cliping
	#that what those lines do
	with fiona.open(add_suffix(watersheds,'buff')) as source:

		meta = source.meta

		for f in source:
			path = os.path.join(process_folder,f['properties']['HUC8'])
			name_out = os.path.join(path,f['properties']['HUC8'] + '_buff.shp')
			if not os.path.exists(path): os.makedirs(path)
			with fiona.open(name_out, 'w', **meta) as sink:
				sink.write(f)

	return add_suffix(watersheds,'buff')

def move_sort( watersheds , NED_grid, process_folder) :
	'''This function has been inspired by 
	http://snorf.net/blog/2014/05/12/using-rtree-spatial-indexing-with-ogr/ 
	it allows the user to download all the NLCD tiles intersecting his area of interest. 
	The function unzip and merge the files into a projected DEM. 
	The output_path defines where the downloading and merging will take place.
	NOTE : AOI and NLCD_grid needs to share the same EPSG'''

	import rtree, wget, zipfile, glob
	import os, shutil
	from osgeo import ogr , gdalconst

	#open the file and set the variables
	ds1 = ogr.Open(NED_grid, gdalconst.GA_ReadOnly)
	ds2 = ogr.Open(watersheds, gdalconst.GA_ReadOnly)
	layer1 = ds1.GetLayer()
	layer2 = ds2.GetLayer()
	dic_index = {}
	index = rtree.index.Index(interleaved=False)
	for fid1 in range(0, layer1.GetFeatureCount()):
		feature1 = layer1.GetFeature(fid1)
		geometry1 = feature1.GetGeometryRef()
		xmin, xmax, ymin, ymax = geometry1.GetEnvelope()
		index.insert(fid1, (xmin, xmax, ymin, ymax))
	distance = 0.02
	for fid2 in range(0, layer2.GetFeatureCount()):
		NLCD_list = []
		feature2 = layer2.GetFeature(fid2)
		HUC_8 = feature2.GetField('HUC8')
		geometry2 = feature2.GetGeometryRef()
		xmin, xmax, ymin, ymax = geometry2.GetEnvelope()
		search_envelope = (xmin-distance, xmax+distance, ymin-distance, ymax+distance)
		for fid1 in list(index.intersection(search_envelope)):
			if fid1 != fid2:
				feature1 = layer1.GetFeature(fid1)
				geometry1 = feature1.GetGeometryRef()
				if geometry1.Distance(geometry2) <= distance:
					print '{} is near {}'.format(fid1, fid2)
					NLCD = feature1.GetField('location')
					NLCD_list.append(NLCD)
		
		dic_index[HUC_8]=NLCD_list

	for k,v in dic_index.iteritems() :
		for file in v:
			shutil.copy(file, os.path.join(process_folder,k))
	return process_folder

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

	meta.update( crs={'init':EPSG}, nodata=nodata_value )

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
		
		meta.update( crs={'init':EPSG}, nodata=nodata_value )

		with rasterio.open(output, "w", **meta) as dst :
			dst.write_band(1,sl_arr.astype(rasterio.float32))

def add_suffix(filename,suffix):
	tmp = filename.split('.')[0] + suffix + '.' + filename.split('.')[1]
	return tmp

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

		with rasterio.open(add_suffix(raster,'sm'), "w", **meta) as dst :
			dst.write_band(1,rast_arr.astype(rasterio.float32))

def CTI(sca,slp,cti) :
	print 'Computing CTI'
	with rasterio.drivers() :
		sca_arr = rasterio.open(sca).read(1)
		slp_arr = rasterio.open(slp).read(1)
		meta = rasterio.open(slp).meta
		sca_arr[sca_arr == nodata_value ] = np.nan
		slp_arr[slp_arr == nodata_value ] = np.nan

		tmp = np.log(np.divide(sca_arr,slp_arr))
		
		tmp[np.isnan(tmp)] = nodata_value
		meta.update( crs={'init':EPSG}, nodata=nodata_value )

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


	base  = os.path.join(folder, 'base_DEM' + HUC + '.tif')
	print 'wraping'

	os.system('gdalwarp -overwrite -s_srs EPSG:4269 -t_srs %s -r bilinear -wm 999 -multi -dstnodata -9999 -q -cutline %s %s %s'% ( EPSG,shp,merged , base) ) 


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
	smoothing_data(os.path.join(process, HUC +'AKPCTR_NoIce_sca.tif'),-1, nodata_value)

	#Run the actual TWI calculation
	CTI(os.path.join(process, HUC +'AKPCTR_NoIce_sca.tif'),os.path.join(process, HUC +'AKPCTR_NoIce_No0slp.tif'),os.path.join(process, HUC + 'CTI_NoIce_full.tif'))


if __name__ == '__main__':
	import pandas as pd
	import os , glob
	import wget
	from pathos import multiprocessing as mp
	import zipfile
	import rasterio
	import numpy as np

	base = '/atlas_scratch/jschroder/NED_2_AK_2'
	if not os.path.exists(base): os.makedirs(base)

	DL_NED2 = os.path.join(base,'img')
	if not os.path.exists(DL_NED2): os.makedirs(DL_NED2)

	process_folder = os.path.join(base,'processing')
	if not os.path.exists(process_folder): os.makedirs(process_folder)
	
	os.chdir('/home/UA/jschroder/src/TauDEM-Develop1/')
	proc = 64
	watersheds = '/workspace/Shared/Users/jschroder/CTI/NED2/shp/HUC8_4269.shp'
	EPSG = 'EPSG:3338'
	#crazy download of all NE
	# file = '/workspace/Shared/Users/jschroder/CTI/NED2/shp/NED_2_list_tiles.csv'
	# df = pd.read_csv(file,index_col=0)

	

	# urls = df.downloadURL
	# pool = mp.Pool( proc )
	# pool.map(dl,urls)
	# pool.close()
	# pool.join()


	# watershed_buffed = dissolve(watersheds , process_folder)
	# tiles = create_index(DL_NED2,os.path.join(base,'shp'))
	# watershed_buffed = add_suffix(watersheds,'buff')
	# print watershed_buffed
	# processing_path = move_sort( watersheds , tiles, process_folder)

	# os.system('rm -rf %s'%(os.path.join(process_folder,'19010500')))

	SEAK=[
	'19010301',
	'19010406',
	'19010302',
	'19010405',
	'19010404',
	'19010206',
	'19010304',
	'19010102',
	'19010103',
	'19010209',
	'19010210',
	'19010212',
	'19010204',
	'19010211',
	'19010303',
	'19010403',
	'19070103',
	'19070101',
	'19070102',
	'19070104',
	'19010208',
	'19010107',
	'19010106',
	'19010500',
	'19010104',
	'19010205',
	'19010105',
	'19010207'
	]


	nodata_value = -9999.0
	folders = [os.path.join(process_folder,i) for i in os.listdir(process_folder) if i not in SEAK]

	# pool = mp.Pool( proc )
	# pool.map(processing,folders)
	# pool.close()
	# pool.join()
	a=0

	for i in folders :
		print '%d on %d are done' %(a,len(folders))
		CTIprocess2(i)
		a =a+ 1

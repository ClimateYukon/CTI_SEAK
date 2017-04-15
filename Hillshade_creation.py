def hillshade(folder):
	in_file = glob.glob(os.path.join(folder,'*base_DEM*'))
	out_file = in_file[0].replace('DEM','Hillshade')
	os.system('gdaldem hillshade %s %s' %( in_file[0] , out_file ))

def merge(listing , title , outpath):

	tiles = ' '.join(map(str,listing))
	merged = os.path.join(outpath, title + '.img')
	# os.system('gdal_merge.py -co COMPRESS=LZW -of GTiff -o  %s %s' % (merged, tiles))
	# os.system('gdal_merge.py -of HFA -o  %s %s' % (merged, tiles))
	os.system('gdal_merge.py -of HFA -n -9999 -o  %s %s' % (merged, tiles))


if __name__ == '__main__':
	import os , glob
	from pathos import multiprocessing as mp
	base = '/big_scratch/jschroder/NED_1_AK/Processing'
	base2 = '/atlas_scratch/jschroder/NED_2_AK/processing'

	result = '/big_scratch/jschroder/NED_1_AK/result'
	if not os.path.exists(result): os.makedirs(result)
	result2 = '/atlas_scratch/jschroder/NED_2_AK/result'
	if not os.path.exists(result2): os.makedirs(result2)


	folders = [os.path.join(base,i) for i in os.listdir(base)  if len(glob.glob(os.path.join(base,i, '*DEM*'))) > 0 ]
	

	proc = 16
	pool = mp.Pool( proc )
	pool.map(hillshade,folders)
	pool.close()
	pool.join()


	# folders2 = [os.path.join(base2,i) for i in os.listdir(base2)  if len(glob.glob(os.path.join(base2,i, '*DEM*'))) > 0]
	
	# pool = mp.Pool( proc )
	# pool.map(hillshade,folders2)
	# pool.close()
	# pool.join()

	# print '###################################################################################################'
	# print 'DONE WITH PART 1'
	# print '###################################################################################################'

	list1 = [glob.glob(os.path.join(base,i, '*base_DEM*'))[0] for i in os.listdir( base ) if len(glob.glob(os.path.join(base,i, '*base_DEM*'))) > 0]
	outpath1 = os.path.join(base , 'merge_hillshade') 
	if not os.path.exists(outpath1): os.makedirs(outpath1)
	merge( list1 , 'NED1_DEM' , outpath1  )

	# list2 = [glob.glob(os.path.join(base2 , i , '*RAW*'))[0] for i in os.listdir( base2 ) if len(glob.glob(os.path.join(base2,i, '*RAW*'))) > 0]
	# outpath2 = os.path.join(base2 , 'merge_hillshade') 
	# if not os.path.exists(outpath2): os.makedirs(outpath2)
	# merge( list2 , 'NED2_hillshade' , outpath2  )

	# print '###################################################################################################'
	# print 'DONE WITH PART 2'
	# print '###################################################################################################'

	# list1 = [glob.glob(os.path.join(base , i , '*CTI*'))[0] for i in os.listdir( base ) if len(glob.glob(os.path.join(base,i, '*CTI*'))) > 0]
	# for i in list1 : os.system('cp %s %s' %(i,result))


	# list2 = [glob.glob(os.path.join(base2 , i , '*CTI*'))[0] for i in os.listdir( base2 ) if len(glob.glob(os.path.join(base2,i, '*CTI*'))) > 0]
	# for i in list2 : os.system('cp %s %s' %(i,result2))


	os.system('gdalwarp -overwrite -s_srs EPSG:4269 -t_srs EPSG:3857 -r bilinear -wm 999 -multi -dstnodata -9999 -q -cutline %s %s %s'% ( EPSG,shp,merged , base) ) 
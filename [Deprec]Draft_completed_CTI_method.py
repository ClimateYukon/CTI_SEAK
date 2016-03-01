#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
## This piece of code has been created as an example for a CTI computation and used as training for Python beginner, as 
## well as a based for possible future project.
## It is not fully functionnal as it is. The library used are at an early developement stage and 
## heavy debugging is needed for any utilisation.
## The libraries were chosen in order to be processed on a powerful machine in a Linux environment.
##
## The documentation is limited as the tool is not intended to be used outside of SNAP.
## Creator : Julien Schroder jschroder@alaska.edu
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# To do :
#           - Clean up code
#           - Fix nodata MESS!!!!
#           - Raise flags
#           - Document
#           - Make it "git-ready"
#           - Depending on what it will be used for, create comand line tool
#           - Make output and input generic
#           - Clean up dependencies
#           - Use name = main to clarify and allow flags
#           - Clean up variables
#           - Would be cool to not rely on IO but can dream about that, or ask pygeo creator

# IN ORDER TO MAKE THIS WORK, YOU NEED TO HAVE TAUDEM INSTALLED AND ATLAS SET UP WITH THOSE LINES :
#                   salloc -p devel -N 1 -w atlas02 --exclusive
#                   srun -p devel -N 1 --pty /bin/bash 

#then for heavy process, it would be smart to use fast_scratch in order to avoid I/O time getting crazy that can be done like that :
#                    cp /workspace/UA/jschroder/CTI_Frances/data/Full_dataset/Akfilled.tif /fast_scratch/jschroder/
#                    cp /workspace/UA/jschroder/CTI_Frances/data/Full_dataset/AKplanefullSEAK.tif /fast_scratch/jschroder/

#After process is done you have to run squeue to identify the JOBID and qdel it




import pygeoprocessing.routing
import pygeoprocessing.geoprocessing
import numpy as np
import rasterio, os, sys, math
import os.path
from osgeo import ogr, gdalconst, osr

base_watershed = 'watershed_proj.shp'
NLCD_grid = "/workspace/jschroder/CTI_Frances/data/Shapefiles/NLCD_grid.shp"

def dissolve(base_watershed):
    ''' This function takes a watershed shapefile and group its feature by HUC_8 in order to reduce the number of DEM to process
    but still have a manageable sized DEM throughout the process, base_watershed is the watershed source and watersheds is the HUC_8 grouped shapefile'''
    base_watershed = 'watershed_proj.shp'
    from shapely.geometry import shape, mapping
    from shapely.ops import unary_union
    import fiona ,sys, os
    import itertools
    from osgeo import ogr
    watersheds = "watersheds.shp"   
    with fiona.open(base_watershed) as input:
        # preserve the schema of the original shapefile, including the crs
        meta = input.meta
        with fiona.open(watersheds, 'w', **meta) as output:
            # groupby clusters consecutive elements of an iterable which have the same key so you must first sort the features by the 'HUC_8' field
            e = sorted(input, key=lambda k: k['properties']['HUC_8'])
            # group by the 'HUC_8' field 
            for key, group in itertools.groupby(e, key=lambda x:x['properties']['HUC_8']):
                properties, geom = zip(*[(feature['properties'],shape(feature['geometry'])) for feature in group])
                # write the feature, computing the unary_union of the elements in the group with the properties of the first element in the group
                output.write({'geometry': mapping(unary_union(geom)), 'properties': properties[0]})
    

    '''from https://gist.github.com/jgomezdans/808194 create buffered watersheds, this part creates buffer around watersheds in order to account for inaccuracy.'''

       
    water_buff = "water_buff.shp"

    ds=ogr.Open(watersheds)
    drv=ds.GetDriver()
    if os.path.exists(water_buff):
        drv.DeleteDataSource(water_buff)
    drv.CopyDataSource(ds,water_buff)
    ds.Destroy()

    ds=ogr.Open(water_buff,1)
    lyr=ds.GetLayer(0)
    for i in range(0,lyr.GetFeatureCount()):
        feat=lyr.GetFeature(i)
        lyr.DeleteFeature(i)
        geom=feat.GetGeometryRef()
        feat.SetGeometry(geom.Buffer(float(500)))
        lyr.CreateFeature(feat)
    ds.Destroy()

    #For later it is good to have the unit (HUC_8) in every folder for further cliping
    #that what those lines do
    with fiona.open(water_buff) as source:

        meta = source.meta

        for f in source:
            path = '%s%s/'%(process_folder,f['properties']['HUC_8'])
            name_out = "%s%s/AOI_buff.shp" %(process_folder,f['properties']['HUC_8'])
            if not os.path.exists(path): os.makedirs(path)
            with fiona.open(name_out, 'w', **meta) as sink:
                sink.write(f)

    with fiona.open(watersheds) as source:

    meta = source.meta

    for f in source:
        process_folder = 'DEM_process/'
        path = '%s%s/'%(process_folder,f['properties']['HUC_8'])
        name_out = "%s%s/AOI.shp" %(process_folder,f['properties']['HUC_8'])
        if not os.path.exists(path): os.makedirs(path)
        with fiona.open(name_out, 'w', **meta) as sink:
            sink.write(f)
    return watersheds
    

def download_dem( watersheds , NLCD_grid) :
    
        
    '''This function has been inspired by 
    http://snorf.net/blog/2014/05/12/using-rtree-spatial-indexing-with-ogr/ 
    it allows the user to download all the NLCD tiles intersecting his area of interest. 
    The function unzip and merge the files into a projected DEM. 
    The output_path defines where the downloading and merging will take place.
    NOTE : AOI and NLCD_grid needs to share the same EPSG'''
    

    import rtree, wget, zipfile, glob
    import os
    from osgeo import ogr , gdalconst

    #open the file and set the variables
    ds1 = ogr.Open(NLCD_grid, gdalconst.GA_ReadOnly)
    ds2 = ogr.Open(watersheds, gdalconst.GA_ReadOnly)
    layer1 = ds1.GetLayer()
    layer2 = ds2.GetLayer()

    index = rtree.index.Index(interleaved=False)
    for fid1 in range(0, layer1.GetFeatureCount()):
        feature1 = layer1.GetFeature(fid1)
        geometry1 = feature1.GetGeometryRef()
        xmin, xmax, ymin, ymax = geometry1.GetEnvelope()
        index.insert(fid1, (xmin, xmax, ymin, ymax))
    distance = 600
    for fid2 in range(0, layer2.GetFeatureCount()):
        NLCD_list = []
        feature2 = layer2.GetFeature(fid2)
        HUC_8 = feature2.GetField(fid2)
        geometry2 = feature2.GetGeometryRef()
        xmin, xmax, ymin, ymax = geometry2.GetEnvelope()
        search_envelope = (xmin-distance, xmax+distance, ymin-distance, ymax+distance)
        for fid1 in list(index.intersection(search_envelope)):
            if fid1 != fid2:
                feature1 = layer1.GetFeature(fid1)
                geometry1 = feature1.GetGeometryRef()
                if geometry1.Distance(geometry2) <= distance:
                    print '{} is near {}'.format(fid1, fid2)
                    NLCD = feature1.GetField('file_id')
                    NLCD_list.append(NLCD)


        print "The following tiles will be downloaded"
        print NLCD_list


        #process the download
        #might want to organize some folders as archive folders and processing folders as
        #for now everything is dumped in the same folder

        log_filename = ( './log.txt' )
        log = open( log_filename, 'w' )
        for i in NLCD_list :
            try:
                if not os.path.isfile("%s%s.zip" %(process_folder , i )) :
                    url = "ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/2/IMG/%s.zip" %i
                    filename = wget.download(url, out= process_folder)
                else :
                    print "The tile exists already"
                zip_ref = zipfile.ZipFile("%s/%s.zip" %(process_folder , i ),'r')
                zip_ref.extractall(process_folder + str(HUC_8))
                zip_ref.close()
            except:
            # the above line could also and most likely more effectively be use like this:
            # except Exception( 'some exception class' ) as e:
            #       log.writelines( e.printerror ) # this is not the correct syntax, but it basically error-logs for you
                log.writelines( i+'\n'+'\n' )
            pass
        log.close()

        full_HUC8_DEM = "%s/full_HUC8_DEM.tif" % (process_folder + str(HUC_8))
        DEM = "%s/DEM.tif" % (process_folder + str(HUC_8))
        tile_list = glob.glob('%s/*.img' %( process_folder + str(HUC_8)))
        tiles = ' '.join(map(str,tile_list))

        #########################################################################################################################################################################
        #########################################################################################################################################################################
        #########################################################################################################################################################################
        #########################################################################################################################################################################
        #########################################################################################################################################################################
        # Above is irrelevent to the actual CTI calculation.    




        #Merge all the tiles together within the HUC8 folder
        os.system('gdal_merge.py -co COMPRESS=NONE -co BIGTIFF=IF_NEEDED -of GTiff -o  %s %s' % (full_HUC8_DEM tiles))
        #Project the mosaic to 26931, will have to check into the nodata values, weird issues there
        os.system('gdalwarp -s_srs EPSG:4269 -t_srs EPSG:26931 -wo NUM_THREADS=ALL_CPUS -srcnodata 0 -dstnodata -9999  -r bilinear -q -cutline %s -crop_to_cutline -of GTiff %s %s' % ('AOI_buff.shp', full_HUC8_DEM , DEM))
       


def CTI_preprocessing(DEM):
    '''This function WILL fill pits (when algorithm will be available) and then calculate : slope, flow accumulation and flow direction'''
    
    #first of first : take the dem and fill the pits - not working yet
    #os.system('mpirun -n 32 pitremove -z /fast_scratch/jschroder/AKplanefullSEAK.tif -fel /fast_scratch/jschroder/Akfilled.tif')
    mpirun -n 32 dinfflowdir -fel AKPCTR_DEM_filled.tif    #pygeoprocessing.routing.fill_pits(init, dem_uri)

    #Memory block the dem
    #As 21st of May this command doesn't work because of a small issues. in order to fix that we need 
    #provide an absolute path to the aligned_dem_uri variable

    aligned_dem_uri = 'aligned_dem_pygeo.tif'

    dem_pixel_size = pygeoprocessing.geoprocessing.get_cell_size_from_uri(DEM)

    pygeoprocessing.geoprocessing.align_dataset_list(
        [DEM], [aligned_dem_uri], ['nearest'], dem_pixel_size,
        'intersection', dataset_to_align_index=0)

    #Calculate the flow direction
    flow_direction_uri = 'flow_direction_pygeo.tif' #output
    pygeoprocessing.routing.flow_direction_d_inf(aligned_dem_uri, flow_direction_uri)

    #Calculate the flow accumulation using Dinf method by Tarboton
    flow_accumulation_uri = 'flow_accumulation_pygeo.tif' #output
    pygeoprocessing.routing.flow_accumulation(flow_direction_uri, aligned_dem_uri, flow_accumulation_uri)

    #Calculate the slope using ArcGIS method
    slope = "slope_perc.tif" #output
    pygeoprocessing.geoprocessing.calculate_slope(DEM , slope)

    #Just noticed that the slope is actually calculated in percent, have to switch it back to degree
    sl = rasterio.open(slope)
    sl_arr = sl.read(1)

    deg_slope = np.arctan(sl_arr/100)
    kwargs = sl.meta  
    with rasterio.open("slope.tif", "w", **kwargs) as dst :
        dst.write_band(1,deg_slope.astype(rasterio.float32))
    deg_slope = slope
    return slope, flow_accumulation_uri



def CTI_calculation(flow_accumulation_uri, slope):
    '''This function calculate the Coumpound Topographic index given a flow accumulation input and a slope.
    it is a translation from an AML script produced by Jeff Evans http://arcscripts.esri.com/details.asp?dbid=11863 
    '''
    with rasterio.drivers() :

        #Open both the slope and flow accumulation raster
        fa = rasterio.open(flow_accumulation_uri)
        fa_arr = fa.read(1)
        
        slope = rasterio.open(slope)
        dem_slope = slope.read(1)
        
    
        #Get the informations about cellsize with geotransform
        #geotransform = dem_slope.GetGeoTransform() 
        #cs = geotransform[1]    #define cellsize
        cs = pygeoprocessing.geoprocessing.get_cell_size_from_uri(flow_accumulation_uri)
        # Needed for numpy conversions
        pi = math.pi
        deg2rad = pi / 180.0

        #Convert the dem slope data from degrees to radians
        tmp1 = dem_slope * deg2rad

        #Reclass value in dem slope in radian from 0 to 0.001
        rcl = np.equal(tmp1, 0)
        np.putmask(tmp1, rcl, .001)

        #Calculate the tangent of dem slope in radian
        tmp3 = np.tan(tmp1) 

        #Calculating the contributing area from the flow acc and the cell size
        #First convert Ca to an array
        tmp4 = (fa_arr + 1) * cs
        rcl = np.equal(tmp4, 0)
        np.putmask(tmp4, rcl, .001)

        #Final step is to take the natural logarithm of the two variable we calculated before
        cti = np.log(tmp4/tmp3)

        #Try to get rid of some NoData weirdness, will have to get back to that once method is found - linked to GDAL iteration in the download part
        cti[cti < -1] = -9999

        #Write to disk
        kwargs = fa.meta  
        with rasterio.open("CTI_pygeoprocessing.tif", "w", **kwargs) as dst :
            dst.write_band(1,cti.astype(rasterio.float32))

#not used yet as pit filling doesnt not exist for Pygeoprocessing, should move to name=main anyway
def Compute_CTI(base_watershed, NLCD_grid):
'''Could use if name = main I guess, next step!'''

    dissolve(base_watershed)
    download_dem( watersheds , NLCD_grid) 
    for HUC in HUC_8 :
        os.chdir(process_folder + HUC)
        CTI_preprocessing(DEM)
        CTI_calculation(flow_accumulation_uri,slope)

##########################################################################################################



#for TAUDEM output no need to use pixel size as we already have the contributing area.

import pygeoprocessing.routing
import pygeoprocessing.geoprocessing
import numpy as np
import rasterio, os, sys, math
import os.path

os.system('mpirun -n 32 pitremove -z /fast_scratch/jschroder/CTI2/AKPCTR_DEM.tif -fel /fast_scratch/jschroder/CTI2/AKPCTR_DEMfel.tif')
#os.system('mpirun -n 32 dinfflowdir -fel /fast_scratch/jschroder/ AKPCTR_DEM_filled.tif')
os.system('mpirun -n 32 dinfflowdir -fel /fast_scratch/jschroder/CTI_input/AKPCTR_DEMfel.tif -slp /fast_scratch/jschroder/CTI_input/AKPCTR_DEMslp.tif -ang /fast_scratch/jschroder/CTI_input/AKPCTR_DEMflw.tif')
os.system('mpirun -n 32 areadinf -ang /fast_scratch/jschroder/CTI_input/AKPCTR_DEMflw.tif -sca /fast_scratch/jschroder/CTI_input/AKPCTR_DEMca.tif')
os.system('mpirun -n 32 twi -sca /fast_scratch/jschroder/CTI_input/AKPCTR_DEMca.tif -slp /fast_scratch/jschroder/CTI_input/AKPCTR_DEMslp.tif -twi /fast_scratch/jschroder/CTI_input/AKPCTR_DEMcti.tif')

slope = '/workspace/Shared/Users/jschroder/CTI_input/AKPCTR_DEMslp.tif'
contributing_area = "/workspace/Shared/Users/jschroder/CTI_input/AKPCTR_DEMca.tif"


with rasterio.drivers() :

    #Open both the slope and contributing area maps
    ca = rasterio.open(contributing_area)
    ca_arr = ca.read(1)
    
    slope = rasterio.open(slope)
    dem_slope = slope.read(1)
    
    # Needed for numpy conversions
    pi = math.pi
    deg2rad = pi / 180.0

    #Convert the dem slope data from degrees to radians
    tmp1 = dem_slope * deg2rad

    #Reclass value in dem slope in radian from 0 to 0.001
    rcl = np.equal(tmp1, 0)
    np.putmask(tmp1, rcl, .001)

    #Calculate the tangent of dem slope in radian
    tmp3 = np.tan(tmp1) 

    #Calculating the contributing area from the flow acc and the cell size
    #First convert Ca to an array

    rcl = np.equal(ca_arr, 0)
    np.putmask(ca_arr, rcl, .001)

    #Final step is to take the natural logarithm of the two variable we calculated before
    cti = np.log(ca_arr/tmp3)

    kwargs = ca.meta  
    with rasterio.open("/workspace/Shared/Users/jschroder/CTI_input/DEM-glaciercti.tif", "w", **kwargs) as dst :
        dst.write_band(1,cti.astype(rasterio.float32))



#Run with just taudem in Atlas
import os
os.system('mpirun -n 32 pitremove -z /fast_scratch/jschroder/CTI_input/AKPCTR_DEM.tif -fel /fast_scratch/jschroder/CTI_input/AKPCTR_DEMfel.tif')
#os.system('mpirun -n 32 dinfflowdir -fel /fast_scratch/jschroder/ AKPCTR_DEM_filled.tif')
os.system('mpirun -n 32 dinfflowdir -fel /fast_scratch/jschroder/CTI_input/AKPCTR_DEMfel.tif -slp /fast_scratch/jschroder/CTI_input/AKPCTR_DEMslp.tif -ang /fast_scratch/jschroder/CTI_input/AKPCTR_DEMflw.tif')
os.system('mpirun -n 32 areadinf -ang /fast_scratch/jschroder/CTI_input/AKPCTR_DEMflw.tif -sca /fast_scratch/jschroder/CTI_input/AKPCTR_DEMca.tif')
os.system('mpirun -n 32 twi -sca /fast_scratch/jschroder/CTI_input/AKPCTR_DEMca.tif -slp /fast_scratch/jschroder/CTI_input/AKPCTR_DEMslp.tif -twi /fast_scratch/jschroder/CTI_input/AKPCTR_DEMcti.tif')

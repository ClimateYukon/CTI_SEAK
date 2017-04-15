import os,rasterio,glob
import numpy as np

full = rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/3338/CTInotmasked_glacierseam.tif' ).read(1)
masked = rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/3338/CTI_NoIce_full.tif' ).read(1)

masked[masked<=0]=np.nan
full[full<=0]=np.nan

dif = full - masked
values = dif[dif!=np.nan]
len(values)
meta =rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/3338/CTI_NoIceSeam_full.tif' ).meta
with rasterio.drivers():
    with rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/3338/calc/npdif3.tif', "w",**meta) as dst :
        dst.write_band(1,dif.astype(rasterio.float32))



import os,rasterio,glob
import numpy as np

full = rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/3338/CTI_full.tif' ).read(1)
masked = rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/3338/CTI_NoIceSeam_full.tif' ).read(1)

masked[masked<=0]=np.nan
full[full<=0]=np.nan

dif = full - masked
values = dif[dif!=np.nan]
len(values)
dif[dif==0]=np.nan
meta =rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/3338/CTI_NoIceSeam_full.tif' ).meta
with rasterio.drivers():
    with rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/3338/calc/npdiffull_glacierseam.tif', "w",**meta) as dst :
        dst.write_band(1,dif.astype(rasterio.float32))


import os,rasterio,glob
import numpy as np
full = rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/new/CTI_Fullext.tif' ).read(1)
masked = rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/new/CTI_NoIce_Fullext.tif' ).read(1)
masked[masked<=0]=np.nan
full[full<=0]=np.nan
dif = full - masked
values = dif[dif!=np.nan]
len(values)
dif[dif==0]=np.nan
meta =rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/new/CTI_NoIce_Fullext.tif' ).meta
with rasterio.drivers():
    with rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/test_difference/new/npdiffull_glacier.tif', "w",**meta) as dst :
        dst.write_band(1,dif.astype(rasterio.float32))


import os,rasterio,glob
import numpy as np
full = rasterio.open('/workspace/Shared/Users/jschroder/TMP/CTI_NoIce_AOIcropped.tif').read(1)
masked = rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/Official_finished/4th_delivery/CTI_NoIce.tif' ).read(1)
masked[masked<=0]=np.nan
full[full<=0]=np.nan
dif = full - masked
values = dif[dif!=np.nan]
len(values)
dif[dif==0]=np.nan
meta =rasterio.open('/workspace/Shared/Users/jschroder/CTI/data/Official_finished/4th_delivery/CTI_NoIce.tif' ).meta
with rasterio.drivers():
    with rasterio.open('/workspace/Shared/Users/jschroder/TMP/npdiffull_glacier.tif', "w",**meta) as dst :
        dst.write_band(1,dif.astype(rasterio.float32))

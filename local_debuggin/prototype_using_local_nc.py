# %%
import numpy as np
import xarray as xr
import netCDF4
from scipy.optimize import minimize
import rioxarray as rio

#%%
test = rio.open_rasterio("/home/dpabon/Nextcloud/work/SIF_downscaling_CDSE/data/SIF_20230701.tif")
test[]
# %%
test.set_coords()
# %%

from datetime import datetime
import rioxarray as rio
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
from rioxarray.merge import merge_datasets

files = sorted(glob.glob("data/2023-07-SIF/*.tif"))

files

datetime.strftime(files[1])

pd.to_datetime(composites_list[0][1][21:29])
# 5 days composite
composites_list = np.array_split(files, len(files) // 5)
composites_list[0][1][21:29]


dummy_dictionary = [[composites_list[i][0][21:29], composites_list[i][-1][21:29]] for i in range(len(composites_list))]

dummy_dictionary[0][1]


for i in range(len(composites_list)):

  rasters_list=[rio.open_rasterio(composites_list[i][d],mask_and_scale=True).expand_dims(dim={"time":[pd.to_datetime(composites_list[i][d][21:29])]}) for d in range(len(composites_list[i]))] 
  test = xr.concat(rasters_list, dim = 'time')
  test_mean = test.mean(dim = "time", skipna=True)

  test_mean.rio.to_raster("data/2023-07-SIF/5_days_composites/SIF_"+dummy_dictionary[i][0]+"_"+dummy_dictionary[i][1]+".tif")


plt.close()
len(composites_list)
for i in range(len(composites_list)):
  for d in range(len(composites_list[i])):
    if d < 10:
      local_d = "0"+str(composites_list[i][d])
    else:
      local_d = str(composites_list[i][d])
    



month_data = rio.open_rasterio("data/2023-07-SIF/SIF_20230701.tif")

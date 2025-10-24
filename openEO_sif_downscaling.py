#%%
import openeo
import os
import numpy as np
import xarray as xr
import zarr
import netCDF4
import pathlib
import pandas as pd
import rioxarray as rio

#%%
connection = openeo.connect('https://openeofed.dataspace.copernicus.eu/')
#%%
connection.list_collections()
#%%
connection.list_processes()

#%%
connection.authenticate_oidc()

#%%

spatial_extent_prototype = { "west": 10.36625050263503,
          "east": 11.550495280885006,
          "south": 50.1971815956687,
          "north": 50.701974766278234}

temporal_extent_prototype = ["2023-07-01", "2023-07-07"]

#%%

# Preparing Sentinel-3 data

cube_LST = connection.load_collection(
  "SENTINEL3_SLSTR_L2_LST",
  spatial_extent=spatial_extent_prototype,
  temporal_extent= temporal_extent_prototype
)

cube_OTCI = connection.load_collection(
  "SENTINEL3_OLCI_L2_LAND",
  spatial_extent=spatial_extent_prototype,
  temporal_extent= temporal_extent_prototype,
  bands = ["OTCI"]
) 

cube_IWV = connection.load_collection(
  "SENTINEL3_OLCI_L2_LAND",
  spatial_extent=spatial_extent_prototype,
  temporal_extent= temporal_extent_prototype,
  bands = ["IWV"]
)

#%%
# Collapsing time dimension

cube_LST_median = cube_LST.reduce_temporal(reducer = "median")

cube_OTCI_median = cube_OTCI.reduce_temporal(reducer = "median")

cube_IWV_median = cube_IWV.reduce_temporal(reducer = "median")

#%%
#cube_IWV_median.download('test.tif')

#%% Loading TROPOMI SIF local
cube_SIF_original_median = connection.load_url(url='https://nextcloud.bgc-jena.mpg.de/s/Bd8ngMmeM86Lr9f', format='GTiff')

#%%
cube_SIF_original = connection.load_collection(
  "SENTINEL_5P_L2",
  spatial_extent=spatial_extent_prototype,
  temporal_extent= temporal_extent_prototype,
  bands = 'CH4'
)

#%%
cube_SIF_original_median = cube_SIF_original.reduce_temporal(reducer = "median")

#%%
result = cube_SIF_original_median.save_result(format = "GTiff")

job = result.execute_batch(
    outputfile="openeo.tif",
    title="SIF",
    description="Testing SIF extraction",
    job_options={"image-name": "python311-staging"}
)
 

# %%
# resample cubes at SIF original resolution

cube_LST_median_low = cube_LST_median.resample_cube_spatial(target = cube_SIF_original_median, method='bilinear')

cube_OTCI_median_low = cube_OTCI_median.resample_cube_spatial(target = cube_SIF_original_median, method='bilinear')

cube_IWV_median_low = cube_IWV_median.resample_cube_spatial(target = cube_SIF_original_median, method='bilinear')

cube_LST_median_low

cube_OTCI_median_low

cube_IWV_median_low

#%%
# Merging all the cubes at low resolution

dataset_SIF_low = cube_SIF_original_median.merge_cubes(  cube_LST_median_low)

dataset_SIF_low = dataset_SIF_low.merge_cubes(cube_OTCI_median_low)

dataset_SIF_low = dataset_SIF_low.merge_cubes(cube_IWV_median_low)

#%%
dataset_SIF_low

#%%
# init parameters
param_ini = np.array([1.0, 2.0, -295.0, 10.0])
param_min = np.array([0.5, 0.1, -310.0, 1.0])
param_max = np.array([1.5, 5.0, -290.0, 50.0])
param_bounds = np.array([param_min, param_max])

my_udf = openeo.UDF(
"""
import xarray
import numpy as np
from scipy.optimize import minimize

def optimize_params_window(SIF_w: xarray.DataArray, VI_w: xarray.DataArray, ET_w: xarray.DataArray, LST_w: xarray.DataArray)  -> xarray.DataArray:

  param_ini = np.array([1.0, 2.0, -295.0, 10.0])
  param_min = np.array([0.5, 0.1, -310.0, 1.0])
  param_max = np.array([1.5, 5.0, -290.0, 50.0])
  param_bounds = np.array([param_min, param_max])
  
  # Nested functions, not elegant but looks like the correct way when defining openEO.UDF

  def vegetation(vi, b1, b2):
    # Ensure inputs are numpy arrays for vectorized operations
    vi_np = np.asarray(vi)
    veg = b2 * np.power(vi_np, b1) 
    return veg

  def water(et, b3, b4):

    et_np = np.asarray(et)
  
    wat = 1.0 / (1.0 + np.exp(b3 * (b4-et_np)))
    return wat

  def temperature(lst, b5, b6):
    
    # Gaussian function, ensure b6 (std dev) is positive
    lst_np = np.asarray(lst)
    temp = np.exp(-0.5 * np.power((lst_np + b5) / b6, 2))
    return temp

  def sif_model(vi, et, lst, params):
    # Calculates SIF based on the model components.
    b1, b2, b3, b4, b5, b6 = params
    sif_pred = (vegetation(vi, b1, b2) *
                water(et, b3, b4) *
                temperature(lst, b5, b6))
    return sif_pred


  def cost_function(params, vi, et, lst, sif_observed):
    
    #Least squares cost function for optimization.
    sif_pred = sif_model(vi, et, lst, params)
    
    # Calculate sum of squares, ignoring NaNs
    cost = np.nansum((sif_pred - sif_observed)**2)
    
    return cost

    
    # Wrapper function to optimize SIF parameters for a single window.
    # Designed to be used with apply_ufunc or iteration.
    # Assumes input arrays (sif_w, etc.) are 2D numpy arrays (window data).
    
    # Flatten the window arrays and filter out NaNs
    # Important: Filter consistently across all variables
    
  mask = ~np.isnan(SIF_w) & ~np.isnan(VI_w) & ~np.isnan(LST_w)
  n_valid = np.sum(mask)

  if n_valid < min_obs:
      # Not enough valid observations in the window
      return np.full(len(param_ini), np.nan, dtype=np.float32)
  

  sif_obs_f = SIF_w[mask]
  vi_f = VI_w[mask]
  et_f = ET_w[mask] 
  lst_f = LST_w[mask]

  # Define bounds for L-BFGS-B
  bounds_scipy = list(zip(*param_bounds)) # [(min1, max1), (min2, max2), ...]

  
  result = minimize(cost_function,
                            x0=param_ini,
                            args=(vi_f, lst_f, sif_obs_f),
                            method='L-BFGS-B',
                            bounds=bounds_scipy,
                            options={'maxiter': 1000, 'ftol': 1e-7, 'gtol': 1e-5}) # Adjust options as needed

  optimized_params = np.array(np.clip(result.x, param_bounds[0], param_bounds[1]))
  
  
  return optimized_params
"""
)
#%%
my_udf
#%%
dataset_SIF_low.execute_batch(
     outputfile="openeo_test.tif",
    title="SIF",
    description="Testing SIF extraction",
    job_options={"image-name": "python311-staging"}
)
# %%
parameters_cube_low = dataset_SIF_low.apply_neighborhood(my_udf,
    size=[
        {"dimension": "x", "value": 7, "unit": "px"},
        {"dimension": "y", "value": 7, "unit": "px"},
    ],
    overlap=[
        {"dimension": "x", "value": 6, "unit": "px"},
        {"dimension": "y", "value": 6, "unit": "px"},
    ],
    )
#%%
job = parameters_cube_low.execute_batch(
     outputfile="openeo_test.tif",
    title="SIF",
    description="Testing SIF extraction",
    job_options={"image-name": "python311-staging"}
)

# %%
parameters_cube_high = parameters_cube_low.resample_cube_spatial(target = cube_LST_median, method="cubic")
# %%

# resampling VI and ET at 1 km to match LST

cube_OTCI_median_1 = cube_OTCI_median.resample_cube_spatial(target = cube_LST_median)

cube_IWV_median_1 = cube_IWV_median.resample_cube_spatial(target = cube_LST_median)


# %%
cube_to_upscale = parameters_cube_high.merge_cubes(cube_LST_median)

cube_to_upscale = cube_to_upscale.merge_cubes(cube_OTCI_median_1) 

cube_to_upscale = cube_to_upscale.merge_cubes(cube_IWV_median_1) 

# %%
cube_to_upscale
# %%
# predicting SIF at 1 km resolution

udf_sif_prediction = openeo.UDF(
    
"""
def sif_downscaling_window(VI_hw, ET_hw LST_hw, PARAMS_hw):

    #Calculates downscaled SIF for a central pixel based on window means.
    #Inputs are 3D numpy arrays (window_x, window_y, [parameters]).
    #Calculate mean parameters within the window (ignore NaNs)
    # Need to handle case where all params in window are NaN
    def vegetation(vi, b1, b2):
      
      # Ensure inputs are numpy arrays for vectorized operations
      vi_np = np.asarray(vi)
      veg = b2 * np.power(vi_np, b1) 
      return veg

    def water(et, b3, b4):

      et_np = np.asarray(et)
    
      wat = 1.0 / (1.0 + np.exp(b3 * (b4-et_np)))
      return wat

    def temperature(lst, b5, b6):
      
      # Gaussian function, ensure b6 (std dev) is positive
      lst_np = np.asarray(lst)
      temp = np.exp(-0.5 * np.power((lst_np + b5) / b6, 2))
      return temp

    def sif_model(vi, et, lst, params):
      # Calculates SIF based on the model components.
      b1, b2, b3, b4, b5, b6 = params
      sif_pred = (vegetation(vi, b1, b2) *
                  water(et, b3, b4) *
                  temperature(lst, b5, b6))
      return sif_pred

    mean_params = np.nanmean(params_hw, axis=(0,1))

    if np.all(np.isnan(mean_params)):
        return np.array([np.nan]) # Cannot calculate if no valid parameters in window)

    # Calculate mean predictors within the window (ignore NaNs)
    mean_ogvi = np.nanmean(ogvi_hw)
    mean_ndwi = np.nanmean(ndwi_hw)
    mean_lst = np.nanmean(lst_hw)

    # Check if any mean predictor is NaN (means all values in window were NaN)
    if np.isnan(mean_ogvi) or np.isnan(mean_lst):
        return np.array([np.nan])

    # Apply SIF model using mean values
    sif_ds = sif_model(mean_ogvi, mean_ndwi, mean_lst, mean_params)
    return np.array([sif_ds])

"""
)
# %%
sif_downscaled = cube_to_upscale.apply_neighborhood(udf_sif_prediction, size=[
        {"dimension": "x", "value": 3, "unit": "px"},
        {"dimension": "y", "value": 3, "unit": "px"},
    ],
    overlap=[
        {"dimension": "x", "value": 2, "unit": "px"},
        {"dimension": "y", "value": 2, "unit": "px"},
    ],)

# %%

#result = sif_downscaled.save_result(format = "GTiff")

job = sif_downscaled.execute_batch(
    outputfile="openeo_test.tif",
    title="SIF",
    description="Testing SIF extraction",
    job_options={"image-name": "python311-staging"}
)


# %%

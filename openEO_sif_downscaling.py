#%%
import openeo
import os
import numpy as np
import openeo.udf
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

temporal_extent_prototype = ["2018-06-25", "2018-06-30"]

#%%

# Testing SIF - Dong Li

cube_SIF = connection.load_stac("https://raw.githubusercontent.com/dpabon/SIF_downscaling_CDSE/refs/heads/main/data/SIF_20180629.json", spatial_extent= spatial_extent_prototype, temporal_extent=temporal_extent_prototype, bands = ["SIF"])

#%%
#cube_SIF.execute_batch(outputfile="openeo_test.tif", title="SIF", description="Testing SIF extraction", job_options={"image-name": "python311-staging"})
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
cube_SIF_original = connection.load_stac("https://raw.githubusercontent.com/dpabon/SIF_downscaling_CDSE/refs/heads/main/data/SIF_20180629.json", bands = ["SIF"],
  spatial_extent=spatial_extent_prototype,
  temporal_extent= temporal_extent_prototype
)

#%%
#cube_SIF_original_median = cube_SIF_original.reduce_temporal(reducer = "median")

cube_SIF_original_median = cube_SIF_original
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
#dataset_SIF_low.execute_batch()

#%%
# init parameters
#param_ini = np.array([1.0, 2.0, -295.0, 10.0])
#param_min = np.array([0.5, 0.1, -310.0, 1.0])
#param_max = np.array([1.5, 5.0, -290.0, 50.0])
#param_bounds = np.array([param_min, param_max])

dataset_SIF_low.execute_batch(
     outputfile="openeo_test.nc",
    title="SIF",
    description="Testing SIF extraction",
    job_options={"image-name": "python311-staging"}
)
#%%

my_udf = openeo.UDF("""
import xarray
import numpy as np
from scipy.optimize import minimize

def optimize_params_window(input_cube: xarray.Dataset, context: dict)  -> xarray.DataArray:
  
    param_ini = np.array(context["param_ini"])
    param_min = np.array(context["param_min"])
    param_max = np.array(context["param_max"])
    param_bounds = np.array([param_min, param_max])
    min_obs = context["min_obs"]
  
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


    # Preparing output xarray.Datarray

    # get the coodinates of the datacube
    size_x = len(input_cube['x'].values)
    size_y = len(input_cube['y'].values)

    center_index_x = size_x // 2
    center_index_y = size_y // 2

    coord_x = input_cube.coords['x'].isel(x = center_index_x)
    coord_y = input_cube.coords['y'].isel(y = center_index_y)
    coord_t = input_cube.coords['t']

    # processing data

    SIF_w = input_cube["SIF"]
    SIF_w_valid = SIF_w.count().to_numpy()

    print('SIF_valid ' + str(SIF_w_valid))
  
    VI_w = input_cube["OTCI"]
    VI_w_valid = VI_w.count().to_numpy()
    print('OTCI_valid ' + str(VI_w_valid))
  
    ET_w = input_cube["IWV"]
    ET_w_valid = ET_w.count().to_numpy()
    print('ET_valid ' + str(ET_w_valid))
  
    LST_w =  input_cube["LST"]
    LST_w_valid = LST_w.count().to_numpy()
    print('LST valid: ' + str(LST_w_valid))

    output_bands = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6']



    if any(local < min_obs for local in [SIF_w_valid, VI_w_valid, ET_w_valid, LST_w_valid]):
        print("conditional activated")
        # Not enough valid observations in the window
        #return np.full(len(param_ini), np.nan, dtype=np.float32)
        return  xarray.DataArray(np.reshape(np.full(len(param_ini), np.nan, dtype=np.float32), shape = (1,1,1,6)),
                                 dims=['t','y', 'x', 'bands'],
                                 coords=dict(t=coord_t, 
                                             x=coord_x,
                                             y=coord_y, 
                                             bands = output_bands))
  
      

    sif_obs_f = SIF_w.to_numpy().flatten()
    vi_f = VI_w.to_numpy().flatten()
    et_f = ET_w.to_numpy().flatten()
    lst_f = LST_w.to_numpy().flatten()

    # Define bounds for L-BFGS-B
    bounds_scipy = list(zip(*param_bounds)) # [(min1, max1), (min2, max2), ...]

  
    result = minimize(cost_function,
                            x0=param_ini,
                            args=(vi_f, et_f, lst_f, sif_obs_f),
                            method='L-BFGS-B',
                            bounds=bounds_scipy,
                            options={'maxiter': 1000, 'ftol': 1e-7, 'gtol': 1e-5}) # Adjust options as needed

    optimized_params = np.array(np.clip(result.x, param_bounds[0], param_bounds[1]))

  # This needs to be transformed into a XarrayDataCube object
  # as we're optimizing a local moving window and saving the results of the optimization in the central moving window the output needs to be a x,y,band, object. Where band correspond to b1 to bx parameters.
  
  
    return xarray.DataArray(np.reshape(optimized_params, (1,1,1,6)),
                                  dims=['t','y', 'x', 'bands'],
                                 coords=dict(t=coord_t, 
                                             x=coord_x,
                                             y=coord_y, 
                                             bands = output_bands))
""",
context={"from_parameter": "context"}
)
#%%
param_ini=[1, 2, 50.0, 0, -295, 10]
param_min=[0.5, 0.1, 0.0, -1, -310, 1]
param_max=[1.5, 5, 500.0, 1, -290, 50]
#%%
my_udf
#%%
parameters_cube_low = dataset_SIF_low.apply_neighborhood(process = my_udf,
    size=[
        {"dimension": "x", "value": 7, "unit": "px"},
        {"dimension": "y", "value": 7, "unit": "px"},
    ],
    overlap=[
        {"dimension": "x", "value": 6, "unit": "px"},
        {"dimension": "y", "value": 6, "unit": "px"},
    ],
     context={"param_ini": param_ini,
              "paramin_min": param_min,
             "param_max": param_max}
    )

#%%
#parameters_cube_low.download("test.nc")
#%%
job = parameters_cube_low.execute_batch(
     outputfile="openeo_sif_parameters.nc",
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

#%%
import numpy as np
import xarray as xr
from scipy.optimize import minimize
import rioxarray as rio
import os
#%%
local_path = os.path.dirname(__file__)
local_path
#%%
local_nc = xr.open_dataset("openeo_sif_low.nc")
local_nc
#%%
local_nc['LST'].plot()

#%%
local_nc['SIF'].plot()

#%%
local_nc['OTCI'].plot()

#%%
local_nc['IWV'].plot()
#%%
SIF_w = local_nc["SIF"]
SIF_w_valid = SIF_w.count().values
SIF_w_valid
#%%
local_nc
# %%
import xarray
def apply_datacube(input_cube: xarray.Dataset, context: dict)  -> xarray.DataArray:
  
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
#%%
param_ini=[1, 2, 50.0, 0, -295, 10]
param_min=[0.5, 0.1, 0.0, -1, -310, 1]
param_max=[1.5, 5, 500.0, 1, -290, 50]

#%%

optimize_params_window(local_nc, context={"param_ini": param_ini,
              "param_min": param_min,
             "param_max": param_max,
             "min_obs": 25})


#%%
# Checking parameters low resolution rename and values

input_low = xr.open_dataset("openeo_sif_low.nc")
input_low
# %%
input_low["SIF"].plot()
# %%
input_low["LST"].plot()
# %%
input_low["OTCI"].plot()
# %%
input_low["IWV"].plot()

# %%
parameters_low = rio.open_rasterio("../data/results_parameters_optim_low_resolution.tif")
parameters_low
# %%
# LST high resolution
lst_high = rio.open_rasterio("../data/lst_high_resolution.tif")
lst_high
# %%
parameters_high = parameters_low.rio.reproject_match(lst_high)
parameters_high
# %%
parameters_high.to_raster("../data/results_paramaters_high_resolution.tif")
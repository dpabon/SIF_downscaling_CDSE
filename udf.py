import xarray
import numpy as np
from scipy.optimize import minimize
from openeo.udf.debug import inspect

def apply_datacube(cube: xarray.DataArray, context: dict)  -> xarray.DataArray:
  input_cube_local = cube.to_dataset(dim="bands")
  inspect(input_cube_local, message="local dataset: ")
  inspect(context["param_ini"], message= "context param_ini")
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
  size_x = len(input_cube_local['x'].values)
  size_y = len(input_cube_local['y'].values)

  center_index_x = size_x // 2
  center_index_y = size_y // 2

  coord_x = input_cube_local.coords['x'].isel(x = center_index_x)
  coord_y = input_cube_local.coords['y'].isel(y = center_index_y)
  coord_t = input_cube_local.coords['t']

  # processing data

  SIF_w = input_cube_local["SIF"]
  SIF_w_valid = SIF_w.count().to_numpy()

  print('SIF_valid ' + str(SIF_w_valid))

  VI_w = input_cube_local["OTCI"]
  VI_w_valid = VI_w.count().to_numpy()
  print('OTCI_valid ' + str(VI_w_valid))

  ET_w = input_cube_local["IWV"]
  ET_w_valid = ET_w.count().to_numpy()
  print('ET_valid ' + str(ET_w_valid))

  LST_w =  input_cube_local["LST"]
  LST_w_valid = LST_w.count().to_numpy()
  print('LST valid: ' + str(LST_w_valid))

  output_bands = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6']



  if any(local < min_obs for local in [SIF_w_valid, VI_w_valid, ET_w_valid, LST_w_valid]):
    print("conditional activated")
    # Not enough valid observations in the window
    #return np.full(len(param_ini), np.nan, dtype=np.float32)
    return  xarray.DataArray(np.reshape(np.full(len(param_ini), np.nan, dtype=np.float32), shape = (1,1,1,6)),dims=['t','y', 'x', 'bands'], coords=dict(t=coord_t, 
 x=coord_x, y=coord_y, bands = output_bands))

    

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
                    
  return xarray.DataArray(np.reshape(optimized_params, (1,1,1,6)),
                                dims=['t','y', 'x', 'bands'],
                                coords=dict(t=coord_t, 
                                            x=coord_x,
                                            y=coord_y, 
                                            bands = output_bands))
# %%

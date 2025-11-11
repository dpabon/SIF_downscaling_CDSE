import xarray
import numpy as np
from scipy.optimize import minimize
from openeo.udf.debug import inspect

def apply_datacube(cube: xarray.DataArray, context: dict)  -> xarray.DataArray:

  input_cube_local = cube.to_dataset(dim="bands") # 
  inspect(input_cube_local, message="local dataset: ")
  inspect(context["param_ini"], message= "context param_ini")

  param_ini = np.array(context["param_ini"])
  param_min = np.array(context["param_min"])
  param_max = np.array(context["param_max"])
  param_bounds = np.array([param_min, param_max])
  min_obs = context["min_obs"]
  window_size_lat = context["window_size_lat"]
  window_size_lon = context["window_size_lon"]

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
  
  # From xarray apply_ufunc
  def optimize_params_window(sif_w, ogvi_w, ndwi_w, lst_w,
                           param_ini, param_bounds, min_obs):
    """
    Wrapper function to optimize SIF parameters for a single window.
    Designed to be used with apply_ufunc or iteration.
    Assumes input arrays (sif_w, etc.) are 2D numpy arrays (window data).
    """
    # Flatten the window arrays and filter out NaNs
    # Important: Filter consistently across all variables
    mask = ~np.isnan(sif_w) & ~np.isnan(ogvi_w) & ~np.isnan(ndwi_w) & ~np.isnan(lst_w)
    n_valid = np.sum(mask)

    if n_valid < min_obs:
        # Not enough valid observations in the window
        return np.full(len(param_ini), np.nan, dtype=np.float32)
    

    sif_obs_f = sif_w[mask]
    vi_f = ogvi_w[mask]
    et_f = ndwi_w[mask] # Using NDWI as proxy for water stress/ET effect
    lst_f = lst_w[mask]

    # Define bounds for L-BFGS-B
    bounds_scipy = list(zip(*param_bounds)) # [(min1, max1), (min2, max2), ...]

   
    result = minimize(cost_function,
                              x0=param_ini,
                              args=(vi_f, et_f, lst_f, sif_obs_f),
                              method='L-BFGS-B',
                              bounds=bounds_scipy,
                              options={'maxiter': 1000, 'ftol': 1e-7, 'gtol': 1e-5}) # Adjust options as needed

    optimized_params = np.array(np.clip(result.x, param_bounds[0], param_bounds[1]))
    #print(optimized_params.dtype)
    #print(optimized_params.size)
    #print(optimized_params.shape)
    return optimized_params


  # Preparing output xarray.Datarray

  # get the coodinates of the datacube
  # processing data

  SIF_w = input_cube_local["SIF"]
  

  VI_w = input_cube_local["OTCI"]
  

  ET_w = input_cube_local["IWV"]
  

  LST_w =  input_cube_local["LST"]

  output_bands = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6']

  parameters_cube = xarray.apply_ufunc(
    optimize_params_window, # Function to apply
    # Input arrays:
    SIF_w.rolling(x=window_size_lat, y=window_size_lon, center=True).construct(x = 'lat_roll', y = 'lon_roll'),

    VI_w.rolling(x=window_size_lat, y=window_size_lon, center=True).construct(x = 'lat_roll', y = 'lon_roll'),

    ET_w.rolling(x=window_size_lat, y=window_size_lon, center=True).construct(x = 'lat_roll', y = 'lon_roll'),

    LST_w.rolling(x=window_size_lat, y=window_size_lon, center=True).construct(x = 'lat_roll', y = 'lon_roll'),

    # Keyword arguments for the function:
    kwargs={'param_ini': param_ini, 'param_bounds': param_bounds, 'min_obs': min_obs},
    # Define input core dimensions (the window dimensions):
    input_core_dims=[['lat_roll', 'lon_roll'], ['lat_roll', 'lon_roll'], ['lat_roll', 'lon_roll'], ['lat_roll', 'lon_roll']],
    # Define output core dimensions (the parameters dimension):
    output_core_dims=[['bands']],
    dask_gufunc_kwargs={'output_sizes': {'bands': 6}},
    # Specify the parameters dimension coordinates:
    dask='parallelized', # Enable Dask parallelization
    output_dtypes=[np.float64], # Specify output data type
    vectorize = True,
    # Add a new dimension 'parameters' to the output
    exclude_dims=set(('lat_roll', 'lon_roll')), # Exclude window dims from output
    ).rename('optimized_parameters')
  
  parameters_cube = parameters_cube.assign_coords({'bands': output_bands})

  return parameters_cube


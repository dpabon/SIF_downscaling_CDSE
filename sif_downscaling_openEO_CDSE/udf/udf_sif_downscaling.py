import xarray
import numpy as np
from scipy.optimize import minimize
from openeo.udf.debug import inspect


def apply_datacube(cube: xarray.DataArray, context: dict) -> xarray.DataArray:
    input_cube_local = cube.to_dataset(dim="bands")

    inspect(input_cube_local, message="local dataset: ")

    window_size_lat = context["window_size_lat"]
    window_size_lon = context["window_size_lon"]

    # get the coodinates of the datacube
    # processing data

    VI_w = input_cube_local["OTCI"]

    ET_w = input_cube_local["IWV"]

    LST_w = input_cube_local["LST"]

    PARAMETERS_w = input_cube_local[["b1", "b2", "b3", "b4", "b5", "b6"]]

    PARAMETERS_w = PARAMETERS_w.to_dataarray(dim="parameters")

    def sif_downscaling_window(vi, et, lst, parameters):
        # Calculates downscaled SIF for a central pixel based on window means.
        # Inputs are 3D numpy arrays (window_x, window_y, [parameters]).
        # Calculate mean parameters within the window (ignore NaNs)
        # Need to handle case where all params in window are NaN
        def vegetation(vi, b1, b2):
            # Ensure inputs are numpy arrays for vectorized operations
            vi_np = np.asarray(vi)
            veg = b2 * np.power(vi_np, b1)
            return veg

        def water(et, b3, b4):
            et_np = np.asarray(et)

            wat = 1.0 / (1.0 + np.exp(b3 * (b4 - et_np)))
            return wat

        def temperature(lst, b5, b6):
            # Gaussian function, ensure b6 (std dev) is positive
            lst_np = np.asarray(lst)
            temp = np.exp(-0.5 * np.power((lst_np + b5) / b6, 2))
            return temp

        def sif_model(vi, et, lst, params):
            # Calculates SIF based on the model components.
            b1, b2, b3, b4, b5, b6 = params
            sif_pred = (
                vegetation(vi, b1, b2) * water(et, b3, b4) * temperature(lst, b5, b6)
            )
            return sif_pred

        mean_params = np.nanmean(parameters, axis=(0, 1))

        if np.all(np.isnan(mean_params)):
            return np.array(
                [np.nan]
            )  # Cannot calculate if no valid parameters in window)

        # Calculate mean predictors within the window (ignore NaNs)
        mean_vi = np.nanmean(vi)
        mean_et = np.nanmean(et)
        mean_lst = np.nanmean(lst)

        # Check if any mean predictor is NaN (means all values in window were NaN)
        if np.isnan(mean_vi) or np.isnan(mean_lst) or np.isnan(mean_et):
            return np.array([np.nan])
        else:
            sif_ds = sif_model(mean_vi, mean_et, mean_lst, mean_params)
            return np.array([sif_ds])

    sif_cube_high = xarray.apply_ufunc(
        sif_downscaling_window,
        # Input arrays with rolling windows constructed
        VI_w.rolling(x=window_size_lat, y=window_size_lon, center=True).construct(
            x="lat_roll", y="lon_roll"
        ),
        ET_w.rolling(x=window_size_lat, y=window_size_lon, center=True).construct(
            x="lat_roll", y="lon_roll"
        ),
        LST_w.rolling(x=window_size_lat, y=window_size_lon, center=True).construct(
            x="lat_roll", y="lon_roll"
        ),
        PARAMETERS_w.rolling(
            x=window_size_lat, y=window_size_lon, center=True
        ).construct(x="lat_roll", y="lon_roll"),
        # Input core dimensions now include the window dims and the parameter dim for the last input
        input_core_dims=[
            ["lat_roll", "lon_roll"],
            ["lat_roll", "lon_roll"],
            ["lat_roll", "lon_roll"],
            ["lat_roll", "lon_roll", "parameters"],
        ],
        # Output is scalar for each pixel (no core dims)
        output_core_dims=[["SIF_downscaled"]],
        dask_gufunc_kwargs={"output_sizes": {"SIF_downscaled": 1}},
        dask="parallelized",
        output_dtypes=[np.float64],
        vectorize=True,
        exclude_dims=set(("lat_roll", "lon_roll")),  # Exclude window dims from output
    )
    inspect(sif_cube_high, message="sif_cube_pre: ")

    sif_cube_high = sif_cube_high.isel(SIF_downscaled=0)

    inspect(sif_cube_high, message="sif_cube_after:")

    output_dataset = sif_cube_high.to_dataset(name="SIF")

    inspect(output_dataset, message="output_dataset pre:")

    # adding dummy_bands
    new_vars = ["dummy_" + str(x) for x in range(8)]
    for var in new_vars:
        output_dataset[var] = xarray.full_like(sif_cube_high, np.nan)

    inspect(output_dataset, message="output_dataset after:")

    return output_dataset.to_dataarray()

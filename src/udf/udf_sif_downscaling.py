import xarray
import numpy as np
from scipy.optimize import minimize
from openeo.udf.debug import inspect


def apply_datacube(cube: xarray.DataArray, context: dict) -> xarray.DataArray:
    input_cube_local = cube.to_dataset(dim="bands")

    inspect(input_cube_local, message="local dataset: ")

    window_size_lat = context["window_size_lat"]
    window_size_lon = context["window_size_lon"]

    window_size_lat_big = (window_size_lat * 2) + 1

    optimal_n = window_size_lat**2

    def get_spiral_indices(n):
        # Create a grid of flat indices (0 to n^2 - 1)
        grid = np.arange(n * n).reshape(n, n)

        # Starting at the center
        r, c = (n - 1) // 2, (n - 1) // 2
        spiral_idx = [grid[r, c]]

        # Movement: Right, Up, Left, Down
        dr = [0, -1, 0, 1]
        dc = [1, 0, -1, 0]

        step_size = 1
        direction = 0  # Start moving Right

        while len(spiral_idx) < n * n:
            # Perform two sides of the spiral for each step_size increment
            for _ in range(2):
                for _ in range(step_size):
                    r += dr[direction]
                    c += dc[direction]
                    # Check boundaries to ensure valid index
                    if 0 <= r < n and 0 <= c < n:
                        spiral_idx.append(grid[r, c])
                direction = (direction + 1) % 4
            step_size += 1

        return np.array(spiral_idx)

    spiral_index = get_spiral_indices(window_size_lat_big)

    # get the coodinates of the datacube
    # processing data

    VI_w = input_cube_local["NIRv"]

    LST_w = input_cube_local["LST"]

    PARAMETERS_w = input_cube_local[["b1", "b2", "b5", "b6"]]

    PARAMETERS_w = PARAMETERS_w.to_dataarray(dim="parameters")

    def sif_downscaling_window(vi, lst, parameters):

        # Calculates downscaled SIF for a central pixel based on window means.
        # Inputs are 3D numpy arrays (window_x, window_y, [parameters]).
        # Calculate mean parameters within the window (ignore NaNs)
        # Need to handle case where all params in window are NaN
        def vegetation(vi, b1, b2):
            # Ensure inputs are numpy arrays for vectorized operations
            vi_np = np.asarray(vi)
            veg = b2 * np.power(vi_np, b1)
            return veg

        def temperature(lst, b5, b6):
            # Gaussian function, ensure b6 (std dev) is positive
            lst_np = np.asarray(lst)
            temp = np.exp(-0.5 * np.power((lst_np + b5) / b6, 2))
            return temp

        def sif_model(vi, lst, params):
            # Calculates SIF based on the model components.
            b1, b2, b5, b6 = params
            sif_pred = vegetation(vi, b1, b2) * temperature(lst, b5, b6)
            return sif_pred

        mean_params = np.nanmean(parameters, axis=(0, 1))

        if np.all(np.isnan(mean_params)) or np.all(np.isnan(vi)) or np.all(np.isnan(lst)):
            return np.array(
                [np.nan]
            )  # Cannot calculate if no valid parameters in window)

        

        # Calculate mean predictors within the window (ignore NaNs)
        mask = ~np.isnan(vi) & ~np.isnan(lst)
        mask = mask.astype(float) 
        mask[mask == 0] = np.nan

        vi = vi * mask

        lst = lst * mask

        vi = vi.flatten()[spiral_index]
        vi = vi[~np.isnan(vi)]

        lst = lst.flatten()[spiral_index]
        lst = lst[~np.isnan(lst)]

        if (len(vi) > optimal_n):
            vi = vi[range(optimal_n)]
            lst = lst[range(optimal_n)]

        # if central pixel is also NaN skip computation
        if np.isnan(vi[0]):
            return np.array([np.nan])

        mean_vi = np.nanmean(vi)
        mean_lst = np.nanmean(lst)

        # Check if any mean predictor is NaN (means all values in window are NaN)
        if np.isnan(np.sum(mean_vi)) or np.isnan(np.sum(mean_lst)):
            return np.array([np.nan])
        else:
            sif_ds = sif_model(mean_vi, mean_lst, mean_params)
            return np.array([sif_ds])

    sif_cube_high = xarray.apply_ufunc(
        sif_downscaling_window,
        # Input arrays with rolling windows constructed
        VI_w.rolling(x=window_size_lat_big, y=window_size_lat_big, center=True).construct(
            x="lat_roll", y="lon_roll"
        ),
        LST_w.rolling(x=window_size_lat_big, y=window_size_lat_big, center=True).construct(
            x="lat_roll", y="lon_roll"
        ),
        PARAMETERS_w.rolling(
            x=window_size_lat_big, y=window_size_lat_big, center=True
        ).construct(x="lat_roll", y="lon_roll"),
        # Input core dimensions now include the window dims and the parameter dim for the last input
        input_core_dims=[
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
        exclude_dims=set(("lat_roll", "lon_roll")),
    )
    inspect(sif_cube_high, message="sif_cube_pre: ")

    sif_cube_high = sif_cube_high.isel(SIF_downscaled=0)

    inspect(sif_cube_high, message="sif_cube_after: ")

    output_dataset = sif_cube_high.to_dataset(name="SIF")

    inspect(output_dataset, message="output_dataset pre: ")

    # adding dummy_bands
    new_vars = ["dummy_" + str(x) for x in range(8)]
    for var in new_vars:
        output_dataset[var] = xarray.full_like(sif_cube_high, np.nan)

    inspect(output_dataset, message="output_dataset after: ")

    return output_dataset.to_dataarray(dim="bands")

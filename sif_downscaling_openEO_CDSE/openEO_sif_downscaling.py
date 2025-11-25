# %%
import openeo
# %%
connection = openeo.connect("https://openeofed.dataspace.copernicus.eu/")
# %%
connection.list_collections()
# %%
connection.list_processes()

# %%
connection.authenticate_oidc()

# %%

spatial_extent_prototype = {
    "west": 0.08431253773658387,
    "east": 8.450284686186247,
    "south": 44.739501978947025,
    "north": 49.062384818441544,
}

temporal_extent_prototype = ["2023-07-01", "2023-07-07"]

temporal_extent_prototype = ["2018-06-20", "2018-06-30"]

# %%

# Testing SIF - Dong Li

cube_SIF = connection.load_stac(
    "https://raw.githubusercontent.com/dpabon/SIF_downscaling_CDSE/refs/heads/main/data/SIF_20180629.json",
    spatial_extent=spatial_extent_prototype,
    temporal_extent=temporal_extent_prototype,
    bands=["SIF"],
)

# %%
# cube_SIF.execute_batch(outputfile="openeo_test.tif", title="SIF", description="Testing SIF extraction", job_options={"image-name": "python311-staging"})
# %%

# Preparing Sentinel-3 data

cube_LST = connection.load_collection(
    "SENTINEL3_SLSTR_L2_LST",
    spatial_extent=spatial_extent_prototype,
    temporal_extent=temporal_extent_prototype,
    bands=["LST"],
)

cube_OTCI = connection.load_collection(
    "SENTINEL3_OLCI_L2_LAND",
    spatial_extent=spatial_extent_prototype,
    temporal_extent=temporal_extent_prototype,
    bands=["OTCI"],
)

cube_IWV = connection.load_collection(
    "SENTINEL3_OLCI_L2_LAND",
    spatial_extent=spatial_extent_prototype,
    temporal_extent=temporal_extent_prototype,
    bands=["IWV"],
)
# %%
cube_LST

# %%
# Collapsing time dimension

cube_LST_median = cube_LST.reduce_temporal(reducer="median")

cube_OTCI_median = cube_OTCI.reduce_temporal(reducer="median")

cube_IWV_median = cube_IWV.reduce_temporal(reducer="median")


# %%
cube_SIF_original = connection.load_stac(
    "https://raw.githubusercontent.com/dpabon/SIF_downscaling_CDSE/refs/heads/main/data/SIF_20180629.json",
    bands=["SIF"],
    spatial_extent=spatial_extent_prototype,
    temporal_extent=temporal_extent_prototype,
)

# %%
# cube_SIF_original_median = cube_SIF_original.reduce_temporal(reducer = "median")

cube_SIF_original_median = cube_SIF_original
# %%
# resample cubes at SIF original resolution

cube_LST_median_low = cube_LST_median.resample_cube_spatial(
    target=cube_SIF_original_median, method="bilinear"
)

cube_OTCI_median_low = cube_OTCI_median.resample_cube_spatial(
    target=cube_SIF_original_median, method="bilinear"
)

cube_IWV_median_low = cube_IWV_median.resample_cube_spatial(
    target=cube_SIF_original_median, method="bilinear"
)

cube_LST_median_low

cube_OTCI_median_low

cube_IWV_median_low

# %%
# Merging all the cubes at low resolution

dataset_SIF_low = cube_SIF_original_median.merge_cubes(cube_LST_median_low)

dataset_SIF_low = dataset_SIF_low.merge_cubes(cube_OTCI_median_low)

dataset_SIF_low = dataset_SIF_low.merge_cubes(cube_IWV_median_low)

# %%
# adding two dummy variables to the datacube to match the number of oputput parameters
dataset_SIF_low = dataset_SIF_low.merge_cubes(
    other=cube_SIF_original_median.rename_labels(dimension="bands", target=["dummy_1"])
)
dataset_SIF_low = dataset_SIF_low.merge_cubes(
    other=cube_SIF_original_median.rename_labels(dimension="bands", target=["dummy_2"])
)

# %%

# checking the output of the merge
"""
dataset_SIF_low.execute_batch(
     outputfile="openeo_sif_low.nc",
    title="SIF",
    description="Testing SIF extraction",
    job_options={"image-name": "python311-staging"}
)
"""
#dataset_SIF_low
# %%
# init parameters
# param_ini = np.array([1.0, 2.0, -295.0, 10.0])
# param_min = np.array([0.5, 0.1, -310.0, 1.0])
# param_max = np.array([1.5, 5.0, -290.0, 50.0])
# param_bounds = np.array([param_min, param_max])


# %%
param_ini = [1, 2, 50.0, 0, -295, 10]
param_min = [0.5, 0.1, 0.0, -1, -310, 1]
param_max = [1.5, 5, 500.0, 1, -290, 50]

# %%

my_udf = openeo.UDF.from_file(
    "udf/udf_parameters_optim_low_res.py",
    context={
        "param_ini": param_ini,
        "param_min": param_min,
        "param_max": param_max,
        "min_obs": 25,
        "window_size_lat": 5,
        "window_size_lon": 5,
    },
)


# %%
my_udf

# %%
# Checking the dimensions of the cube
dataset_SIF_low
# %%
parameters_cube_low = dataset_SIF_low.apply_neighborhood(
    process=my_udf,
    size=[
        {"dimension": "x", "value": 512, "unit": "px"},
        {"dimension": "y", "value": 512, "unit": "px"},
    ],
    overlap=[
        {"dimension": "x", "value": 0, "unit": "px"},
        {"dimension": "y", "value": 0, "unit": "px"},
    ],
)
# %% changing the name of the ouput
output_bands = ["b1", "b2", "b3", "b4", "b5", "b6"]
parameters_cube_low_rename = parameters_cube_low.rename_labels("bands", output_bands)
parameters_cube_low_rename

# %%
# Checking the ouput parameters
parameters_cube_low_rename.execute_batch(
    output_file="openeo_sif_parameters.nc",
    title="SIF_parameters_computation",
    description="Testing SIF extraction",
)

# %%
parameters_cube_high = parameters_cube_low_rename.resample_cube_spatial(
    target=cube_LST_median
)

# %%

# %%
parameters_cube_high = parameters_cube_low.resample_spatial(
    resolution=0.008928571428571, projection=None, method="bilinear"
)
# %%
job = parameters_cube_high.execute_batch(
    outputfile="openeo_sif_parameters_high.nc",
    title="SIF",
    description="Testing SIF extraction",
)

# %%

# resampling VI and ET at 1 km to match LST

cube_OTCI_median_1 = cube_OTCI_median.resample_cube_spatial(target=cube_LST_median)

cube_IWV_median_1 = cube_IWV_median.resample_cube_spatial(target=cube_LST_median)


# %%
cube_to_upscale = parameters_cube_high.merge_cubes(cube_LST_median)

cube_to_upscale = cube_to_upscale.merge_cubes(cube_OTCI_median_1)

cube_to_upscale = cube_to_upscale.merge_cubes(cube_IWV_median_1)

# %%
cube_to_upscale
# %%
# predicting SIF at 1 km resolution

udf_sif_prediction = openeo.UDF.from_file(
    "udf/udf_sif_downscaling.py",
    runtime="python311-staging",
    context={"window_size_lat": 3, "window_size_lon": 3},
)
# %%
sif_downscaled = cube_to_upscale.apply_neighborhood(
    udf_sif_prediction,
    size=[
        {"dimension": "x", "value": 512, "unit": "px"},
        {"dimension": "y", "value": 512, "unit": "px"},
    ],
    overlap=[
        {"dimension": "x", "value": 0, "unit": "px"},
        {"dimension": "y", "value": 0, "unit": "px"},
    ],
)

# %%

# result = sif_downscaled.save_result(format = "GTiff")

job = sif_downscaled.execute_batch(
    outputfile="openeo_sif_downscaled.tif",
    title="SIF_downscaling",
    description="Testing SIF extraction",
    job_options={"image-name": "python311-staging"},
)


# %%

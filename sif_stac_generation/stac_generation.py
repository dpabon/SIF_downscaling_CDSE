# %%
import os
import re
import datetime
import rasterio
from shapely.geometry import box, mapping
import pystac
from pystac.extensions.eo import EOExtension, Band
from pystac.extensions.projection import ProjectionExtension
#%%

# Configuration
DATA_DIR = "./data"
CATALOG_ID = "sif_downscaling_catalog_dong_li"
COLLECTION_ID = "SIF_sample_dong_li"
CATALOG_TITLE = "SIF Downscaling from Dong Li"
COLLECTION_TITLE = "SIF Data from Dong Li"
DESCRIPTION = "Sun Induced Fluorescence (SIF) data provided by Dong Li."

# Define the specific extensions requested
EXTENSIONS = [
    "https://stac-extensions.github.io/eo/v1.1.0/schema.json",
    "https://stac-extensions.github.io/projection/v1.1.0/schema.json"
]
#%%

def get_date_from_filename(filename):
    """
    Attempts to extract a date from the filename.
    Adjust the regex based on the actual file naming convention in the repo.
    Example assumption: 'SIF_20230101.tif'
    """
    match = re.search(r"(\d{8})", filename)
    if match:
        return datetime.datetime.strptime(match.group(1), "%Y%m%d")
    
    # Fallback: use file modification time if no date in name
    # For openEO compliance, accurate temporal extent is crucial.
    return datetime.datetime.now()
#%%

def create_catalog():
    # 1. Create the Main Catalog
    catalog = pystac.Catalog(
        id=CATALOG_ID,
        description=DESCRIPTION,
        title=CATALOG_TITLE
    )

    # 2. Create the Collection
    # openEO requires a Collection to group Items.
    # We initially set dummy extents; these will be updated after processing items.
    spatial_extent = pystac.SpatialExtent([[-180, -90, 180, 90]])
    temporal_extent = pystac.TemporalExtent([[datetime.datetime.now(), None]])
    extent = pystac.Extent(spatial_extent, temporal_extent)

    collection = pystac.Collection(
        id=COLLECTION_ID,
        title=COLLECTION_TITLE,
        description=DESCRIPTION,
        extent=extent,
        license="CC-BY-4.0",  # Adjust as necessary
        providers=[
            pystac.Provider(
                name="SIF Dong Li",
                roles=[pystac.ProviderRole.PRODUCER],
                url="https://github.com/dpabon/SIF_downscaling_CDSE"
            )
        ]
    )
    
    # Add the EO Extension link to the collection (required for summaries)
    collection.stac_extensions.append(EXTENSIONS[0])
    collection.stac_extensions.append(EXTENSIONS[1])

    # 3. Process Files and Create Items
    items = []
    bboxes = []
    dates = []

    # Define Bands (Essential for openEO 'load_collection')
    # Assuming the data is SIF. Adjust common names/descriptions as needed.
    sif_band = Band.create(
        name="SIF",
        common_name="sif",
        description="Solar Induced Fluorescence"
    )
    
    print(f"Scanning {DATA_DIR}...")
    
    if not os.path.exists(DATA_DIR):
        print(f"Error: Directory {DATA_DIR} not found. Please clone the repo first.")
        return

    for filename in os.listdir(DATA_DIR):
        if filename.endswith(".tif") or filename.endswith(".tiff"):
            filepath = os.path.join(DATA_DIR, filename)
            
            # Read metadata using Rasterio
            with rasterio.open(filepath) as src:
                bbox = src.bounds
                transform = src.transform
                shape = src.shape
                crs = src.crs
                # EPSG is required for openEO projection compliance
                epsg = crs.to_epsg() if crs else 4326 
                
                # Geometry
                geom = mapping(box(*bbox))
                
                # DateTime - create start and end datetime
                dt = get_date_from_filename(filename)
                start_dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)
                end_dt = dt.replace(hour=23, minute=0, second=0, microsecond=0)
                dates.append(dt)
                bboxes.append(list(bbox))

                # Create STAC Item
                item = pystac.Item(
                    id=os.path.splitext(filename)[0],
                    geometry=geom,
                    bbox=list(bbox),
                    datetime=None,
                    properties={
                        "start_datetime": start_dt.isoformat() + "Z",
                        "end_datetime": end_dt.isoformat() + "Z"
                    }
                )

                # --- Add Projection Extension ---
                # openEO needs proj:epsg, proj:shape, proj:bbox
                proj_ext = ProjectionExtension.ext(item, add_if_missing=True)
                proj_ext.epsg = epsg
                proj_ext.shape = shape
                proj_ext.transform = list(transform)[0:6]

                # --- Add EO Extension ---
                eo_ext = EOExtension.ext(item, add_if_missing=True)
                eo_ext.bands = [sif_band]

                # Add Asset (The actual file)
                # openEO prefers assets to have roles and types
                item.add_asset(
                    key="data",
                    asset=pystac.Asset(
                        href=f"https://github.com/dpabon/SIF_downscaling_CDSE/raw/refs/heads/main/data/{filename}",
                        media_type=pystac.MediaType.GEOTIFF,
                        roles=["data"],
                        title="SIF Data",
                        # Attach band info to the asset itself as well
                        extra_fields={
                            "eo:bands": [sif_band.to_dict()]
                        }
                    )
                )
                
                items.append(item)

    # 4. Update Collection Extents based on Items
    if items:
        # Calculate full spatial extent
        min_x = min(b[0] for b in bboxes)
        min_y = min(b[1] for b in bboxes)
        max_x = max(b[2] for b in bboxes)
        max_y = max(b[3] for b in bboxes)
        
        # Calculate full temporal extent
        min_date = min(dates)
        max_date = max(dates)

        collection.extent.spatial = pystac.SpatialExtent([[min_x, min_y, max_x, max_y]])
        collection.extent.temporal = pystac.TemporalExtent([[min_date, max_date]])

        # Add Items to Collection
        collection.add_items(items)

        # 5. Add Collection Summaries (Required for openEO discovery)
        # This tells the backend what bands are available across the collection
        summaries = {
            "eo:bands": [sif_band.to_dict()],
            "proj:epsg": [items[0].properties.get('proj:epsg')] # Assuming consistent CRS
        }
        collection.summaries = pystac.Summaries(summaries)

        # Add Collection to Catalog
        catalog.add_child(collection)

        # 6. Normalize and Save
        # This generates the catalog.json, collection.json and item.json files
        catalog.normalize_hrefs("./stac_output")
        catalog.save(catalog_type=pystac.CatalogType.SELF_CONTAINED)

        print("STAC Catalog generated in './stac_output'")
    else:
        print("No GeoTIFF files found to process.")

if __name__ == "__main__":
    create_catalog()
# %%
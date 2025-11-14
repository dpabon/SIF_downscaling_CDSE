#%%
import os
import json
import rasterio
import urllib.request
import pystac

from datetime import datetime, timezone
from shapely.geometry import Polygon, mapping
from tempfile import TemporaryDirectory

#%%

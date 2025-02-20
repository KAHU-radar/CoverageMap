import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import os.path
import sklearn.linear_model
import shapely
import pyproj

import numpy as np
import shapely.geometry
import scipy.spatial
import rasterio
from rasterio.transform import from_bounds
from rasterio.crs import CRS
import numpy as np
import skimage.measure
import scipy.ndimage 
import copy
import rasterio.features
import datetime

import scipy.spatial

xmin, ymin, xmax, ymax = -20037508, -20037508, 20037508, 20037508  # Spherical Mercator bounds
    
ship_classes = {
    "Cargo": {"Beam": [32, 40, 48, 52], "Air Draft": [40, 45, 50, 55]},
    #"Bulk Carriers": {"Beam": [30, 35, 38, 42], "Air Draft": [30, 35, 37, 40]},
    "Tanker": {"Beam": [25, 30, 34, 38], "Air Draft": [25, 28, 30, 32]},
    "Passenger Ship": {"Beam": [35, 38, 42, 45], "Air Draft": [50, 60, 70, 80]},
    "Fishing Vessel": {"Beam": [8, 10, 12, 14, 18, 22], "Air Draft": [8, 10, 12, 15, 20, 25]}
}
ship_classes["Other"] = ship_classes["Tanker"]

ship_beam_length_to_air_draft_regressors = {}
for name, data in ship_classes.items():
    regr = sklearn.linear_model.LinearRegression()
    beam = data["Beam"]
    draft = data["Air Draft"]
    if name != "Fishing Vessel":
        beam = ship_classes["Fishing Vessel"]["Beam"] +  beam
        draft = ship_classes["Fishing Vessel"]["Air Draft"] + draft
    regr.fit([[x] for x in beam], draft)
    ship_beam_length_to_air_draft_regressors[name] = regr


def read_and_filter_data(indir = "data/", outdir="filtered/", outfile="filtered.feather"):
    if os.path.exists(outfile):
        return pd.read_feather(outfile)
    if not os.path.exists(outdir): os.mkdir(outdir)
    files = [f for f in os.listdir(indir) if f.endswith("csv")]
    mmsis = set()
    for filename in sorted(files):
        print(".", end='')
        df = pd.read_csv(os.path.join(indir, filename))
        df.drop_duplicates(["mmsi"], inplace=True)
        new_mmsis = set(df.mmsi) - mmsis
        mmsis.update(new_mmsis)
        df = df.set_index("mmsi").loc[list(new_mmsis)]
        df.reset_index(names="mmsi", inplace=True)
        df.to_csv(os.path.join(outdir, filename))
    dfs = [pd.read_csv(os.path.join(outdir, name)) for name in os.listdir(outdir)] 
    df = pd.concat(dfs)
    df["timestamp_datetime"] = pd.to_datetime(df.timestamp, format='mixed')
    df["ship_type"] = df.ship_type.fillna("Other")
    df.to_feather(outfile)
    return df

def spherical_mercator_project_and_grid(df, resolution=1852 * 60, cache="projected.feather"):
    if os.path.exists(cache):
        return pd.read_feather(cache)
    transformer = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    df['x'], df['y'] = transformer.transform(df['longitude'], df['latitude'])
    df["gridx"] = df.x // resolution
    df["gridy"] = df.y // resolution
    df.to_feather(cache)
    return df

def estimate_air_drafts_and_radar_horizon_from_beam_lengths(df, cache="air_drafts.feather"):
    if os.path.exists(cache):
        return pd.read_feather(cache)
    for ship_type in df.ship_type.unique():
        regr = ship_beam_length_to_air_draft_regressors.get(ship_type, ship_beam_length_to_air_draft_regressors["Other"])
        filt = df.ship_type == ship_type
        lengths = df.loc[filt, ["length"]].fillna(10).values
        lengths = np.where(lengths < 1, 10, lengths)
        df.loc[filt, "air_draft"] = regr.predict(lengths)
    df["air_draft"] = np.where(df["air_draft"] > 100, 100, df["air_draft"])
    df["horizon"] = 2.1 * np.sqrt(df.air_draft) * 1852
    df.to_feather(cache)
    return df


def render_coverage_map(df, resolution, cache="coverage.tiff"):
    if os.path.exists(cache):
        with rasterio.open(cache) as src:
            return src.read(1)

    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x, df.y, crs=3857))
    gdf["geometry"] = gdf.geometry.buffer(gdf.horizon)

    x = np.arange(xmin, xmax, resolution)
    y = np.arange(ymin, ymax, resolution)
    grid = np.zeros((len(y), len(x)), dtype=float)

    # Assume pixels are larger than circles.
    for xidx in range(len(x)):
        start = datetime.datetime.now()
        for yidx in range(len(y)):
            x_val = x[xidx]
            y_val = y[yidx]

            # Non overlapping if A.right < B.left
            xfilt = ~((gdf.x + gdf.horizon < x_val) | (x_val + resolution < gdf.x - df.horizon))
            yfilt = ~((gdf.y + gdf.horizon < y_val) | (y_val + resolution < gdf.y - df.horizon))
            filt = xfilt & yfilt

            if filt.sum() == 0: continue

            pixel = shapely.geometry.box(x_val, y_val, x_val + resolution, y_val + resolution)

            df_overlapping = gdf.loc[filt]

            coverage = df_overlapping.geometry.unary_union.intersection(pixel)

            grid[yidx, xidx] = coverage.area / pixel.area
        end = datetime.datetime.now()
        print("%.2f" % (end - start).total_seconds(), end="; ")
        

    transform = from_bounds(xmin, ymin, xmax, ymax, grid.shape[0], grid.shape[1])
    crs = CRS.from_epsg(3857)

    with rasterio.open(
        "coverage.tiff",
        'w',
        driver='GTiff',
        height=grid.shape[0],
        width=grid.shape[1],
        count=1,
        dtype=grid.dtype,
        crs=crs,
        transform=transform
    ) as dst:
        dst.write(grid, 1)
    return grid

def generate_land_map(height, width):
    coast = gpd.read_file("../coastline50.geojson")
    coast['geometry'] = coast['geometry'].apply(
        lambda geom: geom if geom.is_closed else shapely.geometry.LineString(list(geom.coords) + [geom.coords[0]]))
    crs = CRS.from_epsg(3857)
    coast = coast.to_crs(crs)

    transform = from_bounds(xmin, ymax, xmax, ymin, height, width)

    geoms = shapely.ops.polygonize(coast.geometry.union_all())
    shapes = [(geom, 1) for geom in geoms]
    return rasterio.features.rasterize(shapes, out_shape=(height, width), transform=transform, fill=0, dtype="uint8")    

def scale_to(img, target_size):
    if target_size[0] == img.shape[0] and target_size[1] == img.shape[1]: return img
    scale_y = target_size[0] / img.shape[0]
    scale_x = target_size[1] / img.shape[1]
    
    row_indices = (np.arange(target_size[0]) / scale_y).astype(int)
    col_indices = (np.arange(target_size[1]) / scale_x).astype(int)
    
    return  img[row_indices[:, None], col_indices]

def pyramidize(orig, minsize=10):
    block_size = (2, 2)
    pyramid = [copy.deepcopy(orig)]
    while pyramid[-1].shape[0] > minsize:
        pyramid.append(skimage.measure.block_reduce(pyramid[-1], block_size, np.mean))
    return [scale_to(img, orig.shape) for img in pyramid]

def compose_pyramids(pyramid, landpyramid, cutoff=0.4, cache="composite-coverage.tiff"):
    if os.path.exists(cache):
        with rasterio.open(cache) as src:
            return src.read(1)

    composite = pyramid[-1]
    for idx in range(len(pyramid)-2, -1, -1):
        composite = np.where((composite < cutoff) & (landpyramid[idx+1] == 0), composite, pyramid[idx])

    transform = from_bounds(xmin, ymin, xmax, ymax, composite.shape[0], composite.shape[1])
    crs = CRS.from_epsg(3857)

    with rasterio.open(
        cache,
        'w',
        driver='GTiff',
        height=composite.shape[0],
        width=composite.shape[1],
        count=1,
        dtype=composite.dtype,
        crs=crs,
        transform=transform
    ) as dst:
        dst.write(composite, 1)
        
    return composite


def calculate_vessel_visibility(df):
    tree = scipy.spatial.cKDTree(df[["x", "y"]].values)
    neighbours = tree.query(df[["x", "y"]].values, 100)

    d, i = neighbours
    # Remove first column, as that is the points themselves!
    d = d[:,1:]
    i = i[:,1:]

    isvisible = d < np.column_stack([df.horizon.values] * d.shape[1]) + df.horizon.values[i]

    #isvisible.sum() / (isvisible.shape[0] * isvisible.shape[1])

    return {"distance_to_nearest_vessel": d, "isvisible": isvisible}

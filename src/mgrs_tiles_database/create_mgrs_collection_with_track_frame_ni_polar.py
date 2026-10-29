# Made Mar 4 2026 . for DSWx-NI calval
import argparse
import os

import sqlite3
import numpy as np
import time
from collections import Counter
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
from shapely.wkt import loads as wkt_loads
from shapely.wkb import loads as wkb_loads
from shapely import wkb
from shapely.ops import unary_union
from pyproj import Transformer
from osgeo import ogr, osr, gdal
from rasterio.mask import mask as rasterio_mask
import rasterio
from shapely import wkt
from shapely.ops import transform as shapely_transform
from pyproj import Transformer


ogr.DontUseExceptions()
north_bound_polygon = Polygon(((-180, 84.5, 0),
                               (-180, 89.9, 0),
                               (180, 89.9, 0),
                               (180, 84.5, 0)))

south_bound_polygon = Polygon(((-180, -84, 0),
                               (-180, -89.9, 0),
                               (180, -89.9, 0),
                               (180, -84, 0)))

equator_bound_polygon = Polygon(((-180, -0.5, 0),
                                (-180, 0.5, 0),
                                (180, 0.5, 0),
                                (180, -0.5, 0)))
equator_bound_polygon = Polygon(((-180, -2.5, 0),
                                (-180, 2.5, 0),
                                (180, 2.5, 0),
                                (180, -2.5, 0)))
north_pol_polygon = Polygon(((-180, 70.0, 0),
                            (-180, 85.5, 0),
                            (180, 85.5, 0),
                            (180, 70.0, 0)))
south_pol_polygon = Polygon(((-180, -70.0, 0),
                            (-180, -85.5, 0),
                            (180, -85.5, 0),
                            (180, -70.0, 0)))

north_asc_des_changing_polygon = Polygon(
    ((-180, 76.0, 0),
    (-180, 85.5, 0),
    (180, 85.5, 0),
    (180, 76.0, 0)))
south_asc_des_changing_polygon = Polygon(((-180, -85.0, 0),
                            (-180, -85.5, 0),
                            (180, -85.5, 0),
                            (180, -85.0, 0)))
np2gdal_conversion = {
  "byte": 1,
  "uint8": 1,
  "int8": 1,
  "uint16": 2,
  "int16": 3,
  "uint32": 4,
  "int32": 5,
  "float32": 6,
  "float64": 7,
  "complex64": 10,
  "complex128": 11,
}

def save_raster_gdal(data, output_file, geotransform,
                     projection, scratch_dir='.',
                     datatype='float32'):
    """Save images using Gdal

    Parameters
    ----------
    data: numpy.ndarray
        Data to save into the file
    output_file: str
        full path for filename to save the DSWx-S1 file
    geotransform: gdal
        gdaltransform information
    projection: gdal
        projection object
    scratch_dir: str
        temporary file path to process COG file.
    datatype: str
        Data types to save the file.
    """
    gdal_type = np2gdal_conversion[str(datatype)]
    image_size = data.shape
    #  Set the Pixel Data (Create some boxes)
    # set geotransform
    if data.ndim == 3:
        ndim = image_size[0]
        ny = image_size[1]
        nx = image_size[2]
    elif data.ndim == 2:
        ny = image_size[0]
        nx = image_size[1]
        ndim = 1

    driver = gdal.GetDriverByName("GTiff")
    output_file_path = os.path.join(output_file)
    gdal_ds = driver.Create(output_file_path,
                            nx, ny,
                            ndim, gdal_type)
    gdal_ds.SetGeoTransform(geotransform)
    gdal_ds.SetProjection(projection)

    if data.ndim == 2:
        gdal_ds.GetRasterBand(1).WriteArray(data)
    else:
        for im_ind in range(0, ndim):
            gdal_ds.GetRasterBand(im_ind+1).WriteArray(
                np.squeeze(data[im_ind, :, :]))

    gdal_ds.FlushCache()
    gdal_ds = None
    del gdal_ds  # close the dataset (Python object and pointers)

# determine the path to the world land GPKG file
# LAND_GPKG_FILE = os.path.join(WORKFLOW_SCRIPTS_DIR, 'data',
#                               'GSHHS_l_L1.shp.no_attrs.epsg3413_dissolved.gpkg')
import mgrs_tiles_database

# LAND_GPKG_FILE = "/mnt/aurora-r0/jungkyo/OPERA/burst_db/MGRS_db/GSHHS_l_L1.shp.no_attrs.epsg3413_dissolved.gpkg"
WORKFLOW_SCRIPTS_DIR = os.path.dirname(mgrs_tiles_database.__file__)

LAND_GPKG_FILE = os.path.join(WORKFLOW_SCRIPTS_DIR, 'data',
                              'GSHHS_l_L1.shp.no_attrs.epsg3413_dissolved.gpkg')
LAND_GPKG_FILE = os.path.join('/mnt/aurora-r0/jungkyo/OPERA/burst_db/MGRS_db'
    ,'ocean_land_low_resolution.gpkg')



def get_utm_central_meridian(epsg: int) -> float:
    """
    Return the central meridian of a WGS84 UTM CRS.

    Supports:
        EPSG:32601-32660
        EPSG:32701-32760
    """
    epsg = int(epsg)

    if 32601 <= epsg <= 32660:
        zone = epsg - 32600
    elif 32701 <= epsg <= 32760:
        zone = epsg - 32700
    else:
        raise ValueError(f"EPSG:{epsg} is not a WGS84 UTM CRS.")

    return zone * 6.0 - 183.0


def projected_geometry_bounds_antimeridian_safe(
    geometry,
    dst_epsg: int,
):
    """
    Calculate projected bounds of a WGS84 geometry safely near the
    antimeridian.

    Parameters
    ----------
    geometry : shapely geometry
        Polygon or MultiPolygon in EPSG:4326.
    dst_epsg : int
        Target projected CRS, normally a WGS84 UTM EPSG code.

    Returns
    -------
    xmin, ymin, xmax, ymax : float
        Bounds of the complete transformed geometry.
    """
    if geometry is None or geometry.is_empty:
        raise ValueError("Cannot calculate bounds of an empty geometry.")

    ref_lon = get_utm_central_meridian(dst_epsg)

    def unwrap_coordinates(x, y, z=None):
        # Keep longitude within ±180 degrees of the selected
        # UTM zone's central meridian.
        x_unwrapped = ref_lon + ((np.asarray(x) - ref_lon + 180.0) % 360.0 - 180.0)

        if z is None:
            return x_unwrapped, y

        return x_unwrapped, y, z

    unwrapped_geometry = shapely_transform(
        unwrap_coordinates,
        geometry,
    )

    # force_over=True allows values such as 181 degrees instead of
    # normalizing them back to -179 degrees.
    transformer = Transformer.from_crs(
        "EPSG:4326",
        f"EPSG:{dst_epsg}",
        always_xy=True,
        force_over=True,
    )

    projected_geometry = shapely_transform(
        transformer.transform,
        unwrapped_geometry,
    )

    return projected_geometry.bounds

# Convert 3D geometry to 2D
def _get_land_mask(epsg_cslc: int, geotransform: tuple, shape_mask: tuple):
    '''
    Get the land mask within the CSLC bounding box

    Parameters
    ----------
    epsg_cslc: int
        EPSG code of the CSLC layer
    geotransform: tuple
        Geotransform vector of the CSLC layer
    shape_mask: tuple
        Shape of the raster as numpy array

    Returns
    -------
    mask_land: np.ndarray
        Raster Mask for land area. `1` is land, `0` otherwise
    '''
    # Extract the land polygon
    ds_land = ogr.Open(LAND_GPKG_FILE, 0)
    layer_land = ds_land.GetLayer()
    feature = layer_land.GetNextFeature()
    land_polygon = feature.GetGeometryRef()

    # extract the EPSG of the land polgyon GPKG
    srs_gpkg = layer_land.GetSpatialRef()
    land_epsg = int(srs_gpkg.GetAuthorityCode(None))

    # Compute and create the bounding box
    xmin = geotransform[0]
    ymin = geotransform[3] + geotransform[5] * shape_mask[0]
    xmax = geotransform[0] + geotransform[1] * shape_mask[1]
    ymax = geotransform[3]
    bbox_cslc = ogr.Geometry(ogr.wkbPolygon)
    print('bbox', bbox_cslc)
    ring_cslc = ogr.Geometry(ogr.wkbLinearRing)
    print('ring_cslc', ring_cslc)
    ring_cslc.AddPoint(xmin, ymin)
    ring_cslc.AddPoint(xmax, ymin)
    ring_cslc.AddPoint(xmax, ymax)
    ring_cslc.AddPoint(xmin, ymax)
    ring_cslc.AddPoint(xmin, ymin)
    bbox_cslc.AddGeometry(ring_cslc)

    # Define the SRS for CSLC and land polygon
    srs_cslc = osr.SpatialReference()
    srs_cslc.ImportFromEPSG(epsg_cslc)

    srs_land = osr.SpatialReference()
    srs_land.ImportFromEPSG(land_epsg)

    # Reproject the bounding box (in CSLC EPSG) to land polygon's EPSG
    transformer_cslc_to_land = osr.CoordinateTransformation(srs_cslc, srs_land)
    # transformer_land_to_cslc = osr.CoordinateTransformation(srs_land, srs_cslc)

    bbox_cslc.Transform(transformer_cslc_to_land)
    # bbox_cslc.Transform(transformer_land_to_cslc)
    print('apap')

    # Return a numpy array full of `False` when there is no intersection
    print(epsg_cslc)
    print(land_epsg)
    # print('bbox', bbox_cslc)
    # print('land_polygon', land_polygon)

    if not bbox_cslc.Intersects(land_polygon):
        return np.full(shape_mask, False)
    print('apap2')

    # Compute the intersection and reproject the result back to CSLC's EPSG
    intersection_land = bbox_cslc.Intersection(land_polygon)
    transformer_land_to_cslc = osr.CoordinateTransformation(srs_land, srs_cslc)
    intersection_land.Transform(transformer_land_to_cslc)

    # Build up a vector layer, and add a feature that has `intersection_land`` as geometry
    drv_intersection_polygon = ogr.GetDriverByName('Memory')
    ds_intersection_polygon = drv_intersection_polygon.CreateDataSource(str(time.time_ns()))
    layer_intersection = ds_intersection_polygon.CreateLayer('layer_intersection',
                                                             srs_cslc,
                                                             ogr.wkbPolygon)
    feature_defn = layer_intersection.GetLayerDefn()
    feature = ogr.Feature(feature_defn)
    feature.SetGeometry(intersection_land)
    layer_intersection.CreateFeature(feature)

    # Prepare for output layer for the rasterization
    print('here')
    drv_raster_out = gdal.GetDriverByName('GTIFF')
    rasterized_land = drv_raster_out.Create(f'raster_{xmin}_{xmax}_{ymin}_{ymax}.tif',
                                            shape_mask[1], shape_mask[0],
                                            1, gdal.GDT_Byte)
    rasterized_land.SetGeoTransform(geotransform)
    rasterized_land.SetProjection(srs_cslc.ExportToWkt())

    gdal.RasterizeLayer(rasterized_land, [1], layer_intersection)
    mask_land = rasterized_land.ReadAsArray()
    rasterized_land = None

    save_raster_gdal(mask_land, 'test_landwater2.tif',
                     geotransform=geotransform,
                     projection=srs_cslc.ExportToWkt(),
                     scratch_dir='.',
                     datatype='float32')

    print('--')

    return mask_land


def _get_parser():
    parser = argparse.ArgumentParser(
        description='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-mgrs_db',
                        type=str,
                        dest='mgrs_db',
                        required=True,
                        help='JPL MGRS database sqlite')

    parser.add_argument('-tfdb',
                        type=str,
                        dest='tfdb',
                        required=True,
                        help='NISAR track frame database')

    parser.add_argument('-land_ocean_shapefile',
                        type=str,
                        required=True,
                        dest='land_ocean_shapefile',
                        help='Shapefile for land/ocean mask')

    parser.add_argument('-output',
                        type=str,
                        dest='output',
                        required=True,
                        help='prefix output')

    return parser

def get_label_landcover_esa_10():
    '''Get integer number information what they represent
    ESA 10 m
    https://viewer.esa-worldcover.org/worldcover/
    '''
    label = dict()

    label['Tree Cover'] = 10
    label['Shrubs'] = 20
    label['Grassland'] =30
    label['Crop'] = 40
    label['Urban'] = 50
    label['Bare sparse vegetation'] = 60
    label['Snow and Ice'] = 70
    label['Permanent water bodies'] = 80
    label['Herbaceous wetland'] =90
    label['Mangrove'] = 95
    label['Moss and lichen'] = 100
    label['No_data'] = 0

    return label

# def check_flag_land_ocean_landcover(landcover_data):
#     landcover_label = get_label_landcover_esa_10()
#     pixel_num = landcover_data.size
#     water_pixel = np.sum((landcover_data == landcover_label['No_data']) &
#                          (landcover_data == landcover_label['Permanent water bodies']))
#     if water_pixel == pixel_num:
#         ocean_land_class = 'water'
#     elif water_pixel == 0:
#         ocean_land_class = 'land'
#     else:
#         ocean_land_class = 'water/land'
    
#     return ocean_land_class

def check_flag_land_ocean_landcover(landcover_data):
    print('landcover value', np.nanmean(landcover_data))
    pixel_num = landcover_data.size
    land_pixel = np.sum(landcover_data/255)
    print('check_flag_land_ocean_landcover', int(land_pixel), pixel_num, int(land_pixel)==pixel_num)

    if int(land_pixel) == pixel_num:
        ocean_land_class = 'water'
    elif land_pixel == 0:
        ocean_land_class = 'land'
    else:
        ocean_land_class = 'water/land'

    return ocean_land_class


def utm_to_latlon(utm_epsg, utm_easting, utm_northing):
    """
    Convert UTM coordinates to latitude and longitude using a given EPSG code.

    Parameters:
    -----------
    utm_epsg : int or str
        The EPSG code for the UTM zone in which the coordinates are provided.
    utm_easting : float
        The easting (X) value of the UTM coordinate.
    utm_northing : float
        The northing (Y) value of the UTM coordinate.

    Returns:
    --------
    tuple (float, float)
        A tuple containing the longitude and latitude values, respectively.
    """

    # Create a pyproj transformer from the UTM projection to WGS84 (latitude/longitude)
    transformer = Transformer.from_crs(f"EPSG:{utm_epsg}", "EPSG:4326", always_xy=True)

    # Transform the UTM coordinates to latitude and longitude.
    lon, lat = transformer.transform(utm_easting, utm_northing)
    return lon, lat

def latlon_to_utm(lat, lon, utm_epsg):
    """
    Convert latitude and longitude to UTM coordinates using a given EPSG code.

    Parameters:
    -----------
    lat : float
        Latitude value.
    lon : float
        Longitude value.
    utm_epsg : int or str
        The EPSG code for the desired UTM zone for the output coordinates.

    Returns:
    --------
    tuple (float, float)
        A tuple containing the easting (X) and northing (Y) values of the UTM coordinate.
    """

    # Initialize a transformer object from WGS84 to the specified UTM EPSG.
    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}", always_xy=True)

    # Transform the latitude and longitude to UTM coordinates.
    utm_x, utm_y = transformer.transform(lon, lat)

    return utm_x, utm_y

def transform_to_lower_precision(polygon, decimal_places=4):
    """
    Round polygon coordinates. Works for both 2D and 3D polygons.
    """
    coords_out = []
    for c in list(polygon.exterior.coords):
        if len(c) == 3:
            x, y, z = c
            coords_out.append((round(x, decimal_places),
                               round(y, decimal_places),
                               round(z, decimal_places)))
        else:
            x, y = c
            coords_out.append((round(x, decimal_places),
                               round(y, decimal_places)))
    return Polygon(coords_out)




def reproject_coordinates(x, y, src_epsg: int, dst_epsg: int):
    """
    Reproject coordinates from one EPSG to another.

    Parameters
    ----------
    x : float or array-like
        X coordinate(s). Usually longitude for EPSG:4326, easting for UTM,
        or projected x for polar stereographic.
    y : float or array-like
        Y coordinate(s). Usually latitude for EPSG:4326, northing for UTM,
        or projected y for polar stereographic.
    src_epsg : int
        Source EPSG code.
    dst_epsg : int
        Destination EPSG code.

    Returns
    -------
    x_out, y_out : float or np.ndarray
        Reprojected coordinates.
    """

    transformer = Transformer.from_crs(
        f"EPSG:{src_epsg}",
        f"EPSG:{dst_epsg}",
        always_xy=True
    )

    x_out, y_out = transformer.transform(x, y)

    return x_out, y_out


def get_mgrs_tile_group(mgrs_tile):
    '''Get MGRS tile and return MGRS group id counted referring to equator.
    Starting from equator, MGRSs in row (same latitude) have same mgrs_group_id.
    Also, the MGRS tile collection will consists of two rows of the MGRS tiles.
    The MGRS tiles in adjacent two rows have same group_2row.

    Parameters:
    -----------
    mgrs_tile : str
        MGRS tile representation.

    Returns:
    --------
    int, int
        The MGRS group ID and the secondary group value.

    Notes:
    ------
    Starting from the equator, MGRSs in the same row (same latitude) have the same mgrs_group_id.
    In the northern hemisphere, mgrs_group_id starts from 1 and increases towards the Arctic pole.
    In the southern hemisphere, it starts from -1 and decreases towards the Antarctic pole.
    """
    '''
    # MGRS 100,000-meter square identification letters
    # I and O are omitted in the MGRS 100,000-meter square identification

    # 6 deg x 8 deg
    main_row_letter = 'NPQRSTUVWXCDEFGHJKLM'
    # 100 kmx 100km
    column_letters = 'ABCDEFGHJKLMNPQRSTUVWXYZ'
    row_letters = 'ABCDEFGHJKLMNPQRSTUV'

    # number of subtiles in each character in main_row_letter
    # if mgrs is used.
    # sub_row_num = [9, 10, 10, 11, 10, 10, 10, 10, 10, 15]

    # if kml file is used.
    sub_row_num = [9,  9,  9,  8,  9,  9,  9,  9,  9, 14]

    # Define the number of rows per group
    rows_per_group = 2

    # Extract the UTM zone and 100,000-meter square identification
    # from the MGRS tile
    utm_zone = int(mgrs_tile[:2])
    main_row_zone = str(mgrs_tile[2])
    square_id = mgrs_tile[3:]

    # Get the column and row index
    main_row_index_north_south = main_row_letter.index(main_row_zone)

    # if main_row_zone in 'NPQRSTUVWX'
    if main_row_index_north_south < 10:
        hemisphere = 'North'
        main_row_index = main_row_index_north_south + 1
    else:
        hemisphere = 'South'
        main_row_letter_rev = main_row_letter[::-1]
        main_row_index = main_row_letter_rev.index(main_row_zone) + 1

    row_letters_rep = row_letters * 10

    if hemisphere == 'North':
        if utm_zone % 2 == 1:
            starting_row = row_letters.index('A')
        else:
            starting_row = row_letters.index('F')

        row_rep_start = int(np.sum(sub_row_num[:main_row_index-1]))

        if main_row_index-1 == 0:
            row_rep_start = 0
        row_letters_sub = row_letters_rep[row_rep_start + starting_row:]

    else:
        row_letters_rep_rev = row_letters_rep[::-1]

        if utm_zone % 2 == 1:
            starting_row = row_letters_rep_rev.index('V')
        else:
            starting_row = row_letters_rep_rev.index('E')

        row_rep_start = int(np.sum(sub_row_num[:main_row_index-1]))
        row_letters_sub = row_letters_rep_rev[row_rep_start + starting_row:]

    sub_row_index = row_letters_sub.index(square_id[1])

    # Calculate the group number
    if hemisphere == 'North':
        mgrs_group_id = (row_rep_start) + (sub_row_index) +1
        group_2row = int(np.ceil(mgrs_group_id / rows_per_group))
    else:
        mgrs_group_id = -(row_rep_start) - (sub_row_index) -1
        group_2row = mgrs_group_id // rows_per_group

    return mgrs_group_id, group_2row

def read_polygons_from_tfdb(
        geopakage,
        provided_multipolygonz,
        relative_orbit_number,
        orbit_pass):
    """
    Reads polygons from an SQLite database and returns intersecting polygons
    with a provided multipolygonz.

    Parameters:
    -----------
    sqlite_file : str
        Path to the SQLite database file.
    table_name : str
        Name of the table in the SQLite database.
    provided_multipolygonz : Polygon
        Polygon to check for intersection.
    relative_orbit_number : int
        Relative orbit number to filter the results.
    orbit_pass : str
        Orbit pass direction ('Ascending' or 'Descending').

    Returns:
    --------
    list, list, list, list, list
        Intersecting polygons, and burst IDs for different subswaths.
    """
    # Make sure the SQLite file path is absolute
    absolute_tfdb_geopack = os.path.abspath(geopakage)
    tfdb_gdf = gpd.read_file(absolute_tfdb_geopack)
    
    # Near the equator, the track number changes when the satellite is acending.
    # dual_track_number_flag represents the flag whether we need to search DB with two
    # track numbers.
    dual_track_number_flag = False
    if provided_multipolygonz.intersects(equator_bound_polygon) and \
        orbit_pass == 'Ascending':
        if relative_orbit_number == 1:
            relative_orbit_number_second = 173
            relative_orbit_number_third = 2

        elif relative_orbit_number >= 173:
            relative_orbit_number_second = 1
            relative_orbit_number_third = 172
        else:
            relative_orbit_number_second = relative_orbit_number - 1
            relative_orbit_number_third = relative_orbit_number + 1

        dual_track_number_flag = True
    if provided_multipolygonz.intersects(north_asc_des_changing_polygon) or \
    provided_multipolygonz.intersects(south_asc_des_changing_polygon):
        asc_dec_change_flag = True
    else:
        asc_dec_change_flag = False
    # Prepare the query for filtering the polygons
    if dual_track_number_flag:
        filtered_gdf = tfdb_gdf[((tfdb_gdf['track'] == relative_orbit_number) |
                                 (tfdb_gdf['track'] == relative_orbit_number_second) |
                                 (tfdb_gdf['track'] == relative_orbit_number_third)) &
                                 (tfdb_gdf['passDirection'] == orbit_pass)]

    else:
        if asc_dec_change_flag:
            filtered_gdf = tfdb_gdf[(tfdb_gdf['track'] == relative_orbit_number)]
        else:
            filtered_gdf = tfdb_gdf[(tfdb_gdf['track'] == relative_orbit_number)  &
                                    (tfdb_gdf['passDirection'] == orbit_pass)]
    print('read_polygons_from_tfdb')
    print(filtered_gdf)
    intersecting_polygons = []

    frame_list = []
    tf_list = []
    for index, row in filtered_gdf.iterrows():
        polygon = row.geometry
        if polygon.intersects(provided_multipolygonz):
            frame_list.append(row['frame'])
            intersecting_polygons.append(polygon)
            tf_list.append(f"{row['track']}_{row['frame']}")

            # subswath_name = row[col_id['subswath_name']]
            # if subswath_name in burst_id_lists:
            #     burst_id_lists[subswath_name].append(row[col_id['burst_id']])

    return (intersecting_polygons,
            frame_list,
            tf_list)

def read_polygons_from_jpl_burst_id_db(
        sqlite_file,
        table_name,
        provided_multipolygonz,
        relative_orbit_number,
        orbit_pass):
    """
    Reads polygons from an SQLite database and returns intersecting polygons
    with a provided multipolygonz.

    Parameters:
    -----------
    sqlite_file : str
        Path to the SQLite database file.
    table_name : str
        Name of the table in the SQLite database.
    provided_multipolygonz : Polygon
        Polygon to check for intersection.
    relative_orbit_number : int
        Relative orbit number to filter the results.
    orbit_pass : str
        Orbit pass direction ('Ascending' or 'Descending').

    Returns:
    --------
    list, list, list, list, list
        Intersecting polygons, and burst IDs for different subswaths.
    """
    # Make sure the SQLite file path is absolute
    absolute_sqlite_file = os.path.abspath(sqlite_file)

    # Connect to the SQLite database
    conn_burst = sqlite3.connect(absolute_sqlite_file)
    cur_burst = conn_burst.cursor()

    # Near the equator, the track number changes when the satellite is acending.
    # dual_track_number_flag represents the flag whether we need to search DB with two
    # track numbers.
    dual_track_number_flag = False
    if provided_multipolygonz.intersects(equator_bound_polygon) and \
        orbit_pass == 'Ascending':
        if relative_orbit_number == 1:
            relative_orbit_number_second = 173
        elif relative_orbit_number > 173:
            relative_orbit_number_second = 1
        else:
            relative_orbit_number_second = relative_orbit_number - 1
        dual_track_number_flag = True

    # Prepare the query for filtering the polygons
    if dual_track_number_flag:
        conditions = f"""relative_orbit_number = {relative_orbit_number}
                      AND orbit_pass = '{orbit_pass}'
                      OR relative_orbit_number = {relative_orbit_number_second}"""

    else:
        conditions = f"""relative_orbit_number = {relative_orbit_number}
                      AND orbit_pass = '{orbit_pass}'"""

        relative_orbit_number_second = relative_orbit_number

    query = f"SELECT * FROM {table_name} WHERE {conditions}"
    query_result_all = cur_burst.execute(query)
    result = query_result_all.fetchall()

    col_id = {column_desc[0]: i
              for i, column_desc in enumerate(query_result_all.description)}

    intersecting_polygons = []

    burst_id_lists = {'IW1': [], 'IW2': [], 'IW3': []}
    burst_jpl_list = []

    for row in result:
        polygon = wkb_loads(row[col_id["GEOMETRY"]])
        if polygon.intersects(provided_multipolygonz):
            burst_jpl_list.append(row[col_id['burst_id_jpl']])
            intersecting_polygons.append(polygon)
            subswath_name = row[col_id['subswath_name']]
            if subswath_name in burst_id_lists:
                burst_id_lists[subswath_name].append(row[col_id['burst_id']])

    # Close the SQLite database connection
    cur_burst.close()
    conn_burst.close()
    return (intersecting_polygons,
            burst_id_lists['IW1'],
            burst_id_lists['IW2'],
            burst_id_lists['IW3'],
            burst_jpl_list)

def read_mgrs_polygons_from_sqlite(sqlite_file,
                                   table_name,
                                   north_south='north'):
    """
    Reads polygons from an SQLite database based on their location
    in the northern or southern hemisphere.

    Parameters:
    -----------
    sqlite_file : str
        Path to the SQLite database file.
    table_name : str
        Name of the table in the SQLite database containing polygons.
    north_south : str
        Specifies whether to retrieve polygons from the 'north' or 'south'. Defaults to 'north'.

    Returns:
    --------
    list, list
        List of intersecting polygons and associated MGRS list.
    """
    # Connect to the SQLite database
    conn = sqlite3.connect(os.path.abspath(sqlite_file))
    cur = conn.cursor()

    # Execute a query to get the geometry and name from the 'polygons' table
    cur.execute(f"SELECT mgrs_tile, GEOMETRY FROM {table_name}")

    intersecting_polygons = []
    mgrs_list = []

    # Iterate over the rows in the query result
    for row in cur:
        name = row[0]
        # Parse the geometry from the Well-Known Binary (WKB) format
        geometry = ogr.CreateGeometryFromWkb(bytes(row[1]))
        centroid_y = geometry.Centroid().GetY()

        # Check hemisphere condition
        if (north_south == 'north' and centroid_y >= 0) or \
            (north_south == 'south' and centroid_y < 0):
            mgrs_list.append(name)
            intersecting_polygons.append(geometry)

    # Close the SQLite connection
    conn.close()

    return intersecting_polygons, mgrs_list

def test2(major_epsg, polygon):

    x_min_ext, y_min_ext, x_max_ext, y_max_ext = polygon.bounds
    print('bound', x_min_ext, x_max_ext, y_min_ext, y_max_ext)
    major_epsg = int(major_epsg)

    x_spacing = 5000
    y_spacing = 5000
    x_spacing = 0.05000
    y_spacing = 0.05000

    height = int((y_max_ext - y_min_ext) / y_spacing)
    width = int((x_max_ext - x_min_ext) / x_spacing)

    if (height == 0) or (width == 0):
        x_spacing = 0.00100
        y_spacing = 0.00100
        height = int((y_max_ext - y_min_ext) / y_spacing)
        width = int((x_max_ext - x_min_ext) / x_spacing)

    mask_land = _get_land_mask(
        major_epsg,
        (x_min_ext, x_spacing, 0, y_min_ext, 0, y_spacing),
        (height, width))

    srs_mosaic = osr.SpatialReference()
    srs_mosaic.ImportFromEPSG(major_epsg)
    wkt_projection = srs_mosaic.ExportToWkt()

    land_ocean_mask = check_flag_land_ocean_landcover(mask_land)

    return land_ocean_mask

def chunk_list(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def split_frames_by_latitude(intersecting_polygons, frame_list, tf_list, max_frames=3):
    """
    Sort frames by burst centroid latitude, then chunk into groups of size max_frames.
    Returns list of dicts: {"frames": [...], "tfs": [...], "polys": [...]}
    """
    items = []
    for poly, fr, tf in zip(intersecting_polygons, frame_list, tf_list):
        try:
            lat = poly.centroid.y
        except Exception:
            # fallback
            lat = poly.bounds[1] * 0.5 + poly.bounds[3] * 0.5
        items.append((lat, fr, tf, poly))

    # north->south (descending latitude). If you want south->north, flip reverse.
    items.sort(key=lambda x: x[0], reverse=True)

    groups = []
    for chunk in chunk_list(items, max_frames):
        groups.append({
            "frames": [c[1] for c in chunk],
            "tfs":    [c[2] for c in chunk],
            "polys":  [c[3] for c in chunk],
        })
    return groups

def is_polar_geometry(geom):
    """
    Decide if this is polar-ish. You can tune this threshold.
    Uses your predefined north_pol_polygon / south_pol_polygon (lon/lat).
    """
    try:
        return geom.intersects(north_pol_polygon) or geom.intersects(south_pol_polygon)
    except Exception:
        return False

import numpy as np


def chunk_list(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]


def wrap_lon_relative(lon, ref_lon):
    """
    Wrap lon so it is continuous around ref_lon.

    Output is in [ref_lon - 180, ref_lon + 180).
    """
    lon = ((lon + 180) % 360) - 180
    ref_lon = ((ref_lon + 180) % 360) - 180

    d = lon - ref_lon
    d = ((d + 180) % 360) - 180

    return ref_lon + d


def circular_mean_lon(lons):
    """
    Compute mean longitude safely near the antimeridian.
    """
    lons_rad = np.deg2rad(lons)

    sin_mean = np.mean(np.sin(lons_rad))
    cos_mean = np.mean(np.cos(lons_rad))

    mean_lon = np.rad2deg(np.arctan2(sin_mean, cos_mean))

    return float(mean_lon)

def _min_frame_number(group):
    """
    Return the smallest frame number in a frame group.

    Handles frame values stored as int, float, or numeric string.
    """
    return min(int(fr) for fr in group["frames"])


def split_frames_by_longitude(
    intersecting_polygons,
    frame_list,
    tf_list,
    max_frames=3,
    min_frames=3,
    overlap_frames=1,
    ref_lon=None,
):
    """
    Sort frames by centroid longitude using antimeridian-safe longitude wrapping,
    then split them into overlapping groups.

    Each returned group has:
        - frames/tfs/polys: overlapped frames used for processing context
        - core_frames/core_tfs/core_polys: non-overlapped frames used for unique MGRS assignment
    """

    if max_frames < min_frames:
        raise ValueError(
            f"max_frames must be >= min_frames. "
            f"Got max_frames={max_frames}, min_frames={min_frames}."
        )

    if overlap_frames < 0:
        raise ValueError("overlap_frames must be >= 0.")

    if overlap_frames >= max_frames:
        raise ValueError(
            "overlap_frames must be smaller than max_frames. "
            f"Got overlap_frames={overlap_frames}, max_frames={max_frames}."
        )

    n_items = len(frame_list)

    if not (
        len(intersecting_polygons) == len(frame_list) == len(tf_list)
    ):
        raise ValueError(
            "intersecting_polygons, frame_list, and tf_list must have the same length."
        )

    if n_items == 0:
        return []

    if ref_lon is None:
        lons = []

        for poly in intersecting_polygons:
            try:
                lon = poly.centroid.x
            except Exception:
                b = poly.bounds
                lon = 0.5 * (b[0] + b[2])

            lons.append(lon)

        ref_lon = circular_mean_lon(lons)

    print("splitting collection by longitude")
    print("reference longitude:", ref_lon)

    items = []

    for poly, fr, tf in zip(intersecting_polygons, frame_list, tf_list):
        try:
            lon = poly.centroid.x
        except Exception:
            b = poly.bounds
            lon = 0.5 * (b[0] + b[2])

        lon_u = wrap_lon_relative(lon, ref_lon)

        items.append(
            {
                "lon_u": lon_u,
                "frame": fr,
                "tf": tf,
                "poly": poly,
            }
        )

    items.sort(key=lambda x: x["lon_u"])

    if n_items <= max_frames:
        return [
            {
                "frames": [item["frame"] for item in items],
                "tfs": [item["tf"] for item in items],
                "polys": [item["poly"] for item in items],
                "lon_u": [item["lon_u"] for item in items],

                # no split, so core == full group
                "core_frames": [item["frame"] for item in items],
                "core_tfs": [item["tf"] for item in items],
                "core_polys": [item["poly"] for item in items],
                "core_lon_u": [item["lon_u"] for item in items],
            }
        ]

    step = max_frames - overlap_frames

    groups = []
    ranges = []
    start = 0
    used_ranges = set()

    while start < n_items:
        end = min(start + max_frames, n_items)

        if end == n_items and (end - start) < min_frames:
            start = max(0, n_items - min_frames)
            end = n_items

        range_key = (start, end)
        if range_key in used_ranges:
            break

        used_ranges.add(range_key)
        ranges.append((start, end))

        if end == n_items:
            break

        start += step

    # Define non-overlapping core ranges.
    # Core boundary is halfway between neighboring group starts/ends.
    for i, (start, end) in enumerate(ranges):
        group_items = items[start:end]

        if len(ranges) == 1:
            core_start = start
            core_end = end
        else:
            if i == 0:
                core_start = start
            else:
                prev_start, prev_end = ranges[i - 1]
                core_start = (prev_end + start) // 2

            if i == len(ranges) - 1:
                core_end = end
            else:
                next_start, next_end = ranges[i + 1]
                core_end = (end + next_start) // 2

        # Safety fallback
        if core_end <= core_start:
            core_start = start
            core_end = end

        core_items = items[core_start:core_end]

        groups.append(
            {
                # full overlapped group
                "frames": [item["frame"] for item in group_items],
                "tfs": [item["tf"] for item in group_items],
                "polys": [item["poly"] for item in group_items],
                "lon_u": [item["lon_u"] for item in group_items],

                # non-overlapped core group
                "core_frames": [item["frame"] for item in core_items],
                "core_tfs": [item["tf"] for item in core_items],
                "core_polys": [item["poly"] for item in core_items],
                "core_lon_u": [item["lon_u"] for item in core_items],
            }
        )

    groups.sort(key=_min_frame_number)

    return groups


def assign_mgrs_tiles_to_frame_groups(
    frame_groups,
    target_polygon,
    target_mgrs_list,
    group_list,
    unique_group_number,
):
    """
    Assign each MGRS tile to exactly one frame group.

    Full group polygons may overlap because neighboring groups share frames.
    Therefore, MGRS ownership is decided using core_polys if available.

    Returns
    -------
    assigned : list of dict
        Same length as frame_groups.
        Each dict contains:
            {
                "mgrs_list": [...],
                "mgrs_polys": [...]
            }
    """

    assigned = [
        {
            "mgrs_list": [],
            "mgrs_polys": [],
        }
        for _ in frame_groups
    ]

    # Build group ownership geometries.
    # Use core polygons to avoid overlap-driven duplicate MGRS assignment.
    ownership_geoms = []

    for fg in frame_groups:
        core_polys = fg.get("core_polys", None)

        if core_polys:
            ownership_geoms.append(unary_union(core_polys))
        else:
            ownership_geoms.append(unary_union(fg["polys"]))

    for group_ind, group_number in enumerate(group_list):
        if group_number != unique_group_number:
            continue

        single_mgrs = wkt_loads(target_polygon[group_ind].ExportToWkt())
        single_mgrs_name = target_mgrs_list[group_ind]

        best_group_idx = None
        best_area = 0.0

        for idx, ownership_geom in enumerate(ownership_geoms):
            if not ownership_geom.intersects(single_mgrs):
                continue

            inter_area = ownership_geom.intersection(single_mgrs).area

            if inter_area > best_area:
                best_area = inter_area
                best_group_idx = idx

        # If no core region intersects, fallback to full overlapped group polygons.
        if best_group_idx is None or best_area == 0:
            for idx, fg in enumerate(frame_groups):
                full_geom = unary_union(fg["polys"])

                if not full_geom.intersects(single_mgrs):
                    continue

                inter_area = full_geom.intersection(single_mgrs).area

                if inter_area > best_area:
                    best_area = inter_area
                    best_group_idx = idx

        if best_group_idx is None or best_area == 0:
            continue

        assigned[best_group_idx]["mgrs_list"].append(single_mgrs_name)
        assigned[best_group_idx]["mgrs_polys"].append(single_mgrs)

    return assigned



def run(args):
    '''
    1) Loop starts from track number 1 to track number 173
    2) set processing order Ascending in North Hemisphere -> Descending in North Hemisphere->
    Descending in South Hemisphere -> Ascending in north hemisphere -> next track.
    3) Find MGRS tiles in north-hemisphere first and find sort them using ID.
        two adjacent MGRS tiles in vertical direction have the same ID.

        --------------------
        | 2  | 2  | 2  | 2  |
        --------------------
        | 1  | 1  | 1  | 1  |
        --------------------
        | 1  | 1  | 1  |  1 |
        --------------------    equator
        | -1 | -1 | -1 | -1 |
        --------------------
        | -1 | -1 | -1 | -1 |
        --------------------
    4) First, use the MGRS group ID 1 and define the bounds in y directions.
    5) Read valid bursts near equator and define bounds in x directions.
    6) Find all bursts within x and y bounds
    7) Excludes MGRS tiles from MGRS tile collections not overlapped with all bursts.
    8) Update valid bursts for next loop and re-process for next MGRS tile ID.
    '''

    output_prefix = args.output
    land_ocean_shapefile = args.land_ocean_shapefile

    mgrs_polygon_sqlite_file = args.mgrs_db
    mgrs_collection_db_path = f'{output_prefix}.sqlite'

    if os.path.isfile(mgrs_collection_db_path):
        os.remove(mgrs_collection_db_path)

    land_ocean_flag = False
    if land_ocean_shapefile is not None:
        if os.path.isfile(land_ocean_shapefile):
            land_ocean_flag = True
        else:
            err_str = 'Landcover file is not found'
            raise FileNotFoundError(err_str)

    # Connect to the NISAR track frame database sqlite
    tfdb_geopack = args.tfdb
    tfdb_gdf = gpd.read_file(tfdb_geopack)

    north_polygon, north_mgrs_list = read_mgrs_polygons_from_sqlite(
        mgrs_polygon_sqlite_file,
        table_name='MGRS_tile',
        north_south='north')
    south_polygon, south_mgrs_list = read_mgrs_polygons_from_sqlite(
        mgrs_polygon_sqlite_file,
        table_name='MGRS_tile',
        north_south='south')

    group_list_north = []
    group_name_north = []
    # 'group_list' consists of integer group number.
    for mgrs_cand in north_mgrs_list:
        # MGRS tiles in two adjacent rows are assgined to same group.
        _, group = get_mgrs_tile_group(mgrs_cand)
        group_list_north.append(group)
        group_name_north.append(mgrs_cand[:5])

    # group_list_north,     south_polygon, south_mgrs_list
    group_list_south = []
    group_name_south = []
    for mgrs_cand in south_mgrs_list:
        _, group = get_mgrs_tile_group(mgrs_cand)
        group_list_south.append(group)
        group_name_south.append(mgrs_cand[:5])

    # Connect to the MGRS tile collection database (SQLite)
    mgrs_db_conn = sqlite3.connect(mgrs_collection_db_path)
    mgrs_db_conn.enable_load_extension(True)

    mgrs_db_cursor = mgrs_db_conn.cursor()
    mgrs_db_cursor.execute("SELECT load_extension('mod_spatialite');")
    mgrs_db_cursor.execute("SELECT InitSpatialMetaData(1);")

    # Create the 'mgrs_burst_db' table
    mgrs_db_cursor.execute('''
        CREATE TABLE IF NOT EXISTS mgrs_track_frame_db (
            id INTEGER PRIMARY KEY,
            mgrs_set_id TEXT UNIQUE,
            track_number INTEGER,
            frames TEXT,
            track_frame TEXT,
            number_of_frames INTEGER,
            mgrs_tiles TEXT,
            number_of_mgrs_tiles INTEGER,
            orbit_pass TEXT,
            land_ocean_flag TEXT,
            EPSG INTEGER,
            xmin INTEGER,
            xmax INTEGER,
            ymin INTEGER,
            ymax INTEGER
        );
    ''')

    mgrs_db_cursor.execute('''
        CREATE TABLE IF NOT EXISTS mgrs_frame_check_db (
            id INTEGER PRIMARY KEY,
            mgrs_set_id TEXT,
            track_number INTEGER,
            orbit_pass TEXT,
            unique_group_number INTEGER,
            original_number_of_frames INTEGER,
            original_frames TEXT,
            original_track_frame TEXT,
            split_applied INTEGER,
            reason TEXT
        );
    ''')
    mgrs_db_cursor.execute(
        '''
        SELECT AddGeometryColumn('mgrs_track_frame_db', 'geometry', 4326, 'MULTIPOLYGON', 3)
        ''')
    mgrs_db_cursor.execute(
        "SELECT CreateSpatialIndex('mgrs_track_frame_db', 'geometry');")

    # Loop from track number 1 to track number 173
    for track_number in range(1, 174): #range(8, 9): #range(1, 174):

        mgrs_set_number = 1

        orbit_passes = ['Ascending', 'Descending', 'Descending', 'Ascending']
        group_lists = [group_list_north,
                       group_list_north,
                       group_list_south,
                       group_list_south]
        # orbit_passes = ['Ascending']
        # group_lists = [group_list_north]

        for orbit_pass, group_list in zip(orbit_passes, group_lists):
            # Search Frame with track number and orbit pass
            filtered_gdf = tfdb_gdf[(tfdb_gdf['track'] == track_number) & 
                                    (tfdb_gdf['passDirection'] == orbit_pass)]

            # Find first valid frame
            # for jjjind in range(0, 100):
            for jjjind, row in filtered_gdf.head().iterrows():
                row_start = row
                geometry_valid = row.geometry
                # Bursts in high latitude are not used for MGRS tile collection.
                if north_bound_polygon.contains(geometry_valid) or \
                   north_bound_polygon.intersects(geometry_valid):
                    print('skip', row_start)
                else:
                    break

            # Read burst id
            # burst_id_start = row_start[burst_fields_name['burst_id']]
            burst_utm_zone = int(str(row_start['epsg'])[-2:])

            # Unique MGRS tile id
            unique_group_list = list(set(group_list))

            # mgrs group id increase from south pole to north pole.
            # in descending orbit, the unique group list needs to be revsersed.
            if orbit_pass == 'Descending':
                unique_group_list.reverse()

            # Unique_group_list consists of the unique MGRS tile id starting from
            # equator. Two rows in MGRS tiles have same ID.
            for unique_group_number in unique_group_list:
                # Antarctica is excluded.
                if unique_group_number > -44:
                    print('UNIQUE GROUP ID:', unique_group_number)
                    print('Burst track number:', track_number)
                    print('Orbit pass:', orbit_pass)
                    union_polygon = ogr.Geometry(ogr.wkbPolygon)

                    # Read polygon from MGRS tiles having unique MGRS ID.
                    for group_ind, group_number in enumerate(group_list):

                        # Group number has duplicated IDs in adjacent MGRS tiles.
                        if group_number == unique_group_number:

                            if group_number > 0:
                                polygon_single = north_polygon[group_ind]
                            else:
                                polygon_single = south_polygon[group_ind]
                            union_polygon = union_polygon.Union(polygon_single)

                    mgrs_poly = wkt_loads(union_polygon.ExportToWkt())
                    mgrs_poly0 = mgrs_poly

                    # burst bouding box
                    burst_geometry = row_start.geometry
                    burst_id_in_track = row_start['frame']
                    print(burst_id_in_track)
                    orbit_pass = row_start['passDirection']
                    utm_x_min_frame = row_start['mapTopLeftX']
                    utm_x_max_frame = row_start['mapBottomRightX']
                    utm_y_min_frame = row_start['mapBottomRightY']
                    utm_y_max_frame = row_start['mapTopLeftY']
                    epsg_frame = row_start['epsg']
                    print('row_start')
                    print(row_start)
                    burst_center_lon, burst_center_lat = burst_geometry.centroid.x, burst_geometry.centroid.y
                    print('center ', burst_center_lon, burst_center_lat)
                    burst_hemisphere = 'north' if burst_center_lat >= 0  else 'south'
                    epsg_code = 32600 + burst_utm_zone if burst_hemisphere == 'north' else 32700 + burst_utm_zone

                    intersect_x, intersect_y = latlon_to_utm(burst_center_lat,
                                                             burst_center_lon,
                                                             epsg_code)
                    mean_x_frame = (utm_x_max_frame + utm_x_min_frame) / 2
                    mean_y_frame = (utm_y_max_frame + utm_y_min_frame) / 2
                    if epsg_code != epsg_frame:
                        print(f'NI Frame EPSG {epsg_frame} is not equal to {epsg_code}')

                        intersect_x, intersect_y = reproject_coordinates(
                            mean_x_frame, mean_y_frame, src_epsg=epsg_frame, dst_epsg=epsg_code)
                    else:
                        intersect_x, intersect_y = mean_x_frame, mean_y_frame

                    print('utm coordinate of center', intersect_x, intersect_y, epsg_code)
                    # Intentionally add buffer to get adjacent bursts
                    if orbit_pass == 'Ascending':
                        burst_xmin = intersect_x - 530000
                        burst_xmax = intersect_x + 700000
                    else:
                        burst_xmin = intersect_x - 700000
                        burst_xmax = intersect_x + 530000
                    burst_xmin, _ = utm_to_latlon(epsg_code,
                                                  burst_xmin,
                                                  intersect_y)
                    burst_xmax, _= utm_to_latlon(epsg_code,
                                                 burst_xmax,
                                                 intersect_y)

                    # Check if burst covers antimeridian
                    if (burst_xmin - burst_xmax > 180):
                        print('Antimeridian detected', burst_xmin, burst_xmax)
                        burst_poly_coords = [Polygon(((burst_xmin, -84, 0),
                                                     (burst_xmin, 84.7, 0),
                                                     (180, 84.7, 0),
                                                    (180, -84, 0))),
                                            Polygon(((-180, -84, 0),
                                                    (-180, 80.7, 0),
                                                    (burst_xmax, 80.7, 0),
                                                    (burst_xmax, -84, 0)))]

                    else:
                        burst_poly_coords = [Polygon(((burst_xmin, -84, 0),
                                                     (burst_xmin, 84.7, 0),
                                                     (burst_xmax, 84.7, 0),
                                                     (burst_xmax, -84, 0)))]
                    burst_multi_polygon = MultiPolygon(burst_poly_coords)

                    print("group:", unique_group_number, "apply_gate:", (-30 < unique_group_number < 30))
                    print("burst_xmin/xmax:", burst_xmin, burst_xmax, "span:", burst_xmax - burst_xmin)
                    print("mgrs_poly0 bounds:", mgrs_poly0.bounds)

                    if (unique_group_number < 30) and (
                        unique_group_number > -30):
                        # By finding intersecting areas between polygon computed from bursts and
                        # unique MGRS group, extract the tiles MGRS within bounds in x direction.
                        # y bounds are computed from unique MGRS tile group and x bounds are determined
                        # from bursts.
                        mgrs_poly2 = mgrs_poly.intersection(burst_multi_polygon)
                        print("mgrs_poly bounds before:", mgrs_poly.bounds)
                        print("mgrs_poly bounds after :", mgrs_poly2.bounds)
                        print("changed:", not mgrs_poly2.equals(mgrs_poly), "area ratio:", (mgrs_poly2.area / mgrs_poly.area) if mgrs_poly.area else None)
                        mgrs_poly = mgrs_poly2
                        mgrs_poly = mgrs_poly.intersection(burst_multi_polygon)

                    # Find intersecting bursts from all bursts along the specific track number
                    intersecting_burst, frame_list, tf_list = \
                            read_polygons_from_tfdb(
                                tfdb_geopack,
                                mgrs_poly,
                                relative_orbit_number=track_number,
                                orbit_pass=orbit_pass)
                    print(tf_list)
                    print(intersecting_burst)
                    if len(frame_list) == 1:
                        print('found 1 frame')
                    print('track', track_number)
                    minx_test, miny_test, maxx_test, maxy_test = mgrs_poly.bounds

                    print("min x:", minx_test)
                    print("max x:", maxx_test)
                    print("min y:", miny_test)
                    print("max y:", maxy_test)
                    print('burst')
                    minx_test, miny_test, maxx_test, maxy_test = burst_multi_polygon.bounds

                    print("min x:", minx_test)
                    print("max x:", maxx_test)
                    print("min y:", miny_test)
                    print("max y:", maxy_test)
                    if len(frame_list) == 0 or len(frame_list) == 1:
                        print('here: frame_list is empty')
                        intersecting_burst, frame_list, tf_list = \
                            read_polygons_from_tfdb(
                                tfdb_geopack, 
                                mgrs_poly0,
                                relative_orbit_number=track_number,
                                orbit_pass=orbit_pass)

                    if frame_list:
                        # Decide whether to split (polar + too many frames)
                        do_split = is_polar_geometry(mgrs_poly0) and (len(frame_list) > 3)
                        print(f'splitting collection {do_split} because number of frame {len(frame_list)}')
                        if unique_group_number > 0:
                            target_polygon = north_polygon
                            target_mgrs_list = north_mgrs_list
                        else:
                            target_polygon = south_polygon
                            target_mgrs_list = south_mgrs_list

                        if do_split:
                            print('splitting inputs')
                            print(intersecting_burst)
                            print(frame_list)
                            print(tf_list)
                            print(burst_center_lon)
                            frame_groups = split_frames_by_longitude(
                                intersecting_polygons=intersecting_burst,
                                frame_list=frame_list,
                                tf_list=tf_list,
                                max_frames=3,
                                ref_lon=None  # you already computed this earlier
                                )
                            print('frame_Groups')
                            print(frame_groups)
                        else:
                            frame_groups = [{
                                "frames": frame_list,
                                "tfs": tf_list,
                                "polys": intersecting_burst
                            }]

                        # For each group, select mgrs tiles that intersect ONLY that group's burst union
                        # for group_idx, fg in enumerate(frame_groups, start=1):
                        #     group_frames = fg["frames"]
                        #     group_tfs = fg["tfs"]
                        #     group_polys = fg["polys"]

                        #     union_burst_polygon = unary_union(group_polys)

                        #     final_mgrs_list = []
                        #     final_mgrs_polygon_list = []

                        #     for group_ind, group_number in enumerate(group_list):
                        #         if group_number != unique_group_number:
                        #             continue

                        #         single_mgrs = wkt_loads(target_polygon[group_ind].ExportToWkt())

                        #         if union_burst_polygon.intersects(single_mgrs) or union_burst_polygon.contains(single_mgrs):
                        #             final_mgrs_list.append(target_mgrs_list[group_ind])
                        #             final_mgrs_polygon_list.append(single_mgrs)

                        #     # If somehow nothing matched, skip this subgroup
                        #     if not final_mgrs_list:
                        #         continue

                        #     final_mgrs_polygon = unary_union(final_mgrs_polygon_list).simplify(0.05)
                        #     if final_mgrs_polygon.geom_type == 'Polygon':
                        #         final_mgrs_polygon = MultiPolygon([final_mgrs_polygon])

                        #     print('selected_mgrs : ', final_mgrs_list)
                        #     print('jpl burst id', group_frames)
                        #     print('number of bursts', len(group_frames))

                        #     # IMPORTANT: unique id per subgroup
                        #     # Example: MS_12_34_1, MS_12_34_2, ...
                        #     # (track, mgrs_set_number) identifies the "base", and subgroup splits append _{group_idx}
                        #     mgrs_set_id = f"MS_{track_number}_{mgrs_set_number}"
                        #     print('mgrs_set_id', mgrs_set_id)
                        assigned_mgrs_by_group = assign_mgrs_tiles_to_frame_groups(
                            frame_groups=frame_groups,
                            target_polygon=target_polygon,
                            target_mgrs_list=target_mgrs_list,
                            group_list=group_list,
                            unique_group_number=unique_group_number,
                        )

                        for group_idx, fg in enumerate(frame_groups, start=1):
                            group_frames = fg["frames"]
                            group_tfs = fg["tfs"]
                            group_polys = fg["polys"]

                            final_mgrs_list = assigned_mgrs_by_group[group_idx - 1]["mgrs_list"]
                            final_mgrs_polygon_list = assigned_mgrs_by_group[group_idx - 1]["mgrs_polys"]

                            if not final_mgrs_list:
                                continue

                            final_mgrs_polygon = unary_union(final_mgrs_polygon_list).simplify(0.05)

                            if final_mgrs_polygon.geom_type == "Polygon":
                                final_mgrs_polygon = MultiPolygon([final_mgrs_polygon])

                            print("selected_mgrs:", final_mgrs_list)
                            print("jpl burst id:", group_frames)
                            print("number of bursts:", len(group_frames))

                            # Important: because you may now have multiple subgroups,
                            # include group_idx in the ID.
                            mgrs_set_id = f"MS_{track_number}_{mgrs_set_number}"
                            print("mgrs_set_id", mgrs_set_id)
                            # ---- geometry precision / WKT ----
                            def iter_polygons(geom):
                                if geom is None:
                                    return []
                                if geom.geom_type == "Polygon":
                                    return [geom]
                                if geom.geom_type == "MultiPolygon":
                                    return list(geom.geoms)
                                raise TypeError(f"Unsupported geometry type: {geom.geom_type}")

                            polys = iter_polygons(final_mgrs_polygon)
                            newfinal_mgrs_polygons = [transform_to_lower_precision(p) for p in polys]
                            new_multipolygonz = MultiPolygon(newfinal_mgrs_polygons)

                            geom = new_multipolygonz.wkt  # stored as WKT (2D or 3D depending on inputs)

                            bursts = str(group_frames)
                            number_of_bursts = len(group_frames)
                            mgrs_tiles = str(final_mgrs_list)
                            number_of_mgrs_tiles = len(final_mgrs_list)
                            relative_orbit_number = track_number

                            # ---- EPSG majority selection (same as your logic) ----
                            if burst_hemisphere == 'north':
                                ppp = 32600
                            else:
                                ppp = 32700
                            epsgcodes_list = [int(mgrs_cand[:2]) + ppp for mgrs_cand in final_mgrs_list]
                            major_epsg, _ = Counter(epsgcodes_list).most_common()[0]

                            lon_min_ext, lat_min_ext, lon_max_ext, lat_max_ext = final_mgrs_polygon.bounds
                            # x_min_ext, y_min_ext = latlon_to_utm(lat_min_ext, lon_min_ext, major_epsg)
                            # x_max_ext, y_max_ext = latlon_to_utm(lat_max_ext, lon_max_ext, major_epsg)
                            x_min_ext, y_min_ext, x_max_ext, y_max_ext = (
                                projected_geometry_bounds_antimeridian_safe(
                                    final_mgrs_polygon,
                                    major_epsg,
                                )
                            )
                            unique_epsg_codes = sorted(set(epsgcodes_list))
                            cross_utm_zone = len(unique_epsg_codes) > 1
                            print(
                                "Projected extent:",
                                {
                                    "major_epsg": major_epsg,
                                    "source_epsgs": unique_epsg_codes,
                                    "cross_utm_zone": cross_utm_zone,
                                    "xmin": x_min_ext,
                                    "xmax": x_max_ext,
                                    "ymin": y_min_ext,
                                    "ymax": y_max_ext,
                                    "width_m": x_max_ext - x_min_ext,
                                    "height_m": y_max_ext - y_min_ext,
                                }
                            )
                            # if x_max_ext < x_min_ext:
                            #     x_min_ext, x_max_ext = x_max_ext, x_min_ext
                            # if y_max_ext < y_min_ext:
                            #     y_min_ext, y_max_ext = y_max_ext, y_min_ext

                            # ---- land/ocean flag (keep your existing logic, but use group_frames) ----
                            land_ocean_mask = None
                            if land_ocean_flag:
                                land_ocean_flag_list = []
                                for frame_cand in group_frames:
                                    filtered_gdf2 = tfdb_gdf[(tfdb_gdf['track'] == track_number) &
                                                            (tfdb_gdf['passDirection'] == orbit_pass) &
                                                            (tfdb_gdf['frame'] == frame_cand)]
                                    for _, filt_gdf_cand2 in filtered_gdf2.iterrows():
                                        land_ocean_flag_list.append(filt_gdf_cand2['fractionLand'])

                                land_ocean_sum = np.sum(land_ocean_flag_list) / len(land_ocean_flag_list)
                                if land_ocean_sum == 1:
                                    land_ocean_mask = 'land'
                                elif land_ocean_sum == 0:
                                    land_ocean_mask = 'water'
                                else:
                                    land_ocean_mask = 'water/land'
                                print(f"Land Ocean Mask {land_ocean_mask}")

                            mgrs_db_cursor.execute(
                                """
                                INSERT INTO mgrs_track_frame_db
                                    (mgrs_set_id,
                                    geometry,
                                    frames,
                                    track_frame,
                                    number_of_frames,
                                    mgrs_tiles,
                                    number_of_mgrs_tiles,
                                    track_number,
                                    orbit_pass,
                                    land_ocean_flag,
                                    EPSG,
                                    xmin,
                                    xmax,
                                    ymin,
                                    ymax)
                                VALUES (?, GeomFromText(?, 4326), ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                                """,
                                (mgrs_set_id,
                                geom,
                                bursts, str(group_tfs), number_of_bursts,
                                mgrs_tiles, number_of_mgrs_tiles,
                                relative_orbit_number,
                                orbit_pass, land_ocean_mask, major_epsg,
                                int(x_min_ext), int(x_max_ext),
                                int(y_min_ext), int(y_max_ext))
                            )
                            if len(frame_list) > 3:
                                reason = "More than 3 frames intersected the MGRS collection"
                                if do_split:
                                    reason += "; polar split was applied"
                                else:
                                    reason += "; no split was applied"

                                mgrs_db_cursor.execute(
                                    """
                                    INSERT INTO mgrs_frame_check_db
                                        (mgrs_set_id,
                                        track_number,
                                        orbit_pass,
                                        unique_group_number,
                                        original_number_of_frames,
                                        original_frames,
                                        original_track_frame,
                                        split_applied,
                                        reason)
                                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                                    """,
                                    (mgrs_set_id,
                                    track_number,
                                    orbit_pass,
                                    unique_group_number,
                                    len(frame_list),
                                    str(frame_list),
                                    str(tf_list),
                                    int(do_split),
                                    reason)
                                )
                            # one "mgrs_set_number" per base set (even if split into subgroups)
                            mgrs_set_number += 1

                        print(" ")
                        print(" ")


    # Commit the changes and close the cursor
    mgrs_db_conn.commit()
    mgrs_db_cursor.close()

def main():

    db_parser = _get_parser()

    args = db_parser.parse_args()

    run(args)

if __name__ == '__main__':

    main()
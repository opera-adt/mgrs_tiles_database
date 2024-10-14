import argparse
import os

import sqlite3
import numpy as np
import time
from collections import Counter
from shapely.geometry import Polygon, MultiPolygon
from shapely.wkt import loads as wkt_loads
from shapely.wkb import loads as wkb_loads
from shapely import wkb
from shapely.ops import unary_union
from pyproj import Transformer
from osgeo import ogr, osr, gdal
from rasterio.mask import mask as rasterio_mask
import rasterio


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


# determine the path to the world land GPKG file
# LAND_GPKG_FILE = os.path.join(WORKFLOW_SCRIPTS_DIR, 'data',
#                               'GSHHS_l_L1.shp.no_attrs.epsg3413_dissolved.gpkg')
import mgrs_tiles_database

# LAND_GPKG_FILE = "/mnt/aurora-r0/jungkyo/OPERA/burst_db/MGRS_db/GSHHS_l_L1.shp.no_attrs.epsg3413_dissolved.gpkg"
WORKFLOW_SCRIPTS_DIR = os.path.dirname(mgrs_tiles_database.__file__)

LAND_GPKG_FILE = os.path.join(WORKFLOW_SCRIPTS_DIR, 'data',
                              'GSHHS_l_L1.shp.no_attrs.epsg3413_dissolved.gpkg')
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
    ring_cslc = ogr.Geometry(ogr.wkbLinearRing)
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
    bbox_cslc.Transform(transformer_cslc_to_land)

    # Return a numpy array full of `False` when there is no intersection
    if not bbox_cslc.Intersects(land_polygon):
        return np.full(shape_mask, False)

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

    return mask_land


def _get_parser():
    parser = argparse.ArgumentParser(
        description='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-mgrs_db',
                        type=str,
                        dest='mgrs_db',
                        help='JPL MGRS database sqlite')

    parser.add_argument('-burst_id_db',
                        type=str,
                        dest='burst_id_db',
                        help='JPL burst_id_db sqlite')

    parser.add_argument('-land_ocean_shapefile',
                        type=str,
                        dest='land_ocean_shapefile',
                        help='Shapefile for land/ocean mask')

    parser.add_argument('-output',
                        type=str,
                        dest='output',
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

    pixel_num = landcover_data.size

    land_pixel = np.sum(landcover_data)
    print(land_pixel, pixel_num)
    if land_pixel == pixel_num:
        ocean_land_class = 'land'
    elif land_pixel == 0:
        ocean_land_class = 'water'
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
    transformer = Transformer.from_crs(f"EPSG:{utm_epsg}", "EPSG:4326")

    # Transform the UTM coordinates to latitude and longitude.
    lat, lon = transformer.transform(utm_easting, utm_northing)
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
    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}")

    # Transform the latitude and longitude to UTM coordinates.
    utm_x, utm_y = transformer.transform(lat, lon)

    return utm_x, utm_y

def transform_to_lower_precision(polygon, decimal_places=4):
    """
    Transform the coordinates of a polygon to a lower precision.

    Parameters:
    -----------
    polygon : shapely.geometry.Polygon
        Input polygon whose coordinates need to be rounded.
    decimal_places : int, optional
        Number of decimal places to round each coordinate component to. Default is 4.

    Returns:
    --------
    shapely.geometry.Polygon
        A new polygon with rounded coordinates.

    Note:
    -----
    This function assumes the polygon's coordinates are in 3D (i.e., x, y, z).
    For 2D polygons, the code needs modification.
    """
    # List comprehension to round each component of the coordinate
    coords = [(round(x, decimal_places),
               round(y, decimal_places),
               round(z, decimal_places))
              for x, y, z in list(polygon.exterior.coords)]

    return Polygon(coords)

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
        Orbit pass direction ('ASCENDING' or 'DESCENDING').

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
        orbit_pass == 'ASCENDING':
        if relative_orbit_number == 1:
            relative_orbit_number_second = 175
        elif relative_orbit_number > 175:
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


def run(args):
    '''
    1) Loop starts from track number 1 to track number 175
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

    # Connect to the JPL burst id database sqlite
    absolute_burst_sqlite_file = args.burst_id_db
    conn = sqlite3.connect(absolute_burst_sqlite_file)
    burst_id_cur = conn.cursor()

    query_result_all = burst_id_cur.execute('SELECT * FROM burst_id_map')
    burst_fields_name = {}

    # Read field names from database.
    for i, column_desc in enumerate(query_result_all.description):
        burst_fields_name[column_desc[0]] = i

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
    CREATE TABLE IF NOT EXISTS mgrs_burst_db (
        id INTEGER PRIMARY KEY,
        mgrs_set_id TEXT,
        relative_orbit_number INTEGER,
        bursts TEXT,
        number_of_bursts INTEGER,
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

    mgrs_db_cursor.execute(
        '''
        SELECT AddGeometryColumn('mgrs_burst_db', 'geometry', 4326, 'MULTIPOLYGON', 3)
        ''')
    mgrs_db_cursor.execute(
        "SELECT CreateSpatialIndex('mgrs_burst_db', 'geometry');")

    # Loop from track number 1 to track number 175
    for track_number in range(1, 176):

        mgrs_set_number = 1

        orbit_passes = ['ASCENDING', 'DESCENDING', 'DESCENDING', 'ASCENDING']
        group_lists = [group_list_north,
                       group_list_north,
                       group_list_south,
                       group_list_south]
        orbit_passes = ['ASCENDING']
        group_lists = [group_list_north]

        for orbit_pass, group_list in zip(orbit_passes, group_lists):
            # Search bursts with track number and orbit pass
            query = f"""
                SELECT *
                FROM burst_id_map
                WHERE relative_orbit_number = {track_number}
                    AND orbit_pass = '{orbit_pass}'
                """
            burst_id_cur.execute(query)
            burst_id_query_output = burst_id_cur.fetchall()

            # Find first valid bursts
            for jjjind in range(0, 100):
                row_start = burst_id_query_output[jjjind]
                geometry_valid = wkb_loads(row_start[burst_fields_name['GEOMETRY']])
                # Bursts in high latitude are not used for MGRS tile collection.
                if north_bound_polygon.contains(geometry_valid) or \
                   north_bound_polygon.intersects(geometry_valid):
                    print('skip', row_start[burst_fields_name['burst_id_jpl']])
                else:
                    break

            # Read burst id
            burst_id_start = row_start[burst_fields_name['burst_id']]
            burst_utm_zone = int(str(row_start[burst_fields_name['EPSG']])[-2:])

            # Unique MGRS tile id
            unique_group_list = list(set(group_list))

            # mgrs group id increase from south pole to north pole.
            # in descending orbit, the unique group list needs to be revsersed.
            if orbit_pass == 'DESCENDING':
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
                    burst_geometry = wkb_loads(row_start[burst_fields_name['GEOMETRY']])
                    burst_id_in_track = row_start[burst_fields_name['burst_id']]
                    str_subswath = row_start[burst_fields_name['subswath_name']]
                    burst_id_jpl = row_start[burst_fields_name['burst_id_jpl']]
                    orbit_pass = row_start[burst_fields_name['orbit_pass']]

                    burst_center_lon, burst_center_lat = burst_geometry.centroid.x, burst_geometry.centroid.y

                    burst_hemisphere = 'north' if burst_center_lat >= 0  else 'south'
                    epsg_code = 32600 + burst_utm_zone if burst_hemisphere == 'north' else 32700 + burst_utm_zone

                    intersect_x, intersect_y = latlon_to_utm(burst_center_lat,
                                                             burst_center_lon,
                                                             epsg_code)

                    # Intentionally add buffer to get adjacent bursts
                    if orbit_pass == 'ASCENDING':
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

                    if (unique_group_number < 30) and (
                        unique_group_number > -30):
                        # By finding intersecting areas between polygon computed from bursts and
                        # unique MGRS group, extract the tiles MGRS within bounds in x direction.
                        # y bounds are computed from unique MGRS tile group and x bounds are determined
                        # from bursts.
                        mgrs_poly = mgrs_poly.intersection(burst_multi_polygon)

                    # Find intersecting bursts from all bursts along the specific track number
                    intersecting_burst, burst_ind_iw1, _, _, burst_jpl_id = \
                            read_polygons_from_jpl_burst_id_db(
                                absolute_burst_sqlite_file, 'burst_id_map',
                                mgrs_poly,
                                relative_orbit_number=track_number,
                                orbit_pass=orbit_pass)


                    if len(burst_jpl_id) == 0:
                        print('here: burst_jpl_id is empty')
                        intersecting_burst, burst_ind_iw1, _, _, burst_jpl_id = \
                            read_polygons_from_jpl_burst_id_db(
                                absolute_burst_sqlite_file, 'burst_id_map',
                                mgrs_poly0,
                                relative_orbit_number=track_number,
                                orbit_pass=orbit_pass)

                    union_burst_polygon = unary_union(intersecting_burst)

                    if unique_group_number > 0:
                        target_polygon = north_polygon
                        target_mgrs_list = north_mgrs_list
                    else:
                        target_polygon = south_polygon
                        target_mgrs_list = south_mgrs_list

                    final_mgrs_list = []
                    final_mgrs_polygon_list = []

                    # Find intersecting MGRS tiles with burst
                    for group_ind, group_number in enumerate(group_list):
                        if group_number == unique_group_number:

                            single_mgrs = wkt_loads(target_polygon[group_ind].ExportToWkt())
                            # Excludes the MGRS tiles not overlapped with the polygon of burts
                            if union_burst_polygon.intersects(single_mgrs) or \
                                union_burst_polygon.contains(single_mgrs):

                                final_mgrs_list.append(target_mgrs_list[group_ind])
                                final_mgrs_polygon_list.append(single_mgrs)

                    final_mgrs_polygon = unary_union(final_mgrs_polygon_list).simplify(0.05)

                    if final_mgrs_polygon.geom_type == 'Polygon':
                        final_mgrs_polygon = MultiPolygon([final_mgrs_polygon])

                    print('selected_mgrs : ', final_mgrs_list)
                    print('jpl burst id', burst_jpl_id)
                    print('number of bursts', len(burst_jpl_id))


                    mgrs_set_id = f'MS_{track_number}_{mgrs_set_number}'

                    # transform coordinates to lower precision coordinates
                    newfinal_mgrs_polygons = [transform_to_lower_precision(polygon_sample) for polygon_sample in final_mgrs_polygon]
                    new_multipolygonz = MultiPolygon(newfinal_mgrs_polygons)

                    geom = new_multipolygonz.wkt
                    geom_wkb = wkb.dumps(new_multipolygonz)
                    bursts = str(burst_jpl_id)
                    number_of_bursts = len(burst_jpl_id)
                    mgrs_tiles = str(final_mgrs_list)
                    number_of_mgrs_tiles = len(final_mgrs_list)
                    relative_orbit_number = track_number

                    if burst_hemisphere == 'north':
                        ppp = 32600
                    else:
                        ppp = 32700
                    epsgcodes_list = [ int(mgrs_cand[:2]) + ppp
                                    for mgrs_cand in final_mgrs_list]

                    # Find majority from epsg list and assign it to MGRS collection
                    major_epsg, _ = Counter(epsgcodes_list).most_common()[0]

                    lon_min_ext, lat_min_ext, lon_max_ext, lat_max_ext = \
                        final_mgrs_polygon.bounds
                    x_min_ext, y_min_ext = latlon_to_utm(lat_min_ext,
                                                        lon_min_ext,
                                                        major_epsg)
                    x_max_ext, y_max_ext = latlon_to_utm(lat_max_ext,
                                                        lon_max_ext,
                                                        major_epsg)
                    if x_max_ext < x_min_ext:
                        x_min_ext, x_max_ext = x_max_ext, x_min_ext
                    if y_max_ext < y_min_ext:
                        y_min_ext, y_max_ext = y_max_ext, y_min_ext
                    # Flag for land/ocean
                    land_ocean_mask = 'N/A'
                    if land_ocean_flag:
                        # Read GeoTIFF file

                        # Mask raster data with each polygon and compute mean
                        # out_image, out_transform = rasterio_mask(landcover_src, final_mgrs_polygon, crop=True)
                        x_spacing = 5000
                        y_spacing = 5000
                        height = int((y_max_ext - y_min_ext) / y_spacing)
                        width = int((x_max_ext - x_min_ext) / x_spacing)
                        if (height == 0) or (width == 0):
                            x_spacing = 100
                            y_spacing = 100
                            height = int((y_max_ext - y_min_ext) / y_spacing)
                            width = int((x_max_ext - x_min_ext) / x_spacing)                            
                        print(x_max_ext, x_min_ext)
                        print(height, width)
                        mask_land = _get_land_mask(major_epsg,
                            (x_min_ext, x_spacing, 0, y_min_ext, 0, y_spacing),
                            (height, width))
                        land_ocean_mask = check_flag_land_ocean_landcover(mask_land)
                        print(f"Land Ocean Mask {land_ocean_mask}")
    
                    mgrs_db_cursor.execute(
                        """
                        INSERT INTO mgrs_burst_db
                                        (mgrs_set_id,
                                            geometry,
                                            bursts,
                                            number_of_bursts,
                                            mgrs_tiles,
                                            number_of_mgrs_tiles,
                                            relative_orbit_number,
                                            orbit_pass,
                                            land_ocean_flag,
                                            EPSG,
                                            xmin,
                                            xmax,
                                            ymin,
                                            ymax
                                        ) VALUES (?, GeomFromText(?, 4326), ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (mgrs_set_id, geom,
                        bursts, number_of_bursts,
                        mgrs_tiles, number_of_mgrs_tiles,
                        relative_orbit_number,
                        orbit_pass, land_ocean_mask, major_epsg,
                        int(x_min_ext), int(x_max_ext), int(y_min_ext), int(y_max_ext)
                        ))

                    mgrs_set_number += 1

                    # Update the valid first burst for the next loop
                    for row_sample in burst_id_query_output:
                        burst_id_in_track = row_sample[burst_fields_name['burst_id']]
                        str_subswath = row_sample[burst_fields_name['subswath_name']]
                        if (burst_id_in_track == np.max(burst_ind_iw1) + 1) and str_subswath == 'IW1':
                            burst_utm_zone = int(final_mgrs_list[-1][:2])
                            row_start = row_sample
                            break

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
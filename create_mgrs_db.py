import argparse
import os

from bs4 import BeautifulSoup
from pykml import parser
from osgeo import ogr, osr
from shapely import wkt

def _get_parser():
    parser = argparse.ArgumentParser(
        description='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-kml_path',
                        type=str,
                        dest='kml_path',
                        help='ESA MGRS kml file')

    parser.add_argument('-output',
                        type=str,
                        dest='output',
                        help='prefix output')

    return parser


def find_polygons(element):
    '''
    Recursively search for elements with tags ending in 'Polygon'.

    This function yields elements from an XML-like structure whose tags
    end with 'Polygon'. It traverses the structure in a depth-first manner.

    Parameters
    ----------
    element : xml.etree.ElementTree.Element
        The starting XML element to search from.

    Yields
    ------
    xml.etree.ElementTree.Element
        Elements with tags ending in 'Polygon'.
    '''
    if element.tag.endswith('Polygon'):
        yield element
    for child in element.getchildren():
        yield from find_polygons(child)

def extract_coordinates(polygon):
    '''
    Extracts the coordinates of the outer boundary of a KML polygon.

    This function processes the given KML polygon element, retrieving
    the coordinates of its outer boundary, and returning them as a list
    of tuples.

    Parameters
    ----------
    polygon : xml.etree.ElementTree.Element
        The KML polygon element to extract coordinates from.

    Returns
    -------
    list[tuple[float, float]]
        A list of tuples, where each tuple represents a coordinate pair
        (latitude, longitude).
    '''
    coordinates_text = polygon.outerBoundaryIs.LinearRing.coordinates.text.strip()
    coordinates = [tuple(map(float, coord.split(',')))
                   for coord in coordinates_text.split()]
    return coordinates

def parser_description_html(html):
    '''
    Parses an HTML table to extract values associated with "EPSG"
    and "UTM_WKT".

    The function expects an HTML content with a table containing rows.
    Each row should contain two cells where the first cell might contain
    either "EPSG" or "UTM_WKT" and the second cell contains the corresponding
    value.

    Parameters
    ----------
    html : str
        The HTML content to be parsed.

    Returns
    -------
    tuple
        A tuple containing two values:
            - epsg (str or None): The value associated with "EPSG" in the
              table, or None if not found.
            - utm_wkt (str or None): The value associated with "UTM_WKT"
              in the table, or None if not found.
    '''
    soup = BeautifulSoup(html, 'html.parser')

    # find all rows in the table
    rows = soup.find_all('tr')

    epsg = None
    utm_wkt = None
    for row in rows:
        # find all cells in the row
        cells = row.find_all('td')
        if len(cells) == 2:
            # check if the first cell contains "EPSG"
            if "EPSG" in cells[0].text:
                epsg = cells[1].text.strip()
            # check if the first cell contains "UTM_WKT"
            elif "UTM_WKT" in cells[0].text:
                utm_wkt = cells[1].text.strip()

    return epsg, utm_wkt


def run(args):
    # ESA MGRS tile input will be used to assign the same extents and
    # get necessary tiles.
    # kml_path =
    # 'S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml'

    # Set the paths for input KML and output filenames
    kml_path = args.kml_path
    output_name = args.output

    # Create output data sources for SQLite and Shapefile
    sqlite_driver = ogr.GetDriverByName("SQLite")
    shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")

    output_sqlite_file = f'{output_name}.sqlite'
    output_shapefile = f'{output_name}.shp'

    # Remove existing output files to prevent appending data
    if os.path.isfile(output_sqlite_file):
        os.remove(output_sqlite_file)
    if os.path.isfile(output_shapefile):
        os.remove(output_shapefile)

    # Initialize output data sources
    sqlite_ds = sqlite_driver.CreateDataSource(output_sqlite_file)
    shapefile_ds = shapefile_driver.CreateDataSource(output_shapefile)

    # Set spatial reference (WGS 84)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    # Create output layers with specified options
    sqlite_layer = sqlite_ds.CreateLayer(output_name,
                                         srs,
                                         ogr.wkbPolygon,
                                         options=["SPATIALITE=YES"])

    shapefile_layer = shapefile_ds.CreateLayer(output_name,
                                               srs,
                                               ogr.wkbPolygon)

    # Define fields for layers
    field_names = [output_name,
                   "EPSG",
                   "xmin",
                   "xmax",
                   "ymin",
                   "ymax"]
    field_types = [ogr.OFTString,
                   ogr.OFTInteger64,
                   ogr.OFTInteger64,
                   ogr.OFTInteger64,
                   ogr.OFTInteger64,
                   ogr.OFTInteger64]

    for fname, ftype in zip(field_names, field_types):
        field_defn = ogr.FieldDefn(fname, ftype)
        sqlite_layer.CreateField(field_defn)
        if fname == output_name:  # Only add the name field to the shapefile layer
            shapefile_layer.CreateField(field_defn)

    # Open, read and parse the KML file
    with open(kml_path, 'r') as file:
        kml_content = file.read().encode('utf-8')

    root = parser.fromstring(kml_content)

    placemark_names = []
    mgrs_polygons = []
    epsg_list = []

    # Iterate through all Placemark elements in the KML file
    num_all = len(root.Document.Folder.Placemark)
    for aind, placemark in enumerate(root.Document.Folder.Placemark):
        # Extract necessary details from the placemark
        placemark_name = placemark.name.text
        description_html = placemark.description.text
        epsg, polygon_mgrs = parser_description_html(description_html)
        geom = wkt.loads(polygon_mgrs)

        print(f'processing: {aind + 1:,} / {num_all:,}', end='\r')
        # Get the bounds of the geometry
        xmin, ymin, xmax, ymax = geom.bounds

        placemark_name = placemark.name.text
        placemark_names.append(placemark_name)

        description_html = placemark.description.text
        epsg, polygon_mgrs = parser_description_html(description_html)

        mgrs_polygons.append(polygon_mgrs)
        epsg_list.append(epsg)
        print(f'processing: {aind + 1:,} / {num_all:,}', end='\r')

        geom = wkt.loads(polygon_mgrs)

        # Get the bounds of the geometry
        xmin, ymin, xmax, ymax = geom.bounds

        for polygon_cand in placemark.MultiGeometry.Polygon:
            coordinates = extract_coordinates(polygon_cand)

            # Construct a polygon geometry from coordinates
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for coord in coordinates:
                ring.AddPoint(coord[0], coord[1])
            ring.CloseRings()
            polygon = ogr.Geometry(ogr.wkbPolygon)
            polygon.AddGeometry(ring)

            # Create a feature, set its properties and save it to the layers
            feature = ogr.Feature(sqlite_layer.GetLayerDefn())
            feature.SetGeometry(polygon)
            feature.SetField(output_name, placemark_name)
            feature.SetField("EPSG", epsg)
            feature.SetField("xmin", xmin)
            feature.SetField("xmax", xmax)
            feature.SetField("ymin", ymin)
            feature.SetField("ymax", ymax)
            sqlite_layer.CreateFeature(feature)
            shapefile_layer.CreateFeature(feature)

    # Close the data sources
    sqlite_ds.Destroy()
    shapefile_ds.Destroy()

def main():

    db_parser = _get_parser()

    args = db_parser.parse_args()

    run(args)

if __name__ == '__main__':

    main()

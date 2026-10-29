# issue was found due to the rasterization.
import sqlite3

from mgrs_tiles_database.create_mgrs_collection_db import (_get_land_mask,
                                                           check_flag_land_ocean_landcover)

def test(major_epsg,
         x_min_ext,
         x_max_ext,
         y_min_ext,
         y_max_ext):
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

    return land_ocean_mask


mgrs_collection_db_path = \
    '/mnt/aurora-r0/jungkyo/OPERA/burst_db/MGRS_db/MGRS_tile_collection_v0.2.sqlite'
mgrs_db_conn = sqlite3.connect(mgrs_collection_db_path)
mgrs_db_conn.enable_load_extension(True)

mgrs_db_cursor = mgrs_db_conn.cursor()

query = f"""
    SELECT *
    FROM mgrs_burst_db
    WHERE land_ocean_flag = 'land'
    """
query_result_all = mgrs_db_cursor.execute(query)
mgrs_col_db_name = {}
for i, column_desc in enumerate(query_result_all.description):
    mgrs_col_db_name[column_desc[0]] = i

rows_to_update = mgrs_db_cursor.execute(query).fetchall()


for row in rows_to_update:
    row_id = row[0]
    epsg = row[mgrs_col_db_name['EPSG']]
    xmin = row[mgrs_col_db_name['xmin']]
    xmax = row[mgrs_col_db_name['xmax']]
    ymin = row[mgrs_col_db_name['ymin']]
    ymax  = row[mgrs_col_db_name['ymax']]
    print(row_id)
    print(epsg, xmin, xmax, ymin, ymax)
    # Use the function to get the new value
    new_value = test(epsg, 
                     xmin, 
                     xmax, 
                     ymin,
                     ymax)
    print(new_value)
    # Update the record in the database
    update_query = f"""
        UPDATE mgrs_burst_db
        SET land_ocean_flag = ?
        WHERE id = ?
    """
    mgrs_db_cursor.execute(update_query, (new_value, row_id))
# Commit the changes
mgrs_db_conn.commit()

# Close the database connection
mgrs_db_conn.close()





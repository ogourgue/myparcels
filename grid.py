from parcels import Field

######################################################
# Convert "cell areas" from field property to field. #
######################################################

def cell_areas(field):
    name = 'cell_areas'
    data = field.cell_areas()
    lon = field.grid.lon
    lat = field.grid.lat
    return Field(name = name, data = data, lon = lon, lat = lat)

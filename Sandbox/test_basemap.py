# import urllib.request
# import json
# import matplotlib.pyplot as plt
# from matplotlib import animation
# import cartopy.crs as ccrs
# from cartopy.io.img_tiles import OSM

import cartopy.crs as ccrs
from cartopy.io.img_tiles import OSM
import matplotlib.pyplot as plt


osm_tiles = OSM()

plt.figure(figsize=(16, 16))

# Use the tile's projection for the underlying map.
ax = plt.axes(projection=osm_tiles.crs)

# Specify a region of interest, in this case, Cardiff.

ax.set_extent([-122.335058, -122.222542, 47.426043, 47.472231],
              ccrs.PlateCarree())

# Add the tiles at zoom level 12.
ax.add_image(osm_tiles, 12, interpolation='spline36')

ax.coastlines('10m')

plt.show()


# #DEFINE FIGURE
# fig, ax = plt.subplots()
#
# #SET AXES FOR PLOTTING AREA
#
# # lat_top, lon_top = 47.533953, -122.393251
# # lat_bot, lon_bot = 47.373107, -122.227104
# lat_top, lon_top = 47.472231, -122.335058
# lat_bot, lon_bot = 47.426043, -122.222542
# ax=plt.axes(projection=ccrs.PlateCarree())
# ax.set_ylim(lat_bot, lat_top)
# ax.set_xlim(lon_top, lon_bot)
#
# #ADD OSM BASEMAP
# osm_tiles=OSM()
# ax.add_image(osm_tiles,14, interpolation='spline36', regrid_shape=2000) #Zoom level 13
#
# #PLOT JFK INTL AIRPORT
# # ax.text(-73.778889,40.639722,'JFK Intl',horizontalalignment='right',size='large')
# # ax.plot([-73.778889],[40.639722],'bo') #Plot a point in blue color
#

# plt.show()


# fig = plt.figure(figsize=(10, 10))
#
# imagery = OSM()
#
# ax = plt.axes(projection=imagery.crs, )
# # ax.set_extent((lat_bot, lat_top, lon_top, lon_bot))
# ax.set_extent((47.426043, 47.472231, -122.335058, -122.222542))
# # ax.set_extent(( 153, 153.2, -26.6, -26.4))
#
# # Add the imagery to the map.
# zoom = 12
# ax.add_image(imagery, zoom, interpolation='none')
#
# plt.title('Open Street Map and Cartopy')
# plt.show()
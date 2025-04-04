import earthaccess  
import h5py  
import numpy as np  
import matplotlib.pyplot as plt  
import os
from netCDF4 import Dataset
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean

#Authenticate
earthaccess.login()

#Searching for Cryosat2 Sea Ice Thickness Data for 2024
#Returns a list of DataGranule objects
search_results = earthaccess.search_data(
    short_name="RDEFT4",
    bounding_box=(-163.5, 66.4, -161.1, 67.1), #Zoom in around Kotzebue sound
    temporal=("2023-09-01", "2024-05-31"),
)

#Returns a list of downloaded files
raw_files = earthaccess.download(search_results[0:5], local_path='data/')

#Unpack netCDF files so that we can extract the features we want
data_list = []
for file in os.listdir('data/'):
    file_path = os.path.join('data', file)
    data_list.append(Dataset(file_path))
    print(data_list[0])

thickness_list = []
coords_list = []
for d in range(len(data_list)):
    thickness_list.append(data_list[d].variables['sea_ice_thickness'][:,:].flatten().T)
    coords_list_entry = (data_list[d].variables['lat'][:,:].flatten().T,
                  data_list[d].variables['lon'][:,:].flatten().T)
    coords_list.append(np.vstack(coords_list_entry).T)

ice_thickness = np.hstack(thickness_list)
coords = np.vstack(coords_list)

print(ice_thickness.shape)
print(coords.shape)

df_actual = gpd.GeoDataFrame({"thickness": ice_thickness, "geometry": gpd.points_from_xy(coords[:,1], coords[:,0])})

# PLOTTING #
fig, axes = plt.subplots(1, 1, figsize=(10, 8), constrained_layout=True, squeeze=False)

x_min, x_max = df_actual.total_bounds[0], df_actual.total_bounds[2]
y_min, y_max = df_actual.total_bounds[1], df_actual.total_bounds[3]
zoom_x = 0.9  
zoom_y = 0.8

x_range = (x_max - x_min) * (1 - zoom_x) / 2
y_range = (y_max - y_min) * (1 - zoom_y) / 2
im1 = df_actual.plot(column="thickness", cmap=cmocean.cm.ice, markersize=10, alpha=0.7, legend=False, ax=axes[0])
axes[0].set_title("Ground Truth Ice Thickness")
axes[0].set_xlabel("Longitude")
axes[0].set_ylabel("Latitude")
axes[0].set_xlim(x_min + x_range, x_max - x_range)
axes[0].set_ylim(y_min + y_range, y_max - y_range)

#Add a single colorbar for both subplots
sm = plt.cm.ScalarMappable(cmap=cmocean.cm.ice, norm=norm)
sm._A = [] 
cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.02, pad=0.02)
cbar.set_label("Sea Ice Thickness (m)")

plt.show()
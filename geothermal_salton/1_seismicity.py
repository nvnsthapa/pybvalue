

# %%
import pandas as pd
import sys
import os
sys.path.append('../')
from src import myUtils
from src import plotUtils
#Parent directory
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#Project Directories
data_dir = os.path.join(parent_dir, 'geothermal_salton', 'data')
plot_dir = os.path.join(parent_dir, 'geothermal_salton', 'plots')
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
# %% Read relocated catalog
relocatedCatalogF = os.path.join(data_dir, 'gc_only.gc')
relocatedCat = myUtils.read_cat_relocated(relocatedCatalogF)
print(relocatedCat.keys())
# %% Plot seismicity
# rename column name of relocatedCat to match the plot_seismicity function 
relocatedCat.rename(columns={'lonR': 'longitude', 'latR': 'latitude', 'depR': 'depth', 'mag': 'magnitdue'}, inplace=True)
plotUtils.plot_seismicity(relocatedCat)
# %% 
# Slice catalog data for salton Trough and Imperial Valley area
salton_imperial = relocatedCat[
    (relocatedCat.latitude >= 32.084) & (relocatedCat.latitude <= 33.94) &
    (relocatedCat.longitude >= -117.16) & (relocatedCat.longitude <= -114.376)
]

print(f"Number of Events in Salton Trough: {len(salton_imperial)}")

# %%
plotUtils.plot_seismicity(salton_imperial)
# %%
# Slice catalog data for Salton Sea area only
salton_sea = relocatedCat[
    (relocatedCat.latitude >= 32.99) & (relocatedCat.latitude <= 33.39) &
    (relocatedCat.longitude >= -115.99) & (relocatedCat.longitude <= -115.39)
]
print(f"Number of Events in Salton Sea: {len(salton_sea)}")


# %%
# read geothermal well data from xlsx file
well_data = os.path.join(data_dir, 'Brawely_wellData.xlsx')
well_df = pd.read_excel(well_data, engine='openpyxl')
# %%
# %%
saltonFig = plotUtils.plot_seismicity(salton_sea, show_fig=False)
# Plot well data over saltonFig
saltonFig.plot(
    x=well_df['Longitude'],
    y=well_df['Latitude'],
    style='i0.3s',
    color='red',
    label='Well Location',
)
saltonFig.show()

# This is plotUtils modules for plotting seismic data 

import pygmt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator


def plot_seismicity(catalog_df):
    """
    Plot seismicity using GMT
    Args:
        catalog_df (pd.DataFrame): Catalog of seismicity, with columns longitude, latitude, depth, and magnitude
        filename (str): Filename for saving the plot
    Returns:
        Seismicity plot
    """
    # Define Region for plotting from data
    region = [
        catalog_df.longitude.min(), 
        catalog_df.longitude.max(),
        catalog_df.latitude.min(),
        catalog_df.latitude.max(),
    ]

    topo = "@earth_relief_01s"
    grid_map = pygmt.datasets.load_earth_relief(
        resolution="01s",
        region=region,
    )

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain")
    pygmt.config(FORMAT_GEO_MAP="ddd.xx")
    pygmt.config(MAP_ANNOT_OFFSET="0.2c") 
    pygmt.config(MAP_TICK_LENGTH="0.4c") 
    pygmt.config(MAP_TICK_PEN="0.07c") 

    fig.basemap(region=region, projection="M15i", frame=["af", "WSne"])
    land = fig.grdimage(grid=grid_map, cmap="grayC", shading=True)
    # land = fig.grdimage(grid=grid_map, cmap="grayC", shading=True)
    fig.coast(land=land, borders=["1/1p,black"], water='skyblue', shorelines=True, resolution='f')
    pygmt.makecpt(cmap="seis", series=[0, 30])

    fig.plot(
        x=catalog_df.longitude,
        y=catalog_df.latitude,
        # size=0.02 * (2**catalog_df.magnitude),
        fill=catalog_df.depth,
        cmap=True,
        # style="cc",
        style="c0.1c",
        # pen="red",
        transparency=50,
    )
    fig.colorbar(frame="af+lDepth (km)")
    # Coso geothermal field lat, lon range
    latmin = 35.90
    latmax = 36.17
    lonmin = -118.00
    lonmax = -117.67
    #plot closed polygon around the geothermal field
    fig.plot(
        x=[lonmin, lonmax, lonmax, lonmin, lonmin],
        y=[latmin, latmin, latmax, latmax, latmin],
        pen="1p,black",)
        # color="black",

    # Legend in the upper right
    # fig.legend(spec='quake_mag_sym.txt', position="jTR+o0.1c")
    # fig.savefig(f"Sesimicity_{filename}.pdf", dpi=600)
    
    return fig


def plot_welldata(catalog_df, well_df, faultFile=None):
    """
    Plot seismicity using GMT
    Args:
        catalog_df (pd.DataFrame): Catalog of seismicity, with columns longitude, latitude, depth, and magnitude
        well_df (pd.DataFrame): Well data with columns longitude, latitude, depth, and temperature
    Returns:
        Seismicity plot
    """
    # Define Region for plotting from data
    region = [
        catalog_df.longitude.min(), 
        catalog_df.longitude.max(),
        catalog_df.latitude.min(),
        catalog_df.latitude.max(),
    ]

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain")
    pygmt.config(FORMAT_GEO_MAP="ddd.xx")
    pygmt.config(MAP_ANNOT_OFFSET="0.2c") 
    pygmt.config(MAP_TICK_LENGTH="0.4c") 
    pygmt.config(MAP_TICK_PEN="0.07c") 

    topo = "@earth_relief_01s"
    grid_map = pygmt.datasets.load_earth_relief(
        resolution="01s",
        region=region,
    )

    fig.basemap(region=region, projection="M6i", frame=["af", "WSne"])
    # land = fig.grdimage(grid=grid_map, cmap="grayC", shading=True)
    land = fig.grdimage(grid=grid_map, cmap="grayC", shading=True)
    fig.coast(region = land, borders=["1/1p,black"], water='skyblue', shorelines=True, resolution='f')
    pygmt.makecpt(cmap="seis", series=[catalog_df.depth.min(), 20])

    fig.plot(
        x=catalog_df.longitude,
        y=catalog_df.latitude,
        # size=0.02 * (2**catalog_df.magnitude),
        fill=catalog_df.depth,
        cmap=True,
        # style="cc",
        style="c0.1c",
        # pen="red",
        transparency=50,
    )
    #plot earthquake that are greater than 3.0 magnitude and less than 4.0 magnitude
    # filtered_eq = catalog_df[catalog_df.mag >= 3.5]
    # filtered_eq = filtered_eq[filtered_eq.mag <= 4.0]
    # filtered_eq = filtered_eq[filtered_eq.time <= '2000-01-01']
    # datefilter = catalog_df[catalog_df.mag >= 3.5]
#     print(filtered_eq)


#     fig.plot(x=filtered_eq.longitude, 
#              y=filtered_eq.latitude, 
#             # fill="red",
#              style="c0.3c", 
#              pen="1p,black")
    
#     fig.plot(x=filtered_eq.longitude, 
#              y=filtered_eq.latitude, 
#             fill="red",
#              style="a0.3c", 
#              pen="red")

    inj_well = well_df[well_df['welltype'] == 'I']
    production_well = well_df[well_df['welltype'] == 'P']

    fig.plot(
        x=inj_well.Longitude,
        y=inj_well.Latitude,
        fill="blue",
        style="s0.4s",
        pen="black",
    )

    fig.plot(
    x=production_well.Longitude,
    y=production_well.Latitude,
    fill="green",
    style="s0.4s",
    pen="black",
)
    
    #plot fault data
    fig.plot(
        data = faultFile,

    )
#     fig.plot(
#     x=-113.010829,
#     y=38.3966,
#     fill="white",
#     style="s0.30s",
#     pen="black",
#     label="Milford Town",
# )
    
#     fig.plot(
#     x=-112.853284,
#     y=38.488842,
#     fill="red",
#     style="s0.30s",
#     pen="black",
#     label="Blundell Geothermal Plant",
# )

    # fig.colorbar(frame="af+lDepth (km)")
    # Legend in the upper right
    # fig.legend(spec='quake_mag_sym.txt', position="jTR+o0.1c")
    # fig.savefig(f"Sesimicity_{filename}.pdf", dpi=600)
    
    
    return fig



def plotTemporalSeismicity(mineral_mountains):
    # plt.close('all')
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = 18
    font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 18}

    minorTick = {'which': 'minor', 'direction': 'out', 'length': 4, 'width': 1}
    majorTick = {'which': 'major', 'direction': 'out', 'length': 7, 'width': 3}

    fig = plt.figure(figsize=(14, 6))
    ax = plt.subplot()

    ax.plot(mineral_mountains['time'], mineral_mountains['cumcount'], color='r', linestyle='-', linewidth=2)

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.set_xlabel('Time (Year-Month-Day)')
    ax.set_ylabel('Cumulative Number of Earthquakes')

    Loct_y1 = plt.MultipleLocator(250)
    ax.yaxis.set_minor_locator(Loct_y1)
    ax.tick_params(axis='x', top=True,  **majorTick) 
    ax.tick_params(axis='x', top=True, **minorTick)  
    ax.tick_params(axis='y', right=True, **majorTick)
    ax.tick_params(axis='y', right=True, **minorTick) 

    ax1 = ax.twinx()
    filtered_eq = mineral_mountains[mineral_mountains.mag >= 3.5]

    # ax1.stem(mineral_mountains['time'], mineral_mountains['mag'], linefmt='gray', markerfmt='o', basefmt=' ', bottom=0)
    ax1.scatter(mineral_mountains['time'], mineral_mountains['mag'], color='none', edgecolor= 'k', marker='o', s=10)
    ax1.scatter(filtered_eq['time'], filtered_eq['mag'], color='r', edgecolor= 'r', marker='*', s=50)
    Loct_x1 = plt.MultipleLocator(0.2)
    ax1.yaxis.set_minor_locator(Loct_x1)
    ax1.tick_params(axis='y', right=True, **majorTick)
    ax1.tick_params(axis='y', right=True, **minorTick) 
    ax1.set_ylabel("Magnitude")

    plt.tight_layout()
    ax.set_zorder(ax1.get_zorder()+1)
    ax.patch.set_visible(False)
    plt.show()

    return fig
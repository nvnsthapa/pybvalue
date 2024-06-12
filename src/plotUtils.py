# This is plotUtils modules for plotting seismic data 




def plot_seismicity(catalog_df):
    """
    Plot seismicity using GMT
    Args:
        catalog_df (pd.DataFrame): Catalog of seismicity, with columns longitude, latitude, depth, and magnitude
        filename (str): Filename for saving the plot
    Returns:
        Seismicity plot
    """
    import pygmt
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

    fig.basemap(region=region, projection="M15c", frame=True,)
    land = fig.grdimage(grid=grid_map, cmap="grayC", shading=True)
    fig.coast(land=land, borders=["1/1p,black"], water='skyblue', shorelines=True, resolution='f')
    pygmt.makecpt(cmap="seis", series=[catalog_df.depth.min(), 20])

    fig.plot(
        x=catalog_df.longitude,
        y=catalog_df.latitude,
        # size=0.02 * (2**catalog_df.magnitude),
        fill=catalog_df.depth,
        cmap=True,
        # style="cc",
        style="c0.05c",
        # pen="red",
        transparency=50,
    )
    fig.colorbar(frame="af+lDepth (km)")
    # Legend in the upper right
    # fig.legend(spec='quake_mag_sym.txt', position="jTR+o0.1c")
    # fig.savefig(f"Sesimicity_{filename}.pdf", dpi=600)
    # fig.show()
    return fig

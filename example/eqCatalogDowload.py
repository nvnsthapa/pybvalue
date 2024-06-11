
# %%
import obspy
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime

client = Client("SCEDC")
# %%
starttime=UTCDateTime(1977,1,11)
endtime=UTCDateTime.now()
network = "CI"
station = "PAS"
location = ""
channel = "BH*"
min_lat = 31.0
max_lat = 35.1
min_lon = -117.0
max_lon = -112.9
min_depth = 0
max_depth = 30
min_mag = 0
max_mag = 8


# %%
stations = client.get_stations(network=network, level='channel', filename="stations.xml")
# %%
print(stations)
# %%
stations.plot(projection="local", resolution="i")
# %%
try: 
    events = client.get_events(starttime=starttime, 
                           endtime=endtime,
                            minlatitude=min_lat, maxlatitude=max_lat,
                            minlongitude=min_lon, maxlongitude=max_lon,
                            # mindepth=min_depth, maxdepth=max_depth,
                            # minmagnitude=min_mag, maxmagnitude=max_mag,
                             )
except Exception as e:
    print(f"An error occurred: {e}")

print(events.__str__(print_all=True))
# %%
events.write("events.csv", format="csv")
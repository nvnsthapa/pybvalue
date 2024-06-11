
import obspy
from obspy.clients.fdsn import Client

client = Client("SCEDC")

network = "CI"
station = "PAS"
location = ""
channel = "BH*"

starttime=obspy.UTCDateTime(2016,1,1), 
endtime=obspy.UTCDateTime(2016,2,1)

stations = client.get_stations(network=network, level='channel')
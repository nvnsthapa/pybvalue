#Earthquake catalog download from SCEDC network using obspy 

from pystp import STPClient
from datetime import datetime

client = STPClient('stp3.gps.caltech.edu', 9999) #possible option stp.gps.caltech.edu, stp2.gps.caltech.edu, stp3.gps.caltech.edu
client.connect()

starttime=datetime(1977, 1, 1, 0, 0, 0)
endtime=datetime(2024, 6, 11, 0, 0, 0)

# # Download earthquake catalog from SCEDC network
events = client.get_events(times=[starttime, endtime],
                            mags=[0, 8], 
                            types=['eq'],
                            lats=[31.0, 35.0], lons=[-117.0, -112.0], 
                            depths=[0, 30],
                            # output_file='catalogDownload.txt',
                            # is_xml=True,
                            )

# print(events)
print(events.__str__(print_all=True))
client.disconnect()

import pandas as pd
import numpy as np


def read_cat_SCEDC(catalogFilename):
    """
    Read catalog file from SCEDC txt file
    Args:
        catalogFilename (str): Path to catalog file
        Returns:
        catalog_df (pd.DataFrame): Catalog file with header names as follows:
            date, time, ET, GT, magnitude, M, latitude, longitude, depth, Q, EVID, NPH
    """
    headerName = ["date", "time", "ET", "GT", "magnitude", "M", "latitude", "longitude", "depth", "Q", "EVID", "NPH", "NGRM"]
    catalog_df = pd.read_csv(
        catalogFilename,
        delim_whitespace=True,
        header=0,
        names=headerName,
        comment="#",
        index_col=False,
        skiprows=1,
    )
    print(f"Total Number of Events: {len(catalog_df)}")

    return catalog_df

def read_cat_relocated(relocatedCatalogF):
    """
    Read relocated catalog file by GrowClust
    Args:
        relocatedCatalogF (str): Path to relocated catalog file
    Returns:
        relocatedCat (pd.DataFrame): Relocated catalog file with header names as follows:
            yr, mon, day, hr, min, sec, eID, latR, lonR, depR, mag, qID, cID, nbranch, qnpair, qndiffP, qndiffS, rmsP, rmsS, eh, ez, et, latC, lonC, depC, eType, magType, GrowClust, relocBox
    """
    headerGc = ["yr", "mon", "day", "hr", "min", "sec", "eID", "latR", "lonR", "depR", "mag", "qID", "cID", "nbranch", "qnpair", "qndiffP", "qndiffS", "rmsP", "rmsS", "eh", "ez", "et", "latC", "lonC", "depC", "eType", "magType", "GrowClust", "relocBox"]
    relocatedCat = pd.read_csv(
        relocatedCatalogF,
        delim_whitespace=True,
        header=0,
        names=headerGc,
        index_col=False,
    )
    print(f"Total Number of Events: {len(relocatedCat)}")
    return relocatedCat


def create_magnitude_legend_file(mag_min, mag_max, filename='quake_mag_sym.txt'):
    """
    Create a legend file for GMT plot
    Args:
        mag_min (int): Minimum magnitude
        mag_max (int): Maximum magnitude
        filename (str): Filename for saving the legend file
    Returns:
        txt file e.g.: quake_mag_sym.txt 
    """
    mags = range(mag_min, mag_max)
    mag_sym_size = [0.02 * 2 ** m for m in mags]
    with open(filename, 'w') as sio:
        for mag, ms in zip(mags, mag_sym_size):
            if ms > 4.0:  # need extra space above big symbol
                sio.write('G 10l\n')
            sio.write('S 0.1i c %.1f none 0.50p 0.3i M=%d\n' % (ms, mag))
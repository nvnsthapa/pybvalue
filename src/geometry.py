"""
Sampling geometry.

Everything here works in kilometres on a locally projected plane, never in
degrees. A degree of longitude is ~90 km at Parkfield and ~94 km at the Salton
Sea, so a "radius" expressed in degrees is not a circle on the ground and would
silently distort every map. Coordinates are projected with an azimuthal
equidistant projection centred on the region, which preserves distance from
that centre - the right choice for neighbourhood sampling over a few hundred km.

Two geometries, one interface (`project` / `select`):

    CrossSection - a vertical plane through P1->P2 of a given width. Distances
                   run along the section and in depth; this is what Schorlemmer
                   et al. (2004) map at Parkfield.
    MapRegion    - plain map view, distances are epicentral. The natural
                   product for a basin rather than a single fault strand.
"""

import numpy as np
import pandas as pd
from pyproj import CRS, Transformer


def local_crs(lat0, lon0):
    """Azimuthal equidistant CRS in km, centred on (lat0, lon0)."""
    return CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +units=km +ellps=WGS84 +no_defs"
    )


class CrossSection:
    """
    A vertical cross section from P1 to P2, of finite width.

    Events are described by two coordinates in the section plane - `along`
    (km from P1) and `depth` - plus `perp`, the signed distance out of the
    plane. Sampling radii act within the plane; the width acts across it.

    Args:
        p1, p2 (tuple): (latitude, longitude) endpoints, degrees
        width_km (float): full section width; events are kept within +-width/2
    """

    def __init__(self, p1, p2, width_km=5.0):
        self.p1, self.p2, self.width_km = tuple(p1), tuple(p2), float(width_km)
        mid_lat = (self.p1[0] + self.p2[0]) / 2
        mid_lon = (self.p1[1] + self.p2[1]) / 2
        self.crs = local_crs(mid_lat, mid_lon)
        self._to_xy = Transformer.from_crs("EPSG:4326", self.crs, always_xy=True)

        self._origin = np.array(self._to_xy.transform(self.p1[1], self.p1[0]))
        end = np.array(self._to_xy.transform(self.p2[1], self.p2[0]))
        vector = end - self._origin
        self.length_km = float(np.hypot(*vector))
        if self.length_km == 0:
            raise ValueError("p1 and p2 are the same point")
        self._unit = vector / self.length_km

    @property
    def strike_deg(self):
        """Azimuth of P1->P2, degrees clockwise from north."""
        return float(np.degrees(np.arctan2(self._unit[0], self._unit[1])) % 360)

    def project(self, events):
        """
        Add `along` and `perp` columns (km) to a catalog.

        Args:
            events (pd.DataFrame): needs `latitude`, `longitude`
        Returns:
            pd.DataFrame: a copy with `along` and `perp`
        """
        x, y = self._to_xy.transform(events["longitude"].values,
                                     events["latitude"].values)
        offset = np.column_stack([x - self._origin[0], y - self._origin[1]])
        return events.assign(
            along=offset @ self._unit,
            perp=offset[:, 0] * self._unit[1] - offset[:, 1] * self._unit[0],
        )

    def select(self, events, depth_range=None, trim_ends=True):
        """
        Keep the events inside the section volume.

        Args:
            events (pd.DataFrame): catalog, projected or not
            depth_range (tuple): (min_km, max_km), inclusive; None for no cut
            trim_ends (bool): also require 0 <= along <= length
        Returns:
            pd.DataFrame: projected and filtered
        """
        if "along" not in events.columns or "perp" not in events.columns:
            events = self.project(events)
        keep = events["perp"].abs() <= self.width_km / 2
        if trim_ends:
            keep &= events["along"].between(0, self.length_km)
        if depth_range is not None:
            keep &= events["depth"].between(*depth_range)
        return events[keep]

    def node_grid(self, spacing_km=0.5, depth_range=(0, 16)):
        """
        Regular nodes on the section plane, as used for b-value mapping.

        Args:
            spacing_km (float): node spacing in both along and depth
            depth_range (tuple): (min_km, max_km)
        Returns:
            tuple: (along_1d, depth_1d, along_2d, depth_2d)
        """
        along = np.arange(0, self.length_km + spacing_km, spacing_km)
        depth = np.arange(depth_range[0], depth_range[1] + spacing_km, spacing_km)
        grid_along, grid_depth = np.meshgrid(along, depth)
        return along, depth, grid_along, grid_depth

    def __repr__(self):
        return (f"CrossSection(p1={self.p1}, p2={self.p2}, "
                f"width_km={self.width_km}, length_km={self.length_km:.1f}, "
                f"strike={self.strike_deg:.1f})")


class MapRegion:
    """
    Plain map view over a lat/lon box. Distances are epicentral kilometres in a
    locally projected plane, so a radius is a true circle on the ground.

    Args:
        bbox (tuple): (minlat, maxlat, minlon, maxlon)
    """

    def __init__(self, bbox):
        self.bbox = tuple(bbox)
        mid_lat = (self.bbox[0] + self.bbox[1]) / 2
        mid_lon = (self.bbox[2] + self.bbox[3]) / 2
        self.crs = local_crs(mid_lat, mid_lon)
        self._to_xy = Transformer.from_crs("EPSG:4326", self.crs, always_xy=True)
        self._origin = np.array(self._to_xy.transform(mid_lon, mid_lat))

    def project(self, events):
        """Add `x`, `y` columns (km east/north of the region centre)."""
        x, y = self._to_xy.transform(events["longitude"].values,
                                     events["latitude"].values)
        return events.assign(x=x - self._origin[0], y=y - self._origin[1])

    def select(self, events, depth_range=None):
        """Keep events inside the box (and optionally a depth range)."""
        if "x" not in events.columns:
            events = self.project(events)
        keep = (events["latitude"].between(self.bbox[0], self.bbox[1])
                & events["longitude"].between(self.bbox[2], self.bbox[3]))
        if depth_range is not None:
            keep &= events["depth"].between(*depth_range)
        return events[keep]

    def node_grid(self, spacing_km=0.5):
        """Regular map-view nodes, in the projected km frame."""
        corners_lon = [self.bbox[2], self.bbox[3], self.bbox[3], self.bbox[2]]
        corners_lat = [self.bbox[0], self.bbox[0], self.bbox[1], self.bbox[1]]
        xs, ys = self._to_xy.transform(corners_lon, corners_lat)
        xs = np.asarray(xs) - self._origin[0]
        ys = np.asarray(ys) - self._origin[1]
        x = np.arange(xs.min(), xs.max() + spacing_km, spacing_km)
        y = np.arange(ys.min(), ys.max() + spacing_km, spacing_km)
        grid_x, grid_y = np.meshgrid(x, y)
        return x, y, grid_x, grid_y

    def __repr__(self):
        return f"MapRegion(bbox={self.bbox})"

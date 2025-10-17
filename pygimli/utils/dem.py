#!/usr/bin/env python
"""Digital Elevation Model (DEM) class for interpolating elevations."""
import os.path
import math
import numpy as np


class DEM:
    """Interpolation class for digital elevation models."""

    def __init__(self, demfile, x=None, y=None, **kwargs):
        """Initialize DGM (regular grid) interpolation object.

        Parameters
        ----------
        demfile : str or iterable elevations (list, ndarray)
            digital elevation file:
            * ASC file with lower left corner and spacing in header OR
            * x, y, z list of grid points or irregular points

        x, y : iterable of (unique) x and y positions matching z

        Keyword Arguments
        -----------------
        toLatLon: callable(x, y) [None]
            Custom coordinate translator. If set to None then
            `lambda x_, y_: utm.to_latlon(x_, y_, zone, 'N')` is taken.
        zone: int [32]
            UTM zone to be chosen
        """
        from scipy.interpolate import RegularGridInterpolator

        self.latlon = False
        self.x = x
        self.y = y
        self.fallback = kwargs.pop('fallback', kwargs.pop('z0', None))

        if isinstance(self.fallback, str):
            self.fallback = DEM(self.fallback)

        if isinstance(demfile, (list, tuple)):
            self.__init__(demfile[0], **kwargs)
            for addfile in demfile[1:]:
                self.add(addfile)
            self.dem = RegularGridInterpolator((self.x, self.y),
                                    np.fliplr(self.z.T))
        elif isinstance(demfile, str):
            if demfile[-4:].lower() == '.asc':
                self.loadASC(demfile)
            elif demfile[-4:].lower() == '.hgt':
                self.loadHGT(demfile)
            else:
                self.loadTXT(demfile)
        elif x is not None and y is not None:
            self.z = demfile
        else:
            raise Exception("Either DEM file or z with x and y must be given!")

        if self.latlon:
            self._toLatLon = kwargs.pop('toLatLon', None)
            if self._toLatLon is None:
                import utm
                zone = kwargs.pop('zone', 32)
                self._toLatLon = lambda x_, y_: utm.to_latlon(x_, y_, zone, 'N')


    def __call__(self, x, y=None):
        """Interpolation function."""
        if self.latlon:
            y, x = self._toLatLon(x, y)

        if y is None:
            return self.dem(x)
        else:
            if hasattr(self, 'tri'):
                out = self.dem(x, y)
            else:
                out = self.dem((x, y))

        if isinstance(out, np.ma.MaskedArray):
            out = out.data

        if self.fallback is not None:
            if isinstance(self.fallback, (float, int)):
                out[np.isnan(out)] = self.fallback
            elif isinstance(self.fallback, DEM):
                out[np.isnan(out)] = self.fallback(x, y)

        return out


    def loadTXT(self, demfile):
        """Load column-based DEM."""
        import matplotlib.tri as mtri
        from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator

        xp, yp, zp = np.loadtxt(demfile, unpack=True)
        be = self.fallback is None
        if len(np.unique(xp)) * len(np.unique(yp)) > len(xp):
            self.x = xp
            self.y = yp
            self.z = zp
            self.tri = mtri.Triangulation(self.x, self.y)
            self.dem = mtri.CubicTriInterpolator(self.tri, self.z)
            # self.dem = LinearNDInterpolator(np.column_stack((xp, yp)), zp)
            return

        if np.isclose(yp[0], yp[1], rtol=1e-32, atol=1e-2):
            if np.isclose(xp[0], xp[1], rtol=1e-32, atol=1e-2):
                print('Fatal error! Neither first two x- nor y- coords are '
                      'increasing in the specified xyz file! Aborting...')
                raise SystemExit
            elif xp[1] > xp[0]:
                nx = np.argwhere(np.diff(xp) < 0)[0][0] + 1
                ny = len(xp) // nx
                x = xp[:nx]
                y = yp[::nx]
                zp = zp.reshape((ny, nx))
        else:
            if yp[1] > yp[0]:
                ny = np.argwhere(np.diff(yp) < 0)[0][0] + 1
            else:
                ny = np.argwhere(np.diff(yp) > 0)[0][0] + 1
            nx = len(xp) // ny
            x = xp[::ny]
            y = yp[:ny]
            zp = zp.reshape((nx, ny)).T

        if y[1] < y[0]:
            y.sort()
            zp = np.flipud(zp)
        if x[1] < x[0]:
            x.sort()
            zp = np.fliplr(zp)

        self.z = zp
        self.x = x
        self.y = y
        self.dem = RegularGridInterpolator((self.x, self.y),
                                           self.z.T,
                                           bounds_error=be)


    def loadASC(self, ascfile):
        """Load ASC (DEM matrix with location header) file."""
        from scipy.interpolate import RegularGridInterpolator

        with open(ascfile) as fid:
            header = {}
            sp = []
            nheader = 0
            while len(sp) < 3:
                sp = fid.readline().split()
                if len(sp) < 3:
                    header[sp[0]] = float(sp[1].replace(',', '.'))
                    nheader += 1

        self.z = np.flipud(np.genfromtxt(ascfile, skip_header=nheader))
        if 'NODATA_value' in header:
            self.z[self.z == header['NODATA_value']] = np.nan
        dx = header.pop('cellsize', 1.0)  # just a guess
        self.x = np.arange(header['ncols']) * dx + header['xllcorner']
        self.y = np.arange(header['nrows']) * dx + header['yllcorner']
        be = self.fallback is None
        self.dem = RegularGridInterpolator((self.x, self.y),
                                           self.z.T, bounds_error=be)


    def loadHGT(self, hgtfile):
        """Load ASC (DEM matrix with location header) file."""
        from scipy.interpolate import RegularGridInterpolator
        siz = os.path.getsize(hgtfile)
        samples = int(math.sqrt(siz/2))
        lat = int(hgtfile[-10:-8])
        lon = int(hgtfile[-7:-4])
        self.x = np.linspace(lon,lon+1,samples, endpoint=False)
        self.y = np.linspace(lat,lat+1,samples, endpoint=False)
        be = self.fallback is None
        with open(hgtfile, 'rb') as hgt_data:
            self.z = np.fromfile(hgt_data, np.dtype('>i2'),
                                 samples*samples).reshape((samples, samples))
            self.z[self.z < -32000] = 0

        self.dem = RegularGridInterpolator((self.x, self.y),
                                           np.fliplr(self.z.T), bounds_error=be)
        self.latlon = True


    def add(self, new):
        """Combine two DEM by concatenatation.

        x or y vectors must be equal (e.g. for 1° SRTM models).

        Parameters
        ----------
        new : DEM | str
            DEM instance or string to load
        """
        if isinstance(new, (str, list, tuple)):
            new = DEM(new)
        assert isinstance(new, DEM), "No DEM instance!"
        if np.allclose(self.y, new.y):
            if self.x[0] < new.x[0]:
                self.x = np.concatenate([self.x, new.x])
                self.z = np.hstack([self.z, new.z])
            else:
                self.x = np.concatenate([new.x, self.x])
                self.z = np.hstack([new.z, self.z])
        elif np.allclose(self.x, new.x):
            if self.y[0] < new.y[0]:
                self.y = np.concatenate([self.y, new.y])
                self.z = np.vstack([new.z, self.z])
            else:
                self.y = np.concatenate([new.y, self.y])
                self.z = np.hstack([self.z, new.z])


    def show(self, cmap="terrain", cbar=True, ax=None, **kwargs):
        """Show digital elevation model (i.e. the elevation map).

        Arguments
        ---------
        cmap: str="terrain"
            matplotlib colormap definiton
        cbar: bool=True
            Add colorbar to the plot or not.
        ax: matplotlib axes object
            Add the plot to a given axes object or create a new one.

        Keyword Args
        ------------
        **kwargs, type keyword arguments
            Additional keyword arguments forwarded to mpl plot.
        """
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
        from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator

        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 12)))

        # extract some kwargs for axis setting and colorbar
        orientation = kwargs.pop('orientation', 'vertical')
        xlim = kwargs.pop('xlim', (-9e99, 9e99))
        ylim = kwargs.pop('ylim', (-9e99, 9e99))
        clim = kwargs.pop('clim', (np.min(self.z), np.max(self.z)))
        cmap = kwargs.pop('cmap', 'terrain')
        nl = kwargs.pop("nl", 15)
        if isinstance(self.dem, mtri.TriInterpolator):
            im = ax.tricontourf(self.tri, self.z, cmap=cmap,
                                levels=np.linspace(*clim, nl))
            ax.triplot(self.tri, '-', color='gray', alpha=0.5, lw=0.5)
        elif isinstance(self.dem, LinearNDInterpolator):
            im = ax.tripcolor(self.dem.points[:, 0], self.dem.points[:, 1],
                              self.dem.tri)
        else:
            x = self.dem.grid[0]
            y = self.dem.grid[1]
            ix0 = np.argmin(x < xlim[0])
            ix1 = np.argmax(x > xlim[1]) - 1
            if ix1 < 0:
                ix1 = len(x)
            iy0 = np.argmin(y < ylim[0])
            iy1 = np.argmax(y > ylim[1]) - 1
            if iy1 < 0:
                iy1 = len(y)

            im = ax.pcolormesh(x[ix0:ix1], y[iy0:iy1],
                               self.dem.values[ix0:ix1, iy0:iy1].T,
                               cmap=cmap, **kwargs)

        im.set_clim(clim)
        cb = None
        if cbar:
            # norm = Normalize(vmin=clim[0], vmax=clim[1])
            cb = plt.colorbar(im, ax=ax, orientation=orientation)
            if clim:
                cb.vmin = clim[0]
                cb.vmax = clim[1]

        ax.set_aspect(1.0)
        return ax


if __name__ == '__main__':  # if called directly as a script
    dgm = DEM('dgm5_borkum.asc')
    ax = dgm.show(vmin=0, vmax=10)


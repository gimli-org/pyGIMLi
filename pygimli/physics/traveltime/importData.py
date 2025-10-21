#!/usr/bin/env python
"""Import travel time data from various formats into DataContainerTT."""
import struct
from pathlib import Path
import numpy as np
import pygimli as pg
from .tt import DataContainerTT


def load(fileName, verbose=False, **kwargs):
    """Shortcut to load TravelTime data.

    Import Data and try to assume the file format.

    Parameters
    ----------
    fileName: str

    Returns
    -------
    data: pg.DataContainer
    """
    if fileName.lower().endswith('.gtt'):
        data = importGTT(fileName)
    elif fileName.lower().endswith('.tom'):
        data = readTOMfile(fileName)
    elif "nR" in kwargs or "nS" in kwargs:
        data = importAsciiColumns(fileName, **kwargs)
    else:
        data = DataContainerTT(fileName)

    return data


def importGTT(filename, return_header=False):
    """Import refraction data from Tomo+ GTT data file into DataContainer."""
    header = {}
    with Path(filename).open('rb') as fid:
        block = fid.read(100)
        nshots = struct.unpack(">I", block[:4])[0]
        ngeoph = struct.unpack(">I", block[4:8])[0]
        header['ntrace'] = struct.unpack(">Q", block[8:16])[0]
        header['nchan'] = struct.unpack(">I", block[16:20])[0]
        header['tminmax'] = struct.unpack(">2f", block[20:28])
        header['offsetminmax'] = struct.unpack(">2f", block[28:36])
        header['angle'] = struct.unpack(">f", block[36:40])[0]
        header['origin'] = struct.unpack(">2f", block[40:48])
        header['unit'] = struct.unpack(">I", block[48:52])[0]
        header['shotSpacing'] = struct.unpack(">f", block[52:56])[0]
        header['receiverSpacing'] = struct.unpack(">f", block[56:60])[0]
        SPOS = np.zeros((nshots, 3))
        RPOS = np.zeros((ngeoph*5, 3))
        SHOT, REC, TT, VA = [], [], [], []
        # tmat = np.ones((nshots, ngeoph)) * np.nan
        for i in range(nshots):
            block = fid.read(24)  # shot information
            shotid = struct.unpack(">I", block[:4])[0]
            nci = struct.unpack(">I", block[4:8])[0]  # channels for the shot
            spos = np.array(struct.unpack(">4f", block[8:24]))
            SPOS[i, :] = spos[:3]
            X, Y = [], []
            for _j in range(nci):
                block = fid.read(24)  # receiver information
                recid = struct.unpack(">I", block[:4])[0]
                rpos = np.array(struct.unpack(">4f", block[4:20]))
                RPOS[recid, :] = rpos[:3]
                tt = struct.unpack(">f", block[20:24])[0]
                # print(shotid, recid, tt)
                SHOT.append(i+1)  # shotid)
#                SHOT.append(shotid)
                REC.append(recid)
                TT.append(tt)
                X.append(rpos[0])
                Y.append(rpos[0])
                offset = np.sqrt(np.sum((rpos-spos)**2))
                VA.append(offset/tt)
                # tmat[shotid, recid] = tt

        SHOT = np.array(SHOT, dtype=int) - 1
        REC = np.array(REC, dtype=int) - 1 + len(SPOS)
        pos = np.vstack((SPOS, RPOS[1:max(REC)+1]))
        x, ifwd, irev = np.unique(pos[:, 0],
                                  return_index=True, return_inverse=True)
        data = pg.DataContainer()
        data.registerSensorIndex('s')
        data.registerSensorIndex('g')
        for i in ifwd:
            data.createSensor([pos[i, 0], pos[i, 2], 0])

        data.resize(len(TT))
        data.set('t', np.array(TT))
        data.set('va', np.array(VA))
        data.set('s', irev[SHOT].astype(float))
        data.set('g', irev[REC].astype(float))
        data.markValid(data('t') > 0.)
        data.removeUnusedSensors()
        data.markInvalid(data('g') >= data.sensorCount())
        data.markInvalid(data('s') >= data.sensorCount())
        data.removeInvalid()
        if return_header:
            return data, header
        else:
            print(header)
            return data


def importAsciiColumns(filename, ndig=2, roundto=0,
                       nS=None, nR=None, nT=None, nA=None):
    """Read in columns from ASCII file.

    Parameters
    ----------
    filename : str
        filename holding position and traveltime columns
    nS : [int, int, int]
        columns for shot positions
    nR : [int, int, int]
        columns for shot positions
    """
    if nS is None:
        nS = [0, 1]
    if nR is None:
        nR = [max(nS)+1, max(nS)+2]

    A = np.genfromtxt(filename)
    if nT is None:
        nT = A.shape[1] - 1  # last one

    posS = A[:, nS]
    posR = A[:, nR]
    # if twoD: take only receivers for defining topo
    uR, iR = np.unique(posR[:, 0], return_index=True)
    ux, ix, ii = np.unique(np.vstack([posS[:, 0], posR[:, 0]]),
                           return_index=True, return_inverse=True)
    uz = np.interp(ux, posR[iR, 0], posR[iR, 1])
    data = DataContainerTT()
    for xx, zz in zip(ux, uz, strict=False):
        data.createSensor(pg.Pos(xx, zz))

    ndata = A.shape[0]
    data.resize(ndata)
    data["s"] = ii[:ndata]
    data["g"] = ii[ndata:]
    data["t"] = A[:, nT]
    data["valid"] = 1
    return data

def readTOMfile(filename, ndig=2, roundto=0):
    """Read Reflex tomography (*.TOM) file."""
    t, xT, zT, xR, zR = np.loadtxt(filename, usecols=(0, 2, 3, 4, 5), unpack=1)
    if roundto > 0:
        pT = (np.round(xT/roundto) - np.round(zT/roundto) * 1j) * roundto
        pR = (np.round(xR/roundto) - np.round(zR/roundto) * 1j) * roundto
    else:
        pT = xT.round(ndig) - zT.round(ndig) * 1j
        pR = xR.round(ndig) - zR.round(ndig) * 1j

    pU = np.unique(np.hstack((pT, pR)))
    iT = np.array([np.nonzero(pU == pi)[0][0] for pi in pT], dtype=float)
    iR = np.array([np.nonzero(pU == pi)[0][0] for pi in pR], dtype=float)
    data = pg.DataContainer()
    for pp in pU:
        data.createSensor(pg.RVector3(pp.real, pp.imag))

    for tok in ['s', 'g']:
        data.registerSensorIndex(tok)

    data.resize(len(t))
    data.set('t', t)
    data.set('s', iT)
    data.set('g', iR)
    data.markValid(pg.abs(data('s') - data('g')) > 0)
    return data


class ReadAHL:
    """Class reading seismic refraction format provided by Uppsala University.

    Supply a filename and a delimiter character. The delimiter is used to
    find the header string which in turn is used to label and extract the
    columns.
    """

    def __init__(self, filename, maxsensorid, delimiter='|'):
        self.header = None
        self.delimiter = delimiter
        self.filename = filename
        self.alldata = []  # contains all of the data (unaltered)
        self.skiprows = 0  # how many leading rows are in the file
        self.labels = dict()  # labels and corresponding column
        self.maxsensorid = maxsensorid  # for excluding geophones ( >= )
        self.sensorcols = []  # the columns that contain info about the sensors
        self.shotcols = []  # the columns that contain info about the shot
        self.positionlist = []  # the list of sensor and shot positions
        self.use_xz_only = True

    def __call__(self):
        """Call function."""
        self.load()
        self.convert()

    # def __repr__(self):
    #     """Return string representation."""
    #     s = 'header line: ' + self.header
    #     s += '\nalldata length: {}\nfirst row: '.format(len(self.alldata))
    #     s += str(self.alldata[0, :]) + '\nlast row: '
    #     s += str(self.alldata[-1, :])

    #     return s

    def __repr__(self):
        """Return string representation."""
        in_arg = "'" + self.filename + "'" + \
            ', maxsensorid=' + str(self.maxsensorid) + \
            ', delimiter=' + "'" + self.delimiter + "'"
        return self.__class__.__name__ + "(" + in_arg + ")"

    def _extractlabels(self):
        """Extract the labels from the header line."""
        labels = self.header.split(self.delimiter)
        n = 0
        for l in labels:
            l = l.strip()
            if len(l) > 0:
                self.labels[l] = n  # add labels as keys and colnums as value
                if l.find('DELAY') >= 0 or l.find('REC') >= 0:
                    self.sensorcols.append(n)
                else:
                    self.shotcols.append(n)
                n += 1

    def _extractheader(self):
        """Search for the header and extract it."""
        # open the file and get some info out:
        # number of lines to skip and the header
        with Path(self.filename).open('r') as f:
            # loop until we find the delimiter as the first character of line
            firstchar = ''
            while firstchar is not self.delimiter:
                self.header = f.readline().strip()
                firstchar = self.header[0]
                self.skiprows += 1  # update number of leading rows in the file

        self._extractlabels()  # split the header line into individual labels

    def load(self):
        """Process/read the data file."""
        # this will determine number of leading rows and extract labels etc...
        self._extractheader()

        with Path(self.filename).open('r') as f:
            # skip some rows as determined by _extractheader()
            for _ in range(self.skiprows):
                f.next()

            # build a 2D list where each row is a line from the file
            currentshot = 0
            for line in f:
                # convert to int so we can use
                datarow_str = line.strip().split()
                datarow = [int(d) for d in datarow_str if len(d) > 0]

                # normal case is that we have data without shotid, so add it
                if len(datarow) < len(self.labels):
                    datarow.insert(0, currentshot)
                else:
                    currentshot = datarow[0]

                self.alldata.append(datarow)

        # convert it all to numpy array for easier manipulation later
        self.alldata = np.asarray(self.alldata)

    def convert(self):
        """Extract sensors and build a list for pyGIMLi indexing.

        Also determine relevant column numbers for the data.
        """
        if len(self.alldata) == 0.:
            raise ValueError('Data not read yet!')

        shots_col = self.labels['SHOT_PEG']
        sensors_col = self.labels['REC_PEG']

        # cut out the unwanted sensors
        data = self.alldata[self.alldata[:, shots_col] <= self.maxsensorid, :]
        # ...and shotpositions
        data = data[data[:, sensors_col] <= self.maxsensorid, :]

        shots_uniq = np.unique(data[:, shots_col])
        sensors_uniq = np.unique(data[:, sensors_col])
        all_uniq = np.unique(np.row_stack((shots_uniq[:, np.newaxis],
                                           sensors_uniq[:, np.newaxis])))

        # remap the sensor ids to indices starting from 1 going to N
        sensor_map = zip(all_uniq,
                         np.arange(1, len(all_uniq)+1, dtype=int), strict=False)
        for old, new in sensor_map:
            data[data[:, sensors_col] == old, sensors_col] = new
            data[data[:, shots_col] == old, shots_col] = new

        # now create a list of unique (x, y, z)'s
        final_uniq_list, flat_idx = np.unique(
            data[:, [shots_col, sensors_col]], return_index=True)
        idx = np.unravel_index(flat_idx,
                               data[:, [shots_col, sensors_col]].shape)
        # rows and corresponding columns (column==0 means sensor, 1 means shot)
        # Here we build a list of which columns to pick the positions from
        # TODO: Not sure that -3 is always correct!
        column_list = [np.asarray(self.sensorcols[-3:]) if c == 1
                       else np.asarray(self.shotcols[-3:]) for c in idx[1]]

        unique_xyz = data[idx[0][:, np.newaxis], column_list]

#        mpl.plot(unique_xyz[:, 1], unique_xyz[:, 2], 'b.-')
#        mpl.show(block=False)

        if self.use_xz_only:
            x = unique_xyz[:, 1]
            y = unique_xyz[:, 2]
            xx = np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))
            xx = xx.tolist()
            xx.insert(0, 0.0)
            xx = np.asarray(xx)
            final_pos = np.column_stack((xx[:, np.newaxis],
                                         unique_xyz[:, 0])) / 10.0
        else:
            # TODO: Not sure that [1,2,0] is always the correct order!
            final_pos = unique_xyz[:, [1, 2, 0]] / 10.0

        # change the order of the positions to be in x,y,z order
        # also change to [m] instead of [dm]
        self.save(data=final_pos, desc='pos')

        sgt_cols = [self.labels['SHOT_PEG'], self.labels['REC_PEG'],
                    self.labels['X1:DELAY']]
        sgt = data[:, sgt_cols]
#        valid = sgt[:, -1] > 0
        sgt = sgt.astype(float)
        sgt[:, -1] /= 1000.0  # milliseconds to seconds

#        self.save(data=np.column_stack((sgt, valid)), desc='sgt')
        self.save(data=sgt, desc='sgt')

    def save(self, data, desc):
        """Write the converted data to disk."""
        fname_in = self.filename
        fname_out = fname_in[:fname_in.rfind('.')] + '.sgt'

        mode = 'w' if desc == 'pos' else 'a'

        with Path(fname_out).open(mode) as f:
            if desc == 'pos':
                if self.use_xz_only:
                    f.write(str(data.shape[0]) + ' # No positions\n#x z\n')
                else:
                    f.write(str(data.shape[0]) + ' # No positions\n#x y z\n')

                np.savetxt(f, data, fmt='%.2f')
            elif desc == 'sgt':
                f.write(str(data.shape[0]) + ' # Number of data\n#s g t\n')
                np.savetxt(f, data, fmt='%i %i %.4f')
            else:
                raise ValueError(f'Invalid description of data to be written: {desc}')


if __name__ == '__main__':
    pass

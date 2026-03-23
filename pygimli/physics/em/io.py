#!/usr/bin/env python
"""Import functions for several Em data formats (to be extended)."""
from pathlib import Path
import numpy as np


def readusffile(filename, data=None):
    """Read data from single USF (universal sounding file) file.

    data = readusffile( filename )
    data = readusffile( filename, data ) will append to data
    """
    if not data:
        data = []

    columns = []
    nr = 0
    sounding = {}
    sounding['FILENAME'] = filename
    isdata = False
    with Path(filename).open() as fid:
        for line in fid:
            zeile = line.rstrip('\n').replace(',', '')  # commas useless here
            if zeile:  # anything at all
                if zeile[0] == '/':  # comment-like
                    if zeile[1:4] == 'END':  # end of a sounding
                        if isdata:  # have already read some data
                            sounding['data'] = columns
                            for i, cn in enumerate(sounding['column_names']):
                                sounding[cn] = columns[:, i]

                            data.append(sounding)
                            sounding = {}

                        isdata = not isdata  # turn off data mode
                    elif zeile.find(':') > 0:  # key-value pair
                        key, value = zeile[1:].split(':')
                        sounding[key] = value if isinstance(value, str) else float(value)
                else:
                    if isdata:
                        values = zeile.split()
                        try:
                            for i, v in enumerate(values):
                                columns[nr, i] = float(v)

                            nr += 1
                        except:
                            sounding['column_names'] = values
                            columns = np.zeros((int(sounding['POINTS']),
                                                len(values)))
                            nr = 0

    return data


def readusffiles(filenames):
    """Read all soundings data from a list of usf files.

    Example
    -------
    DATA = readusffiles(filenames)
    """
    DATA = []
    for onefile in filenames:
        DATA = readusffile(onefile, DATA)

    return DATA


def importMaxminData(filename, verbose=False):
    """Import function reading in positions, data, frequencies, geometry."""
    delim = None
    coilspacing = 0.
    freq = []
    x, IP, OP = [], [], []
    with Path(filename).open() as fid:
        i = 0
        for aline in fid:
            if aline.split()[0][0].isdigit():  # number found
                break
            elif aline.upper().find('COIL') > 0:  # [:6] == '/ COIL':
                coilspacing = float(aline.split()[-2])
            elif aline.upper().find('FREQ') > 0:  # [:6] == '/ FREQ':
                freq = np.array([float(aa) for aa in aline[aline.find(
                    ':') + 1:].replace(',', ' ').split() if aa[0].isdigit()])

        fid.close()

        if verbose:
            print("CS=", coilspacing, "F=", freq)
        if aline.find(',') > 0:
            delim = ','

        nf = len(freq)
        if verbose:
            print("delim=", delim, "nf=", nf)

        A = np.loadtxt(filename, skiprows=i, delimiter=delim).T
        x, IP, OP = A[0], A[2:nf * 2 + 2:2].T, A[3:nf * 2 + 2:2].T

    return x, freq, coilspacing, IP, OP


if __name__ == '__main__':
    print("do some tests here")

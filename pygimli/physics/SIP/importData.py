#!/usr/bin/env python
"""Import/Export for SIP data."""

import codecs
from datetime import  datetime
from pathlib import Path
import numpy as np
import re

import pygimli as pg

def load(fileName, verbose=False, **kwargs):
    """Shortcut to load SIP spectral data.

    Import Data and try to assume the file format.

    Parameters
    ----------
    fileName: str

    Returns
    -------
    freqs, amp, phi : np.array
        Frequencies, amplitudes and phases phi in neg. radiant

    """
    firstLine = None
    with codecs.open(fileName, 'r', encoding='iso-8859-15',
                     errors='replace') as fi:
        firstLine = fi.readline()

    f, amp, phi = None, None, None

    fnLow = fileName.lower()

    if 'SIP Fuchs III' in firstLine:
        if verbose:
            pg.info("Reading SIP Fuchs III file")
        f, amp, phi, header = readFuchs3File(fileName,
                                             verbose=verbose, **kwargs)
        phi *= -np.pi/180.
        # print(header) # not used?
    elif 'SIP-Quad' in firstLine:
        if verbose:
            pg.info("Reading SIP Quad file")
        f, amp, phi, header = readFuchs3File(fileName, nfr=9, namp=10, nphi=11,
                                             nk=7, verbose=verbose, **kwargs)
        phi *= -np.pi/180.
    elif 'SIP-Fuchs' in firstLine:
        if verbose:
            pg.info("Reading SIP Fuchs file")
        f, amp, phi, drhoa, dphi = readRadicSIPFuchs(fileName,
                                                     verbose=verbose, **kwargs)
        phi *= -np.pi/180.
    elif fnLow.endswith('.txt') or fnLow.endswith('.csv'):
        f, amp, phi = readTXTSpectrum(fileName)
        amp *= 1.0 # scale it with k if available
    else:
        raise NotImplementedError("Don't know how to read data.")

    return f, amp, phi


def fstring(fri):
    """Format frequency to human-readable (mHz or kHz)."""
    if fri > 1e3:
        fstr = f'{int(np.round(fri/1e3)):d} kHz'
    elif fri < 1.:
        fstr = f'{int(np.round(fri*1e3)):d} mHz'
    elif fri < 10.:
        fstr = f'{fri:3.1f} Hz'
    elif fri < 100.:
        fstr = f'{fri:4.1f} Hz'
    else:
        fstr = f'{int(np.round(fri)):d} Hz'
    return fstr


def readTXTSpectrum(filename, nfr=0, namp=1, nphi=3, sphi=-1):
    """Read spectrum from ZEL device output (txt) data file."""
    f, amp, phi = [], [], []
    with Path(filename).open() as fid:
        lines = fid.readlines()
        for line in lines[1:]:
            snums = line.replace(';', ' ').split()
            if len(snums) > 3:
                f.append(float(snums[nfr]))
                amp.append(float(snums[namp]))
                phi.append(sphi*float(snums[nphi]))
            else:
                break

        return np.asarray(f), np.asarray(amp), np.asarray(phi)


def readFuchs3File(resfile, k=1.0, verbose=False, nfr=11, namp=12, nphi=13, nk=9):
    """Read Fuchs III (SIP spectrum) data file.

    Parameters
    ----------
    k : float
        Overwrite internal geometric factor from device.

    """
    activeBlock = ''
    header = {}
    LINE = []
    dataAct = False
    with codecs.open(resfile, 'r', encoding='iso-8859-15', errors='replace') as fid:
        for line in fid:
            line = line.replace('\r\n', '\n') # correct for carriage return
            if dataAct:
                LINE.append(line)
                if len(line) < 2:
                    f, amp, phi, kIn = [], [], [], []
                    for li in LINE:
                        sline = li.split()
                        if len(sline) > 12:
                            fi = float(sline[nfr])
                            if np.isfinite(fi):
                                f.append(fi)
                                amp.append(float(sline[namp]))
                                phi.append(float(sline[nphi]))
                                kIn.append(float(sline[nk]))

                    if k != 1.0 and verbose is True:
                        pg.info("Geometric value changed to:", k)

                    return np.array(f), np.array(amp)/np.array(kIn) * k, \
                           np.array(phi), header
            elif len(line):
                if line.rfind('Current') >= 0:
                    if dataAct:
                        break
                    else:
                        dataAct = True

                if line[0] == '[':
                    token = line[1:line.rfind(']')].replace(' ', '_')
                    if token[:3] == 'End':
                        header[activeBlock] = np.array(header[activeBlock])
                        activeBlock = ''
                    elif token[:5] == 'Begin':
                        activeBlock = token[6:]
                        header[activeBlock] = []
                    else:
                        value = line[line.rfind(']') + 1:]
                        try:  # direct line information
                            num = float(value) if '.' in value else int(value)
                            header[token] = num
                        except BaseException:
                            # maybe beginning or end of a block
                            pass

                else:
                    if activeBlock:
                        nums = np.array(line.split(), dtype=float)
                        header[activeBlock].append(nums)


def readRadicSIPFuchs(filename, readSecond=False, delLast=True):
    """Read SIP-Fuchs Software rev.: 070903.

    Read Radic instrument res file containing a single spectrum.

    Please note the apparent resistivity value might be scaled with the
    real geometric factor. Default is 1.0.

    Parameters
    ----------
    filename : string

    readSecond: bool [False]
        Read the first data block[default] or read the second that
        consists in the file.

    delLast : bool [True]
        ??

    Returns
    -------
    fr : array [float]
        Measured frequencies

    rhoa : array [float]
        Measured apparent resistivties

    phi : array [float]
        Measured phases

    drhoa : array [float]
        Measured apparent resistivties error

    phi : array [float]
        Measured phase error
    """
    with codecs.open(filename, 'r', encoding='iso-8859-15', errors='replace') as f:
        line = f.readline()
        fr = []
        rhoa = []
        phi = []
        drhoa = []
        dphi = []
        while True:
            line = f.readline()
            if line.rfind('Freq') > -1:
                break
        return
        if readSecond:
            while True:
                if f.readline().rfind('Freq') > -1:
                    break

        while True:
            line = f.readline()
            b = line.split('\t')
            if len(b) < 5:
                break

            fr.append(float(b[0]))
            rhoa.append(float(b[1]))
            phi.append(-float(b[2]) * np.pi / 180.)
            drhoa.append(float(b[3]))
            dphi.append(float(b[4]) * np.pi / 180.)

        f.close()

    if delLast:
        fr.pop(0)
        rhoa.pop(0)
        phi.pop(0)
        drhoa.pop(0)
        dphi.pop(0)

    return np.array(fr), np.array(rhoa), np.array(phi), np.array(drhoa), np.array(dphi)


def toTime(t, d):
    """Convert time format into timestamp 11:08:02, 21/02/2019."""
    tim = [int(_t) for _t in t.split(':')]
    if '/' in d:  # 03/02/1975
        day = [int(_t) for _t in d.split('/')]
        dt = datetime(year=day[2], month=day[1], day=day[0],
                      hour=tim[0], minute=tim[1], second=tim[2])
    elif '.' in d:  # 03.02.1975
        day = [int(_t) for _t in d.split('.')]
        dt = datetime(year=day[2], month=day[1], day=day[0],
                      hour=tim[0], minute=tim[1], second=tim[2])
    else:  # 1975-02-03
        day = [int(_t) for _t in d.split('-')]
        dt = datetime(year=day[0], month=day[1], day=day[2],
                      hour=tim[0], minute=tim[1], second=tim[2])

    return dt.timestamp()


def readSIP256file(resfile, verbose=False):
    """Read SIP256 file (RES format) - mostly used for 2d SIP by pybert.sip.

    Read SIP256 file (RES format) - mostly used for 2d SIP by pybert.sip.

    Parameters
    ----------
    filename: str
        *.RES file (SIP256 raw output file)
    verbose: bool
        do some output [False]

    Returns
    -------
        header - dictionary of measuring setup
        DATA - data AB-list of MN-list of matrices with f, amp, phi, dAmp, dPhi
        AB - list of current injection
        RU - list of remote units

    Examples
    --------
        header, DATA, AB, RU = readSIP256file('myfile.res', True)
    """
    activeBlock = ''
    header = {}
    LINE = []
    dataAct = False

    with codecs.open(resfile, 'r', encoding='iso-8859-15',
                                   errors='replace') as fi:
        content = fi.readlines()

    for line in content:
        if dataAct:
            LINE.append(line)
        elif len(line):
            if line[0] == '[':
                token = line[1:line.rfind(']')].replace(' ', '_')
                # handle early 256D software bug
                if 'FrequencyParameterBegin' in token:
                    token = token.replace('FrequencyParameterBegin',
                                            'Begin_FrequencyParameter')
                if 'FrequencyParameterEnd' in token:
                    token = token.replace('FrequencyParameterEnd',
                                            'End_FrequencyParameter')

                if token.replace(' ', '_') == 'Messdaten_SIP256':
                    dataAct = True
                elif 'Messdaten' in token:
                    # res format changed into SIP256D .. so we are a
                    # little bit more flexible with this.
                    dataAct = True
                elif token[:3] == 'End':
                    header[activeBlock] = np.array(header[activeBlock])
                    activeBlock = ''
                elif token[:5] == 'Begin':
                    activeBlock = token[6:]
                    header[activeBlock] = []
                else:
                    value = line[line.rfind(']') + 1:]
                    try:  # direct line information
                        if '.' in value:
                            num = float(value)
                        else:
                            try:
                                num = int(value)
                            except BaseException:
                                num = 0
                                pass
                        header[token] = num
                    except BaseException as e:
                        # maybe beginning or end of a block
                        print(e)
            else:
                if activeBlock:
                    nums = np.array(line.split(), dtype=float)
                    header[activeBlock].append(nums)

    DATA, dReading, dFreq, AB, RU, ru = [], [], [], [], [], []
    # tMeas = []
    for i, line in enumerate(LINE):
        # print(i, line)
        line = line.replace(' nc ', ' 0 ') # no calibration should 0
        line = line.replace(' c ', ' 1 ') # calibration should 1
        # sline = line.split()
        sline = line.rstrip('\r\n').split()
        if line.find('Reading') == 0:
            rdno = int(sline[1])
            if rdno > 0:
                AB.append((int(sline[4]), int(sline[6])))
            if ru:
                RU.append(ru)
                ru = []
            if rdno > 1 and dReading:
                dReading.append(np.array(dFreq))
                DATA.append(dReading)
                pg.verbose(f'Reading {rdno-1}:{len(dReading)} RUs')
                dReading, dFreq = [], []
        elif line.find('Remote Unit') == 0:
            ru.append(int(sline[2]))
            if dFreq:
                dReading.append(np.array(dFreq))
                dFreq = []
        elif line.find('Freq') >= 0:
            pass
        elif len(sline) > 1 and rdno > 0:  # some data present
            # search for two numbers (with .) without a space inbetween
            # variant 1: do it for every part
            for i, ss in enumerate(sline):
                if re.search(r'\.20[01][0-9]', ss) is None:  # no date
                    fd = re.search(r'\.[0-9-]*\.', ss)
                    if fd:
                        bpos = ss[1:].find('-') + 1 if '-' in ss[1:] else fd.start() + 4

                        # print(ss[:bpos], ss[bpos:])
                        if ss[5:8] != ".20":
                            sline.insert(i, ss[:bpos])
                            sline[i+1] = ss[bpos:]
                        # print(sline)
                    fd = re.search(r'NaN[0-9-]*\.', ss)
                    if fd:
                        bpos = ss.find('-') if '-' in ss[1:] else fd.start() + 3

                        # print(ss[:bpos], ss[bpos:])
                        sline.insert(i, ss[:bpos])
                        sline[i+1] = ss[bpos:]
                        # print(sline)
            # variant 2: do it on whole line
            # cdate = re.search('\.20[01][0-9]', line)
            # if cdate:
            #     n2000 = cdate.start()
            # else:
            #     n2000 = len(line)

            # print(sline)
            # concnums = re.search('\.[0-9-]*\.', line[:n2000])
            # while concnums:
            #     bpos = concnums.span()[0] + 4
            #     line = line[:bpos] + ' ' + line[bpos:]
            #     n2000 += 1
            #     concnums = re.search('\.[0-9-]*\.', line[:n2000])
            #     sline = line.rstrip('\r\n').split()
            #     print(sline)

            # if re.search('[0-9]-', line[:85]):  # missing whitespace before -
            #     sline = re.sub('[0-9]-', '5 -', line).split()
                # not a good idea for dates

            for c in range(7): # this is expensive .. do we really need this?
                if len(sline[c]) > 15:  # too long line / missing space
                    if c == 0:
                        part1 = sline[c][:-15]
                        part2 = sline[c][-15:]   # [10:]
                    else:
                        part1 = sline[c][:-10]
                        part2 = sline[c][-10:]   # [11:]
                    sline = sline[:c] + [part1] + [part2] + sline[c + 1:]

                if sline[c].find('c') >= 0:
                    sline[c] = '1.0'
            try:
                dFreq.append(
                    np.array(sline[:8] + [toTime(sline[-2], sline[-1])],
                             dtype=float))
            except BaseException as e:
                # dFreq.append(np.array(sline[:8], dtype=float))
                print(i, line, sline)
                raise IndexError("Could not convert line!") from e

    dReading.append(np.array(dFreq))
    DATA.append(dReading)
    pg.verbose(f'Reading {rdno}:{len(dReading)} RUs')
    return header, DATA, AB, RU


#Frequency /Hz       RA/Ohmm    PA/�      ERA/%     EPA/�     Cal?     IA/mA     K.-F./m
#     Gains  Time/h:m:s    Date/d.m.y
#20000.00000000        0.4609  -6.72598   0.02234   0.01280    1      20.067        1.00
#     0     11:08:02     21/02/2019
if __name__ == "__main__":
    pass

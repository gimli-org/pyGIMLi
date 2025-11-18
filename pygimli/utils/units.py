#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Collection of frequently used physical quantities with units and
there abbreviations."""

import pygimli as pg
from ..core.config import rc

quants = {
    'rhoa': {
        'name': 'Apparent resistivity',
        'unit': r'$\Omega$m',
        'ger': 'scheinbarer spez. elektr. Widerstand',
        'cMap': 'Spectral_r',
    },
    'res': {
        'name': 'Resistivity',
        'unit': r'$\Omega$m',
        'ger': 'spez. elektr. Widerstand',
        'cMap': 'Spectral_r',
    },
    'ip': {
        'name': 'Phase',
        'unit': 'mrad',
        'ger': 'Phasenwinkel',
        'cMap': 'viridis',
    },
   'ip_n': {
        'name': 'neg. Phase',
        'unit': 'mrad',
        'ger': 'neg. Phasenwinkel',
        'cMap': 'viridis',
    },
    'ipa': {
        'name': 'Apparent phase',
        'unit': 'mrad',
        'ger': 'scheinbarer Phasenwinkel',
        'cMap': 'viridis',
    },
    'ipa_n': {
        'name': 'neg. Apparent phase',
        'unit': 'mrad',
        'ger': 'neg. scheinbarer Phasenwinkel',
        'cMap': 'viridis',
    },
    'u': {
        'name': 'Voltage',
        'unit': 'V',
        'ger': 'Spannung',
    },
    'i': {
        'name': 'Current',
        'unit': 'A',
        'ger': 'Stromstärke',
    },
    'err': {
        'name': 'Error',
        'unit': '-',
        'ger': 'Fehler',
    },
    'e': {
        'name': 'electric field',
        'unit': 'V/m',
        'ger': 'elektrisches Feld',
    },
    'b': {
        'name': 'magnetic flux',
        'unit': 'V/m²',
        'ger': 'magnetische Flussdichte',
    },
    'h': {
        'name': 'magnetic field',
        'unit': 'A/m',
        'ger': 'magnetische Feldstärke',
    },
    'va': {
        'name': 'Apparent velocity',
        'unit': 'm/s',
        'ger': 'Scheingeschwindigkeit',
    },
    'vel': {
        'name': 'Velocity',
        'unit': 'm/s',
        'ger': 'Geschwindigkeit',
    },
    'slo': {
        'name': 'Slowness',
        'unit': 's/m',
        'ger': 'Slowness',
    },
    'por': {
        'name': 'Porosity',
        'unit': None,
        'ger': 'Porosität',
        'cMap': 'pink_r',
    },
}

rc['quants'] = quants


def quantity(name):
    """Return quantity for given name."""
    quantity = None

    if name.lower() not in quants:
        for k, v in quants.items():
            if v['name'].lower() == name.lower():
                quantity = v
                break
    else:
        quantity = quants[name.lower()]
    return quantity


def cmap(name):
    """Return default colormap for physical quantity name."""
    q = quantity(name)
    if q is None:
        pg.warn('No information about quantity name', name)
        return 'viridis'
    return q.get('cMap', 'viridis')


def unit(name, unit='auto'):
    """ Return the name of a physical quantity with its unit.

    TODO
    ----
        * example
        * localization

    Parameters
    ----------
    """
    q = quantity(name)

    if unit == 'auto' and q is None:
        pass
        # fall back if the name is given instead of the abbreviation
        # print(quants)
        # pg.error('Please give abbreviation or full name '
        #          'for the quantity name: {0}'.format(name))
    else:
        if rc['lang'] == 'german' or rc['lang'] == 'de' or rc['lang'] == 'ger':
            name = q['ger']
        else:
            name = q['name']

        if unit == 'auto':
            unit = q['unit']

    if unit is None:
        return '{0}'.format(name)

    if rc['unitStyle'] == 1 or rc['lang'] == 'german' or \
            rc['lang'] == 'de' or rc['lang'] == 'ger':
        return '{0} in {1}'.format(name, unit)
    elif rc['unitStyle'] == 2:
        return '{0} ({1})'.format(name, unit)
    elif rc['unitStyle'] == 3:
        return '{0} [{1}]'.format(name, unit)

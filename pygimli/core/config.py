"""Management of pygimli user configuration file."""

import json
import os
import sys

from .logger import info

# pyGIMLi default configuration
rc = {
    'lang': 'eng', # 'eng', 'german', 'de', 'ger'
    'unitStyle': 2,
    # 1: german style: 'value in unit'
    # 2: default style: 'value (unit)'
    'view3D': 'auto',
    'pyvista.backend': 'client',
    # auto: Use pyvista if installed or set it to 'fallback' to force fallback mode
    'globalCache': True,
    # call pg.wait() before the terminal script ends if there are pending 
    # mpl widgets and your backend this supports
    'waitOnExit': True,
}


def getCPUCount():
    """Return number of processors on multiple platoforms."""
    # Windows
    if os.name == 'nt':
        return int(os.getenv('NUMBER_OF_PROCESSORS'))
    # Linux
    elif sys.platform == 'linux2':
        retv = 0
        with open('/proc/cpuinfo', 'rt') as cpuinfo:
            for line in cpuinfo:
                if line[:9] == 'processor':
                    retv += 1
        return retv

    # Please add similar hacks for MacOSX, Solaris, Irix,
    # FreeBSD, HPUX, etc.
    else:
        raise RuntimeError('unknown platform')


def getConfigPath():
    r"""Return the user space path for configuration files.

    - Windows: 'C:\Documents and Settings\username\Appdata\pygimli'
    - Linux: '~/.config/pygimli' (if not overwritten by $XDG_CONFIG_HOME)
    - Mac: '~/Library/Preferences/pygimli'

    See Also
    --------
    pygimli.getCachePath
    """
    appname = "pygimli"
    system = sys.platform

    if system == "win32":
        path = os.path.join(os.environ['APPDATA'])
    elif system == "darwin":
        path = os.path.expanduser('~/Library/Preferences/')
    else:
        path = os.getenv('XDG_CONFIG_HOME', os.path.expanduser("~/.config"))

    return os.path.join(path, appname)


__configpath = getConfigPath()
__configfile = os.path.join(__configpath, "config.json")

if os.path.isfile(__configfile):
    #info("Loading user configuration file at " + __configfile)
    with open(__configfile) as cfg:
        userrc = json.load(cfg)

    # Check if user config holds all keys and update if necessary
    if len(userrc.keys()) != len(rc.keys()):
        #print(userrc)
        rc.update(userrc)
        #print(rc)
        info("Updating user configuration.")
        # # for key in rc:
        #     if key not in userrc:
        #         info("Updating user configuration with", key, "=", rc[key])
        #         userrc[key] = rc[key]
        with open(__configfile, "w") as cfg:
            json.dump(rc, cfg, indent=4, sort_keys=True)
    rc.update(userrc)
else:
    info("Creating default user configuration file at " + __configfile)
    os.makedirs(os.path.join(__configpath, "pygimli"), exist_ok=True)
    with open(__configfile, "w") as cfg:
        json.dump(rc, cfg, indent=4, sort_keys=True)

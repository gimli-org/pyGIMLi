#!/usr/bin/env python3
"""Management of pygimli user configuration file."""

import json
import os
import sys
from pathlib import Path

from .logger import info, debug

from .core import pgcore

# pyGIMLi default configuration
rc = {
    "lang": "eng",  # 'eng', 'german', 'de', 'ger'
    "unitStyle": 2,
    # 1: german style: 'value in unit'
    # 2: default style: 'value (unit)'
    # auto: Use pyvista if installed or set to 'fallback' to force fallback mode
    "view3D": "auto",
    "pyvista.backend": "client",
    "globalCache": True,
    # call pg.wait() before the terminal script ends if there are pending
    # mpl widgets and your backend this supports
    "waitOnExit": True,
    "matplotlib": "ipympl",
}


def getCPUCount():
    """Return number of processors on multiple platforms."""
    if os.name == "nt":
        return int(os.getenv("NUMBER_OF_PROCESSORS"))
    # Linux
    elif sys.platform == "linux2":
        retv = 0
        with Path("/proc/cpuinfo").open() as cpuinfo:
            for line in cpuinfo:
                if line[:9] == "processor":
                    retv += 1
            else:
                return pgcore.numberOfCPU()
        return retv
    else:
        return pgcore.numberOfCPU()

    # Please add similar hacks for MacOSX, Solaris, Irix,
    # FreeBSD, HPUX, etc.


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
        path = Path(os.environ["APPDATA"])
    elif system == "darwin":
        path = Path("~/Library/Preferences/").expanduser()
    else:
        path = Path(
            os.getenv("XDG_CONFIG_HOME", str(Path("~/.config").expanduser()))
        )

    return str(path / appname)


__configpath = getConfigPath()
__configfile = str(Path(__configpath) / "config.json")


if Path(__configfile).is_file():
    # info("Loading user configuration file at " + __configfile)
    with Path(__configfile).open() as cfg:
        debug("Loading user configuration file at " + __configfile)
        userrc = json.load(cfg)

    # Check if user config holds all keys and update if necessary
    if len(userrc.keys()) != len(rc.keys()):
        # print(userrc)
        rc.update(userrc)
        # print(rc)
        info("Updating user configuration.")
        # # for key in rc:
        #     if key not in userrc:
        #         info("Updating user configuration with", key, "=", rc[key])
        #         userrc[key] = rc[key]
        with Path(__configfile).open("w") as cfg:
            json.dump(rc, cfg, indent=4, sort_keys=True)
    rc.update(userrc)
else:
    info("Creating default user configuration file at " + __configfile)
    (Path(__configpath) / "pygimli").mkdir(parents=True, exist_ok=True)
    with Path(__configfile).open("w") as cfg:
        json.dump(rc, cfg, indent=4, sort_keys=True)

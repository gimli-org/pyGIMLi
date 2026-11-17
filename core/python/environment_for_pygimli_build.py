#!/usr/bin/env python
"""Provide the settings needed to build the pygimli python bindings.

Collects the platform defines and include paths used by
generate_pygimli_code.py when it runs castxml and pyplusplus.
"""

import sys
from pathlib import Path

import distutils.sysconfig

this_module_dir_path = Path(__file__).resolve().parent
project_root = (this_module_dir_path / ".." / ".." / "..").resolve()

complete_path = lambda *args: str(project_root.joinpath(*args))


def _abspath(path):
    """Return path as an absolute, resolved string."""
    return str(Path(path).resolve())


# Get the major Python version without patchlevel
# print("Python version:", distutils.sysconfig.get_python_version())
# print("Python lib:", distutils.sysconfig.get_python_lib(standard_lib=True))


class Settings:
    """Hold the settings used to build the pygimli python bindings."""

    module_name = "_pygimli_"
    gimli_path = "../src"
    gimli_defines = ""
    includes_paths = []

    @staticmethod
    def setup_environment():
        """Add the pyplusplus and pygccxml sources to sys.path."""
        Settings.pygccxml_path = "../../../pygccxml"
        Settings.pyplusplus_path = "../../../pyplusplus"

        if Path("../../../gccxml-bin/bin").exists():
            Settings.gccxml_path = _abspath("../../../gccxml-bin/bin")
        elif Path("../../gccxml-bin/bin").exists():
            Settings.gccxml_path = _abspath("../../gccxml-bin/bin")

        sys.path.append(Settings.pygccxml_path)
        sys.path.append(Settings.pyplusplus_path)


if sys.platform == "linux2" or sys.platform == "linux":
    pass

elif sys.platform == "darwin":
    pass

elif sys.platform == "win32":
    Settings.gimli_defines = "MINGW"
    Settings.python_libs_path = distutils.sysconfig.get_python_lib(
        standard_lib=True
    )

    if Path("../../../boost/include").exists():
        Settings.includes_paths.append(_abspath("../../../boost/include"))

else:
    raise RuntimeError(
        f'There is no configuration for "{sys.platform}" platform.'
    )

Settings.setup_environment()

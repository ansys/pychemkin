# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Platform-specific setup and diagnostics for Chemkin-CFD-API."""

import ctypes
import datetime
import os
from pathlib import Path
import platform
import sys
from typing import Any, cast

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.logger import logger


def __diagnose_windows_dll_dependencies(target_lib: Path) -> list[tuple[str, str]]:
    """Return unresolved direct DLL imports for a Windows shared library."""
    if platform.system() != "Windows":
        return []
    try:
        import pefile
    except ImportError:
        return []

    unresolved: list[tuple[str, str]] = []
    pe = pefile.PE(str(target_lib))
    for entry in getattr(pe, "DIRECTORY_ENTRY_IMPORT", []):
        dep_name = entry.dll.decode(errors="ignore")
        try:
            # getattr returns Any so mypy will not complain about missing WinDLL
            # on some platforms/stubs
            getattr(ctypes, "WinDLL")(dep_name)
        except OSError as exc:
            unresolved.append((dep_name, str(exc)))
        except AttributeError:
            # If WinDLL isn't available on this Python build, record as unresolved
            unresolved.append((dep_name, "WinDLL not available in ctypes"))
    return unresolved


def get_runtime_diagnostics(
    config: tuple[int, str, str, str, list[str]],
    probe_load: bool = False,
) -> dict[str, object]:
    """Collect runtime diagnostics for Chemkin shared library initialization.

    Parameters
    ----------
    config : tuple[int, str, str, str, list[str]]
        Configuration tuple containing Chemkin initialization parameters.
    probe_load : bool, default: False
        Whether to actively probe loading target and dependency DLLs.

    Returns
    -------
    dict[str, object]
        Structured diagnostics with resolved paths and dependency status.

    """
    # unpack the configuration tuple
    ansys_ver = int(config[0])
    ansys_dir = str(config[1])
    target_lib = str(config[3])
    lib_paths = cast(list[str], config[4])
    #
    target_path = Path(target_lib) if target_lib else Path()
    diagnostics: dict[str, object] = {
        "platform": platform.system(),
        "ansys_version": ansys_ver,
        "ansys_dir": ansys_dir,
        "target_library": str(target_path) if target_lib else "",
        "target_exists": bool(target_lib) and target_path.is_file(),
        "configured_library_paths": list(lib_paths),
        "unresolved_direct_imports": [],
        "dependency_probes": [],
        "target_load_probe": None,
    }

    if not target_lib or platform.system() != "Windows":
        if probe_load and target_lib:
            try:
                ctypes.CDLL(target_lib)
                diagnostics["target_load_probe"] = {
                    "success": True,
                    "error": "",
                }
            except OSError as exc:
                diagnostics["target_load_probe"] = {
                    "success": False,
                    "error": str(exc),
                }
        return diagnostics

    unresolved = __diagnose_windows_dll_dependencies(target_path)
    diagnostics["unresolved_direct_imports"] = [name for name, _ in unresolved]

    if probe_load:
        try:
            ctypes.CDLL(target_lib)
            diagnostics["target_load_probe"] = {
                "success": True,
                "error": "",
            }
        except OSError as exc:
            diagnostics["target_load_probe"] = {
                "success": False,
                "error": str(exc),
            }

    dependency_probes: list[dict[str, object]] = []
    for dep_name, dep_error in unresolved:
        dep_result: dict[str, Any] = {
            "name": dep_name,
            "search_path_load_error": dep_error,
            "found_in_configured_paths": [],
            "full_path_probe": None,
        }

        for lib_path in lib_paths:
            candidate = Path(lib_path) / dep_name
            if candidate.is_file():
                dep_result["found_in_configured_paths"].append(str(candidate))

        if probe_load:
            candidate_paths = dep_result["found_in_configured_paths"]
            if candidate_paths:
                full_path = candidate_paths[0]
                try:
                    getattr(ctypes, "WinDLL")(full_path)
                    dep_result["full_path_probe"] = {
                        "success": True,
                        "path": full_path,
                        "error": "",
                    }
                except OSError as exc:
                    dep_result["full_path_probe"] = {
                        "success": False,
                        "path": full_path,
                        "error": str(exc),
                    }
            else:
                dep_result["full_path_probe"] = {
                    "success": False,
                    "path": "",
                    "error": "not found in configured library paths",
                }

        dependency_probes.append(dep_result)

    diagnostics["dependency_probes"] = dependency_probes
    return diagnostics


def __setup_windows(
    min_ver: int, valid_vers: list[int]
) -> tuple[int, tuple[int, str, str, str, list[str]]]:
    """Set up PyChemkin environment on Windows platforms.

    Parameters
    ----------
    min_ver : int
        Minimum version of Ansys to consider.
    valid_vers : list[int]
        List of valid Ansys versions to check.

    Returns
    -------
    tuple[int, tuple[int, str, str, str, list[str]]]
        Status code and configuration tuple.
        Structured diagnostics with resolved paths and dependency status.

    """
    ansys_ver = min_ver
    ansys_dir = ""
    ckbin = ""
    target_lib = ""
    lib_paths: list[str] = []
    ansyshome = Path()
    # set ansys installation directory (Windows)
    for v in valid_vers:
        ansys_ver = v
        if v >= min_ver:
            _ansys_installation = "ANSYS" + str(ansys_ver) + "_DIR"
            _ansys_home = os.environ.get(_ansys_installation, "NA")
            if _ansys_home != "NA":
                ansyshome = Path(_ansys_home).parent
                ansys_dir = str(ansyshome)
                break
        else:
            break

    if ansys_ver >= min_ver:
        if str(ansyshome) == "." or not ansyshome.is_dir():
            # no local Ansys installation
            msg = [
                Color.RED,
                "PyChemkin cannot find the specific Ansys installation:",
                ansys_dir,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            return 1, (ansys_ver, ansys_dir, ckbin, target_lib, lib_paths)
        else:
            # platform label
            plat = "winx64"
            # chemkin bin directory
            ckbin = "chemkin.win64"
            lib_addition = ansyshome / "reaction" / ckbin / "bin"
            lib_paths = [str(lib_addition)]
            # required third-party shared objects
            if ansys_ver <= 252:
                # <= 25R2
                lib_addition = ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat
                lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat
                lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "zlib" / "1.2.13" / plat
                lib_paths.append(str(lib_addition))
            elif ansys_ver == 261:
                # == 26R1
                lib_addition = ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat
                lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat
                lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "zlib" / plat
                lib_paths.append(str(lib_addition))
            else:
                # >= 27R1
                lib_addition = ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat
                lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "IntelMKL" / "2024.2.3" / plat
                if not Path(lib_addition).is_dir():
                    lib_addition = ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat
                lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "zlib" / plat
                lib_paths.append(str(lib_addition))
    else:
        msg = [
            Color.RED,
            "PyChemkin does not support Chemkin versions older than 2025R1.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.critical(this_msg)
        return 2, (ansys_ver, ansys_dir, ckbin, target_lib, lib_paths)
    # set load dll paths
    if sys.platform == "win32":
        for path in lib_paths:
            if Path(path).is_dir():
                os.add_dll_directory(path)
            else:
                msg = [
                    Color.RED,
                    "PyChemkin cannot find the required shared library path:",
                    str(path),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                return 1, (ansys_ver, ansys_dir, ckbin, target_lib, lib_paths)
        # set Chemkin-CFD-API shared object
        my_target = ansyshome / "reaction" / ckbin / "bin" / "KINeticsdll.dll"
        target_lib = str(my_target)
    return 0, (ansys_ver, ansys_dir, ckbin, target_lib, lib_paths)


def __setup_linux(
    min_ver: int, valid_vers: list[int]
) -> tuple[int, tuple[int, str, str, str, list[str]]]:
    """Set up PyChemkin environment on Linux platforms.

    Parameters
    ----------
    min_ver : int
        Minimum version of Ansys to consider.
    valid_vers : list[int]
        List of valid Ansys versions to check.

    Returns
    -------
    tuple[int, tuple[int, str, str, str, list[str]]]
        Status code and configuration tuple.
        Structured diagnostics with resolved paths and dependency status.

    """
    ansys_ver = min_ver
    ansys_dir = ""
    ckbin = ""
    target_lib = ""
    lib_paths: list[str] = []
    #
    ierr = 0
    ansyshome = Path()
    # set ansys installation directory (Linux)
    for v in valid_vers:
        ansys_ver = v
        if v >= min_ver:
            _ansys_installation = "ANSYS" + str(ansys_ver) + "_DIR"
            _ansys_home = os.environ.get(_ansys_installation, "NA")
            if _ansys_home != "NA":
                ansyshome = Path(_ansys_home).parent
                ansys_dir = str(ansyshome)
                break
        else:
            break
    # try using a different method
    if _ansys_home == "NA" and ansys_dir == "":
        # environment variable ANSYSxxx_DIR is NOT defined
        # check local Ansys installation
        _user_home = os.environ.get("HOME", "NA")
        if _user_home != "NA":
            ansyshome = Path(_user_home) / "ansys_inc"
            _ansys_home = str(ansyshome)
            found_home = False
            if ansyshome.is_dir():
                # find all local Ansys installations
                local_versions = [f.name for f in ansyshome.iterdir() if f.is_dir()]
                for v in valid_vers:
                    ansys_ver = v
                    if v >= min_ver:
                        this_version = "v" + str(v)
                        if this_version in local_versions:
                            ansys_dir = _ansys_home + this_version
                            found_home = True
                            break
                    else:
                        ierr = 2
                        break
                if not found_home:
                    ierr = 1
            else:
                # no local Ansys installation
                ierr = 1
        else:
            ierr = 1

    if str(ansyshome) == ".":
        ierr = 1
    # check Ansys version
    if ansys_ver < min_ver:
        ierr = 2

    # check Ansys installation error
    if ierr == 1:
        msg = [
            Color.RED,
            "failed to find local Ansys chemkin installation.\n",
            Color.SPACEx6,
            f"please make sure Ansys v{min_ver} or newer is installed locally\n",
            Color.SPACEx6,
            "otherwise, please set the environment variable",
            '"ANSYSxxx_DIR"\n',
            Color.SPACEx6,
            'with value = "<full path to local Ansys installation>/ANSYS"\n',
            Color.SPACEx6,
            'for example, ANSYS251_DIR = "$HOME/ansys_inc/v251/ANSYS".',
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.critical(this_msg)
        return 1, (ansys_ver, ansys_dir, ckbin, target_lib, lib_paths)
    elif ierr == 2:
        msg = [
            Color.RED,
            "PyChemkin does not support Chemkin versions older than 2025R1.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.critical(this_msg)
        return 2, (ansys_ver, ansys_dir, ckbin, target_lib, lib_paths)
    # platform label
    plat = "linx64"
    # chemkin bin directory
    ckbin = "chemkin.linuxx8664"
    lib_addition = ansyshome / "reaction" / ckbin / "bin"
    lib_paths = [str(lib_addition)]
    # required third-party shared objects
    if ansys_ver <= 252:
        # <= 25R2
        lib_addition = (
            ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat / "lib" / "intel64"
        )
        lib_paths.append(str(lib_addition))
        lib_addition = (
            ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat / "lib" / "intel64"
        )
        lib_paths.append(str(lib_addition))
        lib_addition = ansyshome / "tp" / "zlib" / "1.2.13" / plat / "lib"
        lib_paths.append(str(lib_addition))
    elif ansys_ver == 261:
        # == 26R1
        lib_addition = ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat
        lib_paths.append(str(lib_addition))
        lib_addition = ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat
        lib_paths.append(str(lib_addition))
        lib_addition = ansyshome / "tp" / "zlib" / plat
        lib_paths.append(str(lib_addition))
    else:
        # >= 27R1
        lib_addition = ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat
        lib_paths.append(str(lib_addition))
        lib_addition = ansyshome / "tp" / "IntelMKL" / "2024.2.3" / plat
        if not Path(lib_addition).is_dir():
            lib_addition = ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat
        lib_paths.append(str(lib_addition))
        lib_addition = ansyshome / "tp" / "zlib" / plat
        lib_paths.append(str(lib_addition))
    # set load shared object paths
    combined_path = ":".join(lib_paths)
    if "LD_LIBRARY_PATH" not in os.environ.keys():
        os.environ["LD_LIBRARY_PATH"] = combined_path
    else:
        os.environ["LD_LIBRARY_PATH"] = (
            os.environ["LD_LIBRARY_PATH"] + ":" + combined_path
        )

    if "PATH" not in os.environ.keys():
        os.environ["PATH"] = combined_path
    else:
        os.environ["PATH"] = os.environ["PATH"] + ":" + combined_path
    # set Chemkin-CFD-API shared object
    my_target = ansyshome / "reaction" / ckbin / "bin" / "libKINetics.so"
    target_lib = str(my_target)
    return 0, (ansys_ver, ansys_dir, ckbin, target_lib, lib_paths)


def find_valid_ansys_versions(min_ver: int = 251) -> list[int]:
    """Generate a list of valid ANSYS versions based on the current year.

    Parameters
    ----------
    min_ver : int, default: 251
        Minimum supported ANSYS version.

    Returns
    -------
    list[int]
        List of valid ANSYS versions.

    """
    valid_versions: list[int] = []
    sub_releases = [2, 1]
    # generate possible Ansys versions based on the current year
    this_date = datetime.datetime.now()
    this_year = this_date.year
    # the newest version cannot have release year later than next year
    # get the last two digits of the year: ## of 20##
    # (change to 21## when the 22nd century comes)
    # assemble the release year part
    max_release_year = ((this_year % 100) + 1) * 10
    test_release = max_release_year

    while test_release >= min_ver - (min_ver % 10):
        for r in sub_releases:
            valid_versions.append(test_release + r)
        test_release -= 10

    return valid_versions

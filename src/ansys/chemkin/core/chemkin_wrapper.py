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

"""Chemkin-CFD-API python interfaces."""

import ctypes
from ctypes import cdll
import datetime
import os
from pathlib import Path
import platform
import sys

import numpy as np

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.logger import logger


def __setwindows() -> int:
    """Set up PyChemkin environment on Windows platforms."""
    global _ansys_ver
    global _ansys_dir
    global _ckbin
    global _lib_paths
    global _min_version
    global _target_lib
    global _valid_versions
    ansyshome = Path()
    # set ansys installation directory (Windows)
    for v in _valid_versions:
        _ansys_ver = v
        if v >= _min_version:
            _ansys_installation = "ANSYS" + str(_ansys_ver) + "_DIR"
            _ansys_home = os.environ.get(_ansys_installation, "NA")
            if _ansys_home != "NA":
                ansyshome = Path(_ansys_home).parent
                _ansys_dir = str(ansyshome)
                break
        else:
            break

    if _ansys_ver >= _min_version:
        if str(ansyshome) == "." or not ansyshome.is_dir():
            # no local Ansys installation
            msg = [
                Color.RED,
                "PyChemkin cannot find the specific Ansys installation:",
                str(_ansys_dir),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            return 1
        else:
            plat = "winx64"
            _ckbin = "chemkin.win64"
            # required third-party shared objects
            lib_addition = ansyshome / "reaction" / _ckbin / "bin"
            _lib_paths = [str(lib_addition)]
            if _ansys_ver <= 252:
                # <= 25R2
                lib_addition = ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat
                _lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat
                _lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "zlib" / "1.2.13" / plat
                _lib_paths.append(str(lib_addition))
            else:
                # >= 26R1
                lib_addition = ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat
                _lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat
                _lib_paths.append(str(lib_addition))
                lib_addition = ansyshome / "tp" / "zlib" / plat
                _lib_paths.append(str(lib_addition))
    else:
        msg = [
            Color.RED,
            "PyChemkin does not support Chemkin versions older than 2025R1.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.critical(this_msg)
        return 2
    # set load dll paths
    if sys.platform == "win32":
        for _lib_path in _lib_paths:
            os.add_dll_directory(_lib_path)
        # set Chemkin-CFD-API shared object
        my_target = ansyshome / "reaction" / _ckbin / "bin" / "KINeticsdll.dll"
        _target_lib = str(my_target)
    return 0


def __setlinux() -> int:
    """Set up PyChemkin environment on Linux platforms."""
    global _ansys_ver
    global _ansys_dir
    global _ckbin
    global _lib_paths
    global _min_version
    global _target_lib
    global _valid_versions
    ierr = 0
    ansyshome = Path()
    # set ansys installation directory (Linux)
    for v in _valid_versions:
        _ansys_ver = v
        if v >= _min_version:
            _ansys_installation = "ANSYS" + str(_ansys_ver) + "_DIR"
            _ansys_home = os.environ.get(_ansys_installation, "NA")
            if _ansys_home != "NA":
                ansyshome = Path(_ansys_home).parent
                _ansys_dir = str(ansyshome)
                break
        else:
            break
    # try using a different method
    if _ansys_home == "NA" and _ansys_dir == "":
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
                for v in _valid_versions:
                    _ansys_ver = v
                    if v >= _min_version:
                        this_version = "v" + str(v)
                        if this_version in local_versions:
                            _ansys_dir = _ansys_home + this_version
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
    if _ansys_ver < _min_version:
        ierr = 2

    # check Ansys installation error
    if ierr == 1:
        msg = [
            Color.RED,
            "failed to find local Ansys chemkin installation.\n",
            Color.SPACEx6,
            "please make sure Ansys v251 or newer is installed locally\n",
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
        return 1
    elif ierr == 2:
        msg = [
            Color.RED,
            "PyChemkin does not support Chemkin versions older than 2025R1.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.critical(this_msg)
        return 2
    # required third-party shared objects
    plat = "linx64"
    _ckbin = "chemkin.linuxx8664"
    lib_addition = ansyshome / "reaction" / _ckbin / "bin"
    _lib_paths = [str(lib_addition)]
    if _ansys_ver <= 252:
        # <= 25R2
        lib_addition = (
            ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat / "lib" / "intel64"
        )
        _lib_paths.append(str(lib_addition))
        lib_addition = (
            ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat / "lib" / "intel64"
        )
        _lib_paths.append(str(lib_addition))
        lib_addition = ansyshome / "tp" / "zlib" / "1.2.13" / plat / "lib"
        _lib_paths.append(str(lib_addition))
    else:
        # >= 26R1
        lib_addition = (
            ansyshome / "tp" / "IntelCompiler" / "2023.1.0" / plat / "lib" / "intel64"
        )
        _lib_paths.append(str(lib_addition))
        lib_addition = (
            ansyshome / "tp" / "IntelMKL" / "2023.1.0" / plat / "lib" / "intel64"
        )
        _lib_paths.append(str(lib_addition))
        lib_addition = ansyshome / "tp" / "zlib" / plat / "lib"
        _lib_paths.append(str(lib_addition))
    # set load dll paths
    combined_path = ":".join(_lib_paths)
    if "LD_LIBRARY_PATH" not in os.environ.keys():
        # if os.environ["LD_LIBRARY_PATH"] is None:
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
    my_taget = ansyshome / "reaction" / _ckbin / "bin" / "libKINetics.so"
    _target_lib = str(my_taget)
    return 0


# set ansys version number
_min_version = 251
_valid_versions: list[int] = []
_sub_releases = [2, 1]
_ansys_ver = _min_version
_ansys_dir = ""
_ckbin = ""
status = 0
_target_lib = ""
_lib_paths: list[str] = []
# generate possible Ansys versions based on the current year
this_date = datetime.datetime.now()
_this_year = this_date.year
# the newest version cannot have release year later than next year
# get the last two digits of the year: ## of 20##
# (change to 21## when the 22nd century comes)
# assemble the release year part
_max_release_year = ((_this_year % 100) + 1) * 10
_test_release = _max_release_year
while _test_release >= _min_version - (_min_version % 10):
    for r in _sub_releases:
        _valid_versions.append(_test_release + r)
    _test_release -= 10
# create log
msg = ["minimum Ansys version to run PyChemkin =", str(_min_version)]
this_msg = Color.SPACE.join(msg)
logger.debug(this_msg)
# check platform
if platform.system() == "Windows":
    # set ansys installation directory (Windows)
    status = __setwindows()
elif platform.system() == "Linux":
    # set ansys installation directory (Linux)
    status = __setlinux()
else:
    msg = [
        Color.RED,
        "unsupported platform",
        str(platform.system()),
        "\n",
        "PyChemkin does not support the current os.",
        Color.END,
    ]
    this_msg = Color.SPACE.join(msg)
    logger.critical(this_msg)
    exit()

# check set up status
if status != 0:
    exit()
# load Chemkin-CFD-API shared object
try:
    chemkin = cdll.LoadLibrary(_target_lib)
except OSError:
    inst_dir = Path(_ansys_dir) / "reaction" / _ckbin / "bin"
    msg = [
        Color.RED,
        "error initializing ansys-chemkin.\n",
        Color.SPACEx6,
        "please verify local Chemkin installation at",
        str(inst_dir),
        "\n",
        Color.SPACEx6,
        "run the chemkin set up script",
        "'source chemkin_setup.ksh' in the 'bin' directory.\n",
        Color.SPACEx6,
        "or check for a valid Ansys-chemkin license.",
        Color.END,
    ]
    this_msg = Color.SPACE.join(msg)
    logger.critical(this_msg)
    exit()

# Chemkin-CFD-API
# document: Chemkin-CFD-API Reference Guide (Ansys Help portal)
#
# syntax:
# Specify the return type of the function
# Specify the argument types for the functions
#
# general purpose functions
chemkin.KINSetUnitSystem.restype = ctypes.c_int
chemkin.KINSetUnitSystem.argtypes = [ctypes.POINTER(ctypes.c_int)]
# preprocess
chemkin.KINPreProcess.restype = ctypes.c_int
chemkin.KINPreProcess.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINInitialize.restype = ctypes.c_int
chemkin.KINInitialize.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINFinish.restype = None
chemkin.KINFinish.argtypes = []
chemkin.KINUpdateChemistrySet.restype = ctypes.c_int
chemkin.KINUpdateChemistrySet.argtypes = [
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINSwitchChemistrySet.restype = ctypes.c_int
chemkin.KINSwitchChemistrySet.argtypes = [
    ctypes.POINTER(ctypes.c_int),
]
# size, index, symbols
chemkin.KINGetChemistrySizes.restype = ctypes.c_int
chemkin.KINGetChemistrySizes.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINGetGasSpeciesNames.restype = ctypes.c_int
chemkin.KINGetGasSpeciesNames.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.POINTER(ctypes.c_char)),
]
chemkin.KINGetElementNames.restype = ctypes.c_int
chemkin.KINGetElementNames.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.POINTER(ctypes.c_char)),
]
chemkin.KINGetAtomicWeights.restype = ctypes.c_int
chemkin.KINGetAtomicWeights.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetGasMolecularWeights.restype = ctypes.c_int
chemkin.KINGetGasMolecularWeights.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetGasReactionString.restype = ctypes.c_int
chemkin.KINGetGasReactionString.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_char),
]
chemkin.KINGetReactionStringLength.restype = ctypes.c_int
chemkin.KINGetReactionStringLength.argtypes = [ctypes.POINTER(ctypes.c_int)]
# species information
chemkin.KINGetGasSpecificHeat.restype = ctypes.c_int
chemkin.KINGetGasSpecificHeat.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetGasSpeciesEnthalpy.restype = ctypes.c_int
chemkin.KINGetGasSpeciesEnthalpy.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetGasSpeciesInternalEnergy.restype = ctypes.c_int
chemkin.KINGetGasSpeciesInternalEnergy.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetGasSpeciesComposition.restype = ctypes.c_int
chemkin.KINGetGasSpeciesComposition.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),
]
chemkin.KINGetMassDensity.restype = ctypes.c_int
chemkin.KINGetMassDensity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
]

chemkin.KINGetViscosity.restype = ctypes.c_int
chemkin.KINGetViscosity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetConductivity.restype = ctypes.c_int
chemkin.KINGetConductivity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetDiffusionCoeffs.restype = ctypes.c_int
chemkin.KINGetDiffusionCoeffs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="F_CONTIGUOUS"),
]

chemkin.KINGetGasMixtureSpecificHeat.restype = ctypes.c_int
chemkin.KINGetGasMixtureSpecificHeat.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINGetGasMixtureEnthalpy.restype = ctypes.c_int
chemkin.KINGetGasMixtureEnthalpy.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
]

chemkin.KINGetMixtureViscosity.restype = ctypes.c_int
chemkin.KINGetMixtureViscosity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINGetMixtureConductivity.restype = ctypes.c_int
chemkin.KINGetMixtureConductivity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINGetMixtureDiffusionCoeffs.restype = ctypes.c_int
chemkin.KINGetMixtureDiffusionCoeffs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetOrdinaryDiffusionCoeffs.restype = ctypes.c_int
chemkin.KINGetOrdinaryDiffusionCoeffs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="F_CONTIGUOUS"),
]
chemkin.KINGetThermalDiffusionCoeffs.restype = ctypes.c_int
chemkin.KINGetThermalDiffusionCoeffs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
]
# reaction rate
chemkin.KINGetGasROP.restype = ctypes.c_int
chemkin.KINGetGasROP.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetGasReactionRates.restype = ctypes.c_int
chemkin.KINGetGasReactionRates.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetReactionRateParameters.restype = ctypes.c_int
chemkin.KINGetReactionRateParameters.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINSetAFactorForAReaction.restype = ctypes.c_int
chemkin.KINSetAFactorForAReaction.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
]
# gas-phase equilibrium calculation (limited capabilities)
chemkin.KINCalculateEquil.restype = ctypes.c_int
chemkin.KINCalculateEquil.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINCalculateEquilWithOption.restype = ctypes.c_int
chemkin.KINCalculateEquilWithOption.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINCalculateEqGasWithOption.restype = ctypes.c_int
chemkin.KINCalculateEqGasWithOption.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
# real-gas
chemkin.KINRealGas_SetParameter.restype = ctypes.c_int
chemkin.KINRealGas_SetParameter.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINRealGas_GetEOSMode.restype = ctypes.c_int
chemkin.KINRealGas_GetEOSMode.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_char),
]
chemkin.KINRealGas_SetMixingRule.restype = ctypes.c_int
chemkin.KINRealGas_SetMixingRule.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINRealGas_UseIdealGasLaw.restype = ctypes.c_int
chemkin.KINRealGas_UseIdealGasLaw.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINRealGas_UseCubicEOS.restype = ctypes.c_int
chemkin.KINRealGas_UseCubicEOS.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINRealGas_SetCurrentPressure.restype = ctypes.c_int
chemkin.KINRealGas_SetCurrentPressure.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINRealGas_CheckRealGasStatus.restype = ctypes.c_int
chemkin.KINRealGas_CheckRealGasStatus.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINGetGamma.restype = ctypes.c_int
chemkin.KINGetGamma.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
]
# Batch reactor interfaces
chemkin.KINAll0D_Setup.restype = ctypes.c_int
chemkin.KINAll0D_Setup.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINAll0D_SetupWorkArrays.restype = ctypes.c_int
chemkin.KINAll0D_SetupWorkArrays.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINAll0D_SetupBatchInputs.restype = ctypes.c_int
chemkin.KINAll0D_SetupBatchInputs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_SetupPSRReactorInputs.restype = ctypes.c_int
chemkin.KINAll0D_SetupPSRReactorInputs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_SetupPSRInletInputs.restype = ctypes.c_int
chemkin.KINAll0D_SetupPSRInletInputs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_SetupPFRInputs.restype = ctypes.c_int
chemkin.KINAll0D_SetupPFRInputs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_SetupHCCIInputs.restype = ctypes.c_int
chemkin.KINAll0D_SetupHCCIInputs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_SetupHCCIZoneInputs.restype = ctypes.c_int
chemkin.KINAll0D_SetupHCCIZoneInputs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINAll0D_SetupSIInputs.restype = ctypes.c_int
chemkin.KINAll0D_SetupSIInputs.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINAll0D_Calculate.restype = ctypes.c_int
chemkin.KINAll0D_Calculate.argtypes = [ctypes.POINTER(ctypes.c_int)]
chemkin.KINAll0D_CalculateInput.restype = ctypes.c_int
chemkin.KINAll0D_CalculateInput.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_SetUserKeyword.restype = ctypes.c_int
chemkin.KINAll0D_SetUserKeyword.argtypes = [ctypes.POINTER(ctypes.c_char)]
chemkin.KINAll0D_IntegrateHeatRelease.restype = ctypes.c_int
chemkin.KINAll0D_IntegrateHeatRelease.argtypes = []
chemkin.KINAll0D_SetHeatTransfer.restype = ctypes.c_int
chemkin.KINAll0D_SetHeatTransfer.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINAll0D_SetHeatTransferArea.restype = ctypes.c_int
chemkin.KINAll0D_SetHeatTransferArea.argtypes = [ctypes.POINTER(ctypes.c_double)]
# profile
chemkin.KINAll0D_SetProfilePoints.restype = ctypes.c_int
chemkin.KINAll0D_SetProfilePoints.argtypes = [ctypes.POINTER(ctypes.c_int)]
chemkin.KINAll0D_SetProfileParameter.restype = ctypes.c_int
chemkin.KINAll0D_SetProfileParameter.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_SetProfileKeyword.restype = ctypes.c_int
chemkin.KINAll0D_SetProfileKeyword.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
# batch reactor solver parameters
chemkin.KINAll0D_SetSolverInitialStepTime.restype = ctypes.c_int
chemkin.KINAll0D_SetSolverInitialStepTime.argtypes = [ctypes.POINTER(ctypes.c_double)]
chemkin.KINAll0D_SetSolverMaximumStepTime.restype = ctypes.c_int
chemkin.KINAll0D_SetSolverMaximumStepTime.argtypes = [ctypes.POINTER(ctypes.c_double)]
chemkin.KINAll0D_SetSolverMaximumIteration.restype = ctypes.c_int
chemkin.KINAll0D_SetSolverMaximumIteration.argtypes = [ctypes.POINTER(ctypes.c_int)]
chemkin.KINAll0D_SetRelaxIteration.restype = ctypes.c_int
chemkin.KINAll0D_SetRelaxIteration.argtypes = []
chemkin.KINAll0D_SetMinimumSpeciesBound.restype = ctypes.c_int
chemkin.KINAll0D_SetMinimumSpeciesBound.argtypes = [ctypes.POINTER(ctypes.c_double)]
# get solution
chemkin.KINAll0D_GetSolution.restype = ctypes.c_int
chemkin.KINAll0D_GetSolution.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_GetSolnResponseSize.restype = ctypes.c_int
chemkin.KINAll0D_GetSolnResponseSize.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINAll0D_GetGasSolnResponse.restype = ctypes.c_int
chemkin.KINAll0D_GetGasSolnResponse.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="F_CONTIGUOUS"),
]
chemkin.KINAll0D_GetIgnitionDelay.restype = ctypes.c_int
chemkin.KINAll0D_GetIgnitionDelay.argtypes = [ctypes.POINTER(ctypes.c_double)]
chemkin.KINAll0D_GetHeatRelease.restype = ctypes.c_int
chemkin.KINAll0D_GetHeatRelease.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINAll0D_GetEngineHeatRelease.restype = ctypes.c_int
chemkin.KINAll0D_GetEngineHeatRelease.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINAll0D_GetExitMassFlowRate.restype = ctypes.c_int
chemkin.KINAll0D_GetExitMassFlowRate.argtypes = [ctypes.POINTER(ctypes.c_double)]
# Premix interfaces
chemkin.KINPremix_SetParameter.restype = ctypes.c_int
chemkin.KINPremix_SetParameter.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINPremix_CalculateFlame.restype = ctypes.c_int
chemkin.KINPremix_CalculateFlame.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINPremix_GetSolution.restype = ctypes.c_int
chemkin.KINPremix_GetSolution.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="F_CONTIGUOUS"),
]
chemkin.KINPremix_GetSolutionGridPoints.restype = ctypes.c_int
chemkin.KINPremix_GetSolutionGridPoints.argtypes = [
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINPremix_GetFlameMassFlux.restype = ctypes.c_int
chemkin.KINPremix_GetFlameMassFlux.argtypes = [
    ctypes.POINTER(ctypes.c_double),
]
# Oppdif interfaces
chemkin.KINOppdif_SetInlet.restype = ctypes.c_int
chemkin.KINOppdif_SetParameter.restype = ctypes.c_int
chemkin.KINOppdif_CalculateFlame.restype = ctypes.c_int
chemkin.KINOppdif_GetSolutionGridPoints.restype = ctypes.c_int
chemkin.KINOppdif_GetSolution.restype = ctypes.c_int
chemkin.KINOppdif_SetInlet.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINOppdif_SetParameter.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_double),
]

chemkin.KINOppdif_CalculateFlame.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINOppdif_GetSolutionGridPoints.argtypes = [ctypes.POINTER(ctypes.c_int)]
chemkin.KINOppdif_GetSolution.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),
]  # np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='C')]

chemkin.KINOppdif_GetSolnSpeciesIntegratedROP.restype = ctypes.c_int
chemkin.KINOppdif_GetSolnSpeciesIntegratedROP.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),
]

chemkin.KINGetMassFractionFromMoleFraction.restype = ctypes.c_int
chemkin.KINGetMassFractionFromMoleFraction.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]

chemkin.KINGetMoleFractionFromMassFraction.restype = ctypes.c_int
chemkin.KINGetMoleFractionFromMassFraction.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]

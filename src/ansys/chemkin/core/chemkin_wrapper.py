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
import os
from pathlib import Path
import platform
from typing import cast

import numpy as np

from ansys.chemkin.core._platform_setup import (
    __setup_linux,
    __setup_windows,
    find_valid_ansys_versions,
    get_runtime_diagnostics,
)
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.logger import logger

# main procedure start here:
# The first part is to set up the Chemkin-CFD-API run environment
# for the latest Ansys Chemkin version installed on the local machine,
# including the full paths to the third-party shared objects and
# their dependencies.
# The second part is to load the Chemkin-CFD-API shared object according to
# the platform and Ansys version.
# The third part is to set up the function prototypes.
#
# set ansys version numbers
# the oldest Ansys version that PyChemkin can run with is 25.1 (2025 R1)
_min_version = 251
# the Ansys version number that PyChemkin is currently running with
_ansys_ver = _min_version
# the Ansys installation directory
_ansys_dir = ""
# the Chemkin-CFD-API bin directory
_ckbin = ""
# name of the Chemkin-CFD-API shared object
_target_lib = ""
# path to the Chemkin-CFD-API shared object and its dependencies
_lib_paths: list[str] = []
# pychemkin configuration parameters
pyck_config = (_ansys_ver, _ansys_dir, _ckbin, _target_lib, _lib_paths)
# create log
msg = ["minimum Ansys version to run PyChemkin =", str(_min_version)]
this_msg = Color.SPACE.join(msg)
logger.debug(this_msg)
# find all the valid Ansys versions installed on the local machine
_valid_versions = find_valid_ansys_versions(_min_version)
# check the environment variable "PYCK_CHEMKIN_VER"
# if a specific Ansys version is requested
_requested_ver = int(os.getenv("PYCK_CHEMKIN_VER", "0"))
if _requested_ver in _valid_versions:
    _valid_versions = [_requested_ver]
else:
    if _requested_ver != 0:
        msg = [
            Color.RED,
            "requested Ansys version not found or not valid:",
            str(_requested_ver),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.critical(this_msg)
        exit()
# check os platform
if platform.system() == "Windows":
    # set pychemkin configuration for Windows
    status, pyck_config = __setup_windows(
        min_ver=_min_version, valid_vers=_valid_versions
    )
elif platform.system() == "Linux":
    # set pychemkin configuration for Linux
    status, pyck_config = __setup_linux(
        min_ver=_min_version, valid_vers=_valid_versions
    )
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

# check environment setup status
if status != 0:
    exit()
# unpack the configuration parameters
_ansys_ver = pyck_config[0]
_ansys_dir = pyck_config[1]
_ckbin = pyck_config[2]
_target_lib = pyck_config[3]
_lib_paths = cast(list[str], pyck_config[4])
#
# load Chemkin-CFD-API shared object
if not Path(_target_lib).is_file():
    msg = [
        Color.RED,
        "Chemkin-CFD-API shared object not found: ",
        str(_target_lib),
        Color.END,
    ]
    this_msg = Color.SPACE.join(msg)
    logger.critical(this_msg)
    exit()
else:
    try:
        chemkin = ctypes.CDLL(_target_lib)
        # set an environment variable for the Chemkin-CFD-API to know
        # the calling application is PyChemkin
        os.environ["RD_PY_CHEMKIN"] = "1"
    except OSError as e:
        # error and diagnostic messages
        inst_dir = Path(_ansys_dir) / "reaction" / _ckbin / "bin"
        error_message = e.strerror or str(e)
        msg = [
            Color.RED,
            "error initializing ansys-chemkin.",
            Color.END,
        ]
        logger.critical(Color.SPACE.join(msg))
        # details of the error
        details = [
            Color.SPACEx6,
            "target shared object:",
            str(_target_lib),
            "\n",
            Color.SPACEx6,
            "error message:",
            error_message,
        ]
        if e.filename:
            details.extend(["\n", Color.SPACEx6, "reported file:", str(e.filename)])
        this_msg = Color.SPACE.join(details)
        logger.critical(this_msg)

        if _lib_paths:
            search_paths = [Color.SPACEx6, "configured shared-library search paths:"]
            for lib_path in _lib_paths:
                search_paths.extend(["\n", Color.SPACEx6, "-", lib_path])
            this_msg = Color.SPACE.join(search_paths)
            logger.critical(this_msg)

        diagnostics = get_runtime_diagnostics(config=pyck_config, probe_load=True)
        unresolved_imports = diagnostics.get("dependency_probes", [])
        if unresolved_imports:
            imports_msg = [Color.SPACEx6, "unresolved direct DLL imports detected:"]
            for dep_probe in unresolved_imports:
                dep_name = str(dep_probe.get("name", "<unknown>"))
                dep_error = str(dep_probe.get("search_path_load_error", ""))
                imports_msg.extend(
                    ["\n", Color.SPACEx6, "-", dep_name, "->", dep_error]
                )
                found_paths = dep_probe.get("found_in_configured_paths", [])
                if found_paths:
                    imports_msg.extend(
                        ["\n", Color.SPACEx6, "  located at:", str(found_paths[0])]
                    )
                full_probe = dep_probe.get("full_path_probe")
                if isinstance(full_probe, dict):
                    imports_msg.extend(
                        [
                            "\n",
                            Color.SPACEx6,
                            "  full-path probe:",
                            str(full_probe.get("error", "")),
                        ]
                    )
            this_msg = Color.SPACE.join(imports_msg)
            logger.critical(this_msg)
        elif platform.system() == "Windows":
            msg = [
                Color.YELLOW,
                "direct DLL dependency diagnostics unavailable or inconclusive.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

        msg = [
            Color.RED,
            "please verify local Chemkin installation at",
            str(inst_dir),
            "\n",
            Color.SPACEx6,
            "or check for any missing shared libraries as indicated above.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.critical(this_msg)
        intel_link = "https://www.intel.com/content/www/us/en/developer/articles/tool"
        intel_link += "/compilers-redistributable-libraries-by-version.html"
        msg = [
            Color.YELLOW,
            "For Intel compiler runtime libraries related issues,",
            "please install Intel® oneAPI DPC++/C++",
            "and Fortran Compilers Runtime Libraries",
            "at the link below:\n",
            intel_link,
            Color.END,
        ]
        if platform.system() == "Linux":
            msg = [
                Color.YELLOW,
                "For Linux,",
                "try to run the chemkin set up script",
                "'source chemkin_setup.ksh' in the 'bin' directory.\n",
                Color.SPACEx6,
                "before retrying to initialize ansys-chemkin.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        exit()

# Chemkin-CFD-API prototypes
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
# (2027 R1 feature)
# chemkin.KINExit.restype = ctypes.c_int
# chemkin.KINExit.argtypes = []
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

# surface chemistry
chemkin.KINGetMaterialSizes.restype = ctypes.c_int
chemkin.KINGetMaterialSizes.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINGetMaterialIndex.restype = ctypes.c_int
chemkin.KINGetMaterialIndex.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINGetPhaseSizes.restype = ctypes.c_int
chemkin.KINGetPhaseSizes.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINGetPhaseSpeciesIndices.restype = ctypes.c_int
chemkin.KINGetPhaseSpeciesIndices.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINGetSurfaceSpeciesComposition.restype = ctypes.c_int
chemkin.KINGetSurfaceSpeciesComposition.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),
]
chemkin.KINGetSpeciesOccupancy.restype = ctypes.c_int
chemkin.KINGetSpeciesOccupancy.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
]
chemkin.KINGetSiteDensity.restype = ctypes.c_int
chemkin.KINGetSiteDensity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetAllSpeciesEnthalpy.restype = ctypes.c_int
chemkin.KINGetAllSpeciesEnthalpy.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetSiteMolecularWeights.restype = ctypes.c_int
chemkin.KINGetSiteMolecularWeights.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetBulkDensity.restype = ctypes.c_int
chemkin.KINGetBulkDensity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetBulkMolecularWeights.restype = ctypes.c_int
chemkin.KINGetBulkMolecularWeights.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetBulkSpeciesEnthalpy.restype = ctypes.c_int
chemkin.KINGetBulkSpeciesEnthalpy.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetBulkSpeciesSpecificHeat.restype = ctypes.c_int
chemkin.KINGetBulkSpeciesSpecificHeat.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetBulkSpeciesEntropy.restype = ctypes.c_int
chemkin.KINGetBulkSpeciesEntropy.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetSurfaceReactionRates.restype = ctypes.c_int
chemkin.KINGetSurfaceReactionRates.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetSurfaceProductionRates.restype = ctypes.c_int
chemkin.KINGetSurfaceProductionRates.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINGetSurfaceReactionRateParameters.restype = ctypes.c_int
chemkin.KINGetSurfaceReactionRateParameters.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
]
chemkin.KINSetAFactorForASurfaceReaction.restype = ctypes.c_int
chemkin.KINSetAFactorForASurfaceReaction.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINSetActEnergyForASurfaceReaction.restype = ctypes.c_int
chemkin.KINSetActEnergyForASurfaceReaction.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINGetMaterialNames.restype = ctypes.c_int
chemkin.KINGetMaterialNames.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.POINTER(ctypes.c_char)),
]
chemkin.KINGetPhaseNames.restype = ctypes.c_int
chemkin.KINGetPhaseNames.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.POINTER(ctypes.c_char)),
]
chemkin.KINGetSiteNames.restype = ctypes.c_int
chemkin.KINGetSiteNames.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.POINTER(ctypes.c_char)),
]
chemkin.KINGetBulkNames.restype = ctypes.c_int
chemkin.KINGetBulkNames.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.POINTER(ctypes.c_char)),
]
chemkin.KINGetSurfaceReactionString.restype = ctypes.c_int
chemkin.KINGetSurfaceReactionString.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_char),
]

# mixture information
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
chemkin.KINGetGasMixtureThermicity.restype = ctypes.c_int
chemkin.KINGetGasMixtureThermicity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINGetGasMixtureSoundSpeed.restype = ctypes.c_int
chemkin.KINGetGasMixtureSoundSpeed.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINGetGasPressure.restype = ctypes.c_int
chemkin.KINGetGasPressure.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
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
chemkin.KINRealGas_GetMixtureThermicity.restype = ctypes.c_int
chemkin.KINRealGas_GetMixtureThermicity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
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
chemkin.KINAll0D_SetMaxInletSize.restype = ctypes.c_int
chemkin.KINAll0D_SetMaxInletSize.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]

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
chemkin.KINAll0D_SetTransientSolver.restype = ctypes.c_int
chemkin.KINAll0D_SetTransientSolver.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
]

# get solution
chemkin.KINAll0D_GetSolution.restype = ctypes.c_int
chemkin.KINAll0D_GetSolution.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_GetSurfaceSolution.restype = ctypes.c_int
chemkin.KINAll0D_GetSurfaceSolution.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_GetBulkGrowthRate.restype = ctypes.c_int
chemkin.KINAll0D_GetBulkGrowthRate.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_GetSolution_perPSR.restype = ctypes.c_int
chemkin.KINAll0D_GetSolution_perPSR.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_GetSurfaceSolution_perPSR.restype = ctypes.c_int
chemkin.KINAll0D_GetSurfaceSolution_perPSR.argtypes = [
    ctypes.POINTER(ctypes.c_int),
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
chemkin.KINAll0D_GetSurfaceSolnResponse.restype = ctypes.c_int
chemkin.KINAll0D_GetSurfaceSolnResponse.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="F_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="F_CONTIGUOUS"),
]
chemkin.KINAll0D_GetSurfaceTemperature.restype = ctypes.c_int
chemkin.KINAll0D_GetSurfaceTemperature.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_GetGasHeatReleaseRate.restype = ctypes.c_int
chemkin.KINAll0D_GetGasHeatReleaseRate.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_GetSurfaceHeatReleaseRate.restype = ctypes.c_int
chemkin.KINAll0D_GetSurfaceHeatReleaseRate.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINAll0D_GetSurfaceHeatLossRate.restype = ctypes.c_int
chemkin.KINAll0D_GetSurfaceHeatLossRate.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
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
chemkin.KINAll0D_GetExitMassFlowRate_perPSR.restype = ctypes.c_int
chemkin.KINAll0D_GetExitMassFlowRate_perPSR.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINAll0D_GetEngineIMEP.restype = ctypes.c_int
chemkin.KINAll0D_GetEngineIMEP.argtypes = [ctypes.POINTER(ctypes.c_double)]
chemkin.KINAll0D_GetSIEngineKnockCA.restype = ctypes.c_int
chemkin.KINAll0D_GetSIEngineKnockCA.argtypes = [ctypes.POINTER(ctypes.c_double)]
chemkin.KINAll0D_SuppressOutput.restype = ctypes.c_int
chemkin.KINAll0D_SuppressOutput.argtypes = []

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
chemkin.KINOppdif_SetInlet.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINOppdif_SetStagWallTemp.restype = ctypes.c_int
chemkin.KINOppdif_SetStagWallTemp.argtypes = [
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINOppdif_SetParameter.restype = ctypes.c_int
chemkin.KINOppdif_SetParameter.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINOppdif_CalculateFlame.restype = ctypes.c_int
chemkin.KINOppdif_CalculateFlame.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINOppdif_GetSolutionGridPoints.restype = ctypes.c_int
chemkin.KINOppdif_GetSolutionGridPoints.argtypes = [
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINOppdif_GetSolution.restype = ctypes.c_int
chemkin.KINOppdif_GetSolution.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="F_CONTIGUOUS"),
]
chemkin.KINOppdif_GetGlobalStrainRate.restype = ctypes.c_int
chemkin.KINOppdif_GetGlobalStrainRate.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINOppdif_GetMaxTemperature.restype = ctypes.c_int
chemkin.KINOppdif_GetMaxTemperature.argtypes = [
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINOppdif_GetMixtureFraction.restype = ctypes.c_int
chemkin.KINOppdif_GetMixtureFraction.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINOppdif_GetAxialVelocity.restype = ctypes.c_int
chemkin.KINOppdif_GetAxialVelocity.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINOppdif_GetStrainRates.restype = ctypes.c_int
chemkin.KINOppdif_GetStrainRates.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINOppdif_GetVelocityField.restype = ctypes.c_int
chemkin.KINOppdif_GetVelocityField.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINOppdif_GetSolnSpeciesIntegratedROP.restype = ctypes.c_int
chemkin.KINOppdif_GetSolnSpeciesIntegratedROP.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
# shock tube interfaces
chemkin.KINShock_CalcIncidentShockWithBoundaryLayerCorrection.restype = ctypes.c_int
chemkin.KINShock_CalcIncidentShockWithBoundaryLayerCorrection.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINShock_CalcIncidentShockWithoutBoundaryLayerCorrection.restype = ctypes.c_int
chemkin.KINShock_CalcIncidentShockWithoutBoundaryLayerCorrection.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINShock_CalcZNDWithBoundaryLayerCorrection.restype = ctypes.c_int
chemkin.KINShock_CalcZNDWithBoundaryLayerCorrection.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINShock_CalcZNDWithoutBoundaryLayerCorrection.restype = ctypes.c_int
chemkin.KINShock_CalcZNDWithoutBoundaryLayerCorrection.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINShock_CalcReflectedShock.restype = ctypes.c_int
chemkin.KINShock_CalcReflectedShock.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINShock_SetParameter.restype = ctypes.c_int
chemkin.KINShock_SetParameter.argtypes = [
    ctypes.POINTER(ctypes.c_char),
    ctypes.POINTER(ctypes.c_double),
]
chemkin.KINShock_GetSolution.restype = ctypes.c_int
chemkin.KINShock_GetSolution.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINShock_GetSolnResponseSize.restype = ctypes.c_int
chemkin.KINShock_GetSolnResponseSize.argtypes = [
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINShock_GetGasSolnResponse.restype = ctypes.c_int
chemkin.KINShock_GetGasSolnResponse.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="F_CONTIGUOUS"),
]
chemkin.KINShock_GetInductLengthSize.restype = ctypes.c_int
chemkin.KINShock_GetInductLengthSize.argtypes = [
    ctypes.POINTER(ctypes.c_int),
]
chemkin.KINShock_GetInductLengths.restype = ctypes.c_int
chemkin.KINShock_GetInductLengths.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.double, flags="C_CONTIGUOUS"),
]
chemkin.KINShock_GetGasStates.restype = ctypes.c_int
chemkin.KINShock_GetGasStates.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]

# others
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

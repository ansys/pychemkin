import ctypes
from ctypes import cdll
import logging
import os
import platform
import sys

import numpy as np

logger = logging.getLogger(__name__)

# current_dir = os.getcwd()
# logger.debug("Current working directory =", current_dir)

# set ansys version number
_min_version = 251
_valid_versions = [252, 251, 242]
_ansys_ver = _min_version

# check platform
if platform.system() == "Windows" or platform.system() == "Linux":
    # set ansys installation directory (Windows)
    for v in _valid_versions:
        _ansys_ver = v
        if v >= _min_version:
            _ansys_installation = "ANSYS" + str(_ansys_ver) + "_DIR"
            _ansys_home = os.environ.get(_ansys_installation, "NA")
            print(_ansys_installation + ": " + _ansys_home)
            if _ansys_home != "NA":
                _ansys_dir = os.path.dirname(_ansys_home)
                break
        else:
            break
else:
    print(f"** unsupported platform {platform.system()}")
    print("** PyChemkin does not support the current os")
    exit()

if _ansys_ver >= _min_version:
    if not os.path.isdir(_ansys_dir):
        print(f"** PyChemkin cannot find the specific Ansys installation: {_ansys_dir}")
        exit()
    else:
        plat = "winx64"
        i64 = "/intel64"
        if platform.system() == "Windows":
            plat = "linx64/lib/intel64"
            _lib_paths = [
                os.path.join(_ansys_dir, "tp", "IntelCompiler", "2023.1.0", "win64"),
                os.path.join(_ansys_dir, "tp", "IntelMKL", "2023.1.0", "win64"),
                os.path.join(_ansys_dir, "tp", "zlib", "1.2.13", "win64")
            ]
        else:
            _lib_paths = [
                os.path.join(_ansys_dir, "tp", "IntelCompiler", "2023.1.0", "linx64", "lib", "intel64"),
                os.path.join(_ansys_dir, "tp", "IntelMKL", "2023.1.0", "linx64", "lib", "intel64"),
                os.path.join(_ansys_dir, "tp", "zlib", "1.2.13", "linx64", "lib")
            ]
else:
    print("** PyChemkin does not support Chemkin versions older than 2025R1")
    exit()

if sys.platform == "win32":
    for _lib_path in _lib_paths:
        os.add_dll_directory(_lib_path)
else:
    combined_path = ":".join(_lib_paths)
    if "LD_LIBRARY_PATH" not in os.environ.keys():
    #if os.environ["LD_LIBRARY_PATH"] is None:
        os.environ["LD_LIBRARY_PATH"] = combined_path
    else:
        os.environ["LD_LIBRARY_PATH"] = (
            os.environ["LD_LIBRARY_PATH"] + ":" + combined_path
        )
    if "PATH" not in os.environ.keys():
        os.environ["PATH"] = combined_path
    else:
        os.environ["PATH"] = os.environ["PATH"] + ":" + combined_path

# load KINetics package/module
try:
    target_path = None
    if sys.platform == "win32":
        target_lib = os.path.join(_ansys_dir, "reaction", "chemkin.win64", "bin", "KINeticsdll.dll")
    else:
        target_lib = os.path.join(_ansys_dir, "reaction", "chemkin.linuxx8664", "bin", "libKINetics.so")
    chemkin = cdll.LoadLibrary(target_lib)
except OSError:
    inst_dir = os.path.join(_ansys_dir, "reaction", "chemkin.win64", "bin")
    print("** error initializing ansys-chemkin")
    print("** please check Chemkin installation at " + inst_dir)
    print("** or check for a valid Ansys-chemkin license")
    exit()
except AttributeError:
    inst_dir = os.path.join(_ansys_dir, "reaction", "chemkin.win64", "bin")
    print("** error initializing ansys-chemkin")
    print("** please check Chemkin installation at " + inst_dir)
    print("** or check for a valid Ansys-chemkin license")
    exit()
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
    ctypes.POINTER(ctypes.c_int),
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

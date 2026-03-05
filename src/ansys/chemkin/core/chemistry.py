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

"""Chemkin Chemistry utilities."""

import copy
import ctypes
from ctypes import POINTER, c_char_p, c_double, c_int
from pathlib import Path
from typing import Dict, List, Union

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core import chemkin_wrapper as ck_wrapper
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.constants import P_ATM, R_GAS
from ansys.chemkin.core.info import clear_hints
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.realgaseos import check_realgas_status, set_current_pressure

_symbol_length = 16  # Chemkin element/species symbol length
MAX_SPECIES_LENGTH = _symbol_length + 1  # Chemkin element/species symbol length + 1
LP_c_char = ctypes.POINTER(ctypes.c_char)  # pointer to C type character array
COMPLETE = 0

# string used to identify different chemistry sets in the same project
_chemset_identifiers: List = []
_active_chemistry_set = -10
# verbose mode to turn ON/OFF the print statements that
# do not have the leading '**' character
chemkin_verbose = True
# Chemkin-CFD-API initialization flagfor every Chemistry Set
_CKInitialized: Dict = {}
# == end of global parameters


#
# Chemkin module level methods
#
def verbose() -> bool:
    """Get current Pychemkin verbose mode."""
    """Return the global verbose mode indicating the status (ON/OFF)
    of printing statements that do not have the leading '**' characters.

    Returns
    -------
        mode: boolean, {True, False}, default = True
            the verbose mode

    """
    global chemkin_verbose
    return chemkin_verbose


def set_verbose(onoff: bool):
    """Set Pychemkin verbose mode."""
    """Set the global verbose mode to turn ON(True) or OFF(False) of
    printing statements that do not have the leading '**' characters.

    Parameters
    ----------
        onoff: boolean, {True, False}
            the verbose mode

    """
    global chemkin_verbose
    chemkin_verbose = onoff


def chemkin_version() -> int:
    """Return the Chemkin-CFD-API version number currently in use."""
    """
    Returns
    -------
        version: integer
            Chemkin-CFD-API version number

    """
    return ck_wrapper._ansys_ver


def verify_version(min_version: int) -> bool:
    """Check the  Chemkin-CFD-API vwersion."""
    """Check if the version of Chemkin-CFD-API currently in use meets
    the minimum version required by certain operations.

    Parameters
    ----------
        min_version: integer
            minimum chemkin-CFD-API version required to perform the operation

    Returns
    -------
        status: boolean

    """
    status = chemkin_version() >= min_version
    if not status:
        msg = [
            Color.PURPLE,
            "this operation is NOT supported by the current chemkin version",
            str(chemkin_version()),
            "\n",
            "the minimum chemkin version required for this operation is",
            str(min_version),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
    return status


def chemkin_bin_dir() -> str:
    """Return the local Ansys Chemkin bin directory."""
    """
    Return the local Ansys Chemkin bin directory
    currently used by Pychemkin.

    Returns
    -------
        chemkin_bin_dir: string
            the bin directory of local Ansys Chemkin installation
    """
    local_ansys_dir = Path(ck_wrapper._ansys_dir)
    chemkin_bin_dir = local_ansys_dir / "reaction" / ck_wrapper._ckbin / "bin"
    return str(chemkin_bin_dir)


def ansys_dir() -> str:
    """Return the local Ansys installation."""
    """
    Return the local Ansys installation directory currently in use.

    Returns
    -------
        ansys_dir: string
            the local Ansys installation directory
    """
    local_ansys_dir = ck_wrapper._ansys_dir
    return local_ansys_dir


def done():
    """Release Chemkin license and reset the Chemistry sets."""
    # terminate client (2027 R1 feature)
    # if check_jupyter_notebook():
    # running in Jupyter environment
    # requires addition steps
    #     _ = ck_wrapper.chemkin.KINExit()

    # clean up
    global _active_chemistry_set
    _active_chemistry_set = -10
    global _chemset_identifiers
    _chemset_identifiers.clear()
    global _CKInitialized
    _CKInitialized.clear()
    clear_hints()
    # reset
    global COMPLETE
    COMPLETE = 0
    global chemkin_verbose
    chemkin_verbose = True
    msg = [
        Color.GREEN,
        "done!\n",
        ">>> Chemkin-CFD-API stopped <<<",
        Color.END,
    ]
    this_msg = Color.SPACE.join(msg)
    logger.info(this_msg)


# utilities
def check_chemistryset(chem_index: int) -> bool:
    """Check whether the Chemistry Set is initialized in Chemkin-CFD-API."""
    """
    Check whether the Chemistry Set is initialized in Chemkin-CFD-API.

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set

    Returns
    -------
        status: boolean
            the initialization status of the Chemistry set associated with
            the given Chemistry set index

    """
    global _CKInitialized
    status = _CKInitialized.get(chem_index, False)
    return status


def activate_chemistryset(chem_index: int) -> int:
    """Switch to (re-activate) the Chemistry Set."""
    """Switch to (re-activate) the work spaces of the current Chemistry Set
    when there are multiple Chemistry Sets in the same project.

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set

    Returns
    -------
        error flag: integer

    """
    ierr = ck_wrapper.chemkin.KINSwitchChemistrySet(c_int(chem_index))
    if ierr == 0:
        # mark this chemistry set as active
        global _active_chemistry_set
        _active_chemistry_set = chem_index
    else:
        # failed to activate this chemistry set
        msg = [
            Color.PURPLE,
            "failed to reactivate the Chemistry Set work spaces.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
    return ierr


def force_activate_chemistryset(chem_index: int):
    """Activate the Chemistry Set automatically and silently."""
    """
    Activate the Chemistry Set automatically and silently.

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set

    """
    if not check_active_chemistryset(chem_index):
        # the Chemistry Set is not currently active
        ierr = activate_chemistryset(chem_index)
        if ierr != 0:
            exit()


def chemistryset_new(chem_index: int):
    """Create a new Chemistry Set."""
    """Create a new Chemistry Set initialization flag and
    set the value to False.

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set

    """
    global _CKInitialized
    _CKInitialized[chem_index] = False
    logger.debug("preprocessing done")


def chemistryset_initialized(chem_index: int):
    """Set the Chemistry Set Initialization flag to True."""
    """
    Set the Chemistry Set Initialization flag to True.

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set

    """
    global _CKInitialized
    _CKInitialized[chem_index] = True
    logger.debug(">>> Chemkin-CFD-API initialized <<<")


def check_active_chemistryset(chem_index: int) -> bool:
    """Verify if the chemistry set is currently activated."""
    """
    Verify if the chemistry set is currently activated.

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with a Chemistry Set

    Returns
    -------
        status: boolean
            active status of the Chemistry Set

    """
    global _active_chemistry_set
    return _active_chemistry_set == chem_index


def verify_chemset_surface(chem_set_index: int) -> bool:
    """Check the chemistry set has surface chemistry data."""
    """
    Check the gas chemistry sizes associated with the given chemistry set index.
    Verify the chemistry set contains surface chemistry information/data.

    Parameters
    ----------
        chem_set_index: integer
            chemistry set index

    Returns
    -------
        has_surface: bool
            True = the chemistry set with index chem_set_index
            contains surface chemistry
    """
    # initialization
    chemset_id = c_int(chem_set_index)
    n_elements = c_int(0)
    n_gas_species = c_int(0)
    n_gas_reactions = c_int(0)
    n_materials = c_int(0)
    n_sites = c_int(0)
    n_bulks = c_int(0)
    n_phases = c_int(0)
    n_surface_reactions = c_int(0)
    # get chemistry set sizes
    ierr = ck_wrapper.chemkin.KINGetChemistrySizes(
        chemset_id,
        n_elements,
        n_gas_species,
        n_gas_reactions,
        n_materials,
        n_sites,
        n_bulks,
        n_phases,
        n_surface_reactions,
    )
    if ierr != 0:
        msg = [
            Color.PURPLE,
            "trouble getting size information from the chemistry set:",
            Color.SPACEx6,
            "the error code =",
            str(ierr),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        return False
    else:
        n_surface_species = n_sites.value + n_bulks.value
        if n_materials.value > 0 and n_surface_species > 0:
            return True
        else:
            return False


def verify_chemset_sizes(
    chem_set_index: int,
    num_elements: int = -1,
    num_gas_species: int = -1,
    num_gas_reactions: int = -1,
) -> bool:
    """Verify the gas chemistry sizes."""
    """
    Check the chemistry set given by the index has the gas chemistry sizes
    that match with the sizes given by the parameters of the call.

    Parameters
    ----------
        chem_set_index: integer
            chemistry set index
        num_elements: integer, optional
            number of elements in the gas chemistry
        num_gas_species: integer, optional
            number of gas species in the gas chemistry
        num_gas_reactions: integer, optional
            number of gas reactions in the gas chemistry

    Returns
    -------
        match: bool
            True = the chemistry set with index chem_set_index
            matches the given sizes
    """
    # check inputs
    has_element = num_elements > 0
    has_species = num_gas_species > 0
    has_reaction = num_gas_reactions > 0
    if not (has_element or has_species or has_reaction):
        msg = [
            Color.PURPLE,
            "must provide one of the following chemistry sizes:\n",
            Color.SPACEx6,
            "number of elements, number of gas species,",
            "or number of gas reactions.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.warning(this_msg)
        return False
    # initialization
    chemset_id = c_int(chem_set_index)
    n_elements = c_int(0)
    n_gas_species = c_int(0)
    n_gas_reactions = c_int(0)
    n_materials = c_int(0)
    n_sites = c_int(0)
    n_bulks = c_int(0)
    n_phases = c_int(0)
    n_surface_reactions = c_int(0)
    # get chemistry set sizes
    ierr = ck_wrapper.chemkin.KINGetChemistrySizes(
        chemset_id,
        n_elements,
        n_gas_species,
        n_gas_reactions,
        n_materials,
        n_sites,
        n_bulks,
        n_phases,
        n_surface_reactions,
    )
    if ierr != 0:
        msg = [
            Color.PURPLE,
            "trouble getting size information from the chemistry set:",
            Color.SPACEx6,
            "the error code =",
            str(ierr),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        return False
    else:
        match = False
        if has_element:
            # check number of elements
            if n_elements.value != num_elements:
                msg = [
                    Color.YELLOW,
                    "number of elements does not match\n",
                    Color.SPACEx6,
                    "expected =",
                    str(num_elements),
                    "   chemistry set =",
                    str(n_elements.value),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
                match = False
        if has_species:
            # check number of gas species
            if n_gas_species.value != num_gas_species:
                msg = [
                    Color.YELLOW,
                    "number of gas species does not match\n",
                    Color.SPACEx6,
                    "expected =",
                    str(num_gas_species),
                    "   chemistry set =",
                    str(n_gas_species.value),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
                match = False
        if has_reaction:
            # check number of gas phase reactions
            if n_gas_reactions.value != num_gas_reactions:
                msg = [
                    Color.YELLOW,
                    "number of gas reactions does not match\n",
                    Color.SPACEx6,
                    "expected =",
                    str(num_gas_reactions),
                    "   chemistry set =",
                    str(n_gas_reactions.value),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
                match = False
        return match


def set_temp_array(
    temp: float,
    temp_ele: Union[float, None] = None,
    temp_ion: Union[float, None] = None,
):
    """Set up the temperature array for Chemkin-CFD-API calls."""
    """
    Set up the temperature array to be used in some Chemkin-CFD-API calls.
    The temperature array must be of c_double type of length 3 and consist of
    temperatures of the neutral species, the electron, and the ion species.

    Parameters
    ----------
        temp: double
            temperature of the neutral gas species [K]
        temp_ele: double, optional
            electron temperature [K]
        temp_ion: double, optional
            temperature of the ion/charged species [K].
            Usually ions have the same tempeature as the neutral species

    Returns
    -------
        t_array: 1-D c_double array of length 3
            temperature array containing the temperatures of
            the three different types of species
    """
    # set electron temperature [K]
    if temp_ele is None:
        temp_e = temp
    else:
        temp_e = temp_ele
    # set ion temperature [K]
    if temp_ion is None:
        temp_i = temp
    else:
        temp_i = temp_ion
    # set the temperature used in Chemkin-CFD-API calls
    t = np.ascontiguousarray([temp, temp_e, temp_i], dtype=np.float64)
    t_array = t.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    return t_array


class Chemistry:
    """Define and preprocess Chemkin chemistry set."""

    realgas_cueos = [
        "ideal gas",
        "Van der Waals",
        "Redlich-Kwong",
        "Soave",
        "Aungier",
        "Peng-Robinson",
    ]
    realgas_mixing_rules = ["Van der Waals", "pseudocritical"]

    def __init__(
        self,
        chem: str = "",
        surf: str = "",
        therm: str = "",
        tran: str = "",
        label: str = "",
    ):
        """Create a Chemistry object."""
        """Create a Chemistry object based on given Chemkin mechanism input files,
        thermodynamic data file, and transport data file.

        Parameters
        ----------
            chem: string, optional
                Full path and name of the Chemkin gas-phase mechanism input file
            surf: string, optional
                Full path and name of the Chemkin surface mechanism input file
            therm: string, optional
                Full path and name of the Chemkin thermodynamic data file
            tran: string, optional
                Full path and name of the Chemkin transport data file
            label: string, optional
                label/name of the chemistry set

        """
        # set flags
        self._index_surf = c_int(0)
        self._index_tran = c_int(0)
        # initialization
        self._chemset_index = c_int(-1)  # chemistry set index
        self._num_elements = c_int(0)  # number of elements
        self._num_gas_species = c_int(0)  # number of gas species
        self._num_gas_reactions = c_int(0)  # number of gas-phase reactions
        #
        self._eos = c_int(0)  # equation of state (default 0 = ideal gas)
        self.userealgas = False  # use ideal gas law by default
        self._awt_done = 0
        self._wt_done = 0
        self.ncf_done = 0
        # fake initialization
        self._awt = np.zeros(1, dtype=np.double)
        self._wt = np.zeros(1, dtype=np.double)
        # information about elements
        self._elements: dict[str, int] = {}  # element symbols dictionary
        self.esymbol: list[str] = []
        self._esym_done = 0
        # information about gas-phase species
        self._gas_species: dict[str, int] = {}  # gas species symbols dictionary
        self.ksymbol: list[str] = []
        self._ksym_done = 0
        # chemistry set label
        self.label = " "
        if len(label) > 0:
            self.label = label
        # default linking file names
        self._gas_link = "chem.asc"
        self._surf_link = "surf.asc"
        self._tran_link = "tran.asc"
        # summary file for the preprocessing step
        self._summary_out = "Summary.out"
        # set the mechanism input files names if given
        self.set_file_names(chem, surf, therm, tran)
        # check surface mechanism
        if Path(self._surf_file).is_file():
            self._index_surf = ctypes.c_int(1)
        # check transport data file
        if Path(self._tran_file).is_file():
            self._index_tran = ctypes.c_int(1)
        # surface chemistry properties
        # number of surface materials
        self._num_materials = c_int(0)
        # surface material map {str, int}
        self.material_map: dict[str, int] = {}
        # surface material objects {str, Material object}
        self.materials: dict[str, Material] = {}
        # surface material names [str]
        self.matsymbol: list[str] = []
        self._materialnamedone = 0
        # total number of surface site species
        self._num_max_site_species = c_int(0)
        # total number of bulk species
        self._num_max_bulk_species = c_int(0)
        # number of phases
        self._num_max_phases = c_int(0)
        # total number of surface reactions
        self._num_max_surf_reactions = c_int(0)
        #
        self._material_temperature: list[float] = []
        self._material_area: list[float] = []

    @property
    def chemfile(self) -> str:
        """Get gas-phase mechanism file name of this chemistry set."""
        """
        Get gas-phase mechanism file name of this chemistry set.

        Returns
        -------
            chemfile: string
                Full path and name of the Chemkin gas-phase mechanism input file

        """
        return self._gas_file

    @chemfile.setter
    def chemfile(self, filename: str):
        """Assign the gas-phase mechanism filename."""
        """
        Assign the gas-phase mechanism filename.

        Parameters
        ----------
            filename: string
                name of the gas-phase mechanism file with the full path

        """
        self._gas_file = filename

    @property
    def thermfile(self) -> str:
        """Get thermodynamic data filename of this chemistry set."""
        """
        Get thermodynamic data filename of this chemistry set.

        Returns
        -------
            thermfile: string
                Full path and name of the Chemkin thermodynamic data file

        """
        return self._therm_file

    @thermfile.setter
    def thermfile(self, filename: str):
        """Assign the thermodynamic data filename."""
        """
        Assign the thermodynamic data filename.

        Parameters
        ----------
            filename: string
                name of the thermodynamic data file with the full path

        """
        self._therm_file = filename

    @property
    def tranfile(self) -> str:
        """Get transport data filename of this chemistry set."""
        """
        Get transport data filename of this chemistry set.

        Returns
        -------
            tranfile: string
                Full path and name of the Chemkin thransport data file

        """
        return self._tran_file

    @tranfile.setter
    def tranfile(self, filename: str):
        """Assign the transport data filename."""
        """
        Assign the transport data filename.

        Parameters
        ----------
            filename: string
                name of the transport data file with the full path

        """
        self._tran_file = filename
        if Path(self._tran_file).is_file():
            self._index_tran = ctypes.c_int(1)
        else:
            self._index_tran = c_int(0)
            msg = [
                Color.RED,
                "transport data file",
                self._tran_file,
                "not found.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()

    @property
    def summaryfile(self) -> str:
        """Get the name of the summary file from the preprocessor."""
        """
        Get the name of the summary file from the preprocessor.

        Returns
        -------
            tranfile: string
                Full path and name of the preprocessing summary file

        """
        return self._summary_out

    def preprocess_transportdata(self):
        """Instruct the preprocessor to include the transport data."""
        if self._index_tran.value == 0:
            # send a warning message
            msg = [
                Color.MAGENTA,
                "make sure the gas mechanism contains the 'TRANSPORT ALL' block.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            self._index_tran = ctypes.c_int(1)
        else:
            # send the confirmation message
            if self._index_tran.value == 1:
                msg = [
                    Color.YELLOW,
                    "transport data in file:",
                    self._tran_file,
                    "will be processed.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            else:
                msg = [
                    Color.YELLOW,
                    "transport data in file:",
                    self._gas_file,
                    "will be processed.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)

    @property
    def surffile(self) -> str:
        """Get surface mechanism filename of this chemistry set."""
        """
        Get surface mechanism filename of this chemistry set.

        Returns
        -------
            tranfile: string
                Full path and name of the Chemkin surface mechanism input file

        """
        return self._surf_file

    @surffile.setter
    def surffile(self, filename: str):
        """Assign the surface mechanism filename."""
        """
        Assign the surface mechanism filename.

        Parameters
        ----------
            filename: string
                name of the surface mechanism file with the full path

        """
        self._surf_file = filename
        if Path(self._surf_file).is_file():
            self._index_surf = ctypes.c_int(1)
        else:
            self._index_surf = c_int(0)
            msg = [
                Color.RED,
                "surface mechanism file",
                self._surf_file,
                "not found.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()

    def set_file_names(
        self,
        chem: str = "",
        surf: str = "",
        therm: str = "",
        tran: str = "",
    ):
        """Assign all input files of the chemistry set."""
        """
        Assign all input files of the chemistry set.

        Parameters
        ----------
            chem: string, optional
                name of the gas mechanism file with the full path
            surf: string, optional
                name of the surface mechanism file with the full path
            therm: string, optional
                name of the thermodynamic data file with the full path
            tran: string, optional
                name of the transport data file with the full path

        """
        self._chemset_index = c_int(-1)
        if len(chem) > 1:
            self._gas_file = chem
        else:
            self._gas_file = "chem.inp"
        if len(surf) > 1:
            self._surf_file = surf
            if Path(self._surf_file).is_file():
                self._index_surf = ctypes.c_int(1)
            else:
                self._index_surf = ctypes.c_int(0)
                msg = [
                    Color.RED,
                    "surface mechanism file",
                    self._surf_file,
                    "not found",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()
        else:
            self._surf_file = "surf.inp"
            self._index_surf = c_int(0)
        if len(therm) > 1:
            self._therm_file = therm
        else:
            self._therm_file = "therm.dat"
        if len(tran) > 1:
            self._tran_file = tran
            if Path(self._tran_file).is_file():
                self._index_tran = ctypes.c_int(1)
            else:
                self._index_tran = c_int(0)
                msg = [
                    Color.RED,
                    "transport data file",
                    self._tran_file,
                    "not found.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()
        else:
            self._tran_file = "tran.dat"
            self._index_tran = c_int(0)

    def preprocess(self) -> int:
        """Run Chemkin preprocessor."""
        """
        Run Chemkin preprocessor.

        Returns
        -------
            Error code: integer

        """
        # check minimum set of required files
        if not Path(self._gas_file).is_file():
            msg = [
                Color.RED,
                "gas mechanism file",
                self._gas_file,
                "not found.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        if not Path(self._therm_file).is_file():
            msg = [
                Color.MAGENTA,
                "thermodynamic data file not found/specified.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            msg = [
                Color.YELLOW,
                "make sure the mechanism file",
                self._gas_file,
                "contains the 'THERM ALL' keyword.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

        # verify chemistry set
        # create a new identifier for this chemistry set
        name = self._gas_file + self._therm_file
        if self._index_tran.value == 1:
            name = name + self._tran_file
        if self._index_surf.value == 1:
            name = name + self._surf_file
        identifier = name
        # check if this chemistry set is already processed by this project
        if identifier in _chemset_identifiers:
            # existing chemistry set
            myindex = _chemset_identifiers.index(identifier)
            msg = [
                Color.YELLOW,
                "chemistry set is already processed\n",
                Color.SPACEx6,
                "the chemistry set index is:",
                str(myindex),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        else:
            # new chemistry set
            # add the identifier to the chemistry identifiers list
            _chemset_identifiers.append(identifier)
            # modify linking file names
            myindex = len(_chemset_identifiers)

        if myindex > 1:
            myfilename = "chem_" + str(myindex - 1) + ".asc"
            self._gas_link = myfilename
            myfilename = "surf_" + str(myindex - 1) + ".asc"
            self._surf_link = myfilename
            myfilename = "tran_" + str(myindex - 1) + ".asc"
            self._tran_link = myfilename
            # modify the summary file for the preprocessing step
            myfilename = "Summary_" + str(myindex - 1) + ".out"
            self._summary_out = myfilename

        # run preprocessor
        try:
            self._error_code = ck_wrapper.chemkin.KINPreProcess(
                self._index_surf,
                self._index_tran,
                ctypes.c_char_p(self._gas_file.encode("utf-8")),
                ctypes.c_char_p(self._surf_file.encode("utf-8")),
                ctypes.c_char_p(self._therm_file.encode("utf-8")),
                ctypes.c_char_p(self._tran_file.encode("utf-8")),
                ctypes.c_char_p(self._gas_link.encode("utf-8")),
                ctypes.c_char_p(self._surf_link.encode("utf-8")),
                ctypes.c_char_p(self._tran_link.encode("utf-8")),
                ctypes.c_char_p(self._summary_out.encode("utf-8")),
                self._chemset_index,
            )
        except OSError:
            self._error_code = 1
            return self._error_code

        if self._error_code == 0:
            ierr = ck_wrapper.chemkin.KINGetChemistrySizes(
                self._chemset_index,
                self._num_elements,
                self._num_gas_species,
                self._num_gas_reactions,
                self._num_materials,
                self._num_max_site_species,
                self._num_max_bulk_species,
                self._num_max_phases,
                self._num_max_surf_reactions,
            )

            if ierr != 0:
                # failed to find mechanism sizes
                msg = [
                    Color.RED,
                    "failed to find mechanism parameters\n",
                    Color.SPACEx6,
                    "error code =",
                    str(ierr),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()

            # get species symbols in a dictionary
            self.species_symbols
            # get element symbols
            self.element_symbols
            # get species molecular masses
            self._wt = self.wt
            # get atomic masses
            self._awt = self.awt
            # check real-gas model
            self.verify_realgas_model()
            # setup surface chemistry and materials
            if self._index_surf.value == 1 and self._num_materials.value > 0:
                self.set_surface_chemistry()
            # create a new Chemkin-CFD-API initialization flag for this Chemistry Set
            chemistryset_new(self._chemset_index.value)
            # save the chemkin work spaces for later use (by using active())
            self.save()
        else:
            # fail to preprocess the chemistry files
            details = "\n"
            if self._error_code == 1 or self._error_code == 3203:
                details = "\n" + Color.SPACEx6 + "cannot find a valid license.\n"
            msg = [
                Color.RED,
                "failed to preprocess the chemistry set,",
                "error code =",
                str(self._error_code),
                details,
                Color.SPACEx6,
                "the chemistry set index is:",
                str(self._chemset_index.value),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()

        return self._error_code

    def verify_realgas_model(self):
        """Verify the availability of real-gas data in the mechanism."""
        eos_model = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
        try:
            # check if the mechanism contains the real-gas EOS data
            ierr = ck_wrapper.chemkin.KINRealGas_GetEOSMode(
                self._chemset_index, self._eos, eos_model
            )
            if ierr == 0:
                if self._eos.value > 0 and self._eos.value <= 5:
                    msg = [
                        Color.YELLOW,
                        "real-gas cubic EOS",
                        "'" + str(eos_model.value.decode()) + "'",
                        "is available.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
                    del eos_model
                    return

            self._eos = c_int(0)
        except OSError:
            # mechanism contains no real-gas data
            self._eos = c_int(0)
            if verbose():
                msg = [Color.PURPLE, "accessing the real gas information.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)

        del eos_model
        if verbose():
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

    def verify_transport_data(self) -> bool:
        """Verify the availability of transport property data in the mechanism."""
        """
        Verify the availability of transport property data in the mechanism.

        Returns
        -------
            availability: boolean
                True = the transport property is available

        """
        if self._index_tran.value == 0:
            # no transport data
            return False
        #
        return True

    def verify_surface_mechanism(self) -> bool:
        """Verify the availability of surface chemistry data in the mechanism."""
        """
        Verify the availability of surface chemistry data in the mechanism.

        Returns
        -------
            availability: boolean
                True = the surface chemistry data is available

        """
        if self._index_surf.value == 0:
            # no surface chemistry data
            return False
        #
        return True

    @property
    def species_symbols(self):
        """Get list of gas species symbols."""
        """
        Get list of gas species symbols.

        Returns
        -------
            ksymbol: list of strings
                list of species symbols in the gas-phase mechanism

        """
        global MAX_SPECIES_LENGTH
        if self._ksym_done == 0:
            # recycle existing data
            buff = (LP_c_char * self._num_gas_species.value)()
            for i in range(0, self._num_gas_species.value):
                buff[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp = ctypes.cast(buff, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetGasSpeciesNames(self._chemset_index, pp)
            if ierr == 0:
                self._gas_species.clear()
                for index in range(0, len(buff)):
                    str_val = ctypes.cast(buff[index], c_char_p).value.decode()
                    self._gas_species[str_val] = index
                self._ksym_done == 1
            else:
                # failed to get species symbols
                msg = [Color.PURPLE, "failed to get species symbols.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff, pp

        # convert string type
        mylist = list(self._gas_species.keys())
        self.ksymbol.clear()
        for s in mylist:
            self.ksymbol.append(s)
        del mylist
        return self.ksymbol

    @property
    def element_symbols(self):
        """Get the list of element symbols."""
        """
        Get the list of element symbols.

        Returns
        -------
            esymbol: list of strings
                list of element symbols in the mechanism

        """
        if self._esym_done == 0:
            buff_ele = (LP_c_char * self._num_elements.value)()
            for j in range(0, self._num_elements.value):
                buff_ele[j] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_ele = ctypes.cast(buff_ele, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetElementNames(self._chemset_index, pp_ele)
            if ierr == 0:
                self._elements.clear()
                for index in range(0, len(buff_ele)):
                    ele_val = ctypes.cast(buff_ele[index], c_char_p).value.decode()
                    self._elements[ele_val] = index
                self._esym_done == 1
            else:
                # failed to get element symbols
                msg = [Color.PURPLE, "failed to get element symbols.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_ele, pp_ele

        # convert string type
        my_ele_list = list(self._elements.keys())
        self.esymbol.clear()
        for s_ele in my_ele_list:
            self.esymbol.append(s_ele)
        del my_ele_list
        return self.esymbol

    def get_specindex(self, specname: str) -> int:
        """Get index of the gas species."""
        """
        Get species index corresponding to the given species symbol.
        The species symbols (gas species and the site and the bulk
        species of all surface materials) must be unique.

        Returns
        -------
            specindex: integer
                index of the given species symbols in the gas-phase mechanism

        """
        # gas species
        specindex = self._gas_species.get(specname, -1)
        #
        if specindex < 0:
            if self.verify_surface_mechanism():
                n = 0
                while specindex < 0 and n < self._num_materials.value:
                    # loop over all surface materials
                    mname = self.material_names[n]
                    m = self.materials[mname]
                    # check surface site species
                    k_global, k_local = m.get_site_specindex(specname)
                    if k_local < 0:
                        # check bulk species
                        k_global, k_local = m.get_bulk_specindex(specname)
                        #
                        if k_local < 0:
                            specindex = -1
                        else:
                            specindex = k_global
                    else:
                        specindex = k_global
                    n += 1
        #
        if specindex < 0:
            msg = [Color.PURPLE, "species symbol not found:", specname, Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return specindex

    @property
    def chemid(self) -> int:
        """Get chemistry set index."""
        """
        Get chemistry set index.

        Returns
        -------
            chemidx: integer
                index of the Chemistry set

        """
        if self._chemset_index.value >= 0:
            return self._chemset_index.value
        else:
            return -1

    @property
    def surfchem(self) -> int:
        """Get surface chemistry status."""
        """
        Get surface chemistry status.

        Returns
        -------
            status: integer, {0, 1}
                indicating whether the Chemistry set includes a surface mechanism
                0 = this chemistry set does NOT include a surface chemistry
                1 = this chemistry set includes a  surface chemistry

        """
        return self._index_surf.value

    @property
    def kk(self) -> int:
        """Get number of gas species."""
        """
        Get number of gas species.

        Returns
        -------
            kk: integer
                total number of gas-phase species in the Chemistry set

        """
        return self._num_gas_species.value

    # alias
    number_species = kk

    @property
    def mm(self) -> int:
        """Get number of elements in the chemistry set."""
        """
        Get number of elements in the chemistry set.

        Returns
        -------
            mm: integer
                total number of elements in the Chemistry set

        """
        return self._num_elements.value

    # alias
    number_elements = mm

    @property
    def ii_gas(self) -> int:
        """Get number of gas-phase reactions."""
        """
        Get number of gas-phase reactions.

        Returns
        -------
            ii_gas: integer
                total number of gas-phase reactions in the Chemistry set

        """
        return self._num_gas_reactions.value

    # alias
    number_gas_reactions = ii_gas

    @property
    def awt(self) -> npt.NDArray[np.double]:
        """Compute atomic masses."""
        """
        Compute atomic masses.

        Returns
        -------
            awt: 1-D double array
                masses of the elements in the Chemistry set [g/mole]

        """
        if self._awt_done == 1:
            return self._awt
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        del self._awt  # clear the "original" definition in __init__
        self._awt = np.zeros(self._num_elements.value, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetAtomicWeights(self._chemset_index, self._awt)
        if ierr == 0:
            self._awt_done = 1
        else:
            # failed to find atomic masses
            msg = [Color.PURPLE, "failed to get atomic masses.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return self._awt

    # alias
    atomic_weight = awt

    @property
    def wt(self) -> npt.NDArray[np.double]:
        """Compute gas species molecular masses."""
        """
        Compute gas species molecular masses.

        Returns
        -------
            wt: 1-D double array
                molecular masses of the gas-phase species in the Chemistry set [g/mole]

        """
        if self._wt_done == 1:
            return self._wt
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        del self._wt  # clear the "original" definition in __init__
        self._wt = np.zeros(self._num_gas_species.value, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetGasMolecularWeights(
            self._chemset_index, self._wt
        )
        if ierr == 0:
            self._wt_done = 1
        else:
            # failed to find molecular masses
            msg = [Color.PURPLE, "failed to get species molecular masses.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return self._wt

    # alias
    species_molar_weight = wt

    def species_cp(
        self, temp: float = 0.0, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get species specific heat capacity at constant pressure."""
        """
        Get species specific heat capacity at constant pressure.

        Parameters
        ----------
            temp: double
                Temperature [K]
            pres: double, optional
                Pressure [dynes/cm2] required when real gas model is activated

        Returns
        -------
            cp: 1-D double array
                species specific heat capacities at constant pressure [ergs/mol-K]

        """
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if temp <= 1.0e0:
            msg = [Color.PURPLE, "temperature value is too low.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check real-gas
        if check_realgas_status(self.chemid):
            if pres is None:
                # pressure is not assigned
                msg = [
                    Color.PURPLE,
                    "mixture pressure must be provided",
                    "to evaluate real-gas species properties\n",
                    Color.SPACEx6,
                    "usage: <Chemistry_Obj>.species_cp(temperature, pressure)",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # set current pressure for the real-gas
                set_current_pressure(self.chemid, pres)
        #
        tt = c_double(temp)
        cp = np.zeros(self._num_gas_species.value, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetGasSpecificHeat(self._chemset_index, tt, cp)
        if ierr == 0:
            # convert [ergs/g-K] to [ergs/mol-K]
            # for k in range(len(Cp)):
            #    Cp[k] = Cp[k] * self._wt[k]
            cp *= self._wt
        else:
            # failed to compute specific heats
            msg = [Color.PURPLE, "failed to compute specific heats.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return cp

    def species_cv(
        self, temp: float = 0.0, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get species specific heat capacity at constant volume."""
        """
        Get species specific heat capacity at constant volume (ideal gas only).

        Parameters
        ----------
            temp: double
                Temperature [K]
            pres: double, optional
                Pressure [dynes/cm2] required when real gas model is activated

        Returns
        -------
            cv: 1-D double array
                species specific heat capacities at constant volume [ergs/mol-K]

        """
        if check_realgas_status(self.chemid) and pres is None:
            # pressure is not assigned
            msg = [
                Color.PURPLE,
                "mixture pressure must be provided",
                "to evaluate real-gas species properties\n",
                Color.SPACEx6,
                "usage: <Chemistry_Obj>.species_cv(temperature, pressure)",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        cv = self.species_cp(temp, pres)
        rgas = R_GAS
        # for k in range(len(Cp)):
        #    Cv[k] = Cp[k] - R
        cv -= rgas

        return cv

    def species_h(
        self, temp: float = 0.0, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get species enthalpy."""
        """
        Get species enthalpy.

        Parameters
        ----------
            temp: double
                Temperature [K]
            pres: double, optional
                Pressure [dynes/cm2] required when real gas model is activated

        Returns
        -------
            h: 1-D double array
                species enthalpy [ergs/mol]

        """
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if temp <= 1.0e0:
            msg = [Color.PURPLE, "temperature value is too low.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check real-gas
        if check_realgas_status(self.chemid):
            if pres is None:
                # pressure is not assigned
                msg = [
                    Color.PURPLE,
                    "mixture pressure must be provided",
                    "to evaluate real-gas species properties\n",
                    Color.SPACEx6,
                    "usage: <Chemistry_Obj>.species_h(temperature, pressure)",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # set current pressure for the real-gas
                set_current_pressure(self.chemid, pres)
        tt = c_double(temp)
        h = np.zeros(self._num_gas_species.value, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetGasSpeciesEnthalpy(self._chemset_index, tt, h)
        if ierr == 0:
            # convert [ergs/gm] to [ergs/mol]
            # for k in range(len(H)):
            #    H[k] = H[k], * self._wt[k]
            h *= self._wt
        else:
            # failed to compute enthalpies
            msg = [Color.PURPLE, "failed to compute specific enthalpies.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return h

    def species_u(
        self, temp: float = 0.0, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get species internal energy."""
        """
        Get species internal energy.

        Parameters
        ----------
            temp: double
                Temperature [K]
            pres: double, optional
                Pressure [dynes/cm2] required when real gas model is activated

        Returns
        -------
            u: 1-D double array
                species internal energy [ergs/mol]

        """
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if temp <= 1.0e0:
            msg = [Color.PURPLE, "temperature value is too low.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check real-gas
        if check_realgas_status(self.chemid):
            if pres is None:
                # pressure is not assigned
                msg = [
                    Color.PURPLE,
                    "mixture pressure must be provided",
                    "to evaluate real-gas species properties\n",
                    Color.SPACEx6,
                    "usage: <Chemistry_Obj>.species_u(temperature, pressure)",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # set current pressure for the real-gas
                set_current_pressure(self.chemid, pres)
        tt = c_double(temp)
        u = np.zeros(self._num_gas_species.value, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetGasSpeciesInternalEnergy(
            self._chemset_index, tt, u
        )
        if ierr == 0:
            # convert [ergs/gm] to [ergs/mol]
            # for k in range(len(U)):
            #    U[k] = U[k], * self._wt[k]
            u *= self._wt
        else:
            # failed to compute internal energies
            msg = [
                Color.PURPLE,
                "failed to compute specific internal energies.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return u

    def species_visc(self, temp: float = 0.0) -> npt.NDArray[np.double]:
        """Get species viscosity."""
        """
        Get species viscosity.

        Parameters
        ----------
            temp: double
                Temperature [K]

        Returns
        -------
            visc: 1-D double array
                species viscosity [gm/cm-sec]

        """
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if self._index_tran.value != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if temp <= 1.0e0:
            msg = [Color.PURPLE, "temperature value is too low.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        tt = c_double(temp)
        visc = np.zeros(self._num_gas_species.value, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetViscosity(self._chemset_index, tt, visc)
        if ierr != 0:
            # failed to compute viscosity
            msg = [Color.PURPLE, "failed to compute specific viscosities.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return visc

    def species_cond(self, temp: float = 0.0) -> npt.NDArray[np.double]:
        """Get species conductivity."""
        """
        Get species conductivity.

        Parameters
        ----------
            temp: double
                Temperature [K]

        Returns
        -------
            cond: 1-D double array
                species conductivity [ergs/cm-K-sec]

        """
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if self._index_tran.value != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if temp <= 1.0e0:
            msg = [Color.PURPLE, "temperature value is too low.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        tt = c_double(temp)
        cond = np.zeros(self._num_gas_species.value, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetConductivity(self._chemset_index, tt, cond)
        if ierr != 0:
            # failed to compute conductivities
            msg = [
                Color.PURPLE,
                "failed to compute specific conductivities.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return cond

    def species_diffusioncoeffs(
        self, press: float = 0.0, temp: float = 0.0
    ) -> npt.NDArray[np.double]:
        """Get species diffusion coefficients."""
        """
        Get species diffusion coefficients.

        Parameters
        ----------
            press: double
                Pressure [dynes/cm2]
            temp: double
                Temperature [K]

        Returns
        -------
            diffusioncoeffs: 2-D double array, dimension
            = [number_species, number_species]
                species diffusion coefficients [cm2/sec]

        """
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if self._index_tran.value != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if temp <= 1.0e0:
            msg = [Color.PURPLE, "temperature value is too low.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if press <= 1.0e0:
            msg = [Color.PURPLE, "pressure value is too low.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        pp = c_double(press)
        tt = c_double(temp)
        dim = (self._num_gas_species.value, self._num_gas_species.value)
        diffusioncoeffs = np.zeros(dim, dtype=np.double, order="F")
        ierr = ck_wrapper.chemkin.KINGetDiffusionCoeffs(
            self._chemset_index, pp, tt, diffusioncoeffs
        )
        if ierr != 0:
            # failed to compute diffusion coefficients
            msg = [
                Color.PURPLE,
                "failed to compute specific diffusion coefficients.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return diffusioncoeffs

    def species_composition(self, elemindex: int = -1, specindex: int = -1) -> int:
        """Get elemental composition of a species."""
        """
        Get elemental composition of a species.

        Parameters
        ----------
            elemindex: integer
                index of the element
            specindex: integer
                index of the gas species

        Returns
        -------
            count: integer
                number of the element in the given gas species

        """
        if self.ncf_done == 0:
            # initialize the NCF matrix
            dim = (self._num_elements.value, self._num_gas_species.value)
            self.elementalcomp = np.zeros(dim, dtype=np.int32, order="F")
            # load the NCF matrix
            ierr = ck_wrapper.chemkin.KINGetGasSpeciesComposition(
                self._chemset_index, self.elementalcomp
            )
            if ierr != 0:
                msg = [
                    Color.PURPLE,
                    "failed to compute elemental compositions.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                self.ncf_done = 1

        # check element index
        if elemindex < 0 or elemindex >= self._num_elements.value:
            msg = [Color.PURPLE, "element index is out of bound.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check species index
        if specindex < 0 or specindex >= self._num_gas_species.value:
            msg = [Color.PURPLE, "species index is out of bound.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return self.elementalcomp[elemindex][specindex]

    @property
    def eos(self) -> int:
        """Get the available real-gas EOS model."""
        """
        Get the available real-gas EOS model.

        Returns
        -------
        count: integer
            number of available cubic EOS models in Chemkin

        """
        return self._eos.value

    def use_realgas_cubic_eos(self):
        """Turn ON the real-gas cubic EOS."""
        """Turn ON the real-gas cubic EOS to compute mixture properties
        if the mechanism contains necessary data."""
        if self._eos.value < 1:
            # no real gas EOS data in the mechanism
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return
        # check real-gas EOS status
        iflag = c_int(0)
        ierr = ck_wrapper.chemkin.KINRealGas_UseCubicEOS(self._chemset_index, iflag)
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to turn ON the real-gas EOS model,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if iflag.value == 0:
            msg = [
                Color.YELLOW,
                "real-gas cubic EOS model",
                Chemistry.realgas_cueos[self._eos.value],
                "is turned ON.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = True
        else:
            self.userealgas = False

    def use_idealgas_law(self):
        """Turn on the ideal gas law to compute mixture properties."""
        if self._eos.value < 1:
            # no real gas EOS data in the mechanism
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = False
            return
        # check real-gas EOS status
        iflag = c_int(0)
        ierr = ck_wrapper.chemkin.KINRealGas_UseIdealGasLaw(self._chemset_index, iflag)
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to turn ON ideal gas law,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if iflag.value == 0:
            msg = [Color.YELLOW, "ideal gas law is turned ON.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = False

    def get_reaction_parameters(
        self,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get the Arrhenius reaction rate parameters of all gas-phase reactions."""
        """
        Get the Arrhenius reaction rate parameters of all gas-phase reactions.

        Returns
        -------
            a_factor: 1-D double array
                Arrhenius pre-exponent A-Factor of reaction in [mole-cm3-sec-K]
            beta: 1-D double array
                Arrhenius temperature exponent [-]
            act_energy: 1-D double array
                activation temperature [K]

        """
        reactionsize = self.ii_gas
        # pre-exponent A factor of all gas-phase reactions in the mechanism
        # in cgs units [mole-cm3-sec-K]
        a_factor = np.zeros(shape=reactionsize, dtype=np.double)
        # temperature exponent of all reactions [-]
        beta = np.zeros_like(a_factor, dtype=np.double)
        # activation energy/temperature of all reactions [K]
        act_energy = np.zeros_like(a_factor, dtype=np.double)
        # get the reaction parameters
        ierr = ck_wrapper.chemkin.KINGetReactionRateParameters(
            self._chemset_index, a_factor, beta, act_energy
        )
        if ierr != 0:
            a_factor[:] = 0.0e0
            beta[:] = 0.0e0
            act_energy[:] = 0.0e0
        return a_factor, beta, act_energy

    def set_reaction_afactor(self, reaction_index: int, a_factor: float):
        """Set the Arrhenius A-Factor of the given reaction."""
        """
        (Re)set the Arrhenius A-Factor of the given reaction.

        Parameters
        ----------
            reaction_index: integer
                index of the gas-phase reaction of which the A-Factor to be reset
            a_factor: double
                new A-Factor value in [mole-cm3-sec-K]

        """
        # check inputs
        if reaction_index > self.ii_gas or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.ii_gas) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if a_factor < 0.0e0:
            msg = [Color.PURPLE, "A-Factor must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert the reaction parameters
        ireac = c_int(-reaction_index)  # negative index to "put" A-factor value
        ierr = ck_wrapper.chemkin.KINSetAFactorForAReaction(
            self._chemset_index, ireac, c_double(a_factor)
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to set Arrhenius A-Factor,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_reaction_afactor(self, reaction_index: int) -> float:
        """Get the Arrhenius A-Factor of the given reaction."""
        """
        Get the Arrhenius A-Factor of the given reaction.

        Parameters
        ----------
            reaction_index: integer
                index of the reaction

        Returns
        -------
            a_factor: double
                Arrhenius A-Factor of the given reaction in [mole-cm3-sec-K]

        """
        # initialization
        a_factor = c_double(0.0e0)
        # check inputs
        if reaction_index > self.ii_gas or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.ii_gas) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        # get the A-factor value
        ierr = ck_wrapper.chemkin.KINSetAFactorForAReaction(
            self._chemset_index, ireac, a_factor
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to find Arrhenius A-Factor,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return a_factor.value

    def get_gas_reaction_string(self, reaction_index: int) -> str:
        """Get the reaction string of the given reaction index."""
        """Get the reaction string of the gas-phase reaction specified
        by the reaction index.

        Parameters
        ----------
            reaction_index: integer
                (1-base) gas-phase reaction index

        Returns
        -------
            reactionstring: string
                reaction string of the given reaction

        """
        # initialization
        reactionstring = ""
        if reaction_index > self._num_gas_reactions.value:
            msg = [
                Color.PURPLE,
                "reaction index must be <",
                str(self._num_gas_reactions.value),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        elif reaction_index <= 0:
            msg = [Color.PURPLE, "reaction index must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        i_string_size = c_int(0)
        # get reaction string (might have to be increased to 2048 for 26R1)
        rstring = bytes(" " * 1024, "utf-8")
        ierr = ck_wrapper.chemkin.KINGetGasReactionString(
            self._chemset_index, ireac, i_string_size, rstring
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to find reaction string,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert C string back to string
        # print(rstring.decode()[0:i_string_size.value])  # check
        reactionstring = rstring.decode()[0 : i_string_size.value]
        del rstring
        return reactionstring

    def save(self):
        """Store the work spaces of the current Chemistry Set."""
        """Store the work spaces of the current Chemistry Set
        if new Chemistry Set will be created later in the same project.
        """
        ierr = ck_wrapper.chemkin.KINUpdateChemistrySet(self._chemset_index)
        if ierr == 0:
            # mark this chemistry set as active
            global _active_chemistry_set
            _active_chemistry_set = self._chemset_index.value
            msg = [
                Color.YELLOW,
                "work spaces saved for Chemistry Set",
                self.label,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        else:
            msg = [Color.PURPLE, "saving the Chemistry Set work spaces.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def activate(self):
        """Activate the work spaces of the Chemistry Set."""
        """Switch to (re-activate) the work spaces of the current Chemistry Set
        when there are multiple Chemistry Sets in the same project.
        """
        ierr = activate_chemistryset(self._chemset_index.value)
        if ierr == 0:
            # mark this chemistry set as active
            msg = [
                Color.YELLOW,
                "work spaces saved for Chemistry Set",
                self.label,
                "activated.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        else:
            exit()

    # surface chemistry related
    def set_surface_chemistry(self):
        """Set up the surface chemistry data."""
        """
        Set up the surface chemistry data of the Chemistry Set.
        Get the dimensions of the surface phases, species,
        and reactions. Instantiate the surface Material objects.
        """
        # number of surface materials of the Chemistry Set
        self.numb_mat = self._num_materials.value
        # dimension for variable arrays related to all species (of any material)
        # in the Chemistry Set
        self.max_total_species = (
            self._num_gas_species.value
            + self._num_max_site_species.value
            + self._num_max_bulk_species.value
        )
        # dimension for variable arrays related to all reactions (of any material)
        # in the Chemistry Set
        self.max_total_reactions = (
            self._num_gas_reactions.value + self._num_max_surf_reactions.value
        )
        # set up surface materials
        # initialization
        self.material_map = {}  # surface material map {str, int}
        self.materials = {}  # surface material objects {str, Material object}
        self.matsymbol = []  # surface material names [str]
        self._materialnamedone = 0
        # get material names
        self.material_names
        # set up surface materials
        for i, key in enumerate(self.matsymbol):
            # material index is 0-based
            self.material_map[key] = i
            m = Material(mat_index=i, chem=self, label=key)
            self.materials[key] = copy.deepcopy(m)

    def no_surface_mechanism_declaration(self):
        """Inform users the Chemistry Set does not have surface chemistry."""
        msg = [
            Color.YELLOW,
            "Chemsitry Set",
            self.label,
            "does not contain surface chemistry.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)

    @property
    def number_materials(self) -> int:
        """Get the number of surface materials."""
        """
        Get the number of surface materials in the Chemistry Set.

        Return
        ------
            numb_materials: number of surface materials
        """
        if self.verify_surface_mechanism():
            return self._num_materials.value
        else:
            self.no_surface_mechanism_declaration()
            return 0

    @property
    def max_number_phases(self) -> int:
        """Get the maximum number of phases."""
        """
        Get the maximum number of phases among the surface materials
        of the chemistry set including the gas phase.
        A surface mechanism may contain one or more surface materials, and the
        maximum number of the surface phases (site phases + bulk phases)
        among the materials plus the gas phase will be returned.

        Return
        ------
            max_numb_phases: integer
                the maximum number of phases (including the gas phase)
                of the Chemistry Set
        """
        if self.verify_surface_mechanism():
            return self._num_max_phases.value
        else:
            self.no_surface_mechanism_declaration()
            return 1

    @property
    def max_number_sites(self) -> int:
        """Get the maximum number of surface site species."""
        """
        Get the maximum number of surface site species among the
        surface materials of the surface chemistry.
        A surface mechanism may contain one or more surface materials, and the
        maximum number among the materials will be returned.

        Return
        ------
            max_numb_site_species: integer
                the maximum number of site species of the surface materials
        """
        if self.verify_surface_mechanism():
            return self._num_max_site_species.value
        else:
            self.no_surface_mechanism_declaration()
            return 0

    @property
    def max_number_bulks(self) -> int:
        """Get the maximum number of surface bulk species."""
        """
        Get the maximum number of surface bulk species among the
        surface materials of the surface chemistry.
        A surface mechanism may contain one or more surface materials, and the
        maximum number among the materials will be returned.

        Return
        ------
            max_numb_bulk_species: integer
                the maximum number of bulk species of the surface materials
        """
        if self.verify_surface_mechanism():
            return self._num_max_bulk_species.value
        else:
            self.no_surface_mechanism_declaration()
            return 0

    @property
    def max_total_number_species(self) -> int:
        """Get the maximum number of all species."""
        """
        Get the maximum number of all species (gas, site and bulk) among the surface
        materials of the Chemistry Set. This number
        includes the gas phase species, surface site species, and the bulk species.
        A surface mechanism may contain one or more surface materials, and the
        maximum total species number among the materials will be returned.
        This value can be used to size the arrays holding species
        mole fractions/activities, thermodynamic properties on any material.

        Return
        ------
            max_total_species: integer
                the maximum number of all species among the surface materials
                in the Chemistry Set.
        """
        if self.verify_surface_mechanism():
            return self.max_total_species
        else:
            self.no_surface_mechanism_declaration()
            return self._num_gas_species.value

    @property
    def max_number_surface_reactions(self) -> int:
        """Get the maximum number of surface reactions."""
        """
        Get the maximum number of surface reactions among the
        surface materials of the surface chemistry.
        A surface mechanism may contain one or more surface materials, and the
        maximum number among the materials will be returned. This value can be used
        to size the arrays holding surface reaction rates on any material.

        Return
        ------
            max_numb_surface_rxn: integer
                the maximum number of surface reactions of the surface materials
        """
        if self.verify_surface_mechanism():
            return self._num_max_surf_reactions.value
        else:
            self.no_surface_mechanism_declaration()
            return 0

    @property
    def max_total_number_reactions(self) -> int:
        """Get the maximum number of gas and surface reactions."""
        """
        Get the maximum number of all (gas and surface) reactions among the
        surface materials in the Chemistry Set.
        A surface mechanism may contain one or more surface materials, and the
        maximum number among the materials will be returned. This value can be used
        to size the arrays holding all reaction rates on any material.

        Return
        ------
            max_total_rxn: integer
                the maximum number of surface reactions of the surface materials
        """
        if self.verify_surface_mechanism():
            return self.max_total_reactions
        else:
            self.no_surface_mechanism_declaration()
            return self._num_gas_reactions.value

    @property
    def material_names(self):
        """Get list of surface material names."""
        """
        Get list of surface material names of the Chemistry Set.

        Returns
        -------
            matsymbol: list of strings
                list of surface material names in the surface mechanism
        """
        if not self.verify_surface_mechanism():
            # no surface chemistry data
            self.no_surface_mechanism_declaration()
            return []
        # process surface chemistry data
        if self._materialnamedone == 0:
            # recycle existing data
            buff_m = (LP_c_char * self._num_materials.value)()
            for i in range(0, self._num_materials.value):
                buff_m[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_m = ctypes.cast(buff_m, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetMaterialNames(self._chemset_index, pp_m)
            if ierr == 0:
                self.material_map.clear()
                for index in range(0, len(buff_m)):
                    mat_val = ctypes.cast(buff_m[index], c_char_p).value.decode()
                    self.material_map[mat_val] = index
                self._materialnamedone == 1
            else:
                # failed to get species symbols
                msg = [Color.PURPLE, "failed to get surface material names.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_m
        # convert string type
        mylist = list(self.material_map.keys())
        self.matsymbol.clear()
        for s in mylist:
            self.matsymbol.append(s)
        del mylist
        return self.matsymbol

    def get_material_index(self, mat_name: str) -> int:
        """Get the surface material index."""
        """
        Get the surface material index by its name.

        Parameters
        ----------
            mat_name: str
                surface material name

        Returns
        -------
            index: integer
                index of the named surface material
        """
        index = -1
        if not self.verify_surface_mechanism():
            # no surface chemistry data
            self.no_surface_mechanism_declaration()
        #
        if mat_name in self.matsymbol:
            index = self.material_map.get(mat_name, -1)
        return index


class Material:
    """Surface Chemkin material object."""

    def __init__(self, mat_index: int, chem: Chemistry, label: str = ""):
        """Create a surface material object."""
        """
        Create a Surface Chemkin material object which may contain multiple
        surface site and bulk phases and reactions that take place on
        these surface phases.
        """
        # check
        if not chem.verify_surface_mechanism():
            msg = [
                Color.RED,
                "The Chemistry set does not contain surface chemistry.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # chemistry set index
        self._chemset_index = chem.chemid
        self._material_index = mat_index
        self._dispersed = False
        self._coal = False
        # material name
        if len(label) == 0:
            self.label = "Material_" + str(mat_index)
        else:
            self.label = label
        # dimensions
        chemset_index = c_int(self._chemset_index)
        material_id = c_int(mat_index)
        self.num_element = chem.mm
        self.num_gas_species = chem.kk
        self.num_gas_reactions = chem.ii_gas
        num_sites = c_int(0)
        num_bulks = c_int(0)
        num_phases = c_int(0)
        num_surf_rxns = c_int(0)
        ierr = ck_wrapper.chemkin.KINGetMaterialSizes(
            chemset_index,
            material_id,
            num_sites,
            num_bulks,
            num_phases,
            num_surf_rxns,
        )
        if ierr == 0:
            # number of phases on this material (including the gas phase)
            self.num_phases = num_phases.value
            # number of site phases on this material
            self.num_site_phase = 0
            # number of bulk phases on this material
            self.num_bulk_phase = 0
            # number of site species on this material
            self.num_site_species = num_sites.value
            # number of bulk species on this material
            self.num_bulk_species = num_bulks.value
            # number of surface reactions on this material
            self.num_surf_reactions = num_surf_rxns.value
            # including the gas phase
            self.total_species = chem.kk + self.num_site_species + self.num_bulk_species
            self.total_reactions = self.num_gas_reactions + self.num_surf_reactions
            #
            self.surf_species_map: dict[str, int] = {}
        else:
            # failed to obtain material information from chemkin
            self.num_phases = 0
            self.num_site_phase = 0
            self.num_bulk_phase = 0
            self.num_site_species = 0
            self.num_bulk_species = 0
            self.num_surf_reactions = 0
            self.total_species = 0
            self.total_reactions = 0
            #
            msg = [Color.PURPLE, "failed to process surface material data.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # process additional material data
        ibase0 = c_int(0)  # array index is 1-based
        # mapping of the surface species symbols (for all site and bulk species)
        # and the global index (starting from the gas species) on this material
        # {species_symbol, global_index}
        self.site_species_map: dict[str, int] = {}
        self.bulk_species_map: dict[str, int] = {}
        # get element symbols from the chemistry set
        self.ESymbols = chem.element_symbols
        # get number of phases on this material
        num_site_phase = c_int(0)
        num_bulk_phase = c_int(0)
        nphase = c_int(0)
        ierr = ck_wrapper.chemkin.KINGetPhaseSizes(
            chemset_index,
            material_id,
            nphase,
            num_site_phase,
            num_bulk_phase,
        )
        self.num_site_phase = num_site_phase.value
        self.num_bulk_phase = num_bulk_phase.value
        self.first_site_phase = 0
        self.last_site_phase = 0
        self.first_bulk_phase = 0
        self.last_bulk_phase = 0
        if self.num_site_phase > 0:
            self.first_site_phase = 1
            self.last_site_phase = self.first_site_phase + self.num_site_phase - 1
        if self.num_bulk_phase > 0:
            self.first_bulk_phase = self.last_site_phase + 1
            self.last_bulk_phase = self.first_bulk_phase + self.num_bulk_phase - 1
        # get number of surface/bulk species and their indices on this material
        self.nspecies_of_phase = np.zeros(self.num_phases, dtype=np.int32)
        self.first_species_id_of_phase = np.zeros(self.num_phases, dtype=np.int32)
        self.last_species_id_of_phase = np.zeros(self.num_phases, dtype=np.int32)
        ierr = ck_wrapper.chemkin.KINGetPhaseSpeciesIndices(
            chemset_index,
            material_id,
            self.nspecies_of_phase,
            self.first_species_id_of_phase,
            self.last_species_id_of_phase,
            ibase0,
        )
        # get number of surface reactions on this material
        # get phase names on this material
        self._phasenamesdone = 0
        self.PhaseSymbol: list[str] = []
        self.phase_map: dict[str, int] = {}
        self.phases: dict[str, SurfacePhase] = {}
        phase_type = ""
        # instantiate surface phases
        self.phase_names
        # global index of the site species on this material
        site_global_index = self.num_gas_species
        # global index of the bulk species on this material
        bulk_global_index = self.num_gas_species + self.num_site_species
        # surface species symbols
        # flag indicating the status of species setup of the phase
        self._sitespeciesnamesdone = 0
        self._bulkspeciesnamesdone = 0
        # site species symbols of this material
        self.site_symbols: list[str] = []
        # bulk species symbols of this material
        self.bulk_symbols: list[str] = []
        # species composition flag
        self.ncf_done = 0
        # surface species molecular weight arrays
        self.sitewt_done = 0
        self.bulkwt_done = 0
        self.site_wt = np.zeros(self.num_site_species, dtype=np.double)
        self.bulk_wt = np.zeros(self.num_bulk_species, dtype=np.double)
        #
        if self.num_site_phase > 0:
            # get site species symbols of this material
            self.site_symbols = self.site_species_names
            if len(self.site_symbols) > 0:
                self._sitespeciesnamesdone = 1
                for i, name in enumerate(self.site_symbols):
                    self.site_species_map[name] = site_global_index + i
            # compute site (phase) density [mole/cm2]
            self.site_density = self.get_site_density()
            # get site species molecular weights [g/mole]
            self.site_wt = self.get_site_molar_weights()
        else:
            self.site_density = []
        #
        if self.num_bulk_phase > 0:
            # bulk species symbols of this material
            self.bulk_symbols = self.bulk_species_names
            if len(self.bulk_symbols) > 0:
                self._bulkspeciesnamesdone = 1
                for i, name in enumerate(self.bulk_symbols):
                    self.bulk_species_map[name] = bulk_global_index + i
            # compute bulk (species) density [g/cm3]
            self.bulk_density = self.get_bulk_density()
            # get site species molecular weights [g/mole]
            self.bulk_wt = self.get_bulk_molar_weights()
        else:
            self.bulk_density = []
        # set up surface phases on this materials
        sp_count = 0
        bp_count = 0
        # loop over all phases
        for i, key in enumerate(self.PhaseSymbol):
            # number of species on this phase
            nspec = self.nspecies_of_phase[i]
            # phase index in Chemkin is 1-based
            # phase #0 is the gas phase
            if i == 0:
                # skip the gas phase
                continue
            elif i <= self.last_site_phase:
                # site phase
                phase_type = "site"
                site_global_index += nspec
                sp_count += 1
            else:
                # bulk phase
                phase_type = "bulk"
                bulk_global_index += nspec
                bp_count += 1
            # instantiate the surface phase
            s = SurfacePhase(phase_index=i, phase_type=phase_type, label=key)
            # set the number of species on this surface phase
            s.number_species = self.nspecies_of_phase[i]
            # set (1-based) global index of the first species
            s.first_species_index = self.first_species_id_of_phase[i]
            # set (1-based) global index of the last species
            s.last_species_index = self.last_species_id_of_phase[i]
            if phase_type == "site":
                # set site phase density [mole/cm2]
                s.site_density = self.site_density[i]
            # add to the phases dictionary
            self.phases[key] = copy.deepcopy(s)
            #
            del s

        # additional surface material properties
        self.surftemp = 300.0
        self.activearea = 0.0

    @property
    def material_index(self) -> int:
        """Get the surface material index."""
        return self._material_index

    @property
    def number_phases(self) -> int:
        """Get the number of phases."""
        """
        Get the number of phases (including the gas phase) of the material.

        Returns
        -------
            numb_phases: integer
                number of phases of the material

        """
        return self.num_phases

    @property
    def number_site_phases(self) -> int:
        """Get the number of site phases."""
        """
        Get the number of site phases of the material.

        Returns
        -------
            numb_site_phase: integer
                number of site phases of the material

        """
        return self.num_site_phase

    @property
    def number_bulk_phases(self) -> int:
        """Get the number of bulk phases."""
        """
        Get the number of bulk phases of the material.

        Returns
        -------
            numb_bulk_phase: integer
                number of bulk phases of the material

        """
        return self.num_bulk_phase

    @property
    def number_total_species(self) -> int:
        """Get the number of species in all phases."""
        """
        Get the number of species in all phases (gas, site, and bulk).

        Returns
        -------
            total_species: integer
                total number of species associated with the material
        """
        return self.total_species

    @property
    def number_site_species(self) -> int:
        """Get the total number of site species."""
        """
        Get the total number of site species (of all site phases) of the material.

        Returns
        -------
            numb_site_species: integer
                total number of site species of the material

        """
        return self.num_site_species

    @property
    def number_bulk_species(self) -> int:
        """Get the total number of bulk species."""
        """
        Get the total number of bulk species (of all bulk phases) of the material.

        Returns
        -------
            numb_bulk_species: integer
                total number of bulk species of the material

        """
        return self.num_bulk_species

    @property
    def first_site_phase_index(self) -> int:
        """Get phase index of the first site phase."""
        """
        Get phase index of the first site phase of the material.

        Returns
        -------
            first_site_phase: integer
                phase index (0-based) of the first site phase of the material

        """
        if self.number_site_phases > 0:
            return self.first_site_phase
        else:
            return -1

    @property
    def last_site_phase_index(self) -> int:
        """Get phase index of the last site phase."""
        """
        Get phase index of the last site phase of the material.

        Returns
        -------
            last_site_phase: integer
                phase index (0-based) of the last site phase of the material

        """
        if self.number_site_phases > 0:
            return self.last_site_phase
        else:
            return -1

    @property
    def first_bulk_phase_index(self) -> int:
        """Get phase index of the first bulk phase."""
        """
        Get phase index of the first bulk phase of the material.

        Returns
        -------
            first_bulk_phase: integer
                phase index (0-based) of the first bulk phase of the material

        """
        if self.number_bulk_phases > 0:
            return self.first_bulk_phase
        else:
            return -1

    @property
    def last_bulk_phase_index(self) -> int:
        """Get phase index of the last bulk phase."""
        """
        Get phase index of the last bulk phase of the material.

        Returns
        -------
            last_bulk_phase: integer
                phase index (0-based) of the last bulk phase of the material

        """
        if self.number_bulk_phases > 0:
            return self.last_bulk_phase
        else:
            return -1

    @property
    def number_surface_reactions(self) -> int:
        """Get the number of surface reactions."""
        """
        Get the number of surface reactions of the material.

        Returns
        -------
            ii_surf: integer
                number of surface reactions of the material

        """
        return self.num_surf_reactions

    @property
    def number_total_reactions(self) -> int:
        """Get the number of all reactions."""
        """
        Get the number of all reactions (gas and surface) of the material.

        Returns
        -------
            ii_total: integer
                total number of reactions associated with the material

        """
        return self.total_reactions

    def information(self):
        """List the size information of the material."""
        # surface material information
        print(f"Material Name: {self.label}")
        # gas-phase information
        print("----")
        print(f"  Number of gas species = {self.num_gas_species}")
        print(f"  Number of gas-phase reactions = {self.num_gas_reactions}")
        # surface phases
        print("----")
        print(f"  Number of surface phases = {self.num_phases - 1}")
        # site phase information
        print("----")
        print(f"  Number of Site phases = {self.num_site_phase}")
        for n, p in enumerate(self.PhaseSymbol):
            if n == 0:
                # skip the gas phase
                continue
            if self.phases[p].phase_type == "site":
                print(f"    Site phase: {n} {p} global ID: {self.phase_map[p]}")
        # site species
        print(f"  Number of Site species on phase = {self.num_site_species}")
        if self.num_site_species > 0:
            print("----")
            for k, map in enumerate(self.site_species_map.items()):
                print(
                    f"    Site species index = {k} {map[0]} "
                    f"global ID: {map[1]} WT = {self.site_wt[k]}"
                )
        # bulk phase information
        print("----")
        print(f"  Number of Bulk phases = {self.num_bulk_phase}")
        for n, p in enumerate(self.PhaseSymbol):
            if n == 0:
                # skip the gas phase
                continue
            if self.phases[p].phase_type == "bulk":
                print(f"    Bulk phase: {n} {p} global ID: {self.phase_map[p]}")
        # bulk species
        print(f"  Number of Bulk species = {self.num_bulk_species}")
        if self.num_bulk_species > 0:
            print("----")
            for k, map in enumerate(self.bulk_species_map.items()):
                print(
                    f"    Bulk species index = {k} {map[0]} "
                    f"global ID: {map[1]} WT = {self.bulk_wt[k]}"
                )
        # surface reaction information
        print("----")
        print(f"  Number of surface reactions = {self.num_surf_reactions}")

    def get_site_specindex(self, symbol: str = "") -> tuple[int, int]:
        """Get the local and the global indices of a site species."""
        """
        Get the phase-wise and the global (material-wise) indices of a site species.

        Parameters
        ----------
            symbol: string
                species symbol

        Returns
        -------
            local_index: integer
                (0-based) species index of the phase
            global_index: integer
                (0-based) global index of the matieral (including the gas phase)
        """
        global_index = -1
        local_index = -1
        if self.num_site_species > 0:
            if len(symbol) > 0:
                global_index = self.site_species_map.get(symbol, -1)
                if global_index > 0:
                    first_id = self.first_site_species_index
                    local_index = global_index - first_id
        return global_index, local_index

    def get_bulk_specindex(self, symbol: str = "") -> tuple[int, int]:
        """Get the local and the global indices of a bulk species."""
        """
        Get the phase-wise and the global (material-wise) indices of a bulk species.

        Parameters
        ----------
            symbol: string
                species symbol

        Returns
        -------
            local_index: integer
                (0-based) species index of the phase
            global_index: integer
                (0-based) global index of the matieral (including the gas phase)
        """
        global_index = -1
        local_index = -1
        if self.num_bulk_species > 0:
            if len(symbol) > 0:
                global_index = self.bulk_species_map.get(symbol, -1)
                if global_index > 0:
                    first_id = self.first_bulk_species_index
                    local_index = global_index - first_id
        return global_index, local_index

    def get_surf_specindex(
        self, symbol: str = "", mode: str = "normal"
    ) -> tuple[int, int, str]:
        """Get the local and the global indices and type of a surface species."""
        """
        Get the phase-wise and the global (material-wise) species indices and
        the phase type of the surface species given by the symbol.

        Parameters
        ----------
            symbol: string
                species symbol
            mode: str, {'normal', 'silent'}
                "silent": no print out

        Returns
        -------
            global_index: integer
                (0-based) global index of the matieral (including the gas phase)
            local_index: integer
                (0-based) species index of the phase
            species_type: string, ['site', 'bulk']
                type of the species
        """
        species_type = ""
        # try site phase species
        global_index, local_index = self.get_site_specindex(symbol)
        if global_index > 0:
            species_type = "site"
        else:
            # try bulk phase species
            global_index, local_index = self.get_bulk_specindex(symbol)
            if global_index > 0:
                species_type = "bulk"
            else:
                if mode.lower() == "normal":
                    # failed to find the species
                    msg = [
                        Color.PURPLE,
                        "failed to find species",
                        symbol,
                        "on surface material",
                        self.label,
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
                else:
                    return -1, -1, ""
        return global_index, local_index, species_type

    @property
    def phase_names(self):
        """Get list of phase names."""
        """
        Get list of phase names on the surface material.

        Returns
        -------
            Phasesymbol: list of strings
                list of phase names on the surface material
        """
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        if self._phasenamesdone == 0:
            # recycle existing data
            buff_p = (LP_c_char * self.num_phases)()
            for i in range(0, self.num_phases):
                buff_p[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_p = ctypes.cast(buff_p, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetPhaseNames(chem_id, mat_id, pp_p)
            if ierr == 0:
                self.phase_map.clear()
                for index in range(0, len(buff_p)):
                    str_val = ctypes.cast(buff_p[index], c_char_p).value.decode()
                    self.phase_map[str_val] = index
                self._phasenamesdone == 1
            else:
                # failed to get species symbols
                msg = [Color.PURPLE, "failed to get phase names.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_p
        # convert string type
        mylist = list(self.phase_map.keys())
        self.PhaseSymbol.clear()
        for s in mylist:
            self.PhaseSymbol.append(s)
        del mylist
        return self.PhaseSymbol

    @property
    def site_species_names(self):
        """Set list of site species names."""
        """
        Set list of site species names on the surface phase.

        Returns
        -------
            Phasesymbol: list of strings
                list of site species names on the surface site phase
        """
        if self.num_site_species <= 0:
            msg = [Color.PURPLE, "no site species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        #
        if self._sitespeciesnamesdone == 0:
            buff_s = (LP_c_char * self.num_site_species)()
            for i in range(0, self.num_site_species):
                buff_s[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_s = ctypes.cast(buff_s, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetSiteNames(chem_id, mat_id, pp_s)
            if ierr == 0:
                site_symbols = []
                for index in range(0, len(buff_s)):
                    s_val = ctypes.cast(buff_s[index], c_char_p).value.decode()
                    site_symbols.append(s_val)
            else:
                # failed to get surface site species symbols
                msg = [
                    Color.PURPLE,
                    "failed to get surface site species names.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_s
        else:
            site_symbols = self.site_symbols
        return site_symbols

    @property
    def bulk_species_names(self):
        """Set list of bulk species names."""
        """
        Set list of bulk species names on the surface phase.

        Returns
        -------
            Phasesymbol: list of strings
                list of bulk species names on the surface bulk phase
        """
        if self.num_bulk_species <= 0:
            msg = [Color.PURPLE, "no bulk species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        #
        if self._bulkspeciesnamesdone == 0:
            buff_b = (LP_c_char * self.num_bulk_species)()
            for i in range(0, self.num_bulk_species):
                buff_b[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_b = ctypes.cast(buff_b, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetBulkNames(chem_id, mat_id, pp_b)
            if ierr == 0:
                bulk_symbols = []
                for index in range(0, len(buff_b)):
                    b_val = ctypes.cast(buff_b[index], c_char_p).value.decode()
                    bulk_symbols.append(b_val)
            else:
                # failed to get surface site species symbols
                msg = [
                    Color.PURPLE,
                    "failed to get surface bulk species names.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_b, pp_b
        else:
            bulk_symbols = self.bulk_symbols
        return bulk_symbols

    def get_site_density(self) -> npt.NDArray[np.double]:
        """Get the site phase densities."""
        """
        Get the densities of all site phases of the material.

        Returns
        -------
            density: 1-D double array, dimension = number of phases of the material
                site phase density [mole/cm2]
        """
        if self.num_site_species <= 0:
            msg = [Color.PURPLE, "no site species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        site_density = np.zeros(self.num_phases, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetSiteDensity(chem_id, mat_id, site_density)
        if ierr != 0:
            # failed to get surface site species density
            msg = [Color.PURPLE, "failed to get surface site density.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return site_density

    def get_bulk_density(self) -> npt.NDArray[np.double]:
        """Get bulk species densities."""
        """
        Get the densities of all bulk species of the material.

        Returns
        -------
            density: 1-D double array,
                dimension = number of all bulk species of the material
                bulk species density [g/cm3]
        """
        if self.num_bulk_species <= 0:
            msg = [Color.PURPLE, "no bulk species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        density = np.zeros(self.num_bulk_species, dtype=np.double)
        # "fake" pressure and temperature since density of
        # solid bulk species is ~ constant
        temp = c_double(300.0)
        pres = c_double(P_ATM)
        ierr = ck_wrapper.chemkin.KINGetBulkDensity(
            chem_id, mat_id, pres, temp, density
        )
        if ierr != 0:
            # failed to get surface bulk species density
            msg = [Color.PURPLE, "failed to get surface bulk density.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        return density

    def get_site_molar_weights(self) -> npt.NDArray[np.double]:
        """Get the site species molecular weights."""
        """
        Get the molecular weights of all site species of the material.

        Returns
        -------
            wt: 1-D double array,
                dimension = number of all site species of the material
                molecular weights [g/mole]
        """
        if self.num_site_species == 0:
            wt = np.zeros(1, dtype=np.double)
            msg = [
                Color.MAGENTA,
                "no site species on material",
                self.label,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            return wt
        if self.sitewt_done == 0:
            chem_id = c_int(self._chemset_index)
            mat_id = c_int(self._material_index)
            wt = np.zeros(self.num_site_species, dtype=np.double)
            ierr = ck_wrapper.chemkin.KINGetSiteMolecularWeights(chem_id, mat_id, wt)
            if ierr != 0:
                # failed to get surface site species molecular weights
                msg = [
                    Color.PURPLE,
                    "failed to get site species molecular weights.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            # check for non-positive molar weights (possible for surface chemistry)
            for k, w in enumerate(wt):
                if w > 0.0:
                    continue
                msg = [
                    Color.YELLOW,
                    "site species",
                    self.site_symbols[k],
                    "has non-positive molecular weight.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            self.sitewt_done = 1
        else:
            wt = copy.deepcopy(self.site_wt)
        return wt

    def get_bulk_molar_weights(self) -> npt.NDArray[np.double]:
        """Get the bulk species molecular weights."""
        """
        Get the molecular weights of all bulk species of the material.

        Returns
        -------
            wt: 1-D double array,
                dimension = number of all site species of the material
                molecular weights [g/mole]
        """
        if self.num_bulk_species == 0:
            wt = np.zeros(1, dtype=np.double)
            msg = [
                Color.MAGENTA,
                "no bulk species on material",
                self.label,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            return wt
        if self.bulkwt_done == 0:
            chem_id = c_int(self._chemset_index)
            mat_id = c_int(self._material_index)
            wt = np.zeros(self.num_bulk_species, dtype=np.double)
            ierr = ck_wrapper.chemkin.KINGetBulkMolecularWeights(chem_id, mat_id, wt)
            if ierr != 0:
                # failed to get bulk species molecular weights
                msg = [
                    Color.PURPLE,
                    "failed to get bulk species molecular weights.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            # check for non-positive molar weights (possible for surface chemistry)
            for k, w in enumerate(wt):
                if w > 0.0:
                    continue
                msg = [
                    Color.YELLOW,
                    "bulk species",
                    self.bulk_symbols[k],
                    "has non-positive molecular weight.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            self.bulkwt_done = 1
        else:
            wt = copy.deepcopy(self.bulk_wt)
        return wt

    @property
    def first_site_species_index(self) -> int:
        """Get the global species index of the first site species."""
        """
        Get the global species index of the first site species of the material.

        Returns
        -------
            spec_index: integer
                (0-based) global species index
        """
        if self.num_site_species <= 0:
            return -1
        else:
            return self.site_species_map.get(self.site_symbols[0], -1)

    @property
    def first_bulk_species_index(self) -> int:
        """Get the global species index of the first bulk species."""
        """
        Get the global species index of the first bulk species of the material.

        Returns
        -------
            spec_index: integer
                (0-based) global species index
        """
        if self.num_bulk_species <= 0:
            return -1
        else:
            return self.bulk_species_map.get(self.bulk_symbols[0], -1)

    @property
    def last_site_species_index(self) -> int:
        """Get the global species index of the last site species."""
        """
        Get the global species index of the last site species of the material.

        Returns
        -------
            spec_index: integer
                (0-based) global species index
        """
        if self.num_site_species <= 0:
            return -1
        else:
            id = self.num_site_species - 1
            return self.site_species_map.get(self.site_symbols[id], -1)

    @property
    def last_bulk_species_index(self) -> int:
        """Get the global species index of the last bulk species."""
        """
        Get the global species index of the last bulk species of the material.

        Returns
        -------
            spec_index: integer
                (0-based) global species index
        """
        if self.num_bulk_species <= 0:
            return -1
        else:
            id = self.num_bulk_species - 1
            return self.bulk_species_map.get(self.bulk_symbols[id], -1)

    def get_site_occupancy(self) -> npt.NDArray[np.int32]:
        """Get site species occupancy."""
        """
        Get site species occupancy.

        Returns
        -------
            occupancy: 1-D integer array,
                dimension = total number of site species of the material
                site occupancy of the site species
        """
        if self.num_site_species <= 0:
            msg = [Color.PURPLE, "no site species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        occ = np.zeros(self.total_species, dtype=np.int32)
        ierr = ck_wrapper.chemkin.KINGetSpeciesOccupancy(chem_id, mat_id, occ)
        # get the global index of the first site species
        kstart = self.first_site_species_index
        # get the global index of the first bulk species (= last site species index + 1)
        kstop = kstart + self.num_site_species
        if kstart < 0:
            ierr = 10
        if ierr == 0:
            site_occupancy = np.zeros(self.num_site_species, dtype=np.int32)
            i = 0
            for k in range(kstart, kstop):
                site_occupancy[i] = occ[k]
                i += 1
            del occ
            return site_occupancy
        else:
            # failed to get bulk species molecular weights
            msg = [
                Color.PURPLE,
                "failed to get site species occupancy.\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_site_species_h(
        self, temp: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get site species enthalpy."""
        """
        Get site species enthalpy.

        Parameters
        ----------
            temp: double, optional
                surface material temperature [K]

        Returns
        -------
            h: 1-D double array,
                dimension = total number of site species of the material
                enthalpy of the site species [ergs/mol]
        """
        if self.num_site_species <= 0:
            msg = [Color.PURPLE, "no site species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        h_all = np.zeros(self.total_species, dtype=np.double)
        #
        if temp is None:
            t = max(self.surftemp, 3.0e2)
            tt = set_temp_array(temp=t)
        else:
            tt = set_temp_array(temp=temp)
        ierr = ck_wrapper.chemkin.KINGetAllSpeciesEnthalpy(chem_id, mat_id, tt, h_all)
        if ierr == 0:
            # extract enthalpy values for the site species only [ergs/mol]
            kstart = self.first_site_species_index
            kstop = kstart + self.num_site_species
            h = h_all[kstart:kstop]
            #
            del h_all
            return h
        else:
            # failed to compute enthalpies
            msg = [
                Color.PURPLE,
                "failed to compute site species enthalpies.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_bulk_species_h(
        self, temp: Union[float, None] = None, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get bulk species enthalpy."""
        """
        Get bulk species enthalpy.

        Parameters
        ----------
            temp: double, optional
                surface material temperature [K]
            pres: double, optional
                bulk pressure [dynes/cm2], it might be different from the gas pressure

        Returns
        -------
            h: 1-D double array,
                dimension = total number of bulk species of the material
                enthalpy of the bulk species [ergs/mol]
        """
        if self.num_bulk_species <= 0:
            msg = [Color.PURPLE, "no bulk species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        #
        if temp is None:
            t = max(self.surftemp, 3.0e2)
            tt = set_temp_array(temp=t)
        else:
            tt = set_temp_array(temp=temp)
        #
        if pres is None:
            # bulk pressure is not given
            h_all = np.zeros(self.total_species, dtype=np.double)
            ierr = ck_wrapper.chemkin.KINGetAllSpeciesEnthalpy(
                chem_id, mat_id, tt, h_all
            )
            if ierr == 0:
                # extract enthalpy values for the bulk species only [ergs/mol]
                kstart = self.first_bulk_species_index
                kstop = kstart + self.num_bulk_species
                h = h_all[kstart:kstop]
                #
                del h_all
                return h
            else:
                # failed to compute enthalpies
                msg = [
                    Color.PURPLE,
                    "failed to compute bulk species enthalpies.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        else:
            # bulk pressure is given
            pp = c_double(pres)
            h = np.zeros(self.num_bulk_species, dtype=np.double)
            ierr = ck_wrapper.chemkin.KINGetBulkSpeciesEnthalpy(
                chem_id, mat_id, pp, tt, h
            )
            if ierr == 0:
                # convert [ergs/gm] to [ergs/mol]
                for k in range(len(h)):
                    if self.bulk_wt[k] <= 0.0e0:
                        h[k] = 0.0e0
                    else:
                        h[k] = h[k] * self.bulk_wt[k]
                return h
            else:
                # failed to compute enthalpies
                msg = [
                    Color.PURPLE,
                    "failed to compute bulk species enthalpies.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    def get_bulk_species_cp(
        self, temp: Union[float, None] = None, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get bulk species specific heat capacity."""
        """
        Get bulk species specific heat capacity.

        Parameters
        ----------
            temp: double, optional
                surface material temperature [K]
            pres: double, optional
                bulk pressure [dynes/cm2], it might be different from the gas pressure

        Returns
        -------
            cp: 1-D double array,
                dimension = total number of bulk species of the material
                specific heat capacity of the bulk species [ergs/mol-K]
        """
        if self.num_bulk_species <= 0:
            msg = [Color.PURPLE, "no bulk species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        cp = np.zeros(self.num_bulk_species, dtype=np.double)
        #
        if pres is None:
            pp = c_double(P_ATM)
        else:
            pp = c_double(pres)
        #
        if temp is None:
            t = max(self.surftemp, 3.0e2)
            tt = set_temp_array(temp=t)
        else:
            tt = set_temp_array(temp=temp)
        #
        ierr = ck_wrapper.chemkin.KINGetBulkSpeciesSpecificHeat(
            chem_id, mat_id, pp, tt, cp
        )
        if ierr == 0:
            # convert [ergs/gm-K] to [ergs/mol-K]
            for k in range(len(cp)):
                if self.bulk_wt[k] <= 0.0e0:
                    cp[k] = 0.0e0
                else:
                    cp[k] = cp[k] * self.bulk_wt[k]
        else:
            # failed to compute enthalpies
            msg = [
                Color.PURPLE,
                "failed to compute bulk species specific heat capacity.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return cp

    def material_species_composition(
        self, elemindex: int = -1, specindex: int = -1
    ) -> int:
        """Get elemental composition of a species."""
        """
        Get elemental composition of a species.

        Parameters
        ----------
            elemindex: integer
                (0-based) index of the element
            specindex: integer
                (0-based) index of the species
                (all species of the material including the gas phase species)

        Returns
        -------
            count: integer
                number of the element in the given species
        """
        if self.ncf_done == 0:
            chem_id = c_int(self._chemset_index)
            mat_id = c_int(self._material_index)
            # initialize the NCF matrix
            dim = (self.num_element, self.total_species)
            self.elementalcomp = np.zeros(dim, dtype=np.int32, order="F")
            # load the NCF matrix
            ierr = ck_wrapper.chemkin.KINGetSurfaceSpeciesComposition(
                chem_id, mat_id, self.elementalcomp
            )
            if ierr != 0:
                msg = [
                    Color.PURPLE,
                    "failed to compute elemental compositions.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                self.ncf_done = 1

        # check element index
        if elemindex < 0 or elemindex >= self.num_element:
            msg = [Color.PURPLE, "element index is out of bound.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check species index
        if specindex < 0 or specindex >= self.total_species:
            msg = [Color.PURPLE, "species index is out of bound.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return self.elementalcomp[elemindex][specindex]

    @property
    def surface_temperature(self) -> float:
        """Get surface material temperature."""
        """
        Get the temperature of the material.

        Returns
        -------
            temp: double
                material surface temperature [K]
        """
        return self.surftemp

    @surface_temperature.setter
    def surface_temperature(self, temp: float):
        """Set surface material temperature."""
        """
        Set the surface temperature of the material.

        Parameters
        ----------
            temp: double
                material surface temperature [K]
        """
        if temp > 1.0e2:
            self.surftemp = temp
        else:
            msg = [Color.PURPLE, "surface temperature must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def surface_area(self) -> float:
        """Get surface material area."""
        """
        Get the active surface area of the material.

        Returns
        -------
            area: double
                material surface area [cm2]
        """
        return self.activearea

    @surface_area.setter
    def surface_area(self, area: float):
        """Set surface material area."""
        """
        Set the active surface area of the material.

        Parameters
        ----------
            area: double
                material surface area [cm2]
        """
        if area > 0.0e0:
            self.activearea = area
        else:
            msg = [Color.PURPLE, "surface area must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_surface_reaction_parameters(
        self,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get the Arrhenius reaction rate parameters."""
        """
        Get the Arrhenius reaction rate parameters of all surface reactions
        of the material.

        Returns
        -------
            a_factor: 1-D double array
                Arrhenius pre-exponent A-Factor of reaction in [mole-cm3-sec-K]
            beta: 1-D double array
                Arrhenius temperature exponent [-]
            act_energy: 1-D double array
                activation temperature [K]
        """
        reactionsize = self.num_surf_reactions
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # pre-exponent A factor of all surface reactions in the mechanism
        # in cgs units [mole-cm-sec-K]
        a_factor = np.zeros(shape=reactionsize, dtype=np.double)
        # temperature exponent of all reactions [-]
        beta = np.zeros_like(a_factor, dtype=np.double)
        # activation energy/temperature of all reactions [K]
        act_energy = np.zeros_like(a_factor, dtype=np.double)
        # flag for using A factor as the sticking coefficient
        is_stick = np.zeros_like(a_factor, dtype=np.int32)
        # get the reaction parameters
        ierr = ck_wrapper.chemkin.KINGetSurfaceReactionRateParameters(
            chem_id, mat_id, a_factor, beta, act_energy, is_stick
        )
        if ierr != 0:
            a_factor[:] = 0.0e0
            beta[:] = 0.0e0
            act_energy[:] = 0.0e0
        return a_factor, beta, act_energy

    def set_surface_reaction_afactor(self, reaction_index: int, a_factor: float):
        """Set the Arrhenius A-Factor of the given surface reaction."""
        """
        (Re)set the Arrhenius A-Factor of the given surface reaction of the material.

        Parameters
        ----------
            reaction_index: integer
                (1-based) index of the surface reaction of which
                the A-Factor to be reset
            a_factor: double
                new A-Factor value in [mole-cm-sec-K]
        """
        # check inputs
        if reaction_index > self.num_surf_reactions or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.num_surf_reactions) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if a_factor < 0.0e0:
            msg = [Color.PURPLE, "A-Factor must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        ireac = c_int(-reaction_index)  # negative index to "put" A-factor value
        ierr = ck_wrapper.chemkin.KINSetAFactorForASurfaceReaction(
            chem_id, mat_id, ireac, c_double(a_factor)
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to set Arrhenius A-Factor,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_surface_reaction_afactor(self, reaction_index: int) -> float:
        """Get the Arrhenius A-Factor of the given surface reaction."""
        """
        Get the Arrhenius A-Factor of the given surface reaction of the material.

        Parameters
        ----------
            reaction_index: integer
                (1-based) index of the reaction

        Returns
        -------
            AFactor: double
                Arrhenius A-Factor of the given reaction in [mole-cm3-sec-K]
        """
        # initialization
        a_factor = c_double(0.0e0)
        # check inputs
        if reaction_index > self.num_surf_reactions or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.num_surf_reactions) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        # get the A-factor value
        ierr = ck_wrapper.chemkin.KINSetAFactorForASurfaceReaction(
            chem_id, mat_id, ireac, a_factor
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to find Arrhenius A-Factor,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return a_factor.value

    def set_surface_reaction_act_energy(self, reaction_index: int, act_energy: float):
        """Set the Arrhenius activation energy of the given surface reaction."""
        """
        (Re)set the Arrhenius activation energy of
        the given surface reaction of the material.

        Parameters
        ----------
            reaction_index: integer
                (1-based) index of the surface reaction of which
                the activation energy to be reset
            act_energy: double
                new Aactivation eneergy value in [K]
        """
        # check inputs
        if reaction_index > self.num_surf_reactions or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.num_surf_reactions) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # if act_energy < 0.0e0:
        #     msg = [Color.PURPLE, "activation energy must >= 0.", Color.END]
        #     this_msg = Color.SPACE.join(msg)
        #     logger.error(this_msg)
        #     exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        # use negative index to "put" activation energy value
        ireac = c_int(-reaction_index)
        ierr = ck_wrapper.chemkin.KINSetActEnergyForASurfaceReaction(
            chem_id, mat_id, ireac, c_double(act_energy)
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to set Arrhenius activation energy,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_surface_reaction_act_energy(self, reaction_index: int) -> float:
        """Get the Arrhenius activation energy."""
        """
        Get the Arrhenius activation energy of the given surface
        reaction of the material.

        Parameters
        ----------
            reaction_index: integer
                (1-based) index of the reaction

        Returns
        -------
            act_energy: double
                Arrhenius activation energy of the given reaction in [K]
        """
        # initialization
        act_energy = c_double(0.0e0)
        # check inputs
        if reaction_index > self.num_surf_reactions or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.num_surf_reactions) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        # get the activation energy value
        ierr = ck_wrapper.chemkin.KINSetActEnergyForASurfaceReaction(
            chem_id, mat_id, ireac, act_energy
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to find Arrhenius activation energy,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return act_energy.value

    def get_surface_reaction_string(self, reaction_index: int) -> str:
        """Get the reaction string of the surface reaction."""
        """
        Get the reaction string of the surface reaction
        of the material specified by the reaction index.
        The surface reaction index starts from 1.

        Parameters
        ----------
            reaction_index: integer
                (1-based) surface reaction index

        Returns
        -------
            reactionstring: string
                reaction string of the given surface reaction
        """
        # initialization
        reactionstring = ""
        if reaction_index > self.num_surf_reactions:
            msg = [
                Color.PURPLE,
                "reaction index must be <",
                str(self.num_surf_reactions),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        elif reaction_index <= 0:
            msg = [Color.PURPLE, "reaction index must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        i_string_size = c_int(0)
        # get reaction string (might have to be increased to 2048 for 26R1)
        rstring = bytes(" " * 1024, "utf-8")
        ierr = ck_wrapper.chemkin.KINGetSurfaceReactionString(
            chem_id, mat_id, ireac, i_string_size, rstring
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to find reaction string,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert C string back to string
        # print(rstring.decode()[0:iStringSize.value])
        reactionstring = rstring.decode()[0 : i_string_size.value]
        del rstring
        return reactionstring


class SurfacePhase:
    """Surface phase object."""

    def __init__(self, phase_index: int, phase_type: str, label: str = ""):
        """Create a surface phase object."""
        """The surface phase can be either a site phase or a bulk phase."""
        self.phase_index = phase_index
        if len(label) > 0:
            self.label = label
        else:
            self.label = "phase_" + str(phase_index)
        if phase_type in ["site", "bulk"]:
            self.phase_type = phase_type
        else:
            msg = [
                Color.PURPLE,
                "invalid phase type,",
                "phase type must be either 'site' or 'bulk'",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self._numb_species = 0
        self.first_spec_index = 0
        self.last_spec_index = 0
        self._density = 0.0e0  # site density [mole/cm2] for site phases only

    @property
    def number_species(self) -> int:
        """Get the number of species."""
        """
        Get the number of species on this surface phase.

        Returns
        -------
            numb_species: integer
                number of species
        """
        return self._numb_species

    @number_species.setter
    def number_species(self, count: int):
        """Set the number of species."""
        """
        Set the number of species on this surface phase.

        Parameters
        ----------
            count: integer
                number of species
        """
        if count > 0:
            self._numb_species = count
        else:
            msg = [Color.PURPLE, "number of species must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def site_density(self) -> float:
        """Get surface site density."""
        """
        Get surface site density if this is a surface site phase.

        Returns
        -------
            density: double
                site phase density [mole/cm2]
        """
        return self._density

    @site_density.setter
    def site_density(self, den: float):
        """Set the site density."""
        """
        Set the site density if this is a surface site phase.

        Parameters
        ----------
            den: double
                site phase density [mole/cm2]
        """
        if den > 0.0e0:
            self._density = den
        else:
            msg = [Color.PURPLE, "surface site phase density must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def first_species_index(self) -> int:
        """Get the global index of the first species."""
        """
        Get the global index (1-based) of the first species of the phase.

        Returns
        -------
            global_index: integer
                the (1-based) global index of the first species
        """
        return self.first_spec_index

    @first_species_index.setter
    def first_species_index(self, id: int):
        """Set the global index of the first species."""
        """
        Set the global index (1-based) of the first species of the surface phase.

        Parameters
        ----------
            id: integer
                the (1-based) global index of the first species
        """
        if id > 0:
            self.first_spec_index = id
        else:
            msg = [Color.PURPLE, "global index must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def last_species_index(self) -> int:
        """Get the global index of the last species."""
        """
        Get the global index (1-based) of the last species of the surface phase.

        Returns
        -------
            global_index: integer
                the (1-based) global index of the last species
        """
        return self.last_spec_index

    @last_species_index.setter
    def last_species_index(self, id: int):
        """Set the global index of the last species."""
        """
        Set the global index (1-based) of the last species of the surface phase.

        Parameters
        ----------
            id: integer
                the (1-based) global index of the last species
        """
        if id > 0:
            self.last_spec_index = id
        else:
            msg = [Color.PURPLE, "global index must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

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

"""
    Chemkin Chemistry utilities.
"""

import ctypes
from ctypes import POINTER, c_char_p, c_double, c_int
import os
from typing import Dict, List, Union

from ansys.chemkin import chemkin_wrapper as ck_wrapper
from ansys.chemkin.color import Color
from ansys.chemkin.constants import RGas
from ansys.chemkin.info import clear_hints
from ansys.chemkin.logger import logger
from ansys.chemkin.realgaseos import check_realgas_status, set_current_pressure
import numpy as np
import numpy.typing as npt

_symbol_length = 16  # Chemkin element/species symbol length
MAX_SPECIES_LENGTH = _symbol_length + 1  # Chemkin element/species symbol length + 1
LP_c_char = ctypes.POINTER(ctypes.c_char)  # pointer to C type character array
COMPLETE = 0

_chemset_identifiers: List = (
    []
)  # string used to identify different chemistry sets in the same project
_active_chemistry_set = -10
chemkin_verbose = True  # verbose mode to turn ON/OFF the print statements that do not have the leading '**' characters
_CKInitialized: Dict = {}  # Chemkin-CFD-API initialization flag for every Chemistry Set
# == end of global parameters


#
# Chemkin module level methods
#
def verbose() -> bool:
    """
    return the global verbose mode indicating the status (ON/OFF) of printing statements that do not have the leading '**' characters

    Returns
    -------
        mode: boolean, {True, False}, default = True
            the verbose mode
    """
    global chemkin_verbose
    return chemkin_verbose


def set_verbose(OnOff: bool):
    """
    set the global verbose mode to turn ON(True) or OFF(False) of printing statements that do not have the leading '**' characters

    Parameters
    ----------
        OnOff: boolean, {True, False}
            the verbose mode
    """
    global chemkin_verbose
    chemkin_verbose = OnOff


def done():
    """
    Release Chemkin license and reset the Chemistry sets
    """
    # terminate
    ck_wrapper.chemkin.KINFinish()
    # release
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
    """
    check whether the Chemistry Set is initialized in Chemkin-CFD-API

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set

    Returns
    -------
        status: boolean
            the initialization status of the Chemistry set associated with the given Chemistry set index
    """
    global _CKInitialized
    status = _CKInitialized.get(chem_index, False)
    return status


def activate_chemistryset(chem_index: int) -> int:
    """
    Switch to (re-activate) the work spaces of the current Chemistry Set
    when there are multiple Chemistry Sets in the same project

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set

    Returns
    -------
        error flag: integer
    """
    iErr = ck_wrapper.chemkin.KINSwitchChemistrySet(c_int(chem_index))
    if iErr == 0:
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
    return iErr


def force_activate_chemistryset(chem_index: int):
    """
    activate the Chemistry Set automatically and silently.

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set
    """
    if not check_active_chemistryset(chem_index):
        # the Chemistry Set is not currently active
        iErr = activate_chemistryset(chem_index)
        if iErr != 0:
            exit()


def chemistryset_new(chem_index: int):
    """
    Create a new Chemistry Set initialization flag and set the value to False

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set
    """
    global _CKInitialized
    _CKInitialized[chem_index] = False
    logger.debug("preprocessing done")


def chemistryset_initialized(chem_index: int):
    """
    Set the Chemistry Set Initialization flag to True

    Parameters
    ----------
        chem_index: integer
            chemistry set index associated with the Chemistry Set
    """
    global _CKInitialized
    _CKInitialized[chem_index] = True
    logger.debug(">>> Chemkin-CFD-API initialized <<<")


def check_active_chemistryset(chem_index: int) -> bool:
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


class Chemistry:
    """
    define and preprocess Chemkin chemistry set
    """

    realgas_CuEOS = [
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
        """
        Create a Chemistry object based on given Chemkin mechanism input files, thermodynamic data file,
        and transport data file.

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
        self._num_materials = c_int(0)  # number of materials
        self._num_max_site_species = c_int(0)  # total number of surface site species
        self._num_max_bulk_species = c_int(0)  # total number of bulk species
        self._num_max_phases = c_int(0)  # number of phases
        self._num_max_surf_reactions = c_int(0)  # total number of surface reactions
        self._gas_species: dict[str, int] = {}  # gas species symbols dictionary
        self._elements: dict[str, int] = {}  # element symbols dictionary
        self._EOS = c_int(0)  # equation of state (default 0 = ideal gas)
        self.userealgas = False  # use ideal gas law by default
        self._AWTdone = 0
        self._WTdone = 0
        self._NCFdone = 0
        # fake initialization
        self._AWT = np.zeros(1, dtype=np.double)
        self._WT = np.zeros(1, dtype=np.double)
        self._KSYMdone = 0
        self._ESYMdone = 0
        self.KSymbol: list[str] = []
        self.ESymbol: list[str] = []
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
        if os.path.isfile(self._surf_file):
            self._index_surf = ctypes.c_int(1)
        # check transport data file
        if os.path.isfile(self._tran_file):
            self._index_tran = ctypes.c_int(1)

    @property
    def chemfile(self) -> str:
        """
        Get gas-phase mechanism file name of this chemistry set

        Returns
        -------
            chemfile: string
                Full path and name of the Chemkin gas-phase mechanism input file
        """
        return self._gas_file

    @chemfile.setter
    def chemfile(self, filename: str):
        """
        Assign the gas-phase mechanism filename

        Parameters
        ----------
            filename: string
                name of the gas-phase mechanism file with the full path
        """
        self._gas_file = filename

    @property
    def thermfile(self) -> str:
        """
        Get thermodynamic data filename of this chemistry set

        Returns
        -------
            thermfile: string
                Full path and name of the Chemkin thermodynamic data file
        """
        return self._therm_file

    @thermfile.setter
    def thermfile(self, filename: str):
        """
        Assign the thermodynamic data filename

        Parameters
        ----------
            filename: string
                name of the thermodynamic data file with the full path
        """
        self._therm_file = filename

    @property
    def tranfile(self) -> str:
        """
        Get transport data filename of this chemistry set

        Returns
        -------
            tranfile: string
                Full path and name of the Chemkin thransport data file
        """
        return self._tran_file

    @tranfile.setter
    def tranfile(self, filename: str):
        """
        Assign the transport data filename

        Parameters
        ----------
            filename: string
                name of the transport data file with the full path
        """
        self._tran_file = filename
        if os.path.isfile(self._tran_file):
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
        """
        Get the name of the summary file from the preprocessor

        Returns
        -------
            tranfile: string
                Full path and name of the preprocessing summary file
        """
        return self._summary_out

    def preprocess_transportdata(self):
        """
        Instruct the preprocessor to include the transport data
        """
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
        """
        Get surface mechanism filename of this chemistry set

        Returns
        -------
            tranfile: string
                Full path and name of the Chemkin surface mechanism input file
        """
        return self._surf_file

    @surffile.setter
    def surffile(self, filename: str):
        """
        Assign the surface mechanism filename

        Parameters
        ----------
            filename: string
                name of the surface mechanism file with the full path
        """
        self._surf_file = filename
        if os.path.isfile(self._surf_file):
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
        """
        Assign all input files of the chemistry set

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
            if os.path.isfile(self._surf_file):
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
            if os.path.isfile(self._tran_file):
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
        """
        Run Chemkin preprocessor

        Returns
        -------
            Error code: integer
        """
        # check minimum set of required files
        if not os.path.isfile(self._gas_file):
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
        if not os.path.isfile(self._therm_file):
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
            iErr = ck_wrapper.chemkin.KINGetChemistrySizes(
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

            if iErr != 0:
                # failed to find mechanism sizes
                msg = [
                    Color.RED,
                    "failed to find mechanism parameters\n",
                    Color.SPACEx6,
                    "error code =",
                    str(iErr),
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
            self._WT = self.WT
            # get atomic masses
            self._AWT = self.AWT
            # check real-gas model
            self.verify_realgas_model()
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
        """
        Verify the availability of real-gas data in the mechanism
        """
        EOSModel = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
        try:
            # check if the mechanism contains the real-gas EOS data
            iErr = ck_wrapper.chemkin.KINRealGas_GetEOSMode(
                self._chemset_index, self._EOS, EOSModel
            )
            if iErr == 0:
                if self._EOS.value > 0 and self._EOS.value <= 5:
                    msg = [
                        Color.YELLOW,
                        "real-gas cubic EOS",
                        "'" + str(EOSModel.value.decode()) + "'",
                        "is available.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
                    del EOSModel
                    return

            self._EOS = c_int(0)
        except OSError:
            # mechanism contains no real-gas data
            self._EOS = c_int(0)
            if verbose():
                msg = [Color.PURPLE, "accessing the real gas information.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)

        del EOSModel
        if verbose():
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

    def verify_transport_data(self) -> bool:
        """
        Verify the availability of transport property data in the mechanism

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
        """
        Verify the availability of surface chemistry data in the mechanism

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
        """
        Get list of gas species symbols

        Returns
        -------
            Ksymbol: list of strings
                list of species symbols in the gas-phase mechanism
        """
        global MAX_SPECIES_LENGTH
        if self._KSYMdone == 0:
            # recycle existing data
            buff = (LP_c_char * self._num_gas_species.value)()
            for i in range(0, self._num_gas_species.value):
                buff[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp = ctypes.cast(buff, POINTER(LP_c_char))
            iErr = ck_wrapper.chemkin.KINGetGasSpeciesNames(self._chemset_index, pp)
            if iErr == 0:
                self._gas_species.clear()
                for index in range(0, len(buff)):
                    strVal = ctypes.cast(buff[index], c_char_p).value.decode()
                    self._gas_species[strVal] = index
                self._KSYMdone == 1
            else:
                # failed to get species symbols
                msg = [Color.PURPLE, "failed to get species symbols.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff

        # convert string type
        mylist = list(self._gas_species.keys())
        self.KSymbol.clear()
        for s in mylist:
            self.KSymbol.append(s)
        del mylist
        return self.KSymbol

    @property
    def element_symbols(self):
        """
        Get the list of element symbols

        Returns
        -------
            Esymbol: list of strings
                list of element symbols in the mechanism
        """
        if self._ESYMdone == 0:
            buff_ele = (LP_c_char * self._num_elements.value)()
            for j in range(0, self._num_elements.value):
                buff_ele[j] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_ele = ctypes.cast(buff_ele, POINTER(LP_c_char))
            iErr = ck_wrapper.chemkin.KINGetElementNames(self._chemset_index, pp_ele)
            if iErr == 0:
                self._elements.clear()
                for index in range(0, len(buff_ele)):
                    eleVal = ctypes.cast(buff_ele[index], c_char_p).value.decode()
                    self._elements[eleVal] = index
                self._ESYMdone == 1
            else:
                # failed to get element symbols
                msg = [Color.PURPLE, "failed to get element symbols.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_ele

        # convert string type
        my_ele_list = list(self._elements.keys())
        self.ESymbol.clear()
        for s_ele in my_ele_list:
            self.ESymbol.append(s_ele)
        del my_ele_list
        return self.ESymbol

    def get_specindex(self, specname: str) -> int:
        """
        Get index of the gas species

        Returns
        -------
            specindex: integer
                index of the given species symbols in the gas-phase mechanism
        """
        specindex = self._gas_species.get(specname, -1)
        if specindex <= 0:
            msg = [Color.PURPLE, "species symbol not found:", specname, Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return specindex

    @property
    def chemID(self) -> int:
        """
        Get chemistry set index

        Returns
        -------
            chemIDx: integer
                index of the Chemistry set
        """
        if self._chemset_index.value >= 0:
            return self._chemset_index.value
        else:
            return -1

    @property
    def surfchem(self) -> int:
        """
        Get surface chemistry status

        Returns
        -------
            status: integer, {0, 1}
                indicating whether the Chemistry set includes a surface mechanism
                0 = this chemistry set does NOT include a surface chemistry
                1 = this chemistry set includes a  surface chemistry
        """
        return self._index_surf.value

    @property
    def KK(self) -> int:
        """
        Get number of gas species

        Returns
        -------
            KK: integer
                total number of gas-phase species in the Chemistry set
        """
        return self._num_gas_species.value

    # alias
    number_species = KK

    @property
    def MM(self) -> int:
        """
        Get number of elements in the chemistry set

        Returns
        -------
            MM: integer
                total number of elements in the Chemistry set
        """
        return self._num_elements.value

    # alias
    number_elements = MM

    @property
    def IIGas(self) -> int:
        """
        Get number of gas-phase reactions

        Returns
        -------
            IIGas: integer
                total number of gas-phase reactions in the Chemistry set
        """
        return self._num_gas_reactions.value

    # alias
    number_gas_reactions = IIGas

    @property
    def AWT(self) -> npt.NDArray[np.double]:
        """
        compute atomic masses

        Returns
        -------
            AWT: 1-D double array
                masses of the elements in the Chemistry set [g/mole]
        """
        if self._AWTdone == 1:
            return self._AWT
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        del self._AWT  # clear the "original" definition in __init__
        self._AWT = np.zeros(self._num_elements.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetAtomicWeights(self._chemset_index, self._AWT)
        if iErr == 0:
            self._AWTdone = 1
        else:
            # failed to find atomic masses
            msg = [Color.PURPLE, "failed to get atomic masses.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return self._AWT

    # alias
    atomic_weight = AWT

    @property
    def WT(self) -> npt.NDArray[np.double]:
        """
        compute gas species molecular masses

        Returns
        -------
            WT: 1-D double array
                molecular masses of the gas-phase species in the Chemistry set [g/mole]
        """
        if self._WTdone == 1:
            return self._WT
        if self._chemset_index.value < 0:
            msg = [
                Color.PURPLE,
                "please preprocess the chemistry set first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        del self._WT  # clear the "original" definition in __init__
        self._WT = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasMolecularWeights(
            self._chemset_index, self._WT
        )
        if iErr == 0:
            self._WTdone = 1
        else:
            # failed to find molecular masses
            msg = [Color.PURPLE, "failed to get species molecular masses.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return self._WT

    # alias
    species_molar_weight = WT

    def SpeciesCp(
        self, temp: float = 0.0, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """
        Get species specific heat capacity at constant pressure

        Parameters
        ----------
            temp: double
                Temperature [K]
            pres: double, optional
                Pressure [dynes/cm2] required when real gas model is activated

        Returns
        -------
            Cp: 1-D double array
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
        if check_realgas_status(self.chemID):
            if pres is None:
                # pressure is not assigned
                msg = [
                    Color.PURPLE,
                    "mixture pressure must be provided",
                    "to evaluate real-gas species properties\n",
                    Color.SPACEx6,
                    "usage: <Chemistry_Obj>.SpeciesCp(temperature, pressure)",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # set current pressure for the real-gas
                set_current_pressure(self.chemID, pres)
        #
        TT = c_double(temp)
        Cp = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpecificHeat(self._chemset_index, TT, Cp)
        if iErr == 0:
            # convert [ergs/g-K] to [ergs/mol-K]
            # for k in range(len(Cp)):
            #    Cp[k] = Cp[k] * self._WT[k]
            Cp *= self._WT
        else:
            # failed to compute specific heats
            msg = [Color.PURPLE, "failed to compute specific heats.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return Cp

    def SpeciesCv(
        self, temp: float = 0.0, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """
        Get species specific heat capacity at constant volume (ideal gas only)

        Parameters
        ----------
            temp: double
                Temperature [K]
            pres: double, optional
                Pressure [dynes/cm2] required when real gas model is activated

        Returns
        -------
            Cv: 1-D double array
                species specific heat capacities at constant volume [ergs/mol-K]
        """
        if check_realgas_status(self.chemID) and pres is None:
            # pressure is not assigned
            msg = [
                Color.PURPLE,
                "mixture pressure must be provided",
                "to evaluate real-gas species properties\n",
                Color.SPACEx6,
                "usage: <Chemistry_Obj>.SpeciesCv(temperature, pressure)",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        Cv = self.SpeciesCp(temp, pres)
        R = RGas
        # for k in range(len(Cp)):
        #    Cv[k] = Cp[k] - R
        Cv -= R

        return Cv

    def SpeciesH(
        self, temp: float = 0.0, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """
        Get species enthalpy

        Parameters
        ----------
            temp: double
                Temperature [K]
            pres: double, optional
                Pressure [dynes/cm2] required when real gas model is activated

        Returns
        -------
            H: 1-D double array
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
        if check_realgas_status(self.chemID):
            if pres is None:
                # pressure is not assigned
                msg = [
                    Color.PURPLE,
                    "mixture pressure must be provided",
                    "to evaluate real-gas species properties\n",
                    Color.SPACEx6,
                    "usage: <Chemistry_Obj>.SpeciesH(temperature, pressure)",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # set current pressure for the real-gas
                set_current_pressure(self.chemID, pres)
        TT = c_double(temp)
        H = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpeciesEnthalpy(self._chemset_index, TT, H)
        if iErr == 0:
            # convert [ergs/gm] to [ergs/mol]
            # for k in range(len(H)):
            #    H[k] = H[k], * self._WT[k]
            H *= self._WT
        else:
            # failed to compute enthalpies
            msg = [Color.PURPLE, "failed to compute specific enthalpies.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return H

    def SpeciesU(
        self, temp: float = 0.0, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """
        Get species internal energy

        Parameters
        ----------
            temp: double
                Temperature [K]
            pres: double, optional
                Pressure [dynes/cm2] required when real gas model is activated

        Returns
        -------
            U: 1-D double array
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
        if check_realgas_status(self.chemID):
            if pres is None:
                # pressure is not assigned
                msg = [
                    Color.PURPLE,
                    "mixture pressure must be provided",
                    "to evaluate real-gas species properties\n",
                    Color.SPACEx6,
                    "usage: <Chemistry_Obj>.SpeciesU(temperature, pressure)",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # set current pressure for the real-gas
                set_current_pressure(self.chemID, pres)
        TT = c_double(temp)
        U = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpeciesInternalEnergy(
            self._chemset_index, TT, U
        )
        if iErr == 0:
            # convert [ergs/gm] to [ergs/mol]
            # for k in range(len(U)):
            #    U[k] = U[k], * self._WT[k]
            U *= self._WT
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

        return U

    def SpeciesVisc(self, temp: float = 0.0) -> npt.NDArray[np.double]:
        """
        Get species viscosity

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
        TT = c_double(temp)
        visc = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetViscosity(self._chemset_index, TT, visc)
        if iErr != 0:
            # failed to compute viscosity
            msg = [Color.PURPLE, "failed to compute specific viscosities.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return visc

    def SpeciesCond(self, temp: float = 0.0) -> npt.NDArray[np.double]:
        """
        Get species conductivity

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
        TT = c_double(temp)
        cond = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetConductivity(self._chemset_index, TT, cond)
        if iErr != 0:
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

    def SpeciesDiffusionCoeffs(
        self, press: float = 0.0, temp: float = 0.0
    ) -> npt.NDArray[np.double]:
        """
        Get species diffusion coefficients

        Parameters
        ----------
            press: double
                Pressure [dynes/cm2]
            temp: double
                Temperature [K]

        Returns
        -------
            diffusioncoeffs: 2-D double array, dimension = [number_species, number_species]
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
        PP = c_double(press)
        TT = c_double(temp)
        dim = (self._num_gas_species.value, self._num_gas_species.value)
        diffusioncoeffs = np.zeros(dim, dtype=np.double, order="F")
        iErr = ck_wrapper.chemkin.KINGetDiffusionCoeffs(
            self._chemset_index, PP, TT, diffusioncoeffs
        )
        if iErr != 0:
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

    def SpeciesComposition(self, elemindex: int = -1, specindex: int = -1) -> int:
        """
        Get elemental composition of a species

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
        if self._NCFdone == 0:
            # initialize the NCF matrix
            dim = (self._num_elements.value, self._num_gas_species.value)
            self.elementalcomp = np.zeros(dim, dtype=np.int32, order="F")
            # load the NCF matrix
            iErr = ck_wrapper.chemkin.KINGetGasSpeciesComposition(
                self._chemset_index, self.elementalcomp
            )
            if iErr != 0:
                msg = [
                    Color.PURPLE,
                    "failed to compute elemental compositions.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                self._NCFdone = 1

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
    def EOS(self) -> int:
        """
        Get the available real-gas EOS model

        Returns
        -------
        count: integer
            number of available cubic EOS models in Chemkin
        """
        return self._EOS.value

    def use_realgas_cubicEOS(self):
        """
        Turn ON the real-gas cubic EOS to compute mixture properties if the mechanism contains necessary data
        """
        if self._EOS.value < 1:
            # no real gas EOS data in the mechanism
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return
        # check real-gas EOS status
        iFlag = c_int(0)
        iErr = ck_wrapper.chemkin.KINRealGas_UseCubicEOS(self._chemset_index, iFlag)
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to turn ON the real-gas EOS model,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if iFlag.value == 0:
            msg = [
                Color.YELLOW,
                "real-gas cubic EOS model",
                Chemistry.realgas_CuEOS[self._EOS.value],
                "is turned ON.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = True
        else:
            self.userealgas = False

    def use_idealgas_law(self):
        """
        Turn on the ideal gas law to compute mixture properties
        """
        if self._EOS.value < 1:
            # no real gas EOS data in the mechanism
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = False
            return
        # check real-gas EOS status
        iFlag = c_int(0)
        iErr = ck_wrapper.chemkin.KINRealGas_UseIdealGasLaw(self._chemset_index, iFlag)
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to turn ON ideal gas law,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if iFlag.value == 0:
            msg = [Color.YELLOW, "ideal gas law is turned ON.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = False

    def get_reaction_parameters(
        self,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double], npt.NDArray[np.double]]:
        """
        Get the Arrhenius reaction rate parameters of all gas-phase reactions

        Returns
        -------
            AFactor: 1-D double array
                Arrhenius pre-exponent A-Factor of reaction in [mole-cm3-sec-K]
            TBeta: 1-D double array
                Arrhenius temperature exponent [-]
            AEnergy: 1-D double array
                activation temperature [K]
        """
        reactionsize = self.IIGas
        # pre-exponent A factor of all gas-phase reactions in the mechanism in cgs units [mole-cm3-sec-K]
        AFactor = np.zeros(shape=reactionsize, dtype=np.double)
        # temperature exponent of all reactions [-]
        TBeta = np.zeros_like(AFactor, dtype=np.double)
        # activation energy/temperature of all reactions [K]
        AEnergy = np.zeros_like(AFactor, dtype=np.double)
        # get the reaction parameters
        iErr = ck_wrapper.chemkin.KINGetReactionRateParameters(
            self._chemset_index, AFactor, TBeta, AEnergy
        )
        if iErr != 0:
            AFactor[:] = 0.0e0
            TBeta[:] = 0.0e0
            AEnergy[:] = 0.0e0
        return AFactor, TBeta, AEnergy

    def set_reaction_AFactor(self, reaction_index: int, AFactor: float):
        """
        (Re)set the Arrhenius A-Factor of the given reaction

        Parameters
        ----------
            reaction_index: integer
                index of the gas-phase reaction of which the A-Factor to be reset
            AFctor: double
                new A-Factor value in [mole-cm3-sec-K]
        """
        # check inputs
        if reaction_index > self.IIGas or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.IIGas) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if AFactor < 0.0e0:
            msg = [Color.PURPLE, "A-Factor must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert the reaction parameters
        ireac = c_int(-reaction_index)  # negative index to "put" A-factor value
        iErr = ck_wrapper.chemkin.KINSetAFactorForAReaction(
            self._chemset_index, ireac, c_double(AFactor)
        )
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to set Arrhenius A-Factor,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_reaction_AFactor(self, reaction_index: int) -> float:
        """
        get the Arrhenius A-Factor of the given reaction

        Parameters
        ----------
            reaction_index: integer
                index of the reaction

        Returns
        -------
            AFactor: double
                Arrhenius A-Factor of the given reaction in [mole-cm3-sec-K]
        """
        # initialization
        AFactor = c_double(0.0e0)
        # check inputs
        if reaction_index > self.IIGas or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.IIGas) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        # get the A-factor value
        iErr = ck_wrapper.chemkin.KINSetAFactorForAReaction(
            self._chemset_index, ireac, AFactor
        )
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to find Arrhenius A-Factor,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return AFactor.value

    def get_gas_reaction_string(self, reaction_index: int) -> str:
        """
        Get the reaction string of the gas-phase reaction specified by the reaction index.

        Parameters
        ----------
            reaction_index: integer
                (base-1) gas-phase reaction index

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
        iStringSize = c_int(0)
        # get reaction string
        rstring = bytes(" " * 1024, "utf-8")
        iErr = ck_wrapper.chemkin.KINGetGasReactionString(
            self._chemset_index, ireac, iStringSize, rstring
        )
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to find reaction string,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert C string back to string
        # print(rstring.decode()[0:iStringSize.value])              # check
        reactionstring = rstring.decode()[0 : iStringSize.value]
        del rstring
        return reactionstring

    def save(self):
        """
        Store the work spaces of the current Chemistry Set
        if new Chemistry Set will be created later in the same project
        """
        iErr = ck_wrapper.chemkin.KINUpdateChemistrySet(self._chemset_index)
        if iErr == 0:
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
        """
        Switch to (re-activate) the work spaces of the current Chemistry Set
        when there are multiple Chemistry Sets in the same project
        """
        iErr = activate_chemistryset(self._chemset_index.value)
        if iErr == 0:
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

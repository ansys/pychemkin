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

"""Chemkin general reactor model utilities."""

import copy
import ctypes
from ctypes import c_double, c_int
from typing import Union

import numpy as np
import numpy.typing as npt

from . import chemkin_wrapper
from .chemistry import (
    check_active_chemistryset,
    check_chemistryset,
    chemistryset_initialized,
    verbose,
)
from .color import Color
from .constants import P_ATM
from .inlet import Stream
from .logger import logger
from .mixture import Mixture


#
# Base class for keyword data
#
class Keyword:
    """A Chemkin style keyword."""

    # supported Chemkin keyword data types
    _keyworddatatypes = ["bool", "int", "float", "str"]
    _valuetypes = (bool, int, float, str)
    # required keywords that are given as reactor properties or as mixture properties
    # and will be set by using the KINAll0D_SetupBatchInputs call
    _protectedkeywords = [
        "CONP",
        "CONV",
        "TRAN",
        "STST",
        "TGIV",
        "ENRG",
        "PRES",
        "TEMP",
        "TAU",
        "TIME",
        "XEND",
        "FLRT",
        "VDOT",
        "SCCM",
        "VDOT",
        "DIAM",
        "AREA",
        "REAC",
        "GAS",
        "INIT",
        "XEST",
        "SURF",
        "ACT",
        "TINL",
        "FUEL",
        "OXID",
        "PROD",
        "ASEN",
        "ATLS",
        "RTLS",
        "EPST",
        "EPSS",
    ]
    gasspecieskeywords = ["REAC", "XEST", "FUEL", "OXID"]
    flowratekeywords = ["FLRT", "VDOT", "VEL", "SCCM"]
    profilekeywords = [
        "TPRO",
        "PPRO",
        "VPRO",
        "QPRO",
        "AINT",
        "AEXT",
        "DPRO",
        "FPRO",
        "SCCMPRO",
        "VDOTPRO",
        "VELPRO",
        "TINPRO",
        "AFLO",
    ]
    fourspaces = "    "
    # Under the default API-call mode, important keywords (the _protectedkeywords)
    # are set by direct API calls, the rest of the keywords can be set by
    # keyword input lines (i.e., using the setkeyword method).
    # Under the full-keyword mode, all keywords and their parameters are set
    # by keyword input lines, and specifying those _protectedkeywords via
    # the setkeyword method are required.
    no_fullkeyword = True  # default: API-call mode

    def __init__(self, phrase: str, value: Union[float, bool, str], data_type: str):
        """Initialize the Chemkin keyword."""
        """
        Initialize the Chemkin keyword.

        Parameters
        ----------
            phrase: string
                Chemkin keyword phrase
            value: indicated by the data_type, {int, float, string, bool}
                value assigned to the Chemkin keyword
            data_type: string, {'int', 'float', 'string', or 'bool'}
                data type of value

        """
        self._set = False
        ierr = 0
        # check value data type
        if data_type not in Keyword._keyworddatatypes:
            # the declared data type is not supported
            msg = [
                Color.PURPLE,
                "unsupported data type specified",
                data_type,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            if not isinstance(value, (bool, int, float, str)):
                # value does not match the declared data type
                msg = [
                    Color.PURPLE,
                    "variable has different data type",
                    str(type(value)),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
            ierr = 1
        # block the protected keywords
        if Keyword.no_fullkeyword:
            if phrase.upper() in Keyword._protectedkeywords:
                msg = [
                    Color.PURPLE,
                    "use reactor property setter to assign",
                    phrase,
                    "value\n",
                    Color.SPACEx6,
                    "for example, to set the reactor volume use:",
                    '"MyBatchReactor.volume = 100"',
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr = 2
        if ierr > 0:
            return
        self._key = phrase.upper()  # Chemkin keyword phrase
        self._value = value  # value assigned to the keyword
        self._data_type = data_type  # data type of the values
        self._prefix = ""  # a prefix to the keyword that can be used
        # to comment out/disable the keyword by setting it to '!'
        self._set = True

    @staticmethod
    def setfullkeywords(mode: bool):
        """Require all keywords and their parameters must be specified."""
        """All keywords and their parameters must be specified by
        using the setkeyword method and will be passed to the reactor model
        for further processing.

        Parameters
        ----------
            mode: boolean
                turn the full keyword mode ON/OFF

        """
        if mode:
            # turn ON the full keyword mode (no checking on protected keywords)
            Keyword.no_fullkeyword = False
        else:
            Keyword.no_fullkeyword = True

    def show(self):
        """Display the Chemkin keyword and its parameter value."""
        if self._set:
            if isinstance(self._value, (int, float)):
                msg = [
                    Color.YELLOW,
                    "keyword",
                    "'" + self._key + "':",
                    "value =",
                    str(self._value),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            elif isinstance(self._value, bool):
                if self._value:
                    msg = [
                        Color.YELLOW,
                        "keyword",
                        "'" + self._key + "':",
                        "value = True",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
                else:
                    msg = [
                        Color.YELLOW,
                        "keyword",
                        "'" + self._prefix + self._key + "':",
                        "value = Disabled",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
            else:
                msg = [
                    Color.YELLOW,
                    "keyword",
                    "'" + self._key + "':",
                    "value =",
                    str(self._value),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
        else:
            msg = [
                Color.YELLOW,
                "keyword",
                "'" + self._key + "':",
                "value not set.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

    def resetvalue(self, value: Union[float, bool, str]):
        """Reset the parameter value of an existing keyword."""
        """
        Reset the parameter value of an existing keyword.

        Parameters
        ----------
            value: indicated by the data_type, {int, float, string, bool}
                keyword parameter

        """
        if isinstance(value, Keyword._valuetypes):
            if isinstance(value, bool):
                if value:
                    # true: keep the keyword active
                    self._prefix = ""
                    self._value = True
                else:
                    # false: disable the keyword
                    self._prefix = "!"
                    self._value = False
            else:
                # integer, float, or string parameter
                self._value = value
        else:
            msg = [
                Color.PURPLE,
                "value has a wrong data type",
                type(value),
                "value will not be reset.\n",
                Color.SPACEx6,
                "data type expected by keyword",
                self._key,
                "is",
                self._data_type,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def parametertype(self) -> type:
        """Get parameter type of the keyword."""
        """
        Get parameter type of the keyword.

        Returns
        -------
            parameter_data_type: string, {'int', 'float', 'string', 'bool'}
                parameter data type

        """
        return type(self._value)

    @property
    def value(self) -> Union[int, float, bool, str]:
        """Get parameter value of the keyword."""
        """
        Get parameter value of the keyword.

        Returns
        -------
            parameter value: {integer, floating, string, or boolean}

        """
        # extract the keyword value
        if self._data_type == "bool":
            mode = self._prefix != "!"
            return mode
        else:
            return self._value

    @property
    def keyphrase(self) -> str:
        """Get the phrase of the keyword."""
        """
        Get the phrase of the keyword.

        Returns
        -------
            keyword phrase: string

        """
        return self._key

    @property
    def keyprefix(self) -> bool:
        """Get the status of the keyword."""
        """
        Get the status of the keyword.

        Returns
        -------
            status: boolean
                keyword is ON/OFF

        """
        if self._prefix != "!":
            return True
        else:
            return False

    def getvalue_as_string(self) -> tuple[int, str]:
        """Create the keyword input line for Chemkin applications."""
        """
        Create the keyword input line for Chemkin applications.

        Returns
        -------
            linelength: integer
                line length
            line: string
                keyword value

        """
        # initialization
        line = ""
        linelength = 0
        # assembly the keyword line
        if self._data_type == "bool":
            # boolean keyword (active or disabled by '!')
            line = self._prefix + self._key
        else:
            # integer, double, or string parameter
            line = self._prefix + self._key + Keyword.fourspaces + str(self._value)

        linelength = len(line)
        return linelength, line


#
# This keyword type is used to distinguish keywords that act as
# on/off switches by their presence
#
class BooleanKeyword(Keyword):
    """Chemkin boolean keyword."""

    def __init__(self, phrase: str):
        """Set up a Chemkin keyword with a boolean parameter."""
        """
        Set up a Chemkin keyword with a boolean parameter or
        with no parameter.

        Parameters
        ----------
            phrase: string
                Chemkin keyword phrase

        """
        value = True
        super().__init__(phrase, value, "bool")


#
# This keyword type is used to hold integer keyword types
# (not sure if there actually are any of these)
#
class IntegerKeyword(Keyword):
    """A Chemkin integer keyword."""

    def __init__(self, phrase: str, value: int = 0):
        """Set up a Chemkin keyword with an integer parameter."""
        """
        Set up a Chemkin keyword with an integer parameter.

        Parameters
        ----------
            phrase: string
                Chemkin keyword phrase
            value: integer
                parameter value

        """
        super().__init__(phrase, value, "int")


#
# This keyword type is used to hold real keyword types
#
class RealKeyword(Keyword):
    """A Chemkin real keyword."""

    def __init__(self, phrase: str, value: float = 0.0e0):
        """Set up a Chemkin keyword with a real number parameter."""
        """Set up a Chemkin keyword with a real number
        (floating number) parameter.

        Parameters
        ----------
            phrase: string
                Chemkin keyword phrase
            value: double
                parameter value

        """
        super().__init__(phrase, value, "float")


#
# This keyword type is used to hold string keyword types
#
class StringKeyword(Keyword):
    """A Chemkin string keyword."""

    def __init__(self, phrase: str, value: str = ""):
        """Set up a Chemkin keyword with a string parameter."""
        """Set up a Chemkin keyword with a string parameter.

        Parameters
        ----------
            phrase: string
                Chemkin keyword phrase
            value: string
                parameter value

        """
        if len(value) <= 0:
            msg = [Color.PURPLE, "no string parameter given", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        super().__init__(phrase, value, "str")


class Profile:
    """Chemkin profile keyword class."""

    def __init__(self, key: str, x: npt.NDArray[np.double], y: npt.NDArray[np.double]):
        """Create a profile object."""
        """Create a profile object.

        Parameters
        ----------
            key: string
                profile keyword
            x: 1-D double array
                position of the profile data points
            y: 1-D double array
                variable value of the profile data

        """
        # initialization
        self._profilekeyword = ""
        self._status = 0
        # check
        if key.upper() in Keyword.profilekeywords:
            self._profilekeyword = key.upper()
        else:
            msg = [
                Color.PURPLE,
                "profile is not available under the reactor model",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            self._status = -1
            return
        # profile data sizes
        xsize = len(x)
        ysize = len(y)
        if xsize == ysize:
            self._size = xsize
            # independent variable (time, location, grid, ...)
            if isinstance(x, np.double):
                self._pos = copy.deepcopy(x)
            else:
                self._pos = np.array(x, dtype=np.double)
            # dependent variable value at the corresponding position
            if isinstance(y, np.double):
                self._val = copy.deepcopy(y)
            else:
                self._val = np.array(y, dtype=np.double)
        else:
            msg = [
                Color.PURPLE,
                "the number of positions does not match the number of values\n",
                Color.SPACEx6,
                "number of positions =",
                str(xsize),
                "\n",
                Color.SPACEx6,
                "number of values    =",
                str(ysize),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            self._status = -2

    @property
    def size(self) -> int:
        """Get number of data points in the profile."""
        """
        Returns
        -------
            size: integer
                number of profile data points

        """
        return self._size

    @property
    def status(self) -> int:
        """Get the validity of the profile object."""
        """
        Returns
        -------
            status: integer
                profile status

        """
        return self._status

    @property
    def pos(self) -> npt.NDArray[np.double]:
        """Get position values of profiles data."""
        """
        Returns
        -------
            pos: 1-D double array
                position [sec, cm]

        """
        return self._pos

    @property
    def value(self) -> npt.NDArray[np.double]:
        """Get variable values of profile data."""
        """
        Returns
        -------
            val: 1-D double array
                variable value

        """
        return self._val

    @property
    def profilekey(self) -> str:
        """Get profile keyword."""
        """
        Returns
        -------
            profilekeyword: string
                keyword associated with the profile data

        """
        return self._profilekeyword

    def show(self):
        """Show the profile data."""
        print(f"profile size: {self._size:d}")
        print(f" position           {self._profilekeyword:s}  ")
        for i in range(self._size):
            print(f"{self._pos[i]:f}         {self._val[i]}")

    def resetprofile(
        self, size: int, x: npt.NDArray[np.double], y: npt.NDArray[np.double]
    ):
        """Reset the profile data."""
        """
        Parameters
        ----------
            size: integer
                number of points of the new profile data
            x: 1-D double array
                position of the new profile data points
            y:1-D double array
                variable value of the new profile data

        """
        # check array size
        if size == self._size:
            # new profile has the same size
            self._pos[:] = x[:]
            self._val[:] = y[:]
        else:
            # new profile has different size
            self._size = size
            # resize the arrays
            self._pos.resize(size, refcheck=False)
            self._val.resize(size, refcheck=False)
        # fill the arrays with new values
        self._pos[:] = x[:]
        self._val[:] = y[:]

    def getprofile_as_string_list(self) -> tuple[int, list[str]]:
        """Create the keyword input lines as a list for Chemkin applications."""
        """Create the keyword input lines as a list for Chemkin applications.

        Returns
        -------
            size: integer
                number of profile lines
            line: list of strings
                list of profile related keywords

        """
        # initialization
        line = []
        # special treatment for pressure profile
        factor = 1.0e0
        if self._profilekeyword == "PPRO":
            if Keyword.no_fullkeyword:
                # use API calls: pressure profile units = dynes/cm2
                pass
            else:
                # use Full Keywords: pressure units = atm
                factor = P_ATM
        # assembly the profile keyword lines
        for i in range(self._size):
            thisline = ""
            thisline = (
                self._profilekeyword
                + Keyword.fourspaces
                + str(self._pos[i])
                + Keyword.fourspaces
                + str(self._val[i] / factor)
            )
            line.append(thisline)
        return self._size, line


#
# Framework and generic base classes for running Chemkin reactor models,
# defining methods to set chemistry, process keywords, and run
#
class ReactorModel:
    """Serve as a generic Chemkin reactor model framework."""

    def __init__(self, reactor_condition: Stream, label: str):
        """Initialize the basic parameters of Chemkin reactor model."""
        """Initialize the basic parameters of Chemkin reactor model.

        Parameters
        ----------
            reactor_condition: Chemistry or Mixture object
                mixture containing the initial/estimate reactor pressure, temperature,
                and gas composition
            label: string
                reactor label/name

        """
        # check mixture
        if isinstance(reactor_condition, (Mixture, Stream)):
            # if a Mixture/Stream object is passed in , verify the Mixture/Stream
            ierr = reactor_condition.validate()
            if ierr != 0:
                msg = [Color.PURPLE, "the mixture is not fully defined.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            # check if the chemstry set is active
            if not check_active_chemistryset(reactor_condition.chemid):
                msg = [
                    Color.PURPLE,
                    "the Chemistry Set associated with the"
                    "Mixture is not currently active.\n",
                    Color.SPACEx6,
                    "activate Chemistry Set using the 'active()' method.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            # chemistry set index
            self._chemset_index = ctypes.c_int(reactor_condition.chemid)
            # mixture
            self.reactormixture = copy.deepcopy(reactor_condition)
            # mixture gas species symbols
            self._specieslist = reactor_condition._specieslist  # gas species symbols
            # mixture temperature [K]
            self._temperature = ctypes.c_double(reactor_condition.temperature)
            # mixture pressure [dynes/cm2]
            self._pressure = ctypes.c_double(reactor_condition.pressure)
            self.numbspecies = self.reactormixture._kk
        # elif isinstance(reactor_condition, Chemistry):
        # if a Chemistry Set object is passed in (for flame models)
        # chemistry set index
        # self._chemset_index = ctypes.c_int(reactor_condition.chemid)
        # gas species symbols
        # self._specieslist = reactor_condition.species_symbols  # gas species symbols
        # self.numbspecies = reactor_condition.KK
        else:
            msg = [
                Color.PURPLE,
                "the first argument must be either",
                "a Chemistry object or a Mixture object.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # initialization
        self.label = label
        # number of required input
        self._numb_requiredinput = 0
        self._inputcheck: list[str] = []
        # gas reaction rate multiplier
        self._gasratemultiplier = 1.0e0
        # write text output file
        self._TextOut = True
        # FORTRAN file unit of the text output file
        self._mylout = c_int(154)
        # write XML solution file
        self._XMLOut = True
        # number of keywords used
        self._numbkeywords = 0
        # list of keyword phrases used for easy searching
        self._keyword_index: list[str] = []
        # list of keyword objects defined
        self._keyword_list: list[Keyword] = []
        # list of keyword lines
        # (each line is a string consists of:
        # '<keyword> <parameter>',
        # i.e., _keyword_index + _keyword_parameters)
        self._keyword_lines: list[str] = []
        # number of keyword lines
        self._numblines = 0
        # length of each keyword line
        self._linelength: list[int] = []
        # number of profile assigned
        self._numbprofiles = 0
        # list of profile keywords used for easy searching
        self._profiles_index: list[str] = []
        # list of profile objects defined
        self._profiles_list: list[Profile] = []
        # simulation run status
        #  -100 = not yet run
        #     0 = run success
        # other = run failed
        self.runstatus = -100
        # raw solution data structure
        self._solution_tags: list[str] = [
            "time",
            "distance",
            "temperature",
            "pressure",
            "volume",
            "velocity",
            "flowrate",
        ]
        self._speciesmode = "mass"
        self._numbsolutionpoints = 0
        self._solution_rawarray: dict[str, npt.ArrayLike] = {}
        self._numbsolutionmixtures = 0
        self._solution_mixturearray: list[Mixture] = []
        # initialize Chemkin-CFD-API
        if not check_chemistryset(self._chemset_index.value):
            # need to initialize Chemkin-CFD-API
            msg = [
                Color.YELLOW,
                "initializing Chemkin ...",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            ierr = chemkin_wrapper.chemkin.KINInitialize(self._chemset_index, c_int(0))
            if ierr == 0:
                chemistryset_initialized(self._chemset_index.value)
            else:
                msg = [
                    Color.RED,
                    "Chemkin-CFD-API initialization failed;",
                    "code =",
                    str(ierr),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()

    def usefullkeywords(self, mode: bool):
        """Specify all necessary keywords explicitly."""
        """Specify all necessary keywords explicitly.

        Parameters
        ----------
            mode: boolean, default = False
                turn full keyword mode ON/OFF

        """
        Keyword.setfullkeywords(mode)
        if mode:
            msg = [
                Color.YELLOW,
                "reactor",
                self.label,
                "will be run with full keyword input mode",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

    def __findkeywordslot(self, key: str) -> tuple[int, bool]:
        """Find the proper index in the global keyword list."""
        """Find the proper index in the global keyword list to add
        a new keyword or to modify the keyword parameter.

        Parameters
        ----------
            key: string
                Chemkin keyword

        Returns
        -------
            index: integer
                location of the keyword in the global keyword list
            status: boolean
                whether this is a new keyword

        """
        # check existing keyword
        if self._numbkeywords == 0:
            return 0, True
        else:
            if key in self._keyword_index:
                return self._keyword_index.index(key), False
            else:
                # new keyword
                return self._numbkeywords, True

    def setkeyword(self, key: str, value: Union[bool, float, str]):
        """Set a Chemkin keyword and its parameter."""
        """Set a Chemkin keyword and its parameter.

        Parameters
        ----------
            key: string
                Chemkin keyword phrase
            value: integer, double, string, or boolean depending on the keyword
                value associated with the keyword phrase

        """
        # find the keyword
        i, newkey = self.__findkeywordslot(key.upper())
        # add the keyword to the keywords list
        if newkey:
            # a new keyword
            if isinstance(value, str):
                # value is a string
                self._keyword_list.append(StringKeyword(key.upper(), value))
                self._keyword_index.append(key.upper())
            elif isinstance(value, bool):
                # value is a boolean value
                if value:
                    # set the keyword only if the value is True
                    self._keyword_list.append(BooleanKeyword(key.upper()))
                    self._keyword_index.append(key.upper())
                else:
                    # remove the count
                    self._numbkeywords -= 1
            elif isinstance(value, int):
                # value is an integer
                self._keyword_list.append(IntegerKeyword(key.upper(), value))
                self._keyword_index.append(key.upper())
            elif isinstance(value, float):
                # value is a real number
                self._keyword_list.append(RealKeyword(key.upper(), value))
                self._keyword_index.append(key.upper())
            else:
                msg = [Color.PURPLE, "invalid keyword value data type.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            self._numbkeywords += 1
        else:
            # an existing keyword, just update its value
            if isinstance(value, (str, bool)):
                self._keyword_list[i].resetvalue(value)
            elif isinstance(value, (int, float)):
                self._keyword_list[i].resetvalue(value)
            else:
                msg = [Color.PURPLE, "invalid keyword value data type.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    def removekeyword(self, key: str):
        """Remove an existing Chemkin keyword and its parameter."""
        """Remove an existing Chemkin keyword and its parameter.

        Parameters
        ----------
            key: string
                Chemkin keyword phrase

        """
        # find the keyword
        i, newkey = self.__findkeywordslot(key.upper())
        if newkey:
            msg = [Color.YELLOW, "keyword", key, "not found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        else:
            # remove keyword from the keyword index and the keyword list
            if self._keyword_list[i].keyphrase != key.upper():
                msg = [
                    Color.YELLOW,
                    "keyword index error.\n",
                    Color.SPACEx6,
                    "expected keyword",
                    key.upper(),
                    "   actual keyword",
                    self._keyword_list[i].keyphrase,
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.warning(this_msg)
                exit()
            # remove key from the keyword list and index
            del self._keyword_list[i]
            self._keyword_index.remove(key.upper())
            self._numbkeywords -= 1

    def showkeywordinputlines(self):
        """List all currently-defined keywords and their parameters line by line."""
        # header
        print("** INPUT KEYWORDS: \n")
        print("=" * 40)
        # display the keyword and the parameters line by line
        for k in self._keyword_list:
            n, line = k.getvalue_as_string()
            print(f"{line[:n]:s}")
        print("=" * 40)

    def createkeywordinputlines(self) -> tuple[int, int]:
        """Create keyword input lines for Chemkin applications."""
        """
        Remove an existing Chemkin keyword and its parameter.
        one keyword per line: <keyword>     <parameter>

        Returns
        -------
            Error code: integer
            number of lines: integer
                number of keyword lines to be added to the inputs

        """
        # initialization
        self._numblines = 0
        self._linelength.clear()
        self._keyword_lines.clear()
        # create the keyword lines from the keyword objects in the keyword list
        for k in self._keyword_list:
            n, line = k.getvalue_as_string()
            self._linelength.append(n)
            self._keyword_lines.append(line)
            self._numblines += 1
        # print the entire keyword input block
        if verbose() and self._numblines > 0:
            print("** INPUT KEYWORDS:")
            # print(f'number of keyword input lines:
            # {self._numblines:d} == {self._numbkeywords:d} \n')
            print("=" * 40)
            for line in self._keyword_lines:
                print(line)
            print("=" * 40)
        ierr = self._numbkeywords - self._numblines
        return ierr, self._numblines

    def showkeywordinputlines_with_tag(self, tag: str = ""):
        """List all currently-defined keywords."""
        """List all currently-defined keywords, their parameters, and an
        extra tag string line by line.

        Parameters
        ----------
            tag: string
                additional tag for the keywords, for example, the reactor index

        """
        # header
        print("** INPUT KEYWORDS: \n")
        print("=" * 40)
        # display the keyword and the parameters line by line
        for k in self._keyword_list:
            n, line = k.getvalue_as_string()
            print(f"{line[:n]:s}    {tag}")
        print("=" * 40)

    def createkeywordinputlines_with_tag(self, tag: str = "") -> tuple[int, int]:
        """Create keyword input lines for Chemkin applications."""
        """
        Create keyword input lines for Chemkin applications.
        one keyword per line: <keyword>     <parameter>    <tag>

        Parameters
        ----------
            tag: string
                additional tag for the keywords, for example, the reactor index

        Returns
        -------
            Error code: integer
            number of lines: integer
                number of keyword lines to be added to the inputs

        """
        # initialization
        self._numblines = 0
        self._linelength.clear()
        self._keyword_lines.clear()
        # create the keyword lines from the keyword objects in the keyword list
        for k in self._keyword_list:
            n, line = k.getvalue_as_string()
            # append the tag to the end of the line
            line = line + Keyword.fourspaces + tag
            # re-calculate the line length
            n = len(line)
            self._linelength.append(n)
            self._keyword_lines.append(line)
            self._numblines += 1
        # print the entire keyword input block
        if verbose() and self._numblines > 0:
            print("** INPUT KEYWORDS:")
            # print(f'number of keyword input lines:
            # {self._numblines:d} == {self._numbkeywords:d} \n')
            print("=" * 40)
            for line in self._keyword_lines:
                print(line)
            print("=" * 40)
        ierr = self._numbkeywords - self._numblines
        return ierr, self._numblines

    def __findprofileslot(self, key: str) -> tuple[int, bool]:
        """Find the proper index in the global profile list."""
        """Find the proper index in the global profile list either to add
        a new profile or to modify the existing profile parameter.

        Parameters
        ----------
            key: string
                Chemkin profile keyword

        Returns
        -------
            index: integer
                location of the keyword in the global keyword list
            status: boolean
                whether this is a new keyword

        """
        # check existing keyword
        if self._numbprofiles == 0:
            return 0, True
        else:
            if key in self._profiles_index:
                return self._profiles_index.index(key), False
            else:
                # new keyword
                return self._numbprofiles, True

    def setprofile(
        self, key: str, x: npt.NDArray[np.double], y: npt.NDArray[np.double]
    ) -> int:
        """Set a Chemkin profile and its parameter."""
        """Set a Chemkin profile and its parameter.

        Parameters
        ----------
            key: string
                Chemkin profile keyword phrase
            x: 1-D double array
                position values of the profile data
            y: 1-D double array
                variable values of the profile data

        Returns
        -------
            Error code: integer

        """
        #
        ierr = 0
        # find the keyword
        i, newprofile = self.__findprofileslot(key.upper())
        # add the profile to the profiles index list
        if newprofile:
            # a new profile
            self._profiles_list.append(Profile(key.upper(), x, y))
            status = self._profiles_list[i].status
            if status == 0:
                self._profiles_index.append(key.upper())
                self._numbprofiles += 1
            else:
                msg = [
                    Color.PURPLE,
                    "failed to create the profile",
                    '"' + key + '"\n',
                    Color.SPACEx6,
                    "error code =",
                    str(status),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr = status
        else:
            # an existing keyword, just update its value
            xsize = len(x)
            ysize = len(y)
            if xsize == ysize:
                self._profiles_list[i].resetprofile(xsize, x, y)
            else:
                msg = [
                    Color.PURPLE,
                    "the number of positions does not match the number of values\n",
                    Color.SPACEx6,
                    "number of position data = ",
                    str(xsize),
                    "\n",
                    Color.SPACEx6,
                    "number of value data   =",
                    str(ysize),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr = 1
        return ierr

    def createprofileinputlines(self) -> tuple[int, int, list[str]]:
        """Create profile keyword input lines for Chemkin applications."""
        """
        Create profile keyword input lines for Chemkin applications.
        one keyword per line: <profile keyword>     <position>  <value>

        Returns
        -------
            Error code: integer
            numblines: integer
                total number of profile keyword lines
            keyword_lines: list of string lists, [[string, ...], [string, ...], ...]
                string list containing lists of profile keywords

        """
        # initialization
        numblines: int = 0
        numbprofiles = 0
        keyword_lines: list[str] = []
        # create the keyword lines from the keyword objects in the profile list
        for p in self._profiles_list:
            n, lines = p.getprofile_as_string_list()
            keyword_lines.extend(lines)
            numblines += n
            numbprofiles += 1
            # print the entire keyword input block per profile
            if verbose():
                print("** PROFILE KEYWORDS:")
                print(f"{n:d} keyword input lines in {p._profilekeyword} profile\n")
                print("=" * 40)
                for line in lines:
                    print(line)
                print("=" * 40)
        # lines: list of strings of a profile ['VPRO x1 v1', 'VPRO x2 v2', ...]
        # keyword_lines: list of lines:  [['VPRO x1 v1', 'VPRO x2 v2', ...],
        # ['PPRO x1 p1', 'PPRO x2 p2', ..] , ... ]
        ierr: int = numbprofiles - self._numbprofiles
        return ierr, numblines, keyword_lines

    def createspeciesinputlines(
        self,
        solvertype: int,
        threshold: float = 1.0e-12,
        molefrac: npt.NDArray[np.double] = None,
    ) -> tuple[int, list[str]]:
        """Create keyword input lines."""
        """Create keyword input lines for initial/estimated
        species mole fraction inside the batch reactor.

        Parameters
        ----------
            solvertype: integer
                solver type of the reactor model
            threshold: double
                minimum species mole fraction value to be
                included in the species keyword
            molefrac: 1-D double array
                species composition in mole fractions

        Returns
        -------
            numb_lines: integer
                Number of keyword lines
            lines: list of strings
                list of keyword line strings

        """
        # initial(transient)/estimate(steady-state) composition keyword
        # depends on the solver type
        key = Keyword.gasspecieskeywords[solvertype - 1]
        ksym = self._specieslist
        lines = []
        numb_lines = 0
        for i in range(len(molefrac)):
            if molefrac[i] > threshold:
                thisline = (
                    key
                    + Keyword.fourspaces
                    + ksym[i].rstrip()
                    + Keyword.fourspaces
                    + str(molefrac[i])
                )
                lines.append(thisline)
                numb_lines += 1
        return numb_lines, lines

    def createspeciesinputlineswithaddon(
        self,
        key: str = "XEST",
        threshold: float = 1.0e-12,
        molefrac: npt.NDArray[np.double] = None,
        addon: str = "",
    ) -> tuple[int, list[str]]:
        """Create keyword input lines."""
        """Create keyword input lines for initial/estimated
        species mole fraction inside the batch reactor.

        Parameters
        ----------
            key: string
                Chemkin reactor keyword for species value
            threshold: douoble
                minimum species mole fraction value to be included in
                the species keyword
            molefrac: 1-D double array
                species composition in mole fractions
            addon: string
                add-on string to the species input, usually the
                reactor/zone number

        Returns
        -------
            numb_lines: integer
                Number of keyword lines
            lines: list of keyword line strings

        """
        # must use estimate composition keyword 'XEST'
        # (the 'REAC' keyword does not accept reactor/zone number)
        ksym = self._specieslist
        ksize = len(molefrac)
        if ksize != len(ksym):
            msg = [
                Color.PURPLE,
                "species mole fraction array has size",
                str(ksize),
                "but",
                str(len(ksym)),
                "is expected.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        lines = []
        numb_lines = 0
        for i in range(ksize):
            if molefrac[i] > threshold:
                thisline = (
                    key.rstrip()
                    + Keyword.fourspaces
                    + ksym[i].rstrip()
                    + Keyword.fourspaces
                    + str(molefrac[i])
                    + Keyword.fourspaces
                    + addon.rstrip()
                )
                lines.append(thisline)
                numb_lines += 1
        return numb_lines, lines

    def chemid(self) -> int:
        """Get chemistry set index."""
        """
        Get chemistry set index.

        Returns
        -------
            chemid: integer
                chemistry set index

        """
        return self._chemset_index.value

    @property
    def temperature(self) -> float:
        """Get reactor initial temperature."""
        """
        Get reactor initial temperature.

        Returns
        -------
            temperature: double
                reactor temperature [K]

        """
        return self.reactormixture.temperature

    @temperature.setter
    def temperature(self, t: float):
        """(Re)set reactor temperature."""
        """
        (Re)set reactor temperature.

        Parameters
        ----------
            t: double
                temperature [K]

        """
        if t <= 1.0e1:
            msg = [Color.PURPLE, "invalid temperature value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self._temperature = c_double(t)
        self.reactormixture.temperature = t

    @property
    def pressure(self) -> float:
        """Get reactor pressure."""
        """
        Returns
        -------
            pressure: double
                reactor pressure [dynes/cm2]

        """
        return self.reactormixture.pressure

    @pressure.setter
    def pressure(self, p: float):
        """(Re)set reactor pressure."""
        """
        Parameters
        ----------
            p: double
                pressure [dynes/cm2]

        """
        if p <= 0.0e0:
            msg = [Color.PURPLE, "invalid pressure value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self._pressure = c_double(p)
        self.reactormixture.pressure = p

    @property
    def massfraction(self) -> float:
        """Get the gas species mass fractions."""
        """Get the initial/guessed/estimate gas species mass fractions
        inside the reactor.

        Returns
        -------
            reactormixture: 1-D double array
                mixture mass fraction [-]

        """
        return self.reactormixture.y

    @massfraction.setter
    def massfraction(self, recipe: list[tuple[str, float]]):
        """Reset the initial/guessed/estimate gas species mass fractions."""
        """(Re)set the initial/guessed/estimate gas species mass fractions
        inside the reactor.

        Parameters
        ----------
            recipe: list of tuples, [(species_symbol, fraction), ... ]
                non-zero mixture composition corresponding to
                the given mole/mass fraction array

        """
        self.reactormixture.y(recipe)

    @property
    def molefraction(self) -> npt.NDArray[np.double]:
        """Get the gas species mole fractions."""
        """Get the initial/guessed/estimate gas species mole fractions
        inside the reactor.

        Returns
        -------
            X: 1-D double array
                mixture mole fraction

        """
        return self.reactormixture.x

    @molefraction.setter
    def molefraction(self, recipe: list[tuple[str, float]]):
        """Reset the initial/guessed/estimate gas species mole fractions."""
        """(Re)set the initial/guessed/estimate gas species mole fractions
        inside the reactor.

        Parameters
        ----------
            recipe: list of tuples, [(species_symbol, fraction), ... ]
                non-zero mixture composition corresponding to the given
                mole/mass fraction array

        """
        self.reactormixture.x(recipe)

    @property
    def concentration(self) -> npt.NDArray[np.double]:
        """Get the gas species molar concentrations."""
        """Get the initial/guessed/estimate gas species molar concentrations
        inside the reactor.

        Returns
        -------
            concentration: 1-D double array
                mixture molar concentration [mole/cm3]

        """
        return self.reactormixture.concentration

    def set_molefractions(self, molefrac: npt.NDArray[np.double]):
        """Reset the reactor species mole fractions."""
        """(Re)set the reactor initial/guessed species mole fractions.

        Parameters
        ----------
            molefrac: 1-D double array, dimension = number of gas species

        """
        self.reactormixture.x = molefrac

    def set_massfractions(self, massfrac: npt.NDArray[np.double]):
        """Reset the reactor species mass fractions."""
        """
        (Re)set the reactor initial/guessed species mass fractions.

        Parameters
        ----------
            molefrac: 1-D double array, dimension = number of gas species

        """
        self.reactormixture.y = massfrac

    def list_composition(self, mode: str, option: str = " ", bound: float = 0.0e0):
        """List the gas mixture composition inside the reactor."""
        """List the gas mixture composition inside the reactor.

        Parameters
        ----------
            mode: string, {'mass', 'mole'}, default = 'mole'
                flag specifies the fractions returned are 'mass' or 'mole' fractions
            option: string, {'all', ' '}, default = ' '
                flag specifies to list 'all' species or just the species with
                non-zero fraction
            bound: double, default = 0.0
                minimum fraction value for the species to be printed

        """
        self.reactormixture.list_composition(mode=mode, option=option, bound=bound)

    @property
    def gasratemultiplier(self) -> float:
        """Get the value of the gas-phase reaction rate multiplier."""
        """Get the value of the gas-phase reaction rate multiplier.

        Returns
        -------
            rate_factor: double
                gas-phase reaction rate multiplier

        """
        return self._gasratemultiplier

    @gasratemultiplier.setter
    def gasratemultiplier(self, value: float = 1.0e0):
        """Set the value of the gas-phase reaction rate multiplier."""
        """Set the value of the gas-phase reaction rate multiplier (optional).

        Parameters
        ----------
            value: double, default = 1.0
                gas-phase reaction rate multiplier

        """
        if value < 0.0:
            msg = [Color.PURPLE, "reaction rate multiplier must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._gasratemultiplier = value
            self.setkeyword(key="GFAC", value=value)

    @property
    def std_output(self) -> bool:
        """Get text output status."""
        """Get text output status.

        Returns
        -------
            status: boolean
                text output ON=True/OFF=False

        """
        return self._TextOut

    @std_output.setter
    def std_output(self, mode: bool):
        """Set text output status."""
        """"Set text output status (optional).

        Parameters
        ----------
            mode: boolean, default = True: always write to the text output file
                turn ON/turn OFF

        """
        off = not mode
        self.setkeyword(key="NO_SDOUTPUT_WRITE", value=off)
        self._TextOut = mode

    @property
    def xml_output(self) -> bool:
        """Get XML solution output status."""
        """Get XML solution output status.

        Returns
        -------
            status: boolean
                XML solution output ON=True/OFF=False

        """
        return self._XMLOut

    @xml_output.setter
    def xml_output(self, mode: bool):
        """Set XML solution output status."""
        """Set XML solution output status (optional).

        Parameters
        ----------
            mode: boolean, default = True: always create the XML solution file
                turn ON/turn OFF the XML solution output

        """
        off = not mode
        self.setkeyword(key="NO_XMLOUTPUT_WRITE", value=off)
        self._XMLOut = mode

    def setsensitivityanalysis(
        self,
        mode: bool = True,
        absolute_tolerance: Union[float, None] = None,
        relative_tolerance: Union[float, None] = None,
        temperature_threshold: Union[float, None] = None,
        species_threshold: Union[float, None] = None,
    ):
        """Switch ON/OFF A-factor sensitivity analysis."""
        """Switch ON/OFF A-factor sensitivity analysis.

        Parameters
        ----------
            mode: boolean
                turn A-factor sensitivity ON/OFF
            absolute_tolerance: double
                absolute tolerance of the sensitivity parameters
            relative_tolerance: double
                relative tolerance of the sensitivity parameters
            temperature_threshold: double
                threshold normalized temperature sensitivity parameter value
                to print out to the text output file
            species_threshold: double
                threshold normalized species sensitivity parameter value
                to print out to the text output file

        """
        if "ASEN" in self._keyword_index:
            # already defined
            i = self._keyword_index.index("ASEN")
            if mode:
                # reactivate the keyword if it is disabled
                if self._keyword_list[i]._prefix == "!":
                    self._keyword_list[i]._prefix = ""
                # set tolerances if given
                if absolute_tolerance is not None:
                    self.setkeyword(key="ATLS", value=absolute_tolerance)
                if relative_tolerance is not None:
                    self.setkeyword(key="RTLS", value=relative_tolerance)
                # reset the thresholds
                if temperature_threshold is not None:
                    self.setkeyword(key="EPST", value=temperature_threshold)
                if species_threshold is not None:
                    self.setkeyword(key="EPSS", value=species_threshold)
            else:
                # disable the keyword
                if self._keyword_list[i]._prefix != "!":
                    self._keyword_list[i]._prefix = "!"
        else:
            # not defined
            if mode:
                # enable sensitivity analysis
                self.setkeyword(key="ASEN", value=mode)
                # set sensitivity analysis related parameters
                if absolute_tolerance is not None:
                    self.setkeyword(key="ATLS", value=absolute_tolerance)
                if relative_tolerance is not None:
                    self.setkeyword(key="RTLS", value=relative_tolerance)
                if temperature_threshold is not None:
                    self.setkeyword(key="EPST", value=temperature_threshold)
                if species_threshold is not None:
                    self.setkeyword(key="EPSS", value=species_threshold)
            else:
                # do nothing
                pass

    def set_rop_analysis(self, mode=True, threshold=None):
        """Switch ON/OFF the ROP (Rate Of Production) analysis."""
        """Switch ON/OFF the ROP (Rate Of Production) analysis.

        Parameters
        ----------
            mode: boolean, default = False
                turn ROP ON/OFF
            threshold: double
                threshold ROP value to print out to the text output file

        """
        if "AROP" in self._keyword_index:
            # already defined
            i = self._keyword_index.index("AROP")
            if mode:
                # reactivate the keyword if it is disabled
                if self._keyword_list[i]._prefix == "!":
                    self._keyword_list[i]._prefix = ""
                # reset the threshold
                if threshold is not None:
                    self.setkeyword(key="EPSR", value=threshold)
            else:
                # disable the keyword
                if self._keyword_list[i]._prefix != "!":
                    self._keyword_list[i]._prefix = "!"
        else:
            # not defined
            if mode:
                # enable ROP analysis
                self.setkeyword(key="AROP", value=mode)
                if threshold is not None:
                    self.setkeyword(key="EPSR", value=threshold)
            else:
                # do nothing
                pass

    @property
    def realgas(self) -> bool:
        """Get the real gas EOS status."""
        """Get the real gas EOS status.

        Returns
        -------
            status: boolean
                status of the real-gas EOS model
                True: real gas EOS is turned ON

        """
        if "RLGAS" in self._keyword_index:
            # already defined
            i = self._keyword_index.index("RLGAS")
            if self._keyword_list[i]._prefix == "!":
                # commented out
                return False
            else:
                # is turned ON
                return True
        else:
            # has not been turned ON
            return False

    def userealgas_eos(self, mode: bool):
        """Set the option to turn ON/OFF the real gas model."""
        """Set the option to turn ON/OFF the real gas model
        for cubic EOS enabled gas-phase mechanism.

        Parameters
        ----------
            mode: boolean
                turn the Chemkin real-gas cubic EOS model ON/OFF

        """
        # turn ON/OFF the real gas EOS
        self.setkeyword(key="RLGAS", value=mode)
        # reset the real gas flag
        if not mode:
            # switch to the ideal gas law
            ierr = chemkin_wrapper.chemkin.KINRealGas_UseIdealGasLaw(
                self._chemset_index, c_int(0)
            )
            if ierr != 0:
                msg = [
                    Color.PURPLE,
                    "failed to turn OFF the real-gas EOS model,",
                    "error code =",
                    str(ierr),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    def setrealgasmixingmodel(self, model: int):
        """Set the real gas mixing rule/model."""
        """Set the real gas mixing rule/model
        for cubic EOS enabled gas-phase mechanism.

        Parameters
        ----------
            model: integer, {0, 1}
                Chemkin real-gas mixing rule method
                0 = Van der Waals
                1 = pseudocritical

        """
        # set the real gas mixing model
        _mixingmodels = ["Van der Waals", "pseudocritical"]
        if model in [0, 1]:
            msg = [
                Color.YELLOW,
                "the",
                _mixingmodels[model],
                "mixing rule is used.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.setkeyword(key="RLMIX", value=model)
        else:
            msg = [
                Color.PURPLE,
                "the real-gas mixing rule model index",
                str(model),
                "is invalid\n",
                Color.SPACEx6,
                "set model = 0 to use the",
                _mixingmodels[0],
                "mixing model\n",
                Color.SPACEx6,
                "set model = 1 to use the",
                _mixingmodels[1],
                "mixing model",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def setrunstatus(self, code: int):
        """Set the simulation run status."""
        """Set the simulation run status.

        Parameters
        ----------
            run_status: integer
                error code

        """
        self.runstatus = code

    def getrunstatus(self, mode: str = "silent") -> int:
        """Get the reactor model simulation status."""
        """Get the reactor model simulation status.

        Parameters
        ----------
            mode: string {'verbose', 'silent'}, default = 'silent'
                option for additional print information

        Returns
        -------
            run_status: integer
                error code: 0=success; -100=not run; other=failed

        """
        if mode.lower() == "verbose":
            if self.runstatus == -100:
                msg = [Color.YELLOW, "simulation yet to be run.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            elif self.runstatus == 0:
                msg = [Color.GREEN, "simulation run successfully.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            else:
                msg = [
                    Color.PURPLE,
                    "simulation failed with code",
                    str(self.runstatus),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)

        return self.runstatus

    def __process_keywords(self) -> int:
        """Serve as a dummy Chemkin reactor keyword processing method."""
        """
        Serve as a dummy Chemkin reactor keyword processing method
        to be overridden by child classes.

        Returns
        -------
            error code: integer

        """
        # a shell method to be overridden by child classes
        ierr = 0
        return ierr

    def __run_model(self) -> int:
        """Run reactor model simulation."""
        """
        Serve as a dummy simulation execution procedures
        to a specific Chemkin reactor model.
        It is intended to be overridden by child classes.

        Returns
        -------
            error code: integer

        """
        # a shell method to be overridden by child classes
        ierr = 0
        return ierr

    def run(self) -> int:
        """Serve as a generic Chemkin run reactor model method."""
        """
        Serve as a generic Chemkin run reactor model method
        to be overridden by child classes.

        Returns
        -------
            error code: integer

        """
        # a shell method to be overridden by child classes
        logger.debug("Running " + str(self.__class__.__name__) + " " + self.label)
        # keyword processing
        logger.debug("Processing keywords")
        ret_val = (
            self.__process_keywords()
        )  # each reactor model subclass to perform its own keyword processing
        logger.debug("Processing keywords complete")
        # run reactor model
        logger.debug("Running model")
        ret_val = self.__run_model()
        logger.debug("Running model complete, status = " + str(ret_val))

        return ret_val

    def setsolutionspeciesfracmode(self, mode: str = "mass"):
        """Set the type of species fractions in the solution."""
        """Set the type of species fractions in the solution.

        Parameters
        ----------
            mode: string {'mass', 'mole'}
                species fraction type to be returned by the post-processor

        """
        if mode.lower() in ["mole", "mass"]:
            self._speciesmode = mode.lower()
        else:
            # wrong mode value
            msg = [
                Color.PURPLE,
                "invalid species fraction mode,",
                'use mode = "mass" or mode = "mole"',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def getrawsolutionstatus(self) -> bool:
        """Get the status of the post-process."""
        """Get the status of the post-process.

        Returns
        -------
            status: boolean
                True = raw solution is ready,
                False = raw solution is yet to be processed

        """
        status = False
        if self._numbsolutionpoints > 0:
            status = True
        return status

    def getmixturesolutionstatus(self) -> bool:
        """Get the status of the post-process."""
        """Get the status of the post-process.

        Returns
        -------
            status: boolean
                True = solution mixtures is ready,
                False = solution mixtures are yet to be processed

        """
        status = False
        if len(self._solution_mixturearray) > 0:
            status = True
        return status

    def get_solution_size(self) -> tuple[int, int]:
        """Get the number of reactors and the number of solution points."""
        """Get the number of reactors and the number of solution points.

        Returns
        -------
            nreactor: integer
                number of reactors
            npoints: integer
                number of solution points

        """
        return 1, self._numbsolutionpoints

    def getnumbersolutionpoints(self) -> int:
        """Get  the number of solution points per reactor."""
        """Get  the number of solution points per reactor.

        Returns
        -------
            npoints: integer
                number of solution points

        """
        return self._numbsolutionpoints

    def parsespeciessolutiondata(self, frac: npt.NDArray[np.double]):
        """Parse the 2-D species fraction solution data."""
        """
        Parse the species fraction solution data that are stored
        in a 2D array (number_species x numb_solution).

        Parameters
        ----------
            frac: 2-D double array, dimension = [number_species, numb_solution]
                species fractions of each solution point

        """
        # create a temporary array to hold the solution data of one species
        y = np.zeros(self._numbsolutionpoints, dtype=np.double)
        #
        for k in range(self.numbspecies):
            y[:] = frac[k, :]
            # add to the raw solution data
            self._solution_rawarray[self._specieslist[k].rstrip()] = copy.deepcopy(y)
            y[:] = 0.0e0
        # clean up
        del y

    def process_solution(self):
        """Post-process solution."""
        """Post-process solution to extract the raw solution variable data
        to be overridden by child classes.
        """
        # a shell method to be overridden by child classes
        pass

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
from typing import Any, Union

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core import chemkin_wrapper
from ansys.chemkin.core.chemistry import (
    Chemistry,
    check_active_chemistryset,
    check_chemistryset,
    chemistryset_initialized,
    verbose,
)
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.inlet import Stream
from ansys.chemkin.core.keyword import (
    BooleanKeyword,
    IntegerKeyword,
    Keyword,
    RealKeyword,
    StringKeyword,
)
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.mixture import Mixture
from ansys.chemkin.core.profile import Profile
from ansys.chemkin.core.surface_chemistry_controller import SurfaceChemistryController
from ansys.chemkin.core.utilities import (
    critical_and_exit,
    error_and_exit,
    log_error_message,
    warning_and_exit,
)
from ansys.chemkin.core.validation import (
    validate_file_exists,
    validate_minimum_value,
    validate_species_array_size,
    validate_temperature,
)


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
                error_and_exit(msg)
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
                error_and_exit(msg)
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
            self.num_gas_reactions = self.reactormixture._ii_gas
            # RGP table
            self._use_rgp_table = False
        else:
            msg = [
                Color.PURPLE,
                "the first argument must be either",
                "a Chemistry object or a Mixture object.",
                Color.END,
            ]
            error_and_exit(msg)
        # initialization
        self.label = label
        # number of required input
        self._numb_requiredinput = 0
        self._inputcheck: list[str] = []
        # gas reaction rate multiplier
        self._gasratemultiplier = 1.0e0
        # write text output file
        self._textout = True
        # FORTRAN file unit of the text output file
        self._mylout = ctypes.c_int(154)
        # write XML solution file
        self._xmlout = True
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
            "thermicity",
        ]
        self._speciesmode = "mass"
        self._numbsolutionpoints = 0
        self._solution_rawarray: dict[str, npt.ArrayLike] = {}
        self._numbsolutionmixtures = 0
        self._solution_mixturearray: list[Mixture] = []
        # single point solution variables for transient models
        self._solution_parameters: dict[str, Union[int, float]] = {}
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
            ierr = chemkin_wrapper.chemkin.KINInitialize(
                self._chemset_index, ctypes.c_int(0)
            )
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
                critical_and_exit(msg)
        # surface chemistry
        self.has_surface_chemistry = reactor_condition._has_surface_chemistry
        self.numbmaterials = 0
        self.surface_chemistry: Any = None
        self._surface_controller = SurfaceChemistryController(self)
        if self.has_surface_chemistry:
            # activate surface chemistry
            self.activate_surface_chemistry(reactor_condition._chem_set, mode="silent")

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
        i, newkey = self.__findkeywordslot(key)
        # add the keyword to the keywords list
        if newkey:
            # a new keyword
            if isinstance(value, str):
                # value is a string
                self._keyword_list.append(StringKeyword(key, value))
                self._keyword_index.append(key)
            elif isinstance(value, bool):
                # value is a boolean value
                if value:
                    # set the keyword only if the value is True
                    self._keyword_list.append(BooleanKeyword(key))
                    self._keyword_index.append(key)
                else:
                    # remove the count
                    self._numbkeywords -= 1
            elif isinstance(value, int):
                # value is an integer
                self._keyword_list.append(IntegerKeyword(key, value))
                self._keyword_index.append(key)
            elif isinstance(value, float):
                # value is a real number
                self._keyword_list.append(RealKeyword(key, value))
                self._keyword_index.append(key)
            else:
                msg = [Color.PURPLE, "invalid keyword value data type.", Color.END]
                error_and_exit(msg)
            self._numbkeywords += 1
        else:
            # an existing keyword, just update its value
            if isinstance(value, (str, bool)):
                self._keyword_list[i].resetvalue(value)
            elif isinstance(value, (int, float)):
                self._keyword_list[i].resetvalue(value)
            else:
                msg = [Color.PURPLE, "invalid keyword value data type.", Color.END]
                error_and_exit(msg)

    def removekeyword(self, key: str):
        """Remove an existing Chemkin keyword and its parameter."""
        """Remove an existing Chemkin keyword and its parameter.

        Parameters
        ----------
            key: string
                Chemkin keyword phrase

        """
        # find the keyword
        i, newkey = self.__findkeywordslot(key)
        if newkey:
            msg = [Color.YELLOW, "keyword", key, "not found.", Color.END]
            warning_and_exit(msg)
        else:
            # remove keyword from the keyword index and the keyword list
            if self._keyword_list[i].keyphrase != key:
                msg = [
                    Color.YELLOW,
                    "keyword index error.\n",
                    Color.SPACEx6,
                    "expected keyword",
                    key,
                    "   actual keyword",
                    self._keyword_list[i].keyphrase,
                    Color.END,
                ]
                warning_and_exit(msg)
            # remove key from the keyword list and index
            del self._keyword_list[i]
            self._keyword_index.remove(key)
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
        self,
        key: str,
        x: npt.NDArray[np.double],
        y: npt.NDArray[np.double],
        label: bool = False,
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
            label: boolean, optional
                True: key includes a reactor/species/inlet name

        Returns
        -------
            Error code: integer

        """
        #
        ierr = 0
        # find the keyword
        i, newprofile = self.__findprofileslot(key)
        # add the profile to the profiles index list
        if newprofile:
            # a new profile
            self._profiles_list.append(Profile(key, x, y, label))
            status = self._profiles_list[i].status
            if status == 0:
                self._profiles_index.append(key)
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

    def get_profile_value(self, key: str, x: float) -> float:
        """Get the y value at given x from the profile data."""
        """
        Return the value of the y variable at the given x position
        from the piecewise linear profile data by interpolation.

        Parameters
        ----------
            key: string
                name of the profile data
            x: double
                the position value

        Returns
        -------
            value: double
                the interpolated y variable value at the given x position
        """
        # find the keyword
        i, newprofile = self.__findprofileslot(key)
        if newprofile or i < 0:
            # profile data does not exist
            msg = [
                Color.PURPLE,
                "cannot find the profile data",
                key,
                ", please verify the name of the profile data set.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 0.0
        # get the profile data
        this_x = self._profiles_list[i].pos
        this_val = self._profiles_list[i].value
        # check x range
        if x < np.min(this_x) or x > np.max(this_x):
            msg = [
                Color.PURPLE,
                "the given x value",
                str(x),
                "is out of the profile data range of [",
                str(np.min(this_x)),
                ",",
                str(np.max(this_x)),
                "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 0.0
        values = np.interp(x, this_x, this_val)
        return values[0]

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
            show = verbose()
            show = False
            if show:
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

    def _require_molefrac(
        self, molefrac: Union[None, npt.NDArray[np.double]]
    ) -> npt.NDArray[np.double]:
        """Return mole fractions when provided; terminate otherwise."""
        if molefrac is None:
            msg = [Color.PURPLE, "species composition is not provided.", Color.END]
            log_error_message(msg)
            exit()
        return molefrac

    def createspeciesinputlines(
        self,
        solvertype: int,
        threshold: float = 1.0e-12,
        molefrac: Union[None, npt.NDArray[np.double]] = None,
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
        molefrac = self._require_molefrac(molefrac)
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
        molefrac: Union[None, npt.NDArray[np.double]] = None,
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
        molefrac = self._require_molefrac(molefrac)
        ksym = self._specieslist
        ksize = len(molefrac)
        validate_species_array_size(
            expected=len(ksym),
            actual=ksize,
            context="mole fraction",
        )

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

    def create_inletspeciesinputlines(
        self,
        inlet_name: str,
        threshold: float = 1.0e-12,
        molefrac: Union[None, npt.NDArray[np.double]] = None,
    ) -> tuple[int, list[str]]:
        """Create keyword input lines."""
        """Create keyword input lines for initial/estimated
        species mole fraction inside the batch reactor.

        Parameters
        ----------
            inlet_name: string
                external inlet name
            threshold: douoble
                minimum species mole fraction value to be included in
                the species keyword
            molefrac: 1-D double array
                species composition in mole fractions

        Returns
        -------
            numb_lines: integer
                Number of keyword lines
            lines: list of keyword line strings

        """
        # must use inlet composition keyword 'REAC'
        molefrac = self._require_molefrac(molefrac)
        ksym = self._specieslist
        ksize = len(molefrac)
        validate_species_array_size(
            expected=len(ksym),
            actual=ksize,
            context="mole fraction",
        )

        lines = []
        numb_lines = 0
        if len(inlet_name.rstrip()) == 0:
            name_tag = ""
        else:
            name_tag = inlet_name.rstrip() + Keyword.fourspaces
        for i in range(ksize):
            if molefrac[i] > threshold:
                thisline = (
                    "REAC"
                    + Keyword.fourspaces
                    + name_tag
                    + ksym[i].rstrip()
                    + Keyword.fourspaces
                    + str(molefrac[i])
                )
                lines.append(thisline)
                numb_lines += 1
        return numb_lines, lines

    def clear_all_keywords(self):
        """Delete all existing keyword components."""
        # number of keywords used
        self._numbkeywords = 0
        # list of keyword phrases used for easy searching
        self._keyword_index.clear()
        # list of keyword objects defined
        self._keyword_list.clear()
        # list of keyword lines
        # (each line is a string consists of: '<keyword> <parameter>',
        # i.e., _keyword_index + _keyword_parameters)
        self._keyword_lines.clear()
        # number of keyword lines
        self._numblines = 0
        # length of each keyword line
        self._linelength.clear()

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
        """Set reactor temperature."""
        """
        (Re)set reactor temperature.

        Parameters
        ----------
            t: double
                temperature [K]

        """
        validate_temperature(t)
        self._temperature = ctypes.c_double(t)
        self.reactormixture.temperature = t

    @property
    def pressure(self) -> float:
        """Get reactor pressure."""
        """
        Get reactor pressure.

        Returns
        -------
            pressure: double
                reactor pressure [dynes/cm2]

        """
        return self.reactormixture.pressure

    @pressure.setter
    def pressure(self, p: float):
        """Set reactor pressure."""
        """
        (Re)set reactor pressure.

        Parameters
        ----------
            p: double
                pressure [dynes/cm2]

        """
        validate_minimum_value(
            value=p,
            minimum=0.0e0,
            message="invalid pressure value.",
        )
        self._pressure = ctypes.c_double(p)
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
        """Set the initial/guessed/estimate gas species mass fractions."""
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
        """Set the initial/guessed/estimate gas species mole fractions."""
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
        """Set the reactor species mole fractions."""
        """(Re)set the reactor initial/guessed species mole fractions.

        Parameters
        ----------
            molefrac: 1-D double array, dimension = number of gas species

        """
        self.reactormixture.x = molefrac

    def set_massfractions(self, massfrac: npt.NDArray[np.double]):
        """Set the reactor species mass fractions."""
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
        validate_minimum_value(
            value=value,
            minimum=0.0,
            message="reaction rate multiplier must >= 0.",
        )
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
        return self._textout

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
        self._textout = mode

    @property
    def xml_output(self) -> bool:
        """Get XML solution output status."""
        """Get XML solution output status.

        Returns
        -------
            status: boolean
                XML solution output ON=True/OFF=False

        """
        return self._xmlout

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
        self._xmlout = mode

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
                self._chemset_index, ctypes.c_int(0)
            )
            if ierr != 0:
                msg = [
                    Color.PURPLE,
                    "failed to turn OFF the real-gas EOS model,",
                    "error code =",
                    str(ierr),
                    Color.END,
                ]
                error_and_exit(msg)

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
            error_and_exit(msg)

    def use_rgp_table(self, data_file: str = ""):
        """Set the option to turn ON/OFF the use of the RGP table."""
        """Specify the RGP table file and turn ON the use of the RGP table.

        Parameters
        ----------
            data_file: string
                The RGP table file name with the full path
                if it is not in the current working directory.

        """
        # check file existence
        validate_file_exists(data_file, context_label="RGP table file")
        # check the mechanism to make sure that there is no gas-phase reaction
        # and there is only one gas species.
        if self.numbspecies > 1 or self.num_gas_reactions > 0:
            msg = [
                Color.PURPLE,
                "the RGP table can only be used for a mechanism with one gas species",
                "and no gas-phase reaction.",
                Color.END,
            ]
            error_and_exit(msg)
        # check the cubic EOS real-gas model status
        if self.realgas:
            msg = [
                Color.YELLOW,
                "the RGP table cannot be used with the real-gas EOS model.",
                "the real-gas EOS model is turned OFF.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas_eos(mode=False)
        # set the keyword
        self._use_rgp_table = True
        self.setkeyword(key="RGPTABLE", value=data_file)

    def disable_rgp_table(self):
        """Turn OFF the use of the RGP table."""
        if self._use_rgp_table:
            # reset RGP table flag
            self._use_rgp_table = False
            # remove the keyword
            self.removekeyword(key="RGPTABLE")

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
            error_and_exit(msg)

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
        """Post-process solution to extract the raw solution variable data."""
        # a shell method to be overridden by child classes
        pass

    def activate_surface_chemistry(self, chem: Chemistry, mode: str = "normal"):
        """Activate the surface chemistry."""
        self._surface_controller.activate_surface_chemistry(chem=chem, mode=mode)

    def verify_surface_chemistry(self) -> bool:
        """Verify whether the reactor model is surfac chemistry capable."""
        return self._surface_controller.verify_surface_chemistry()

    def no_surface_mechanism_declaration(self):
        """Inform users surface chemistry is not available for this reactor."""
        self._surface_controller.no_surface_mechanism_declaration()

    def no_surface_material(self, mat_name: str):
        """Send surface material not found message."""
        self._surface_controller.no_surface_material(mat_name=mat_name)

    @property
    def get_numb_material(self) -> int:
        """Get the number of surface material."""
        return self._surface_controller.get_numb_material

    def get_material_names(self) -> list[str]:
        """Get all surface material names."""
        return self._surface_controller.get_material_names()

    def get_site_species_names(self) -> list[list[str]]:
        """Get site species names on the surface materials."""
        return self._surface_controller.get_site_species_names()

    def get_bulk_species_names(self) -> list[list[str]]:
        """Get bulk species names on the surface materials."""
        return self._surface_controller.get_bulk_species_names()

    def get_total_site_species(self) -> int:
        """Get total number of site species from all materials."""
        return self._surface_controller.get_total_site_species()

    def get_total_bulk_species(self) -> int:
        """Get total number of bulk species from all materials."""
        return self._surface_controller.get_total_bulk_species()

    @property
    def surface_ratemultiplier(self) -> float:
        """Get the value of the surface reaction rate multiplier."""
        return self._surface_controller.surface_ratemultiplier

    @surface_ratemultiplier.setter
    def surface_ratemultiplier(self, value: float = 1.0e0):
        """Set the value of the surface reaction rate multiplier."""
        self._surface_controller.surface_ratemultiplier = value

    def check_material_area_fraction(self, mat_name: str) -> int:
        """Verify the surface area fraction is given."""
        return self._surface_controller.check_material_area_fraction(mat_name=mat_name)

    def get_material_area_fraction(self, mat_name: str) -> float:
        """Get the surface area fraction."""
        return self._surface_controller.get_material_area_fraction(mat_name=mat_name)

    def set_material_area_fraction(self, mat_name: str, fraction: float):
        """Set the value of the surface area fraction."""
        self._surface_controller.set_material_area_fraction(
            mat_name=mat_name, fraction=fraction
        )

    def check_material_temperature(self, mat_name: str) -> int:
        """Verify the surface temperature is given."""
        return self._surface_controller.check_material_temperature(mat_name=mat_name)

    def get_material_temperature(self, mat_name: str) -> float:
        """Get the surface temperature."""
        return self._surface_controller.get_material_temperature(mat_name=mat_name)

    def set_material_temperature(self, mat_name: str, temp: float):
        """Set the value of the surface temperature."""
        self._surface_controller.set_material_temperature(mat_name=mat_name, temp=temp)

    def get_all_surface_temperature(self) -> npt.NDArray[np.double]:
        """Get the temperatures of all materials in the mechanism."""
        return self._surface_controller.get_all_surface_temperature()

    def get_site_fraction(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the surface site fractions."""
        return self._surface_controller.get_site_fraction(mat_name=mat_name)

    def set_site_fraction(self, mat_name: str, recipe: list[tuple[str, float]]):
        """Set the surface site fractions."""
        self._surface_controller.set_site_fraction(mat_name=mat_name, recipe=recipe)

    def get_all_site_fractions(self) -> npt.NDArray[np.double]:
        """Get site fractions of all materials in the mechanism."""
        return self._surface_controller.get_all_site_fractions()

    def get_bulk_activity(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the bulk species activities."""
        return self._surface_controller.get_bulk_activity(mat_name=mat_name)

    def set_bulk_activity(self, mat_name: str, recipe: list[tuple[str, float]]):
        """Set the bulk species activities."""
        self._surface_controller.set_bulk_activity(mat_name=mat_name, recipe=recipe)

    def get_all_bulk_growth_rates(self) -> npt.NDArray[np.double]:
        """Get bulk growth rates of all materials in the mechanism."""
        return self._surface_controller.get_all_bulk_growth_rates()

    def set_init_surface_coverage(
        self,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Set up the initial surface coverage arrays."""
        return self._surface_controller.set_init_surface_coverage()

    def set_surface_chemistry_keywords(self):
        """Set up surface chemistry related keywords."""
        self._surface_controller.set_surface_chemistry_keywords()

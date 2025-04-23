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
Steady state, 1-D burner stabilized premixed flame models.
"""

import copy
from ctypes import c_double, c_int
from typing import Union

from ansys.chemkin import chemkin_wrapper
from ansys.chemkin.chemistry import (
    check_chemistryset,
    chemistryset_initialized,
    force_activate_chemistryset,
)
from ansys.chemkin.color import Color as Color
from ansys.chemkin.flame import Flame
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
from ansys.chemkin.mixture import interpolate_mixtures
from ansys.chemkin.reactormodel import Keyword
from ansys.chemkin.utilities import find_interpolate_parameters
import numpy as np
import numpy.typing as npt


class PremixedFlame(Flame):
    def __init__(self, inlet: Stream, label: Union[str, None] = None):
        # check reactor Mixture object
        if not isinstance(inlet, Stream):
            # wrong argument type
            msg = [Color.RED, "the first argument must be a Mixture object.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # set label
        if label is None:
            self.label = "premixedflame"
        else:
            self.label = label
        # set flow area to unity for easy conversion from mass flow rate to mass flux in the flame models
        if not inlet._haveflowarea:
            inlet.flowarea = 1.0  # [cm2]
        # initialization
        super().__init__(fuelstream=inlet, label=label)
        # mass flow rate [g/sec]
        # constant
        self._final_mass_flow_rate = -1.0

    def set_inlet(self, extinlet: Stream):
        """
        Add an external inlet to the reactor

        Parameters
        ----------
            extinlet: Stream object
                external inlet to the open reactor
        """
        # There is only ONE inlet allowed for the premixed flame models.
        msg = [
            Color.MAGENTA,
            "Premixed flame models do NOT allow the second inlet stream.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()

    def unburnt_temperature(self, temperature: float):
        """
        Set the unburnt fuel-oxidizer gas temperature for the flame speed calculation.
        By default, the inlet Stream temperature will be used.

        Parameters
        ----------
            temperature: double
                unburnt gas temperature [K]
        """
        if temperature <= 200.0:
            msg = [Color.PURPLE, "invalid temperature value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # set unburnt temperature
        self.temperature = temperature
        self.setkeyword("TUNB", value=temperature)

    def lump_diffusion_imbalance(self, mode: bool = True):
        """
        Lamp the "mass flux imbalance" due to species transport to the last species.
        The net diffusion flux at any interface should be zero. Use the lumping option to
        assign all mass imbalance to the last gas species of the mechanism by forcing
        its mass fraction to be 1 - (sum of all other species mass fractions). By default,
        the correction velocity formulism is used to distribute the mass flux imbalance evenly
        to all species.

        Parameters
        ----------
            mode: boolean {True, False}
                ON/OFF
        """
        # activate the lumping option to conserve mass
        self.setkeyword("TRCE", value=mode)

    def set_profilekeywords(self) -> int:
        """
        Create profile keywords for Chemkin flame applications

        one keyword per line: <profile keyword>     <position>  <value>

        Returns
        -------
            Error code: integer
        """
        # initialization
        tag = "TPRO"
        numblines = 0
        # create the keyword lines from the keyword objects in the profile list
        if tag in self._profiles_index:
            profile_ID = self._profiles_index.index(tag)
            T_profile = self._profiles_list[profile_ID]
            npoints = T_profile.size
            #
            positions = T_profile.pos
            y = T_profile.value
            # loop over all data points
            for x in positions:
                this_key = ""
                this_key = tag + " " + str(x) + Keyword.fourspaces + str(y[numblines])
                self.setkeyword(this_key, True)
                numblines += 1
            # check error
            iErr = numblines - npoints
            return iErr
        else:
            # no temperature profile found
            msg = [Color.PURPLE, "no temperature profile found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return -1

    def use_TPRO_grids(self, mode: bool = True):
        """
        Use the position values of the temperature profile data as
        the initial grid points to start the simulation

        Parameters
        ----------
            mode: boolean {True, False}
                ON/OFF
        """
        # use the TPRO grids
        self.setkeyword("USE_TPRO_GRID", value=mode)

    def set_gridkeywords(self) -> int:
        """
        Create 1-D grid profile keywords for Chemkin flame applications

        one keyword per line: <profile keyword>     <position>

        Returns
        -------
            Error code: integer
        """
        # initialization
        tag = "GRID"
        numblines = 0
        # create the keyword lines from the keyword objects in the profile list
        npoints = self.numb_grid_profile
        # loop over all data points
        for x in self.grid_profile:
            self.setkeyword(tag, x)
            numblines += 1
        # check error
        iErr = numblines - npoints
        return iErr

    def __run_model(self) -> int:
        """
        Run the reactor model after the keywords are processed

        Returns
        -------
            Error code: integer
        """
        # estimated reactor mass fraction
        Y_init = self.reactormixture.Y
        # run the premixed flame simulation
        iErr = chemkin_wrapper.chemkin.KINPremix_CalculateFlame(
            self._myLOUT,
            self._chemset_index,
            self._pressure,
            self._temperature,
            Y_init,
            c_double(self.starting_x),
            c_double(self.ending_x),
        )

        return iErr

    def __process_keywords(self):
        """
        Process input keywords for the reactor model

        Returns
        -------
            Error code: integer
        """
        iErr = 0
        iErrc = 0
        # set_verbose(True)
        # set inlet mass flux (estimated)
        # this keyword is optional
        if self.mass_flow_rate > 0.0:
            self.setkeyword("FLRT", self.mass_flow_rate)
        # set inlet/unburnt gas temperature
        if self._flamemode == 0:
            # freely propagating flame
            if "TUNB" not in self._keyword_index:
                self.setkeyword("TUNB", self.temperature)
        elif self._flamemode == 1:
            # burner stabilized flame with given temperature profile
            # check temperature profile
            if not self.temp_profile_set:
                # no temperature profile data
                iErr = 1
                msg = [
                    Color.PURPLE,
                    "no temperature profile data given for",
                    "the premixed burner stabilized flame model",
                    "with given temperature profile.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return iErr
        else:
            # burner stabilized flame and solve the energy equation
            if "TUNB" not in self._keyword_index:
                self.setkeyword("TUNB", self.temperature)
        # prepare mesh keywords
        self.set_mesh_keywords()
        if iErr == 0:
            # set additional keywords
            self.set_SSsolver_keywords()
            # set profile keywords
            iErrkey = 0
            if self._numbprofiles > 0:
                iErrc = self.set_profilekeywords()
                iErrkey += iErrc
            # prepare mesh keywords
            iErrc = self.set_mesh_keywords()
            iErrkey += iErrc
        # set keywords
        if iErr + iErrkey == 0:
            # pass all the keywords to the flame model
            for k in self._keyword_list:
                this_key = bytes(k.keyphrase, "utf-8")  # Chemkin keyword phrase
                this_value = c_double(k.value)  # value assigned to the keyword
                this_type = k.parametertype()  # data type of the values
                #
                if k.keyprefix:
                    # active keyword:
                    if this_type is bool:
                        # boolean type value: just assign the keyword value to 0.0
                        this_value = c_double(0.0)
                    elif this_type is str:
                        # string type value: just assign the keyword value to 0.0
                        this_value = c_double(0.0)
                    # set the keyword
                    iErrc = chemkin_wrapper.chemkin.KINPremix_SetParameter(
                        this_key, this_value
                    )
                    if iErrc == 2:
                        # keyword is not available
                        msg = [
                            Color.PURPLE,
                            "keyword,",
                            k.keyphrase,
                            "is not available through PyChemkin.",
                            Color.END,
                        ]
                        this_msg = Color.SPACE.join(msg)
                        logger.error(this_msg)
                        iErr += iErrc
                    elif iErrc != 0:
                        msg = [
                            Color.PURPLE,
                            "failed to process keyword,",
                            k.keyphrase,
                            "error code =",
                            str(iErrc),
                            Color.END,
                        ]
                        this_msg = Color.SPACE.join(msg)
                        logger.error(this_msg)
                        iErr += iErrc
        #
        self.showkeywordinputlines()
        #
        return iErr + iErrkey

    def run(self) -> int:
        """
        Chemkin run premixed flame model method

        Returns
        -------
            Error code: integer
        """
        #
        # activate the Chemistry set associated with the Reactor instance
        force_activate_chemistryset(self._chemset_index.value)
        #
        # get ready to run the reactor model
        # initialize Chemkin-CFD-API
        msg = [
            Color.YELLOW,
            "running model",
            self.__class__.__name__,
            self.label,
            "...\n",
            Color.SPACEx6,
            "initialization =",
            str(check_chemistryset(self._chemset_index.value)),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        if not check_chemistryset(self._chemset_index.value):
            # Chemkin-CFD-API is not initialized: reinitialize Chemkin-CFD-API
            msg = [Color.YELLOW, "initializing Chemkin ...", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            retVal = chemkin_wrapper.chemkin.KINInitialize(
                self._chemset_index, self._solvertype
            )
            if retVal != 0:
                msg = [
                    Color.RED,
                    "Chemkin-CFD-API initialization failed;",
                    "code =",
                    str(retVal),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()
            else:
                chemistryset_initialized(self._chemset_index.value)

        # output initialization
        logger.debug("clearing output ...")

        # keyword processing
        msg = [
            Color.YELLOW,
            "processing and generating keyword inputs ...",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        retVal = (
            self.__process_keywords()
        )  # each reactor model subclass to perform its own keyword processing
        if retVal != 0:
            msg = [
                Color.RED,
                "generating the keyword inputs,",
                "error code =",
                str(retVal),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            return retVal
        logger.debug("processing keywords complete")

        # run reactor model
        msg = [Color.YELLOW, "running premixed flame simulation ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        # use API calls
        retVal = self.__run_model()
        # update run status
        self.setrunstatus(code=retVal)
        msg = ["simulation completed,", "status =", str(retVal), Color.END]
        if retVal == 0:
            msg.insert(0, Color.GREEN)
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        else:
            msg.insert(0, Color.RED)
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)

        return retVal

    def continuation(self) -> int:
        """
        Perform a continuation run after the original flame simulation is
        completed successfully.

        Returns
        -------
            Error code: integer
        """
        # check if the model is already run once
        status = self.getrunstatus(mode="silent")
        if status == -100:
            msg = [Color.MAGENTA, "please run the reactor simulation first.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        elif status != 0:
            msg = [
                Color.PURPLE,
                "simulation was failed.\n",
                Color.SPACEx6,
                "please correct the error(s) and rerun the reactor simulation.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # insert the continuation keyword
        key_continue = bytes("CNTN", "utf-8")
        this_value = c_double(0.0)
        iErr = chemkin_wrapper.chemkin.KINPremix_SetParameter(key_continue, this_value)
        status += iErr
        if status == 0:
            msg = [
                Color.YELLOW,
                "continuation run starting...",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            # run the model
            iErr = self.run()
            status += iErr
        #
        return status

    def get_solution_size(self) -> int:
        """
        Get the number of solution points

        Returns
        -------
            npoints: integer
                number of solution points
        """
        # check run completion
        status = self.getrunstatus(mode="silent")
        if status == -100:
            msg = [Color.MAGENTA, "please run the reactor simulation first.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        elif status != 0:
            msg = [
                Color.PURPLE,
                "simulation was failed.\n",
                Color.SPACEx6,
                "please correct the error(s) and rerun the reactor simulation.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # number of time points in the solution
        npoints = c_int(0)
        # get solution size of the batch reactor
        iErr = chemkin_wrapper.chemkin.KINPremix_GetSolutionGridPoints(npoints)
        if iErr == 0 and npoints.value > 2:
            # return the solution sizes
            self._numbsolutionpoints = (
                npoints.value
            )  # number of time points in the solution profile
            return self._numbsolutionpoints
        else:
            # fail to get solution sizes
            msg = [
                Color.PURPLE,
                "failed to get the solution size,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def process_solution(self):
        """
        Post-process solution to extract the raw solution variable data
        """
        # check existing raw data
        if self.getrawsolutionstatus():
            msg = [
                Color.YELLOW,
                "the solution has been processed before,",
                "any existing solution data will be deleted from the memory.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

        # reset raw and mixture solution parameters
        self._numbsolutionpoints = 0
        self._solution_rawarray.clear()
        self._solution_mixturearray.clear()
        # get solution sizes
        npoints = self.get_solution_size()
        # check values
        if npoints <= 2:
            msg = [
                Color.PURPLE,
                "solution size error(s).\n",
                Color.SPACEx6,
                "number of solution points =",
                str(npoints),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._numbsolutionpoints = npoints
        # create arrays to hold the raw solution data
        pos = np.zeros(self._numbsolutionpoints, dtype=np.double)
        temp = np.zeros_like(pos, dtype=np.double)
        # create a species mass fraction array to hold the solution species fraction profiles
        frac = np.zeros(
            (
                self.numbspecies,
                self._numbsolutionpoints,
            ),
            dtype=np.double,
            order="F",
        )
        msg = [Color.YELLOW, "post-processing raw solution data ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        # create a species mass fraction array to hold the steady-state solution
        frac = np.zeros(
            (
                self.numbspecies,
                self._numbsolutionpoints,
            ),
            dtype=np.double,
            order="F",
        )
        # get raw solution data
        npoint = c_int(npoints)
        nspecies = c_int(self.reactormixture.KK)
        iErr = chemkin_wrapper.chemkin.KINPremix_GetSolution(
            npoint, nspecies, pos, temp, frac
        )
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to fetch the raw solution data from memory,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # get the flame mass flux [g/sec-cm2]
        massflux = c_double(0.0)
        iErr = chemkin_wrapper.chemkin.KINPremix_GetFlameMassFlux(massflux)
        if iErr == 0:
            self._final_mass_flow_rate = max(0.0, massflux.value)
        else:
            msg = [
                Color.RED,
                "failed to get the flame mass flux,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # store the raw solution data in a dictionary
        # grid
        self._solution_rawarray["distance"] = copy.deepcopy(pos)
        # temperature
        self._solution_rawarray["temperature"] = copy.deepcopy(temp)
        # species mass fractions
        self.parsespeciessolutiondata(frac)
        # create solution mixture
        iErr = self.create_solution_streams(frac)
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "forming solution streams",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # compute laniar flame speed
        if self._flamemode == 0:
            inlet_density = self._solution_mixturearray[0].RHO
            self.flamespeed = self._final_mass_flow_rate / inlet_density
        # clean up
        del pos, temp, frac

    def get_solution_variable_profile(self, varname: str) -> npt.NDArray[np.double]:
        """
        Get the profile of the solution variable specified

        Parameters
        ----------
            varname: string
                name of the solution variable

        Returns
        -------
            solution value profile: 1D double array
        """
        if not self.getrawsolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'getsolution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return 1
        # check variable name
        vname = varname.rstrip()
        if vname.lower() in ["distance", "temperature"]:
            # is a property variable?
            vname = vname.lower()
        else:
            if vname not in self._specieslist:
                # is not a species?
                msg = [
                    Color.PURPLE,
                    "variable name",
                    vname,
                    "is not part of the known solution variable.\n",
                    Color.SPACEx6,
                    "and has to be derived from other variable(s).",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

        # create variable arrays to hold the solution profile
        var = np.zeros(self._numbsolutionpoints, dtype=np.double)
        # get variable profile from the raw solution data
        var = self._solution_rawarray.get(vname)
        return var

    def create_solution_streams(self, specfrac: npt.NDArray[np.double]) -> int:
        """
        Create a list of Streams that represent the gas mixture at a solution point

        Parameters
        ----------
            specfrac: 2D double array, dimensions = [number_species, numb_solution_point]
                species fractions of all time points [species_fraction, time_point]

        Returns
        -------
            iError: integer
                 error code
        """
        if not self.getrawsolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'getsolution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return 1
        # create a temporary Stream object to hold the mixture properties at current solution point
        sstream = copy.deepcopy(self.reactormixture)
        # create variable arrays to hold the solution profile
        species = []
        # create a species fraction array to hold the solution species fraction profiles
        frac = np.zeros(self.numbspecies, dtype=np.double)
        # get solution variable profile from the raw solution arrays
        temp = self.get_solution_variable_profile("temperature")
        # loop over all species
        for sp in self._specieslist:
            species.append(self.get_solution_variable_profile(sp))
        # loop over all solution points
        for i in range(self._numbsolutionpoints):
            # get stream properties at the current solution point
            # pressure [dynes/cm2]
            sstream.pressure = self.pressure
            # temperature [K]
            sstream.temperature = temp[i]
            # stream mass flow rate [g/sec]
            sstream.mass_flowrate = self._final_mass_flow_rate
            # species composition
            for k in range(self.numbspecies):
                frac[k] = specfrac[k, i]
            # set mixture composition
            if self._speciesmode == "mass":
                # mass fractions
                sstream.Y = frac
            else:
                # mole fractions
                sstream.X = frac
            # add to the solution stream list
            self._solution_mixturearray.append(copy.deepcopy(sstream))
        # clean up
        species.clear()
        del temp, frac, species, sstream
        return 0

    def get_solution_stream(self, x: float) -> Stream:
        """
        Get the Stream representing the solution state at the given location

        Parameters
        ----------
            x: double
                grid point value [cm]

        Returns
        -------
            mixturetarget: Stream object
                a Stream representing the gas properties in the flame domain at the specific location
        """
        # check status
        if not self.getmixturesolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'process_solution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            exit()
        # get the time point array
        posarray = self.get_solution_variable_profile("distance")
        # find the interpolation parameters
        ileft, ratio = find_interpolate_parameters(x, posarray)
        # find the mixture
        if ratio == 0.0e0:
            # get the mixtures
            mixtureleft = copy.deepcopy(self._solution_mixturearray[ileft])
            return mixtureleft
        elif ratio == 1.0e0:
            # get the mixtures
            mixtureright = copy.deepcopy(self._solution_mixturearray[ileft + 1])
            return mixtureright
        else:
            # get the mixtures
            mixtureleft = copy.deepcopy(self._solution_mixturearray[ileft])
            mixtureright = copy.deepcopy(self._solution_mixturearray[ileft + 1])
            # interpolate the mixture properties
            mixturetarget = interpolate_mixtures(mixtureleft, mixtureright, ratio)
            # set mass flow rate
            mixturetarget.mass_flowrate = mixtureleft.mass_flowrate
            # clean up
            del mixtureleft, mixtureright
            #
            return mixturetarget

    def get_solution_stream_at_grid(self, grid_index: int) -> Stream:
        """
        Get the Stream representing the solution state at the given solution point index

        Parameters
        ----------
            grid_index: integer
                0-base grid point index

        Returns
        -------
            mixturetarget: Stream object
                a Stream representing the gas properties at the specific time
        """
        # check status
        if not self.getmixturesolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'process_solution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            exit()
        # check index
        if grid_index > self._numbsolutionpoints - 1:
            msg = [
                Color.PURPLE,
                "the given time point index:",
                str(grid_index),
                "> the maximum number of grid points:",
                str(self._numbsolutionpoints - 1),
                "\n",
                Color.SPACEx6,
                "the solution grid point index is 0-based.\n",
                Color.SPACEx6,
                "[ 0 ->",
                str(self._numbsolutionpoints - 1),
                "]",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get the mixture
        mixturetarget = copy.deepcopy(self._solution_mixturearray[grid_index])
        return mixturetarget


class BurnedStabilized_GivenTemperature(PremixedFlame):
    def __init__(self, inlet: Stream, label: Union[str, None] = None):
        #
        if label is None:
            label = "Premixed Burner"
        # initialization
        super().__init__(inlet=inlet, label=label)
        # energy equation
        self._energytype = self.EnergyTypes.get("TGiven", c_int(2))
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(161)
        # flame model mode: burner stabilized
        self._flamemode = 1
        # set required keywords for this model
        self.setkeyword("BURN", True)
        # use the user profile
        self.setkeyword("TGIV", True)


class BurnedStabilized_EnergyEquation(PremixedFlame):
    def __init__(self, inlet: Stream, label: Union[str, None] = None):
        # initialization
        if label is None:
            label = "Premixed Burner"
        super().__init__(inlet=inlet, label=label)
        # energy equation
        self._energytype = self.EnergyTypes.get("ENERGY", c_int(1))
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(162)
        # flame model mode: burner stabilized
        self._flamemode = 2
        # set required keywords for this model
        self.setkeyword("BURN", True)
        # solve the energy equation
        self.setkeyword("ENRG", True)

    def skip_fix_T_solution(self, mode: bool = True):
        """
        Skip the step of finding the intermediate solution with fixed temperature

        Parameters
        ----------
            mode: boolean {True, False}
                ON/OFF
        """
        # skip the fixed temperature solution
        self.setkeyword("NOFT", value=mode)

    def automatic_temperature_profile_estimate(self, mode: bool = True):
        """
        Let the premixed flame model to construct an estimated temperature profile
        based on the equilibrium state to start the calculation

        Parameters
        ----------
            mode: boolean {True, False}
                ON/OFF
        """
        # use the automatic temperature profile estimate function
        self.setkeyword("TPROF", value=mode)


class FreelyPropagating(PremixedFlame):
    def __init__(self, inlet: Stream, label: Union[str, None] = None):
        # initialization
        if label is None:
            label = "Premixed Propagating"
        super().__init__(inlet=inlet, label=label)
        # energy equation
        self._energytype = self.EnergyTypes.get("ENERGY", c_int(1))
        # flame model mode: freely propagating
        self._flamemode = 0
        # set required keywords for this model
        self.setkeyword("FREE", True)
        # solve the energy equation
        self.setkeyword("ENRG", True)
        # laminar flame speed [cm/sec]
        self.flamespeed = -1.0

    def skip_fix_T_solution(self, mode: bool = True):
        """
        Skip the step of finding the intermediate solution with fixed temperature

        Parameters
        ----------
            mode: boolean {True, False}
                ON/OFF
        """
        # skip the fixed temperature solution
        self.setkeyword("NOFT", value=mode)

    def automatic_temperature_profile_estimate(self, mode: bool = True):
        """
        Let the premixed flame model to construct an estimated temperature profile
        based on the equilibrium state to start the calculation

        Parameters
        ----------
            mode: boolean {True, False}
                ON/OFF
        """
        # use the automatic temperature profile estimate function
        self.setkeyword("TPROF", value=mode)
        if "TFIX" in self._keyword_index:
            msg = [
                Color.MAGENTA,
                "auto temperature profile option is ON,",
                "the pinned temperature is ignired.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            # remove the pinned temperature
            self.removekeyword("TFIX")

    def pinned_temperature(self, temperature: float = 400.0):
        """
        Pin the assigned temperature to the mesh to anchor the freely propagating flame.
        This temperature should be slightly higher than the inlet/unburnt gas temperature
        and less than the ignition temperature of the inlet gas mixture.

        Parameters
        ----------
            temperature: double, default = 400.0 [K]
                pinned gas temperature [K]
        """
        if temperature <= self.temperature:
            msg = [Color.PURPLE, "invalid temperature value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check
        if "TPROF" in self._keyword_index:
            msg = [
                Color.MAGENTA,
                "auto temperature profile option is ON,",
                "the pinned temperature is ignired.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        else:
            # set the pinned temperature
            self.setkeyword("TFIX", value=temperature)

    def get_flame_speed(self) -> float:
        """
        Get the computed laminar flame speed

        Returns
        -------
            flame_speed: double
                laminar flame speed [cm/sec]
        """
        # check solution
        if not self.getrawsolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'getsolution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return 0.0
        # return the computed laminar flame speed
        return self.flamespeed

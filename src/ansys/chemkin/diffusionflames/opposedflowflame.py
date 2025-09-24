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
Steady state, 1-D opposed-flow flame model.
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
from ansys.chemkin.chemistry import verify_version
from ansys.chemkin.flame import Flame
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
from ansys.chemkin.mixture import interpolate_mixtures
from ansys.chemkin.reactormodel import Keyword
from ansys.chemkin.utilities import find_interpolate_parameters
import numpy as np
import numpy.typing as npt


class OpposedFlame(Flame):
    def __init__(self, fuel_stream: Stream, label: Union[str, None] = None):
        """
        Axisymmetric/Cylindrical opposed-flow flame model

        Parameters
        ----------
            fuel_stream: Stream object
                the inlet stream on the "FUEL" side
            label: string
                reactor name
        """
        # check minimum version requirement = 2026 R1
        if not verify_version(261):
            exit()
        # check reactor Mixture object
        if not isinstance(fuel_stream, Stream):
            # wrong argument type
            msg = [Color.RED, "the first argument must be a Stream object.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # set label
        if label is None:
            self.label = "oppdifflame"
        else:
            self.label = label
        # set flow area to unity for easy conversion from mass flow rate to mass flux in the flame models
        if not fuel_stream._haveflowarea:
            fuel_stream.flowarea = 1.0  # [cm2]
        # initialization
        super().__init__(fuelstream=fuel_stream, label=label)
        # define coordinates
        self._geomkey = "AXIS"
        # define inlets
        # fuel-side
        self.fuelstream: Stream = self.reactormixture
        if self.fuelstream.label is None:
            self.fuelstream.label = "FUEL"
        # oxidizer-side
        self.oxidstream: Stream = None
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(163)
        # raw solution data structure
        self._solution_tags = [
            "distance",
            "temperature",
            "axial_velocity",
            "radial_velocity_gradient",
            "mixture_fraction",
        ]

    def set_oxidizer_inlet(self, oxid_stream: Stream):
        """
        Set the properties of the oxidizer side inlet

        Parameters
        ----------
            oxid_stream: Stream object
                oxidizer-side inlet stream
        """
        # There is only ONE additional inlet allowed for the opposed-flow flame model.
        if self.oxidstream is None:
            # check Inlet
            if not isinstance(oxid_stream, Stream):
                # wrong argument type
                msg = [Color.RED, "the argument must be a Stream object", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()
            else:
                # oxidizer stream is not set
                # clone the stream
                self.oxidstream = copy.deepcopy(oxid_stream)
                if self.oxidstream.label is None:
                    self.oxidstream.label = "OXIDIZER"
                # set flow area to unity for easy conversion from mass flow rate to mass flux in the flame models
                if not self.oxidstream._haveflowarea:
                    self.oxidstream.flowarea = 1.0  # [cm2]
        else:
            # the exidizer stream already exists
            msg = [
                Color.MAGENTA,
                "the oxidizer-side stream is already defined.\n",
                "opposed-flow flame model does NOT allow more than TWO inlet streams.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

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
        the initial grid points to start the simulation???

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
            this_key = ""
            this_key = tag + " " + str(x)
            self.setkeyword(this_key, True)
            numblines += 1
        # check error
        iErr = numblines - npoints
        return iErr

    def set_max_flame_temperature(self, max_temp: float):
        """
        Set the maximum temperature value in the plateau-like
        guessed tempeature profile in the gap space between the
        inlet nozzles when the default "plateau" profile is used

        Parameters
        ----------
            max_temp: double scalar, [K]
                estimated maximum flame temperature value
        """
        tag = "TPRO"
        if tag in self._profiles_index:
            # user temperature profile is provided
            msg =[
                Color.PURPLE,
                "a user temperature profile is already given,",
                "this maximum temperature setting is ignored.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        #
        if max_temp <= 200.0:
            # max temperature is too small
            msg =[Color.PURPLE, "max temperature value must > 200K.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            # set maximum flame temperature value
            self.setkeyword("TMAX", max_temp)

    def set_inlet_keywords(self,inlet: Stream) -> int:
        """
        Set the properties of the inlet stream

        Parameters
        ----------
            inlet: Stream object
                the stream object representing the inlet

        Returns
        -------
            error code: integer
        """
        iErr = 0
        # inlet label
        c_tag = bytes(inlet.label, "utf-8")
        # inlet temperature [K]
        temp = c_double(inlet.temperature)
        # inlet axial velocity [cm/sec]
        axial_vel = c_double(inlet.velocity)
        # flow type {0: velocity, 1: mass flux}
        flowtype = c_int(0)
        # number of gas species
        nspec = c_int(self.numbspecies)
        # inlet gas mass fractions
        Y_inlet = inlet.Y
        iErr = chemkin_wrapper.chemkin.KINOppdif_SetInlet(
            c_tag,
            nspec,
            temp,
            Y_inlet,
            axial_vel,
            flowtype,
        )
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to set up inlet",
                inlet.label,
                ", error code = [",
                str(iErr),
                "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        return iErr


    def __run_model(self) -> int:
        """
        Run the reactor model after the keywords are processed

        Returns
        -------
            Error code: integer
        """
        # run the opposed-flow flame simulation
        iErr = chemkin_wrapper.chemkin.KINOppdif_CalculateFlame(
            self._myLOUT,
            self._chemset_index,
            self._pressure,
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
        iErrkey = 0
        # set_verbose(True)
        # geometry (AXIS, PLAN)
        self.setkeyword(self._geomkey, 0.0)
        # energy equation: always solve the energy equation
        self.setkeyword("ENRG", 0.0)
        # set inlets
        # fuel stream
        iErr = self.set_inlet_keywords(inlet=self.fuelstream)
        if self.fuelstream.velocity_gradient > 0.0:
            tag = "AINL " + self.fuelstream.label
            self.setkeyword(tag, self.fuelstream.velocity_gradient)
        # oxidizer stream
        iErr = self.set_inlet_keywords(inlet=self.oxidstream)
        if self.oxidstream.velocity_gradient > 0.0:
            tag = "AINL " + self.oxidstream.label
            self.setkeyword(tag, self.oxidstream.velocity_gradient)
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
                    iErrc = chemkin_wrapper.chemkin.KINOppdif_SetParameter(this_key, this_value)
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
            msg = [Color.MAGENTA, "please run the flame simulation first.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        elif status != 0:
            msg = [
                Color.PURPLE,
                "simulation was failed.\n",
                Color.SPACEx6,
                "please correct the error(s) and rerun the flame simulation.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # insert the continuation keyword
        key_continue = bytes("CNTN", "utf-8") 
        this_value = c_double(0.0)
        iErr = chemkin_wrapper.chemkin.KINOppdif_SetParameter(key_continue, this_value)
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
            msg = [Color.MAGENTA, "please run the flame simulation first.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        elif status != 0:
            msg = [
                Color.PURPLE,
                "simulation was failed.\n",
                Color.SPACEx6,
                "please correct the error(s) and rerun the flame simulation.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # number of time points in the solution
        npoints = c_int(0)
        # get solution size of the opposed-flow flame solution
        iErr = chemkin_wrapper.chemkin.KINOppdif_GetSolutionGridPoints(npoints)
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
        iErr = chemkin_wrapper.chemkin.KINOppdif_GetSolution(
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
        # get the flow field
        npoint = c_int(npoints)
        # axial velocity [cm/sec]
        axial_vel = np.zeros_like(pos, dtype=np.double)
        # radial velocity gradient [1/sec]
        radial_vel = np.zeros_like(pos, dtype=np.double)
        iErr = chemkin_wrapper.chemkin.KINOppdif_GetVelocityField(
            npoint, axial_vel, radial_vel
        )
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to get the velocities,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # get mixture fraction [-]
        npoint = c_int(npoints)
        mix_frac = np.zeros_like(pos, dtype=np.double)
        iErr = chemkin_wrapper.chemkin.KINOppdif_GetMixtureFraction(
            npoint, mix_frac,
        )
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to get the mixture fractions,",
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
        # axial velocity
        self._solution_rawarray["axial_velocity"] = copy.deepcopy(axial_vel)
        # radial velocity gradient
        self._solution_rawarray["radial_velocity_gradient"] = copy.deepcopy(radial_vel)
        # mixture fraction
        self._solution_rawarray["mixture_fraction"] = copy.deepcopy(mix_frac)
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
        # clean up
        del pos, temp, axial_vel, radial_vel, mix_frac, frac

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
        if vname.lower() in self._solution_tags:
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
        # get the 
        axial_vel = self.get_solution_variable_profile("axial_velocity")
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
            # compute gas density [g/cm3]
            den = sstream.RHO
            # stream mass flux [g/cm2-sec]
            sstream.mass_flowrate = den * abs(axial_vel[i])
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


class OpposedFlame_Planar(OpposedFlame):
    def __init__(self, fuel_stream: Stream, label: Union[str, None] = None):
        """
        Planar opposed-flow flame model

        Parameters
        ----------
            fuel_stream: Stream object
                the inlet stream on the "FUEL" side
            label: string
                reactor name
        """
        # initialization
        OpposedFlame.__init__(fuel_stream=fuel_stream, label=label)
        # define coordinates
        self._geomkey = "PLAN"
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(164)

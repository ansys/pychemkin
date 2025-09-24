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
    Chemkin shock tube reactor model.
"""

import copy
from ctypes import c_double, c_int

from ansys.chemkin import chemkin_wrapper
from ansys.chemkin.chemistry import (
    check_chemistryset,
    chemistryset_initialized,
    force_activate_chemistryset,
    verify_version,
)
from ansys.chemkin.color import Color as Color
from ansys.chemkin.constants import Patm
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
from ansys.chemkin.mixture import equilibrium, interpolate_mixtures
from ansys.chemkin.reactormodel import ReactorModel as reactor
from ansys.chemkin.batchreactors.batchreactor import calculate_effective_activation_energy
from ansys.chemkin.utilities import find_interpolate_parameters, interpolate_point
import numpy as np
import numpy.typing as npt


class ShockTubeReactors(reactor):
    """
    Generic Chemkin 0-D transient shock tube reactor models
    """
    def __init__(self, reactor_condition: Stream, label: str):
        """
        Initialize a generic Shock Tube Reactor object

        Parameters
        ----------
            reactor_condition: Stream object
                a stream representing the initial gas properties relative to the shock wave
            label: string, optional
                reactor name
        """
        # initialize the base module
        super().__init__(reactor_condition, label=label)
        #
        # reactor parameters (required)
        self._volume = c_double(0.0e0)
        self._end_time = c_double(0.0e0)
        self._start_time = c_double(0.0e0)
        # boundary layer correction
        self.BL_correction = 0
        # shock tube diameter [cm]
        self.reactor_diameter = c_double(0.0e0)
        # initial gas viscocity at 300 [K]
        self.init_visc = c_double(0.0e0)
        # inlet property (pressure, temperature, and velocity) location
        # 1: before the incident shock
        # 2: behind the incident shock (for incident shock application)
        # or behind the reflected shock (for reflected shock application)
        self.iType = c_int(1)
        # flag for the reflected shock simulation with initial gas behind the reflected shock
        self.PFR_like = False
        # solver parameters
        self._absolute_tolerance = 1.0e-12
        self._relative_tolerance = 1.0e-6
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(170)
        # check required inputs
        self._numb_requiredinput = 0
        self._requiredlist: list[str] = []
        self._inputcheck: list[str] = []
        # raw solution data structure
        self._solution_tags.append(["density"])
        # single point solution variables
        self.soln_parameters_list = ["P", "T", "RHO", "V", "A", "MA", "CDET"]

    @property
    def time(self) -> float:
        """
        Get simulation end time (required) [sec]

        Returns
        -------
            endtime: double
                simulation duration or simulation end time [sec]
        """
        return self._end_time.value

    @time.setter
    def time(self, value: float = 0.0e0):
        """
        Set simulation end time (required)

        Parameters
        ----------
            value: double, default = 0.0
                simulation end time [sec]
        """
        if value <= 0.0e0:
            msg = [Color.PURPLE, "simulation end time must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._end_time = c_double(value)

    @property
    def starttime(self) -> float:
        """
        Get simulation start time [sec]

        Returns
        -------
            endtime: double
                simulation duration or simulation start time [sec]
        """
        return self._start_time.value

    @starttime.setter
    def starttime(self, value: float = 0.0e0):
        """
        Set simulation start time

        Parameters
        ----------
            value: double, default = 0.0
                simulation start time [sec]
        """
        self._start_time = c_double(value)

    @property
    def diameter(self) -> float:
        """
        Get shock tube diameter

        Returns
        -------
            diam: double
                Reactor diameter [cm]
        """
        return self.reactor_diameter.value

    @diameter.setter
    def diameter(self, diam: float):
        """
        Set the shock tube diameter for boundary layer correction

        Parameters
        ----------
            diam: double
                reactor diameter [cm]
        """
        if diam <= 0.0:
            msg = [Color.PURPLE, "reactor diameter must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self.BL_correction += 1
            self.reactor_diameter = c_double(diam)
            # set flow area at the inlet
            self.reactormixture._haveflowarea = True
            area = np.pi * diam * diam / 4.0
            self.reactormixture._flowarea = area
            self.reactorflowarea = area

    @property
    def tolerances(self) -> tuple:
        """
        Get transient solver tolerances

        Returns
        -------
            tolerances: tuple, [absolute_tolerance, relative_tolerance]
                absolute_tolerance: double
                    absolute tolerance
                relative_tolerance: double
                    relative tolerance
        """
        return (self._absolute_tolerance, self._relative_tolerance)

    @tolerances.setter
    def tolerances(self, tolerances: tuple[float, float]):
        """
        Set transient solver tolerances

        Parameters
        ----------
            tolerances: tuple, [absolute_tolerance, relative_tolerance]
                absolute_tolerance: double
                    absolute tolerance
                relative_tolerance: double
                    relative tolerance
        """
        # set tolerances
        if tolerances is not None:
            # set absolute tolerance
            self._absolute_tolerance = max(tolerances[0], 1.0e-20)
            # set keywords
            self.setkeyword(key="ATOL", value=self._absolute_tolerance)
            # set relative tolerance
            self._relative_tolerance = max(tolerances[1], 1.0e-12)
            # set keywords
            self.setkeyword(key="RTOL", value=self._relative_tolerance)

    def set_inlet_viscosity(self, visc: float):
        """
        Set the initial gas mixture viscocity for boundary layer correction

        Parameters
        ----------
            visc: double, default = 0.0
                mixture viscosity [g/cm-sec] or [Poise]
        """
        if visc <= 0.0:
            msg = [Color.PURPLE, "gas mixture viscosity must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            # set the initial gas viscosity
            self.BL_correction += 1
            self.init_visc = c_double(visc)

    def set_solver_max_timestep_size(self, size: float):
        """
        Set the maximum time step size allowed by the solver

        Parameters
        ----------
            size: double, default = 1/100 of the simulation duration
                step size [sec] or [cm]
        """
        if size > 0.0e0:
            self.setkeyword(key="STPT", value=size)
        else:
            msg = [Color.PURPLE, "solver timestep size must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_solution_printing_timestep_size(self, delta_time: float):
        """
        Set the timestep size between printing the solution data to the text output file

        Parameters
        ----------
            delta_time: double, default = 1/100 of the simulation duration
                timestep size between printing solution data [sec]
        """
        if delta_time > 0.0e0:
            self.setkeyword(key="DELT", value=delta_time)
        else:
            msg = [Color.PURPLE, "solution printing timestep size must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def get_solution_size(self) -> int:
        """
        Get solution size

        Returns
        -------
            npoints: integer
                number of time points in the solution profile
        """
        # check run completion
        status = self.getrunstatus(mode="silent")
        if status == -100:
            msg = [Color.MAGENTA, "please run the reactor simultion first.", Color.END]
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
        iErr = chemkin_wrapper.chemkin.KINShock_GetSolnResponseSize(npoints)
        if iErr == 0:
            # return the solution sizes
            self._numbsolutionpoints = (
                npoints.value
            )  # number of time points in the solution profile
            return self._numbsolutionpoints
        else:
            #
            return npoints.value

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

        msg = [Color.YELLOW, "post-processing raw solution data ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        # reset raw and mixture solution parameters
        self._numbsolutionpoints = 0
        self._solution_rawarray.clear()
        self._solution_mixturearray.clear()
        # get solution sizes
        npoints = self.get_solution_size()
        # check values
        if npoints == 0:
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
        time = np.zeros(self._numbsolutionpoints, dtype=np.double)
        pres = np.zeros_like(time, dtype=np.double)
        temp = np.zeros_like(time, dtype=np.double)
        vel = np.zeros_like(time, dtype=np.double)
        density = np.zeros_like(time, dtype=np.double)
        distance = np.zeros_like(time, dtype=np.double)
        # create a species mass fraction array to hold the solution species fraction profiles
        frac = np.zeros(
            (
                self.numbspecies,
                self._numbsolutionpoints,
            ),
            dtype=np.double,
            order="F",
        )
        # get raw solution data
        icnspec = c_int(self.numbspecies)
        icnpoints = c_int(npoints)
        iErr = chemkin_wrapper.chemkin.KINShock_GetGasSolnResponse(
            icnpoints, icnspec, time, temp, pres, vel, density, distance, frac
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
        # store the raw solution data in a dictionary
        # time
        self._solution_rawarray["time"] = copy.deepcopy(time)
        # temperature
        self._solution_rawarray["temperature"] = copy.deepcopy(temp)
        # pressure
        self._solution_rawarray["pressure"] = copy.deepcopy(pres)
        # velocity [cm/sec]
        self._solution_rawarray["velocity"] = copy.deepcopy(vel)
        # density [g/cm3]
        self._solution_rawarray["density"] = copy.deepcopy(density)
        # distance [cm]
        self._solution_rawarray["distance"] = copy.deepcopy(distance)
        # species mass fractions
        self.parsespeciessolutiondata(frac)
        # create solution mixture
        iErr = self.create_solution_streams(frac)
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "forming solution mixtures",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # calculate total termicity [1/sec]
        thermicity = np.zeros_like(time, dtype=np.double)
        for i, m in enumerate(self._solution_mixturearray):
            thermicity[i] = m.thermicity()
        # total thermicity
        self._solution_rawarray["thermicity"] = copy.deepcopy(thermicity)
        # set up single point solution parameters
        # gas state proeprties at the three reference locations
        # 1 = before incident shock
        # 2 = after incident shock
        # 3 = after reflected shock
        p_state = c_double(0.0)
        t_state = c_double(0.0)
        rho_state = c_double(0.0)
        vel_state = c_double(0.0)
        sound_state = c_double(0.0)
        mach_state = c_double(0.0)
        cdetwave = c_double(0.0)
        for i in range(3):
            location = c_int(i + 1)
            iErr = chemkin_wrapper.chemkin.KINShock_GetGasStates(
                location,
                p_state,  # pressure [atm]
                t_state,  # gas temperature [K]
                rho_state,  # gas density [g/cm3]
                vel_state,  # local gas velocity [cm/sec]
                sound_state,  # local speed of sound [cm/sec]
                mach_state,  # Mach number [-]
                cdetwave,  # detonation wave speed [cm/sec]
            )
            if iErr == 0:
                # convert pressure from [atm] to [dynes/cm2]
                p_state = c_double(p_state.value * Patm)
                soln_values = [
                    p_state,
                    t_state,
                    rho_state,
                    vel_state,
                    sound_state,
                    mach_state,
                    cdetwave,
                ]
                state_tag = str(i + 1)
                for j, v in enumerate(self.soln_parameters_list):
                    tag = v + state_tag
                    self._solution_parameters[tag] = soln_values[j].value
            else:
                msg = [
                    Color.RED,
                    "failed to fetch the signle point solution data from memory,",
                    "error code =",
                    str(iErr),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()
        # maximum thermicity [1/sec] and induction length [cm]
        # induction length where the maximum thermicity is located
        nindlength = self.get_inductionlength_size()
        if nindlength > 0:
            sigma = np.zeros(nindlength, dtype=np.double)
            indlen = np.zeros_like(sigma, dtype=np.double)
            sigma, indlen = self.get_induction_lengths(nindlength)
            if len(sigma) > 1:
                sigma_max = max(sigma)
                self._solution_parameters["max_thermicity"] = sigma_max
                self._solution_parameters["induction_length"] = indlen[sigma.index(sigma_max)]
            else:
                self._solution_parameters["max_thermicity"] = sigma[0]
                self._solution_parameters["induction_length"] = indlen[0]
        else:
            msg = [
                Color.MAGENTA,
                "No induction length found.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
        # clean up
        del time, pres, temp, vel, density, distance, thermicity, frac

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
        Create a list of Streams that represent the gas inside the reactor at a solution point

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
        # create a temporary Mixture object to hold the mixture properties at current solution point
        sstream = copy.deepcopy(self.reactormixture)
        # create variable arrays to hold the solution profile
        species = []
        # create a species fraction array to hold the solution species fraction profiles
        frac = np.zeros(self.numbspecies, dtype=np.double)
        # get solution variable profile from the raw solution arrays
        pres = self.get_solution_variable_profile("pressure")
        temp = self.get_solution_variable_profile("temperature")
        vel = self.get_solution_variable_profile("velocity")
        # for the case both the initial gas and the simulation domain are entirely after the reflected shock 
        if self.PFR_like:
            # set the velocities to a constant value of 1.0 [cm/sec]
            vel[:] = 1.0
        # loop over all species
        for sp in self._specieslist:
            species.append(self.get_solution_variable_profile(sp))
        # loop over all solution points
        for i in range(self._numbsolutionpoints):
            # get mixture properties at the current solution point
            # pressure [dynes/cm2]
            sstream.pressure = pres[i]
            # temperature [K]
            sstream.temperature = temp[i]
            # mixture volume [cm3]
            sstream.volume = 1.0
            # stream flow area [cm2]
            sstream.flowarea = 1.0
            # stream velocity [cm/sec]
            sstream.velocity = vel[i]
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
            # add to the solution mixture list
            self._solution_mixturearray.append(copy.deepcopy(sstream))
        # clean up
        species.clear()
        del pres, temp, vel, frac, species, sstream
        return 0

    def get_solution_stream(self, time: float) -> Stream:
        """
        Get the stream representing the solution state inside the reactor at the given time

        Parameters
        ----------
            time: double
                time point value [sec]

        Returns
        -------
            mixturetarget: Stream object
                a Stream representing the gas properties in the reactor at the specific time
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
        timearray = self.get_solution_variable_profile("time")
        # find the interpolation parameters
        ileft, ratio = find_interpolate_parameters(time, timearray)
        # find the stream
        if ratio == 0.0e0:
            # get the streams
            mixtureleft = copy.deepcopy(self._solution_mixturearray[ileft])
            return mixtureleft
        elif ratio == 1.0e0:
            # get the streams
            mixtureright = copy.deepcopy(self._solution_mixturearray[ileft + 1])
            return mixtureright
        else:
            # get the streams
            mixtureleft = copy.deepcopy(self._solution_mixturearray[ileft])
            mixtureright = copy.deepcopy(self._solution_mixturearray[ileft + 1])
            # interpolate the mixture properties
            mixturetarget = interpolate_mixtures(mixtureleft, mixtureright, ratio)
            # set velocity
            ratiom = 1.0e0 - ratio
            mixturetarget.velocity = ratiom * mixtureleft.velocity + ratio * mixtureright.velocity
            # clean up
            del mixtureleft, mixtureright
            #
            return mixturetarget

    def get_solution_stream_at_index(self, solution_index: int) -> Stream:
        """
        Get the stream representing the solution state inside the reactor at the given solution point index

        Parameters
        ----------
            solution_index: integer
                0-base solution time point index

        Returns
        -------
            mixturetarget: Stream Object
                a Stream representing the gas properties in the reactor at the specific time
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
        if solution_index > self._numbsolutionpoints - 1:
            msg = [
                Color.PURPLE,
                "the given time point index:",
                str(solution_index),
                "> the maximum number of time points:",
                str(self._numbsolutionpoints - 1),
                "\n",
                Color.SPACEx6,
                "the solution time point index is 0-based.\n",
                Color.SPACEx6,
                "[ 0 ->",
                str(self._numbsolutionpoints - 1),
                "]",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get the stream
        mixturetarget = copy.deepcopy(self._solution_mixturearray[solution_index])
        return mixturetarget

    def get_inductionlength_size(self) -> int:
        """
        Get number of induction length data points

        Returns
        -------
            npoints: integer
                number of induction length data points
        """
        # check run completion
        status = self.getrunstatus(mode="silent")
        if status == -100:
            msg = [Color.MAGENTA, "please run the reactor simultion first.", Color.END]
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
        iErr = chemkin_wrapper.chemkin.KINShock_GetInductLengthSize(npoints)
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to access the induction length solution data.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if npoints.value == 0:
            msg = [
                Color.YELLOW,
                "no induction length data available.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        return npoints.value

    def get_induction_lengths(self, npoints: int) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """
        Get induction lengths and local maximum total thermicities

        Parameters
        ----------
            npoints: integer
                number of induction length data points

        Returns
        -------
            induct_length: double array, dimension = npoints
                induction length [cm]
            sigmamax: double array, dimension = npoints
                local maximum total thermicity [1/sec]
        """
        # check run completion
        status = self.getrunstatus(mode="silent")
        if status == -100:
            msg = [Color.MAGENTA, "please run the reactor simultion first.", Color.END]
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
        npts = c_int(npoints)
        # create data lists
        induct_length = np.zeros (npoints, dtype=np.double)
        sigmamax = np.zeros_like (induct_length, dtype=np.double)
        # get solution size of the batch reactor
        iErr = chemkin_wrapper.chemkin.KINShock_GetInductLengths(npts, induct_length, sigmamax)
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to process the induction length solution data.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return induct_length, sigmamax

    def get_single_point_solution(self, varname: str) -> float:
        """
        Get the value of single point solution variable specified.

        Parameters
        ----------
            varname: string
                single point solution variable name

        Returns
        -------
            value: double
                value of the solution variable
        """
        # check run completion
        status = self.getrunstatus(mode="silent")
        if status == -100:
            msg = [Color.MAGENTA, "please run the reactor simultion first.", Color.END]
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
        # get the variable value
        this_value = self._solution_parameters.get(varname, None)
        if this_value is None:
            msg = [Color.MAGENTA, varname, "is not recognized.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        else:
            return this_value

class IncidentShock(ShockTubeReactors):
    """
    Chemkin incident shock reactor models
    """
    def __init__(self, mixture_condition: Stream, location: int, label: str):
        """
        Initialize a generic Shock Tube Reactor object

        Parameters
        ----------
            mixture_condition: Stream object
                a stream representing the initial gas mixture and velocity
                at the location speficied by the "location" parameter
            location: integer, {1, 2}
                location of the initial gas mixture relative to the incident shock; 1 = before, 2 = behind
            label: string, optional
                reactor name
        """
        # check minimum version requirement = 2026 R1
        if not verify_version(261):
            exit()
        # initialize the base module
        super().__init__(reactor_condition=mixture_condition, label=label)
        # initial gas location code
        if location == 2:
            self.iType = c_int(2)

    def __process_keywords(self) -> int:
        """
        Process input keywords for the shock tube model

        Returns
        -------
            error code: integer
        """
        iErr = 0
        iErrc = 0
        # set_verbose(True)
        # set keywords
        # pass all the keywords to the shock tube model
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
                iErrc = chemkin_wrapper.chemkin.KINShock_SetParameter(this_key, this_value)
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
        return iErr

    def __run_model(self) -> int:
        """
        Run the incident shock model after the keywords are processed

        Returns
        -------
            error code: integer
        """
        # run the simulation without keyword inputs
        shock_velocity = c_double(self.reactormixture.velocity)
        temp = c_double(self.reactormixture.temperature)
        pres = c_double(self.reactormixture.pressure)
        # inlet mass fraction
        y_init = self.reactormixture.Y
        #
        if self.BL_correction == 0:
            # without bounday layer correction
            iErr = chemkin_wrapper.chemkin.KINShock_CalcIncidentShockWithoutBoundaryLayerCorrection(
                self._myLOUT,
                self._chemset_index,
                shock_velocity,
                self.iType,
                temp,
                pres,
                y_init,
                self._start_time,
                self._end_time,
            )
        else:
            if self.BL_correction < 2:
                # missing parameter
                missed_parameter = ""
                if self.reactor_diameter.value <= 0.0e0:
                    missed_parameter = "tube diameter"
                if self.init_visc.value <= 0.0e0:
                    missed_parameter = "gas viscosity"
                msg = [
                    Color.PURPLE,
                    "missing parameter for boundary layer correction:",
                    missed_parameter,
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # with boundary layer correction
                #
                iErr = chemkin_wrapper.chemkin.KINShock_CalcIncidentShockWithBoundaryLayerCorrection(
                    self._myLOUT,
                    self._chemset_index,
                    shock_velocity,
                    self.iType,
                    temp,
                    pres,
                    y_init,
                    self._start_time,
                    self._end_time,
                    self.reactor_diameter,
                    self.init_visc,
                )

        return iErr

    def run(self) -> int:
        """
        Generic Chemkin run shock tube model method

        Returns
        -------
            error code: integer
        """
        # activate the Chemistry set associated with the Reactor instance
        force_activate_chemistryset(self._chemset_index.value)
        #
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
                self._chemset_index, c_int(0)
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
        # use API calls
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
        logger.debug("Processing keywords complete")

        # run reactor model
        msg = [Color.YELLOW, "running reactor simulation ...", Color.END]
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


class ReflectedShock(ShockTubeReactors):
    """
    Chemkin reflected shock reactor models
    """
    def __init__(self, mixture_condition: Stream, location: int, label: str):
        """
        Initialize a generic Shock Tube Reactor object

        Parameters
        ----------
            mixture_condition: Stream object
                a stream representing the gas mixture and velocity
                at the location speficied by the "location" parameter
            location: integer, {1, 2}
                location of the initial gas mixture; 1 = before incident shock, 2 = behind reflected shock
            label: string, optional
                reactor name
        """
        # check minimum version requirement = 2026 R1
        if not verify_version(261):
            exit()
        # initialize the base module
        super().__init__(reactor_condition=mixture_condition, label=label)
        # initial gas location code
        if location == 2:
            self.iType = c_int(2)
            # flag for the reflected shock simulation with initial gas behind the reflected shock
            self.PFR_like = True
            if mixture_condition._flowratemode < 0:
                # arbitrarily set the reflected shock velocity if it is not given
                self.reactormixture.velocity = 1.0
            if not mixture_condition._haveflowarea:
                # arbitrarily set the flow area if it is not given
                self.reactormixture._flowarea = 1.0
                self.reactorflowarea = 1.0
        # incident shock velocity [cm/sec]
        self.incidentshock_velocity = c_double(0.0e0)

    def __process_keywords(self) -> int:
        """
        Process input keywords for the shock tube model

        Returns
        -------
            error code: integer
        """
        iErr = 0
        iErrc = 0
        # set_verbose(True)
        # set keywords
        # pass all the keywords to the shock tube model
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
                iErrc = chemkin_wrapper.chemkin.KINShock_SetParameter(this_key, this_value)
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
        return iErr

    def __run_model(self) -> int:
        """
        Run the reflected shock model after the keywords are processed

        Returns
        -------
            error code: integer
        """
        # run the simulation without keyword inputs
        if self.iType.value == 1:
            # gas mixture is before the incident shock
            inc_shock_velocity = c_double(self.reactormixture.velocity)
            ref_shock_velocity = c_double(0.0e0)
        else:
            # gas mixture is behind the reflected shock
            inc_shock_velocity = c_double(0.0e0)
            ref_shock_velocity = c_double(self.reactormixture.velocity)
        temp = c_double(self.reactormixture.temperature)
        pres = c_double(self.reactormixture.pressure)
        # inlet mass fraction
        y_init = self.reactormixture.Y
        #
        iErr = chemkin_wrapper.chemkin.KINShock_CalcReflectedShock(
            self._myLOUT,
            self._chemset_index,
            inc_shock_velocity,
            ref_shock_velocity,
            self.iType,
            temp,
            pres,
            y_init,
            self._start_time,
            self._end_time,
        )

        return iErr

    def run(self) -> int:
        """
        Generic Chemkin run shock tube model method

        Returns
        -------
            error code: integer
        """
        # activate the Chemistry set associated with the Reactor instance
        force_activate_chemistryset(self._chemset_index.value)
        #
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
                self._chemset_index, c_int(0)
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
        # use API calls
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
        logger.debug("Processing keywords complete")

        # run reactor model
        msg = [Color.YELLOW, "running reactor simulation ...", Color.END]
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


class ZNDCalculator(ShockTubeReactors):
    """
    Chemkin Zeldovich-von Neumann-Doering (ZND) model for the detonation wave behind the incident shock
    """
    def __init__(self, mixture_condition: Stream, label: str):
        """
        Initialize a generic Shock Tube Reactor object

        Parameters
        ----------
            mixture_condition: Stream object
                a mixture representing the initial gas properties
                at the location speficied by the "location" parameter
            label: string, optional
                reactor name
        """
        # check minimum version requirement = 2026 R1
        if not verify_version(261):
            exit()
        # check velocity
        if mixture_condition._flowratemode == -1:
            # velocity is not set
            mixture_condition.velocity = 1.0e0
        # initialize the base module
        super().__init__(reactor_condition=mixture_condition, label=label)

    def __process_keywords(self) -> int:
        """
        Process input keywords for the shock tube model

        Returns
        -------
            error code: integer
        """
        iErr = 0
        iErrc = 0
        # set_verbose(True)
        # set keywords
        # pass all the keywords to the shock tube model
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
                iErrc = chemkin_wrapper.chemkin.KINShock_SetParameter(this_key, this_value)
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
        return iErr

    def __run_ZND_model(self) -> int:
        """
        Run the incident shock ZND model after the keywords are processed

        Returns
        -------
            error code: integer
        """
        # run the simulation without keyword inputs
        temp = c_double(self.reactormixture.temperature)
        pres = c_double(self.reactormixture.pressure)
        # inlet mass fraction
        y_init = self.reactormixture.Y
        #
        if self.BL_correction == 0:
            # without bounday layer correction
            iErr = chemkin_wrapper.chemkin.KINShock_CalcZNDWithoutBoundaryLayerCorrection(
                self._myLOUT,
                self._chemset_index,
                temp,
                pres,
                y_init,
                self._start_time,
                self._end_time,
            )
        else:
            if self.BL_correction < 2:
                # missing parameter
                missed_parameter = ""
                if self.reactor_diameter.value <= 0.0e0:
                    missed_parameter = "tube diameter"
                if self.init_visc.value <= 0.0e0:
                    missed_parameter = "gas viscosity"
                msg = [
                    Color.PURPLE,
                    "missing parameter for boundary layer correction:",
                    missed_parameter,
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # with boundary layer correction
                iErr = chemkin_wrapper.chemkin.KINShock_CalcZNDWithBoundaryLayerCorrection(
                    self._myLOUT,
                    self._chemset_index,
                    temp,
                    pres,
                    y_init,
                    self._start_time,
                    self._end_time,
                    self.reactor_diameter,
                    self.init_visc,
                )

        return iErr

    def run(self) -> int:
        """
        Generic Chemkin run shock tube model method

        Returns
        -------
            error code: integer
        """
        # activate the Chemistry set associated with the Reactor instance
        force_activate_chemistryset(self._chemset_index.value)
        #
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
                self._chemset_index, c_int(0)
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
        # use API calls
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
        logger.debug("Processing keywords complete")

        # run reactor model
        msg = [Color.YELLOW, "running reactor simulation ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        # use API calls
        retVal = self.__run_ZND_model()
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

    def calculate_Chi_ng(self, duration: float = 1.0e-3) -> float:
        """
        Calculate the Chi parameter

        Ng, H.D., Higgins, A.J., Kiyanda, C.B., Radulecu, M.I., Lee, J.H.S., Bates, K.R.,
        and Nikiforakis, N., Combustion Theory and Modeling, 9:159-170 (2005)

        Parameters
        ----------
            duration: double, default = 1.0e-3, optional
            reactor simulation time [sec]

        Returns
        -------
            chi: double
                the Ng Chi parameter
        """
        # gas mixture condition before the incident shock
        location1_mixture = copy.deepcopy(self.reactormixture)
        # pressure behind the incident shock [dynes/cm2]
        pres2 = self.get_single_point_solution("P2")
        # gas temperature behind the incident shock [K]
        temp2 = self.get_single_point_solution("T2")
        # gas velocity behind the incident shock [cm/sec]
        vel2 = self.get_single_point_solution("V2")
        # change the intial condition from 'before incident shock' to 'behind incident shock'
        location1_mixture.pressure = pres2
        location1_mixture.temperature = temp2
        # induction length where the maximum thermicity is located
        induct_length = self.get_single_point_solution("induction_length")
        sigma_max = self.get_single_point_solution("max_thermicity")
        # calculate the reduce activation energy
        reduced_EA, act_energy = calculate_effective_activation_energy(
            location1_mixture,
            duration=duration,
            model="CONV",
        )
        # reaction zone length
        if np.isclose(sigma_max, 0.00e0, atol=1.0e-8): 
            msg = [
                Color.PURPLE,
                "Maximum thermicity ~ 0.0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            react_length = vel2 / sigma_max
            # Ng's Chi parameter
            chi = reduced_EA * induct_length / react_length
            # clean up
            del location1_mixture
            return max(chi, 0.0e0)

    def calculate_cell_width_ng(self, duration: float = 1.0e-3) -> tuple[float, float]:
        """
        Calculate the detonation cell width based on the correlation by Ng et al.

        Ng, H.D., Ju, Y., and Lee, J.H.S., International Journal of Hydrogen Energy,
        32(1):93-99 (2007)

        Parameters
        ----------
            duration: double, default = 1.0e-3, optional
            reactor simulation time [sec]

        Returns
        -------
            width: double
                cell size [cm]
            chi: double
                the Ng Chi parameter
        """
        # Chi_ng
        chi = self.calculate_Chi_ng(duration)
        if chi <= 0.0e0:
            msg = [
                Color.PURPLE,
                "The Ng Chi parameter",
                str(chi),
                "<= 0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # induction length where the maximum thermicity is located
        induct_length = self.get_single_point_solution("induction_length")
        # Ng cell size model parameters
        a0 = 3.0465860763763e1
        a = [8.955438805808153e1, -1.30792822369483e2, 4.202450507117405e1]
        b = [-2.92912838385e-2, 1.0263250730647101e-5, -1.031921244571857e-9]
        # calculate cell width
        terma = 0.0e0
        termb = 0.0e0
        m = chi
        for i, q in enumerate(a):
            terma += q / m
            termb += b[i] * m
            m *= chi
        width = a0 + terma + termb
        width *= induct_length
        return width, chi

    def calculate_cell_width_gavrikov(self, limiting_species: str) -> float:
        """
        Calculate the detonation cell width based on the correlation by Ng et al.

        Gavrikov, A.I., Efimenko, A.A., and Dorofeev, S.B., Combustion and Flame,
        120:19-33 (2000)

        Parameters
        ----------
            limiting_species: string
                symbol of the limiting species
        Returns
        -------
            width: double
                cell size [cm]
        """
        # gas mixture condition before the incident shock
        location1_mixture = copy.deepcopy(self.reactormixture)
        # gas temperature before the incident shock [K]
        temp1 = self.get_single_point_solution("T1")
        # gas temperature behind the incident shock [K]
        temp2 = self.get_single_point_solution("T2")
        # pressure behind the incident shock [dynes/cm2]
        pres2 = self.get_single_point_solution("P2")
        # change the intial condition from 'before incident shock' to 'behind incident shock'
        location1_mixture.pressure = pres2
        location1_mixture.temperature = temp2
        # calculate the reduce activation energy
        reduced_EA, act_energy = calculate_effective_activation_energy(location1_mixture, model="CONV")
        if reduced_EA <= 0.0e0:
            msg = [
                Color.PURPLE,
                "The effective reduced activation energy",
                str(reduced_EA),
                "<= 0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # calculate the equilibrium mixture with const internal energy and specific volume
        equil_mixture = equilibrium(location1_mixture, opt=7)
        # index of the limiting species
        limit_index = equil_mixture.get_specindex(limiting_species)
        # the mole fraction profile of the limiting species
        molfractions: float = []
        for m in self._solution_mixturearray:
            molfractions.append(m.X[limit_index])
        # the initial mole fraction of the limiting species
        molfrac_init = location1_mixture.X[limit_index]
        # the equilibrium mole fraction of the limiting species
        # it should correspond to the minimum mole fraction for any fuel component
        molfrac_final = equil_mixture.X[limit_index]
        # 50% consumption mole fraction of the limiting species
        molfrac_50 = (molfrac_init + molfrac_final) * 0.5
        # find the shock tube distance where the limiting species mole fraction equals the 50% consumption value
        distance = self.get_solution_variable_profile("distance")
        # molfractions = self.get_solution_variable_profile(limiting_species)
        index_soln, dist_50 = interpolate_point(molfrac_50, molfractions, distance)

        # model constants
        a = -7.843787493e-3
        b = 1.777662961e-1
        c = 2.371845901e-2
        d = 1.477047968
        e = 1.545112957e-1
        f = 1.547021569e-2
        g = -1.446582357
        h = 8.730494354
        i = 4.599907939
        j = 7.443410379
        k = 4.058325462e-1
        exponent = 1.453392165
        #
        r1 = reduced_EA
        r2 = temp2 / temp1
        r3 = r1 ** exponent
        term = r2 * (a * r2 - b)
        term += r1 * (c * r1 - d + r2 * (e - f * r2 ))
        term += g * np.log(r2)
        term += h * np.log(r1)
        term += r2 * (i / r1 - k * r2 / r3)
        term -= j
        width = 1.0e1 ** term
        width *= dist_50
        # clean up
        del location1_mixture, distance, molfractions
        return width

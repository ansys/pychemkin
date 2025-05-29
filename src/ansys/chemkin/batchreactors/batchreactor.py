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
    Chemkin closed homogeneous reactor model.
"""

import copy
from ctypes import c_double, c_int
from typing import Union

from ansys.chemkin import chemkin_wrapper
from ansys.chemkin.chemistry import (
    check_chemistryset,
    chemistryset_initialized,
    force_activate_chemistryset,
    set_verbose,
    verbose,
)
from ansys.chemkin.color import Color as Color
from ansys.chemkin.constants import Patm
from ansys.chemkin.info import show_ignition_definitions
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
from ansys.chemkin.mixture import Mixture, interpolate_mixtures
from ansys.chemkin.reactormodel import Keyword
from ansys.chemkin.reactormodel import ReactorModel as reactor
from ansys.chemkin.utilities import find_interpolate_parameters
import numpy as np
import numpy.typing as npt

# import numpy.typing as npt


class BatchReactors(reactor):
    """
    Generic Chemkin 0-D transient closed homogeneous reactor model
    """

    # set possible types in batch reactors
    ReactorTypes: dict[str, int] = {
        "Batch": int(1),
        "PSR": int(2),
        "PFR": int(3),
        "HCCI": int(4),
        "SI": int(5),
        "DI": int(6),
    }
    SolverTypes: dict = {"Transient": int(1), "SteadyState": int(2)}
    EnergyTypes: dict = {"ENERGY": int(1), "GivenT": int(2)}
    ProblemTypes: dict = {"CONP": int(1), "CONV": int(2), "ICEN": int(3)}

    def __init__(self, reactor_condition: Stream, label: str):
        """
        Initialize a generic Batch Reactor object

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside the batch reactor
            label: string, optional
                reactor name
        """
        # initialize the base module
        super().__init__(reactor_condition, label)
        #
        # reactor parameters (required)
        self._volume = c_double(0.0e0)
        self._endtime = c_double(0.0e0)
        self._reactivearea = c_double(0.0e0)
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        # solver parameters
        self._absolute_tolerance = 1.0e-12
        self._relative_tolerance = 1.0e-6
        # check required inputs
        self._numb_requiredinput = 0
        self._requiredlist: list[str] = []
        self._inputcheck: list[str] = []
        # default number of reactors
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # default energy type
        self._reactortype = c_int(self.ReactorTypes.get("Batch", 1))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("CONP", 1))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # profile points
        self._profilesize = int(0)

    @property
    def volume(self) -> float:
        """
        Get reactor volume

        Returns
        -------
            volume: double
                reactor volume [cm3]
        """
        return self._volume.value

    @volume.setter
    def volume(self, value: float):
        """
        Set reactor volume (required)

        Parameters
        ----------
            value: double, default = 0.0
                reactor volume [cm3]
        """
        if value > 0.0e0:
            # set reactor volume
            self._volume = c_double(value)
            # set initial mixture volume
            self.reactormixture.volume = value
            # set volume keyword (not set by the setup calls)
            self.setkeyword(key="VOL", value=value)
        else:
            print(Color.PURPLE + "** reactor volume must be > 0", end=Color.END)

    @property
    def area(self) -> float:
        """
        Get reactive surface area

        Returns
        -------
            area: double
                surface area [cm2]
        """
        return self._reactivearea.value

    @area.setter
    def area(self, value: float = 0.0e0):
        """
        Set reactive surface area

        Parameters
        ----------
            value: double, default = 0.0
                surface area [cm2]
        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "reactor active surface area must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        else:
            self._reactivearea = c_double(value)
            # set internal surface area keyword for PFR
            if (
                self._reactortype.value == self.ReactorTypes.get("PFR", 3)
                or Keyword.noFullKeyword
            ):
                self.setkeyword(key="AREA", value=value)

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

    @property
    def force_nonnegative(self) -> bool:
        """
        Get the status of the forcing non-negative option of the transient solver

        Returns
        -------
            mode: boolean
                status of the non-negative solver option
        """
        if "NNEG" in self._keyword_index:
            # defined: find index
            i = self._keyword_index.index("NNEG")
            return self._keyword_list[i].value
        else:
            # not defined: return default value
            return False

    @force_nonnegative.setter
    def force_nonnegative(self, mode: bool = False):
        """
        Set the forcing non-negative solution option

        Parameters
        ----------
            mode: boolean, default = False
                turn the option ON/OFF
        """
        # set keyword
        self.setkeyword(key="NNEG", value=mode)

    def set_solver_initial_timestep_size(self, size: float):
        """
        Set the initial time step size to be used by the solver

        Parameters
        ----------
            size: double, default = determined by the solver
                step size [sec] or [cm]
        """
        if size > 0.0e0:
            self.setkeyword(key="HO", value=size)
        else:
            msg = [Color.PURPLE, "solver timestep size must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

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

    @property
    def timestep_for_saving_solution(self) -> float:
        """
        Get the timestep size between saving the solution data

        Returns
        -------
            delta_time: double
                timestep size between saving solution data [sec]
        """
        if "DTSV" in self._keyword_index:
            # defined: find index
            i = self._keyword_index.index("DTSV")
            return self._keyword_list[i].value
        else:
            # return default value (100th of the end time)
            if self._endtime.value > 0.0e0:
                return self._endtime.value / 1.0e2
            else:
                # not defined yet
                msg = [
                    Color.MAGENTA,
                    "solution saving timestep is not defined",
                    'because "end time" has not been set,',
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.warning(this_msg)
                return 0.0

    @timestep_for_saving_solution.setter
    def timestep_for_saving_solution(self, delta_time: float):
        """
        Set the timestep size between saving the solution data

        Parameters
        ----------
            delta_time: double, default = 1/100 of the simulation duration
                timestep size between saving solution data [sec]
        """
        if delta_time > 0.0e0:
            self.setkeyword(key="DTSV", value=delta_time)
        else:
            msg = [Color.PURPLE, "solution saving timestep size must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    @property
    def timestep_for_printing_solution(self) -> float:
        """
        Get the timestep size between printing the solution data to the text output file

        Returns
        -------
            delta_time: double
                timestep size between printing solution data [sec]
        """
        if "DELT" in self._keyword_index:
            # defined: find index
            i = self._keyword_index.index("DELT")
            return self._keyword_list[i].value
        else:
            # return default value (100th of the end time)
            if self._endtime.value > 0.0e0:
                return self._endtime.value / 1.0e2
            else:
                # not defined yet
                msg = [
                    Color.MAGENTA,
                    "solution printing timestep is not defined",
                    'because "end time" has not been set,',
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.warning(this_msg)
                return 0.0

    @timestep_for_printing_solution.setter
    def timestep_for_printing_solution(self, delta_time: float):
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

    def adaptive_solution_saving(
        self,
        mode: bool,
        value_change: Union[float, None] = None,
        target: Union[str, None] = None,
        steps: Union[int, None] = None,
    ):
        """
        Set up adaptive solution data saving

        Parameters
        ----------
            mode: boolean
                switch adaptive solution saving ON/OFF
            value_change: double, optional
                change in solution variable value between saving additional solution data
            target: string, optional
                the target variable that is used by the value_change option
            steps: integer, optional
                number of solver time steps between saving additional solution data
        """
        # turn ON/OFF the adaptive solution saving option
        self.setkeyword(key="ADAP", value=mode)
        self.setkeyword(key="NADAP", value=not mode)
        if not mode:
            # turn OFF the adaptive solution saving
            return
        # set options
        if steps is not None:
            # use number of solver time steps option:
            if steps <= 0.0:
                # non-positive number of steps
                msg = [
                    Color.PURPLE,
                    "the number of steps per adaptive solution saving must > 0.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
            else:
                # set parameters
                self.setkeyword(key="ADAP", value=True)
                self.setkeyword(key="ASTEPS", value=int(steps))
        elif value_change is not None:
            # use change in the solution variable value option:
            if target is None:
                # target variable is not given
                msg = [
                    Color.PURPLE,
                    "a reference variable is required for value-change adaptive saving.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
            elif not isinstance(target, str):
                # not given as string
                msg = [
                    Color.PURPLE,
                    "a reference variable is assigned as a string, e.g., 'OH'.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
            elif value_change <= 0.0:
                # non-positive change value
                msg = [
                    Color.PURPLE,
                    "the value change per adaptive solution saving must > 0.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
            else:
                # set parameters
                self.setkeyword(key="ADAP", value=True)
                self.setkeyword(key="AVAR", value=target)
                self.setkeyword(key="AVALUE", value=value_change)
        else:
            # use error
            msg = [
                Color.PURPLE,
                "need to specify either the number of steps\n",
                Color.SPACEx6,
                "or the 'change value' + 'target variable' pair.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_ignition_delay(
        self,
        method: str = "T_inflection",
        val: float = 0.0,
        target: str = "",
    ):
        """
        Set ignition detection criterion

        Parameters
        ----------
            method: string
                ignition definition/detection method
            val: double, optional
                temperature or temperature rise value associated with the ignition detection method specified
            target: string, optional
                target species symbol if the 'Species_peak' method is used
        """
        if isinstance(method, str):
            # ignition detection method assigned
            if method == "T_inflection":
                # use inflection points in the temperature profile
                self.setkeyword(key="TIFP", value=True)
            elif method == "T_rise":
                # use temperature rise
                if val <= 0.0:
                    msg = [
                        Color.PURPLE,
                        "temperature rise value must > 0.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                else:
                    self.setkeyword(key="DTIGN", value=val)
            elif method == "T_ignition":
                # use temperature value
                if val <= 0.0:
                    msg = [
                        Color.PURPLE,
                        "ignition temperature value must > 0.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                else:
                    self.setkeyword(key="TLIM", value=val)
            elif method == "Species_peak":
                # use species peak location
                if target not in self._specieslist:
                    # no species given
                    msg = [
                        Color.PURPLE,
                        "target species is assigned as a string, e.g., 'OH'.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                else:
                    self.setkeyword(key="KLIM", value=target)
            else:
                # incorrect ignition detection method given
                msg = [
                    Color.PURPLE,
                    "ignition definition",
                    method,
                    "is not recognized.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                show_ignition_definitions()
        else:
            # use error
            show_ignition_definitions()

    def stop_after_ignition(self):
        """
        Set the option to stop the simulation after ignition is detected
        """
        # stop the simulation after ignition is detected
        self.setkeyword(key="IGN_STOP", value=True)

    def get_ignition_delay(self) -> float:
        """
        Get the predicted ignition delay time from the transient reactor simulation

        Returns
        -------
            ignition_delay_time: double
                ignition delay time [msec] or [CA]
        """
        # initialization
        ignitiondelaytime = c_double(0.0e0)
        # check run status
        status = self.getrunstatus(mode="silent")
        if status == -100:
            msg = [
                Color.YELLOW,
                "simulation has yet to be run.\n",
                Color.SPACEx6,
                "please run the reactor simulation first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return ignitiondelaytime.value
        elif status != 0:
            msg = [
                Color.YELLOW,
                "simulation was failed.\n",
                Color.SPACEx6,
                "please correct the error(s) and rerun the reactor simulation.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return ignitiondelaytime.value

        # get the ignition delay time (batch reactor model [sec], engine model [CA])
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetIgnitionDelay(ignitiondelaytime)
        if iErr != 0:
            msg = [
                Color.MAGENTA,
                "potential bad ignition delay time value,",
                "error code =",
                str(iErr),
                "\n",
                Color.SPACEx6,
                "please check the text output and revisit the reactor/solver settings.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            return ignitiondelaytime.value
        # check reactor model
        if self._reactortype.value == self.ReactorTypes.get("Batch", 1):
            # check ignition delay time value
            if ignitiondelaytime.value <= 0.0:
                msg = [
                    Color.MAGENTA,
                    "potential bad ignition delay time value.\n",
                    Color.SPACEx6,
                    "please check the text output and revisit the reactor/solver settings.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.warning(this_msg)
                return ignitiondelaytime.value * 1.0e3
            else:
                # convert ignition  delay time from [sec] to [msec]
                return ignitiondelaytime.value * 1.0e3
        elif self._reactortype.value in [
            self.ReactorTypes.get("HCCI", 4),
            self.ReactorTypes.get("SI", 5),
            self.ReactorTypes.get("DI", 6),
        ]:
            # engine models
            msg = ["cylinder-averaged ignition delay time in CA"]
            Color.ckprint("info", msg)
            return ignitiondelaytime.value
        elif self._reactortype.value == self.ReactorTypes.get("PFR", 3):
            msg = ["PFR ignition delay in [cm]"]
            Color.ckprint("info", msg)
            # check ignition distance value
            if ignitiondelaytime.value <= 0.0:
                msg = [
                    Color.MAGENTA,
                    "potential bad ignition distance value.\n",
                    Color.SPACEx6,
                    "please check the text output and revisit the reactor/solver settings.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.warning(this_msg)
                return ignitiondelaytime.value
            else:
                # return ignition distance [cm]
                return ignitiondelaytime.value
        else:
            return ignitiondelaytime.value

    def set_volume_profile(
        self, x: npt.NDArray[np.double], vol: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor volume profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            vol: 1D double array
                volume value of the profile data [cm3]

        Returns
        -------
            error code: integer
        """
        if (
            self.ProblemTypes.get(self._problemtype.value) == "CONP"
            and self.EnergyTypes.get(self._energytype.value) == "GivenT"
        ):
            msg = [
                Color.PURPLE,
                "cannot constrain volume of a Given-Pressure Fixed-Temperature",
                "batch reactor",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 10
        else:
            keyword = "VPRO"
            iErr = self.setprofile(key=keyword, x=x, y=vol)
            return iErr

    def set_pressure_profile(
        self, x: npt.NDArray[np.double], pres: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor pressure profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            pres: 1D double array
                pressure value of the profile data [dynes/cm2]

        Returns
        -------
            error code: integer
        """
        if (
            self.ProblemTypes.get(self._problemtype.value) == "CONV"
            and self.EnergyTypes.get(self._energytype.value) == "GivenT"
        ):
            msg = [
                Color.PURPLE,
                "cannot constrain pressure of a Given-Volume Fixed-Temperature",
                "batch reactor",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 10
        else:
            keyword = "PPRO"
            iErr = self.setprofile(key=keyword, x=x, y=pres)
            return iErr

    def set_surfacearea_profile(
        self, x: npt.NDArray[np.double], area: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor reactive surface area profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            area: 1D double array
                reactive surface area value of the profile data [cm2]

        Returns
        -------
            error code: integer
        """
        keyword = "AINT"
        iErr = self.setprofile(key=keyword, x=x, y=area)
        return iErr

    def set_reactortype_keywords(self):
        """
        Set reactor type keywords under the Full-Keywords mode
        """
        # keyword headers
        # set solver types
        if self._solvertype.value == self.SolverTypes.get("Transient", 1):
            self.setkeyword(key="TRAN", value=True)
        else:
            self.setkeyword(key="STST", value=True)
        # set reactor related keywords
        if self._reactortype.value == self.ReactorTypes.get("Batch", 1):
            # batch reactors
            # set problem type
            if self._problemtype.value == self.ProblemTypes.get("CONP", 1):
                self.setkeyword(key="CONP", value=True)
            else:
                self.setkeyword(key="CONV", value=True)
            # set energy equation
            if self._energytype.value == self.EnergyTypes.get("ENERGY", 1):
                self.setkeyword(key="ENRG", value=True)
            else:
                self.setkeyword(key="TGIV", value=True)
        elif self._reactortype.value == self.ReactorTypes.get("PSR", 2):
            # PSR
            # set energy equation
            if self._energytype.value == self.EnergyTypes.get("ENERGY", 1):
                self.setkeyword(key="ENRG", value=True)
            else:
                self.setkeyword(key="TGIV", value=True)
        elif self._reactortype.value == self.ReactorTypes.get("PFR", 3):
            # PFR
            self.setkeyword(key="PLUG", value=True)
            # set energy equation
            if self._energytype.value == self.EnergyTypes.get("ENERGY", 1):
                self.setkeyword(key="ENRG", value=True)
            else:
                self.setkeyword(key="TGIV", value=True)
        else:
            # IC engine models
            self.setkeyword(key="ICEN", value=True)
            self.setkeyword(key="TRAN", value=True)
            # set energy equation
            self.setkeyword(key="ENRG", value=True)

    def set_reactorcondition_keywords(self):
        """
        Set reactor initial/estimated condition keywords under the Full-Keywords mode
        """
        self.setkeyword(key="PRES", value=self._pressure.value / Patm)
        self.setkeyword(key="TEMP", value=self._temperature.value)
        self.setkeyword(key="TIME", value=self._endtime.value)
        # initial mole fraction
        nspecieslines, species_lines = self.createspeciesinputlines(
            self._solvertype.value, threshold=1.0e-12, molefrac=self.reactormixture.X
        )
        for line in species_lines:
            self.setkeyword(key=line, value=True)

    def validate_inputs(self) -> int:
        """
        Check the required inputs before running the reactor simulation

        Returns
        -------
            error code: integer
        """
        iErr = 0
        # required inputs:
        if self._numb_requiredinput <= 0:
            # no required input
            return iErr
        else:
            if len(self._inputcheck) < self._numb_requiredinput:
                msg = [Color.PURPLE, "some required inputs are missing.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
            # verify required inputs one by one
            for k in self._requiredlist:
                if k not in self._inputcheck:
                    iErr += 1
                    msg = [Color.PURPLE, "missing required input:", k, Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)

            return iErr

    def __process_keywords_withFullInputs(self) -> int:
        """
        Process input keywords for the batch reactor model under the Full-Keyword mode

        Returns
        -------
            Error code: integer
        """
        iErr = 0
        set_verbose(True)
        # verify required inputs
        iErr = self.validate_inputs()
        if iErr != 0:
            msg = [Color.PURPLE, "missing required input keywords.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return iErr
        # re-size work arrays if profile is used
        if self._numbprofiles > 0:
            # find total profile data points
            numbprofilepoints = 0
            for p in self._profiles_list:
                numbprofilepoints += p.size
            if numbprofilepoints != self._profilesize:
                # re-size work arrays
                self._profilesize = numbprofilepoints
                ipoints = c_int(numbprofilepoints)
                iErrc = chemkin_wrapper.chemkin.KINAll0D_SetProfilePoints(ipoints)
                # setup reactor model working arrays
                if iErrc == 0:
                    iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                        self._myLOUT, self._chemset_index
                    )
                iErr += iErrc
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "profile data generation error, error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return iErr
        # prepare initial conditions
        # initial mass fraction
        Y_init = self.reactormixture.Y
        # surface sites (not applicable)
        Site_init = np.zeros(1, dtype=np.double)
        # bulk activities (not applicable)
        Bulk_init = np.zeros_like(Site_init, dtype=np.double)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("Batch"):
            iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupBatchInputs(
                self._chemset_index,
                self._endtime,
                self._temperature,
                self._pressure,
                self._volume,
                self._heat_loss_rate,
                self._reactivearea,
                Y_init,
                Site_init,
                Bulk_init,
            )
            iErr += iErrc
            if iErrc != 0:
                logger.error("failed to set up basic reactor keywords")
                return iErrc

        # set reactor type
        self.set_reactortype_keywords()
        # reactor initial/estimated condition
        self.set_reactorcondition_keywords()
        if iErr == 0 and self._numbprofiles > 0:
            # get keyword lines of all profiles
            iErrProf, nproflines, prof_lines = self.createprofileinputlines()
            iErr += iErrProf
            if iErrProf == 0:
                # set the profile keywords
                for pkey in prof_lines:
                    for line in pkey:
                        self.setkeyword(key=line, value=True)
            else:
                msg = [
                    Color.PURPLE,
                    "failed to set up profile keywords,",
                    "error code =",
                    str(iErrProf),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return iErrProf
        # solve integrated heat release rate due to chemical reactions
        if self.EnergyTypes.get("ENERGY") == self._energytype.value:
            self.setkeyword(key="QRGEQ", value=True)
        # add the END keyword
        self.setkeyword(key="END", value=True)
        # create input lines from additional user-specified keywords
        iErr, nlines = self.createkeywordinputlines()
        if iErr == 0:
            if verbose():
                msg = [
                    Color.YELLOW,
                    str(nlines),
                    "input lines are added.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
        else:
            msg = [
                Color.PURPLE,
                "failed to create additional keyword lines.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

        return iErr

    def __run_model_withFullInputs(self) -> int:
        """
        Run the batch reactor model after the keywords are processed under the Full-Keyword mode
        All keywords must be assigned

        Returns
        -------
            error code: integer
        """
        # get information about the keyword inputs
        # convert number of keyword lines
        nlines = c_int(self._numblines)
        # combine the keyword lines into one single string
        lines = "".join(self._keyword_lines)
        # convert string to byte
        longline = bytes(lines, "utf-8")
        # convert line lengths array
        linelength = np.zeros(shape=self._numblines, dtype=np.int32)
        linelength[:] = self._linelength[:]
        # run the simulation with keyword inputs
        iErr = chemkin_wrapper.chemkin.KINAll0D_CalculateInput(
            self._myLOUT, self._chemset_index, longline, nlines, linelength
        )
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to set up reactor keywords in Full mode,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

        return iErr

    def __process_keywords(self) -> int:
        """
        Process input keywords for the batch reactor model

        Returns
        -------
            Error code: integer
        """
        iErr = 0
        iErrc = 0
        iErrKey = 0
        iErrInputs = 0
        # set_verbose(True)
        # verify required inputs
        iErr = self.validate_inputs()
        if iErr != 0:
            msg = [Color.PURPLE, "missing required input keywords.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return iErr
        # re-size work arrays if profile is used
        if self._numbprofiles > 0:
            # find total profile data points
            numbprofilepoints = 0
            for p in self._profiles_list:
                numbprofilepoints += p.size
            if numbprofilepoints != self._profilesize:
                # re-size work arrays
                self._profilesize = numbprofilepoints
                ipoints = c_int(numbprofilepoints)
                iErrc = chemkin_wrapper.chemkin.KINAll0D_SetProfilePoints(ipoints)
                # setup reactor model working arrays
                if iErrc == 0:
                    iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                        self._myLOUT, self._chemset_index
                    )
                iErr += iErrc
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "profile data generation error, error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return iErr
        # prepare initial conditions
        # initial mass fraction
        Y_init = self.reactormixture.Y
        # surface sites (not applicable)
        Site_init = np.zeros(1, dtype=np.double)
        # bulk activities (not applicable)
        Bulk_init = np.zeros_like(Site_init, dtype=np.double)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("Batch"):
            iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupBatchInputs(
                self._chemset_index,
                self._endtime,
                self._temperature,
                self._pressure,
                self._volume,
                self._heat_loss_rate,
                self._reactivearea,
                Y_init,
                Site_init,
                Bulk_init,
            )
            iErr += iErrc
            if iErrc != 0:
                msg = [
                    Color.PURPLE,
                    "failed to set up basic reactor keywords,",
                    "error code =",
                    str(iErrc),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return iErrc
            # heat transfer (use additional keywords)
            # solver parameters (use additional keywords)
            # output controls (use additional keywords)
            # ROP (use additional keywords)
            # sensitivity (use additional keywords)
            # ignition delay (use additional keywords)
            # solve integrated heat release rate due to chemical reactions
            if self.EnergyTypes.get("ENERGY") == self._energytype.value:
                iErrc = chemkin_wrapper.chemkin.KINAll0D_IntegrateHeatRelease()
                iErr += iErrc
                if iErrc != 0:
                    msg = [
                        Color.PURPLE,
                        "failed to set up heat release keyword,",
                        "error code =",
                        str(iErrc),
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    return iErrc
        else:
            pass
        if iErr == 0 and self._numbprofiles > 0:
            for p in self._profiles_list:
                key = bytes(p.profilekey, "utf-8")
                npoints = c_int(p.size)
                x = p.pos
                y = p.value
                iErrProf = chemkin_wrapper.chemkin.KINAll0D_SetProfileParameter(
                    key, npoints, x, y
                )
                iErr += iErrProf
            if iErrProf != 0:
                msg = [
                    Color.PURPLE,
                    "failed to set up profile keywords,",
                    "error code =",
                    str(iErrProf),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return iErrProf
        if iErr == 0:
            # set additional keywords
            # create input lines from additional user-specified keywords
            iErrInputs, nlines = self.createkeywordinputlines()
            if iErrInputs == 0:
                # process additional keywords in _keyword_index and _keyword_lines
                for s in self._keyword_lines:
                    # convert string to byte
                    line = bytes(s, "utf-8")
                    # set additional keyword one by one
                    iErrKey = chemkin_wrapper.chemkin.KINAll0D_SetUserKeyword(line)
                    iErrInputs += iErrKey
                if iErrInputs == 0:
                    if verbose():
                        msg = [
                            Color.YELLOW,
                            str(nlines),
                            "additional input lines are added.",
                            Color.END,
                        ]
                        this_msg = Color.SPACE.join(msg)
                        logger.info(this_msg)
                else:
                    msg = [
                        Color.PURPLE,
                        "failed to create additional input lines,",
                        "error code =",
                        str(iErrInputs),
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
            else:
                msg = [
                    Color.PURPLE,
                    "failed to process additional keywords, error code =",
                    str(iErrInputs),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
        #
        iErr = iErr + iErrInputs + iErrKey

        return iErr

    def __run_model(self) -> int:
        """
        Run the batch reactor model after the keywords are processed

        Returns
        -------
            error code: integer
        """
        # run the simulation without keyword inputs
        iErr = chemkin_wrapper.chemkin.KINAll0D_Calculate(self._chemset_index)
        return iErr

    def run(self) -> int:
        """
        Generic Chemkin run reactor model method

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
        if Keyword.noFullKeyword:
            # use API calls
            retVal = (
                self.__process_keywords()
            )  # each reactor model subclass to perform its own keyword processing
        else:
            # use full keywords
            retVal = self.__process_keywords_withFullInputs()
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
        if Keyword.noFullKeyword:
            # use API calls
            retVal = self.__run_model()
        else:
            # use full keywords
            retVal = self.__run_model_withFullInputs()
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

    def get_solution_size(self) -> tuple[int, int]:
        """
        Get the number of reactors and the number of solution points

        Returns
        -------
            nreactor: integer
                number of reactors
            npoints: integer
                number of solution points
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
        # number of reactor
        nreac = c_int(0)
        # number of time points in the solution
        npoints = c_int(0)
        # get solution size of the batch reactor
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetSolnResponseSize(nreac, npoints)
        nreactors = nreac.value
        if iErr == 0 and nreactors == self._nreactors:
            # return the solution sizes
            self._numbsolutionpoints = (
                npoints.value
            )  # number of time points in the solution profile
            return self._nreactors, self._numbsolutionpoints
        elif self._nreactors == nreactors:
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
        else:
            # incorrect number of reactor (batch reactor is single reactor)
            msg = [
                Color.PURPLE,
                "incorrect number of reactor.\n",
                Color.SPACEx6,
                "the reactor model expects",
                str(self._nreactors),
                "reactors\n",
                Color.SPACEx6,
                str(nreactors),
                "found in the solution.",
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

        msg = [Color.YELLOW, "post-processing raw solution data ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        # reset raw and mixture solution parameters
        self._numbsolutionpoints = 0
        self._solution_rawarray.clear()
        self._solution_mixturearray.clear()
        # get solution sizes
        nreac, npoints = self.get_solution_size()
        # check values
        if npoints == 0 or nreac != self._nreactors:
            msg = [
                Color.PURPLE,
                "solution size error(s).\n",
                Color.SPACEx6,
                "number of solution points =",
                str(npoints),
                "\n",
                Color.SPACEx6,
                "number of reactor(s)/zone(s) =",
                str(nreac),
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
        vol = np.zeros_like(time, dtype=np.double)
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
        icreac = c_int(nreac)
        icnpts = c_int(npoints)
        icnspec = c_int(self.numbspecies)
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetGasSolnResponse(
            icreac, icnpts, icnspec, time, temp, pres, vol, frac
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
        # volume
        self._solution_rawarray["volume"] = copy.deepcopy(vol)
        # species mass fractions
        self.parsespeciessolutiondata(frac)
        # create solution mixture
        iErr = self.create_solution_mixtures(frac)
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
        # clean up
        del time, pres, temp, vol, frac

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

    def create_solution_mixtures(self, specfrac: npt.NDArray[np.double]) -> int:
        """
        Create a list of Mixtures that represent the gas inside the reactor at a solution point

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
        smixture = copy.deepcopy(self.reactormixture)
        # create variable arrays to hold the solution profile
        species = []
        # create a species fraction array to hold the solution species fraction profiles
        frac = np.zeros(self.numbspecies, dtype=np.double)
        # get solution variable profile from the raw solution arrays
        pres = self.get_solution_variable_profile("pressure")
        temp = self.get_solution_variable_profile("temperature")
        vol = self.get_solution_variable_profile("volume")
        # loop over all species
        for sp in self._specieslist:
            species.append(self.get_solution_variable_profile(sp))
        # loop over all solution points
        for i in range(self._numbsolutionpoints):
            # get mixture properties at the current solution point
            # pressure [dynes/cm2]
            smixture.pressure = pres[i]
            # temperature [K]
            smixture.temperature = temp[i]
            # mixture volume [cm3]
            smixture.volume = vol[i]
            # species composition
            for k in range(self.numbspecies):
                frac[k] = specfrac[k, i]
            # set mixture composition
            if self._speciesmode == "mass":
                # mass fractions
                smixture.Y = frac
            else:
                # mole fractions
                smixture.X = frac
            # add to the solution mixture list
            self._solution_mixturearray.append(copy.deepcopy(smixture))
        # clean up
        species.clear()
        del pres, temp, vol, frac, species, smixture
        return 0

    def get_solution_mixture(self, time: float) -> Mixture:
        """
        Get the mixture representing the solution state inside the reactor at a given time

        Parameters
        ----------
            time: double
                time point value [sec]

        Returns
        -------
            mixturetarget: Mixture object
                a Mixture representing the gas properties in the reactor at the specific time
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
            # clean up
            del mixtureleft, mixtureright
            #
            return mixturetarget

    def get_solution_mixture_at_index(self, solution_index: int) -> Mixture:
        """
        Get the mixture representing the solution state inside the reactor at the given solution point index

        Parameters
        ----------
            solution_index: integer
                0-base solution time point index

        Returns
        -------
            mixturetarget: Mixture object
                a Mixture representing the gas properties in the reactor at the specific time
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
        # get the mixture
        mixturetarget = copy.deepcopy(self._solution_mixturearray[solution_index])
        return mixturetarget


class GivenPressureBatchReactor_FixedTemperature(BatchReactors):
    """
    Chemkin 0-D transient closed homogeneous reactor model
    with given reactor pressure (CONP) and reactor temperature (TGIV)
    """

    def __init__(self, reactor_condition: Stream, label: str = "CONPT"):
        """
        Initialize a Given Pressure Batch Reactor object with given reactor temperature

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside the batch reactor
            label: string, optional
                reactor name
        """
        # initialize the base module
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("Batch", 1))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("CONP", 1))
        self._energytype = c_int(self.EnergyTypes.get("GivenT", 2))
        # defaults for all closed homogeneous reactor models
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        # solver parameters
        self._absolute_tolerance = 1.0e-12
        self._relative_tolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]
        # set up basic batch reactor parameters
        iErr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if iErr == 0:
            # setup reactor model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to initialize the reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # if full-keyword mode is turned ON
        if not Keyword.noFullKeyword:
            # populate the reactor setup keywords
            self.set_reactortype_keywords()

    @property
    def time(self) -> float:
        """
        Get simulation end time (required) [sec]

        Returns
        -------
            endtime: double
                simulation duration or simulation end time [sec]
        """
        return self._endtime.value

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
            self._inputcheck.append("TIME")
            self._endtime = c_double(value)

    def set_temperature_profile(
        self, x: npt.NDArray[np.double], temp: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor temperature profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            temp: 1D double array
                temperature value of the profile data [K]

        Returns
        -------
            error code: integer
        """
        keyword = "TPRO"
        iErr = self.setprofile(key=keyword, x=x, y=temp)
        return iErr


class GivenPressureBatchReactor_EnergyConservation(BatchReactors):
    """
    Chemkin 0-D transient closed homogeneous reactor model
    with given reactor pressure (CONP) and
    solving the energy equation (ENRG)
    """

    def __init__(self, reactor_condition: Stream, label: str = "CONP"):
        """
        Initialize a Given Pressure Batch Reactor object that solves the Energy Equation

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside the batch reactor
            label: string, optional
                reactor name
        """
        # initialize the base module
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("Batch", 1))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("CONP", 1))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # defaults for all closed homogeneous reactor models
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        self._heat_transfer_coefficient = 0.0e0
        self._ambient_temperature = 3.0e2
        self._heat_transfer_area = 0.0e0
        # solver parameters
        self._absolute_tolerance = 1.0e-12
        self._relative_tolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]
        # set up basic batch reactor parameters
        iErr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if iErr == 0:
            # setup reactor model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to initialize the batch reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # if full-keyword mode is turned ON
        if not Keyword.noFullKeyword:
            # populate the reactor setup keywords
            self.set_reactortype_keywords()

    @property
    def time(self) -> float:
        """
        Get simulation end time (required)

        Returns
        -------
            endtime: double
                simulation duration or simulation end time [sec]
        """
        return self._endtime.value

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
            self._inputcheck.append("TIME")
            self._endtime = c_double(value)

    @property
    def heat_loss_rate(self) -> float:
        """
        Get heat loss rate from the reactor to the surroundings

        Returns
        -------
            heat_loss_rate: double
                heat loss rate [cal/sec]
        """
        return self._heat_loss_rate.value

    @heat_loss_rate.setter
    def heat_loss_rate(self, value: float):
        """
        Set the heat loss rate from the reactor to the surroundings

        Parameters
        ----------
            value: double, default = 0.0
                heat loss rate [cal/sec]
        """
        self._heat_loss_rate = c_double(value)
        if not Keyword.noFullKeyword:
            self.setkeyword(key="QLOS", value=value)

    @property
    def heat_transfer_coefficient(self) -> float:
        """
        Get heat transfer coefficient between the reactor and the surroundings

        Returns
        -------
            heat_transfer_coefficient: double
                heat transfer coefficient [cal/cm2-K-sec]
        """
        return self._heat_transfer_coefficient

    @heat_transfer_coefficient.setter
    def heat_transfer_coefficient(self, value: float = 0.0e0):
        """
        Set heat transfer coefficient between the reactor and the surroundings

        Parameters
        ----------
            value: double, optional, default = 0.0
                heat transfer coefficient [cal/cm2-K-sec]
        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer coefficient must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._heat_transfer_coefficient = value
            # set the corresponding keyword
            self.setkeyword(key="HTC", value=value)

    @property
    def ambient_temperature(self) -> float:
        """
        Get ambient temperature

        Returns
        -------
            temperature: double
                ambient temperature [K]
        """
        return self._ambient_temperature

    @ambient_temperature.setter
    def ambient_temperature(self, value: float = 0.0e0):
        """
        Set ambient temperature

        Parameters
        ----------
            value: double, optional, default = 300.0
                ambient temperature [K]
        """
        if value <= 0.0e0:
            msg = [Color.PURPLE, "ambient temperature must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._ambient_temperature = value
            # set the corresponding keyword
            self.setkeyword(key="TAMB", value=value)

    @property
    def heat_transfer_area(self) -> float:
        """
        Get heat transfer area between the reactor and the surroundings

        Returns
        -------
            area: double, optional, default = 0.0
                heat transfer area [cm2]
        """
        return self._heat_transfer_area

    @heat_transfer_area.setter
    def heat_transfer_area(self, value: float = 0.0e0):
        """
        Set heat transfer area between the reactor and the surroundings

        Parameters
        ----------
            value: double, optional, default = 0.0
                heat transfer area [cm2]
        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer area must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._heat_transfer_area = value
            # set the corresponding keyword
            self.setkeyword(key="AREAQ", value=value)

    def set_heat_transfer_area_profile(
        self, x: npt.NDArray[np.double], area: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor heat transfer area profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            area: 1D double array
                heat transfer area value of the profile data [cm2]

        Returns
        -------
            error code: integer
        """
        if self.EnergyTypes.get(self._energytype.value) == "GivenT":
            msg = [
                Color.PURPLE,
                "cannot specify heat transfer area to",
                "a Fixed-Temperature batch reactor",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 10
        else:
            keyword = "AEXT"
            iErr = self.setprofile(key=keyword, x=x, y=area)
            return iErr

    def set_heat_loss_profile(
        self, x: npt.NDArray[np.double], Qloss: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor heat loss rate profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            Qloss: 1D double array
                heat loss rate value of the profile data [cal/sec]

        Returns
        -------
            error code: integer
        """
        if self.EnergyTypes.get(self._energytype.value) == "GivenT":
            msg = [
                Color.PURPLE,
                "cannot specify heat loss rate to",
                "a Fixed-Temperature batch reactor",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 10
        else:
            keyword = "QPRO"
            iErr = self.setprofile(key=keyword, x=x, y=Qloss)
            return iErr


class GivenVolumeBatchReactor_FixedTemperature(BatchReactors):
    """
    Chemkin 0-D transient closed homogeneous reactor model
    with given reactor volume (CONV) and reactor temperature (TGIV)
    """

    def __init__(self, reactor_condition: Stream, label: str = "CONVT"):
        """
        Initialize a Given Volume Batch Reactor object with given reactor temperature

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside the batch reactor
            label: string, optional
                reactor name
        """
        # initialize the base module
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("Batch", 1))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("CONV", 2))
        self._energytype = c_int(self.EnergyTypes.get("GivenT", 2))
        # defaults for all closed homogeneous reactor models
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        # solver parameters
        self._absolute_tolerance = 1.0e-12
        self._relative_tolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]
        # set up basic batch reactor parameters
        iErr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if iErr == 0:
            # setup reactor model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to initialize the reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # if full-keyword mode is turned ON
        if not Keyword.noFullKeyword:
            # populate the reactor setup keywords
            self.set_reactortype_keywords()

    @property
    def time(self) -> float:
        """
        Get simulation end time (required)

        Returns
        -------
            endtime: double
                simulation duration or simulation end time [sec]
        """
        return self._endtime.value

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
            self._inputcheck.append("TIME")
            self._endtime = c_double(value)

    def set_temperature_profile(
        self, x: npt.NDArray[np.double], temp: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor temperature profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            temp: 1D double array
                temperature value of the profile data [K]

        Returns
        -------
            error code: integer
        """
        keyword = "TPRO"
        iErr = self.setprofile(key=keyword, x=x, y=temp)
        return iErr


class GivenVolumeBatchReactor_EnergyConservation(BatchReactors):
    """
    Chemkin 0-D transient closed homogeneous reactor model
    with given reactor volume (CONV) and
    solving the energy equation (ENRG)
    """

    def __init__(self, reactor_condition: Stream, label: str = "CONV"):
        """
        Initialize a Given Volume Batch Reactor object that solves the Energy Equation

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside the batch reactor
            label: string, optional
                reactor name
        """
        # initialize the base module
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("Batch", 1))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("CONV", 2))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # defaults for all closed homogeneous reactor models
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        self._heat_transfer_coefficient = 0.0e0
        self._ambient_temperature = 3.0e2
        self._heat_transfer_area = 0.0e0
        # solver parameters
        self._absolute_tolerance = 1.0e-12
        self._relative_tolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]
        # set up basic batch reactor parameters
        iErr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if iErr == 0:
            # setup reactor model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to initialize the reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # if full-keyword mode is turned ON
        if not Keyword.noFullKeyword:
            # populate the reactor setup keywords
            self.set_reactortype_keywords()

    @property
    def time(self) -> float:
        """
        Get simulation end time (required)

        Returns
        -------
            endtime: double
                simulation duration or simulation end time [sec]
        """
        return self._endtime.value

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
            self._inputcheck.append("TIME")
            self._endtime = c_double(value)

    @property
    def heat_loss_rate(self) -> float:
        """
        Get heat loss rate from the reactor to the surroundings

        Returns
        -------
            heat_loss_rate: double
                heat loss rate [cal/sec]
        """
        return self._heat_loss_rate.value

    @heat_loss_rate.setter
    def heat_loss_rate(self, value: float):
        """
        Set the heat loss rate from the reactor to the surroundings (required)

        Parameters
        ----------
            value: double, default = 0.0
                heat loss rate [cal/sec]
        """
        self._heat_loss_rate = c_double(value)
        if not Keyword.noFullKeyword:
            self.setkeyword(key="QLOS", value=value)

    @property
    def heat_transfer_coefficient(self) -> float:
        """
        Get heat transfer coefficient between the reactor and the surroundings

        Returns
        -------
            heat_transfer_coefficient: double
                heat transfer coefficient [cal/cm2-K-sec]
        """
        return self._heat_transfer_coefficient

    @heat_transfer_coefficient.setter
    def heat_transfer_coefficient(self, value: float = 0.0e0):
        """
        Set heat transfer coefficient between the reactor and the surroundings

        Parameters
        ----------
            value: double, optional, default = 0.0
                heat transfer coefficient [cal/cm2-K-sec]
        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer coefficient must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._heat_transfer_coefficient = value
            # set the corresponding keyword
            self.setkeyword(key="HTC", value=value)

    @property
    def ambient_temperature(self) -> float:
        """
        Get ambient temperature

        Returns
        -------
            temperature: double
                ambient temperature [K]
        """
        return self._ambient_temperature

    @ambient_temperature.setter
    def ambient_temperature(self, value: float = 0.0e0):
        """
        Set ambient temperature

        Parameters
        ----------
            value: double, optional, default = 300.0
                ambient temperature [K]
        """
        if value <= 0.0e0:
            msg = [Color.PURPLE, "ambient temperature must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._ambient_temperature = value
            # set the corresponding keyword
            self.setkeyword(key="TAMB", value=value)

    @property
    def heat_transfer_area(self) -> float:
        """
        Get heat transfer area between the reactor and the surroundings

        Returns
        -------
            area: double, optional, default = 0.0
                heat transfer area [cm2]
        """
        return self._heat_transfer_area

    @heat_transfer_area.setter
    def heat_transfer_area(self, value: float = 0.0e0):
        """
        Set heat transfer area between the reactor and the surroundings

        Parameters
        ----------
            value: double, optional, default = 0.0
                heat transfer area [cm2]
        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer area must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._heat_transfer_area = value
            # set the corresponding keyword
            self.setkeyword(key="AREAQ", value=value)

    def set_heat_transfer_area_profile(
        self, x: npt.NDArray[np.double], area: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor heat transfer area profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            area: 1D double array
                heat transfer area value of the profile data [cm2]

        Returns
        -------
            error code: integer
        """
        if self.EnergyTypes.get(self._energytype.value) == "GivenT":
            msg = [
                Color.PURPLE,
                "cannot specify heat transfer area to",
                "a Fixed-Temperature batch reactor",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 10
        else:
            keyword = "AEXT"
            iErr = self.setprofile(key=keyword, x=x, y=area)
            return iErr

    def set_heat_loss_profile(
        self, x: npt.NDArray[np.double], Qloss: npt.NDArray[np.double]
    ) -> int:
        """
        Specify reactor heat loss rate profile

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            Qloss: 1D double array
                heat loss rate value of the profile data [cal/sec]

        Returns
        -------
            error code: integer
        """
        if self.EnergyTypes.get(self._energytype.value) == "GivenT":
            msg = [
                Color.PURPLE,
                "cannot specify heat loss rate to",
                "a Fixed-Temperature batch reactor",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 10
        else:
            keyword = "QPRO"
            iErr = self.setprofile(key=keyword, x=x, y=Qloss)
            return iErr

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

"""Chemkin closed homogeneous reactor model."""

import copy
from ctypes import c_double, c_int
from typing import Union

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core import chemkin_wrapper
from ansys.chemkin.core.chemistry import (
    check_chemistryset,
    chemistryset_initialized,
    force_activate_chemistryset,
    set_verbose,
    verbose,
)
from ansys.chemkin.core.color import Color as Color
from ansys.chemkin.core.constants import P_ATM, R_GAS_CAL
from ansys.chemkin.core.info import show_ignition_definitions
from ansys.chemkin.core.inlet import Stream
from ansys.chemkin.core.mixture import Mixture, interpolate_mixtures
from ansys.chemkin.core.reactormodel import Keyword, ReactorModel as Reactor
from ansys.chemkin.core.utilities import (
    critical_and_exit,
    error_and_exit,
    find_interpolate_parameters,
    log_critical_message,
    log_error_message,
    log_info_message,
    log_warning_message,
)
from ansys.chemkin.core.validation import validate_minimum_value


class BatchReactors(Reactor):
    """Generic Chemkin 0-D transient closed homogeneous reactor model."""

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
        """Initialize a generic Batch Reactor object."""
        """
        Initialize a generic Batch Reactor object.

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside
                the batch reactor
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
        self._ninlets = np.zeros(1, dtype=np.int32)
        self._nzones = c_int(0)
        # default energy type
        self._reactortype = c_int(self.ReactorTypes.get("Batch", 1))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("CONP", 1))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # profile points
        self._profilesize = int(0)
        # text output control flag
        self.suppress_output = False
        # report ignition delay flag
        self.report_ignition_delay: bool = False

    @property
    def time(self) -> float:
        """Get simulation end time."""
        """
        Get simulation end time (required).

        Returns
        -------
            endtime: double
                simulation duration or simulation end time [sec]

        """
        return self._endtime.value

    @time.setter
    def time(self, value: float = 0.0e0):
        """Set simulation end time."""
        """
        Set simulation end time (required).

        Parameters
        ----------
            value: double, default = 0.0
                simulation end time [sec]

        """
        validate_minimum_value(
            value=value,
            minimum=np.finfo(float).eps,
            message="simulation end time must > 0.",
        )
        self._inputcheck.append("TIME")
        self._endtime = c_double(value)

    @property
    def volume(self) -> float:
        """Get reactor volume."""
        """
        Get reactor volume.

        Returns
        -------
            volume: double
                reactor volume [cm3]

        """
        return self._volume.value

    @volume.setter
    def volume(self, value: float):
        """Set reactor volume."""
        """
        Set reactor volume (required).

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
        """Get reactive surface area."""
        """
        Get reactive surface area.

        Returns
        -------
            area: double
                surface area [cm2]

        """
        return self._reactivearea.value

    @area.setter
    def area(self, value: float = 0.0e0):
        """Set reactive surface area."""
        """
        Set reactive surface area.

        Parameters
        ----------
            value: double, default = 0.0
                surface area [cm2]

        """
        validate_minimum_value(
            value=value,
            minimum=0.0,
            message="reactor active surface area must >= 0.",
        )
        self._reactivearea = c_double(value)
        # set internal surface area keyword for PFR
        if (
            self._reactortype.value == self.ReactorTypes.get("PFR", 3)
            or Keyword.no_fullkeyword
        ):
            self.setkeyword(key="AREA", value=value)

    @property
    def tolerances(self) -> tuple:
        """Get transient solver tolerances."""
        """
        Get transient solver tolerances.

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
        """Set transient solver tolerances."""
        """
        Set transient solver tolerances.

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
        """Get the status of the forcing non-negative option."""
        """Get the status of the forcing non-negative option
        of the transient solver.

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
        """Set the forcing non-negative solution option."""
        """
        Set the forcing non-negative solution option.

        Parameters
        ----------
            mode: boolean, default = False
                turn the option ON/OFF

        """
        # set keyword
        self.setkeyword(key="NNEG", value=mode)

    def set_legacy_option(self, option: bool):
        """Turn ON or OFF the legacy transient solver."""
        """Turn ON (True) or OFF (False) the legacy transient solver.

        Parameters
        ----------
            option: boolean
                turn On the legacy solver

        """
        if isinstance(option, bool):
            if option:
                self.setkeyword(key="USE_LEGACY_TECHNIQUE", value="4X")
        else:
            msg = [Color.PURPLE, "parameter must be either True or False.", Color.END]
            log_error_message(msg)

    def set_solver_initial_timestep_size(self, size: float):
        """Set the initial time step size."""
        """
        Set the initial time step size to be used by the solver.

        Parameters
        ----------
            size: double, default = determined by the solver
                step size [sec] or [cm]

        """
        if size > 0.0e0:
            self.setkeyword(key="HO", value=size)
        else:
            msg = [Color.PURPLE, "solver timestep size must > 0.", Color.END]
            log_error_message(msg)

    def set_solver_max_timestep_size(self, size: float):
        """Set the maximum time step size allowed."""
        """Set the maximum time step size allowed by the solver.

        Parameters
        ----------
            size: double, default = 1/100 of the simulation duration
                step size [sec] or [cm]

        """
        if size > 0.0e0:
            self.setkeyword(key="STPT", value=size)
        else:
            msg = [Color.PURPLE, "solver timestep size must > 0.", Color.END]
            log_error_message(msg)

    @property
    def timestep_for_saving_solution(self) -> float:
        """Get the timestep size between saving the solution data."""
        """
        Get the timestep size between saving the solution data.

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
                log_warning_message(msg)
                return 0.0

    @timestep_for_saving_solution.setter
    def timestep_for_saving_solution(self, delta_time: float):
        """Set the timestep size between saving the solution data."""
        """
        Set the timestep size between saving the solution data.

        Parameters
        ----------
            delta_time: double, default = 1/100 of the simulation duration
                timestep size between saving solution data [sec]

        """
        if delta_time > 0.0e0:
            self.setkeyword(key="DTSV", value=delta_time)
        else:
            msg = [Color.PURPLE, "solution saving timestep size must > 0.", Color.END]
            log_error_message(msg)

    @property
    def timestep_for_printing_solution(self) -> float:
        """Get the timestep size between printing."""
        """Get the timestep size between printing the solution data to
        the text output file.

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
                log_warning_message(msg)
                return 0.0

    @timestep_for_printing_solution.setter
    def timestep_for_printing_solution(self, delta_time: float):
        """Set the timestep size between printing."""
        """Set the timestep size between printing the solution data to
        the text output file.

        Parameters
        ----------
            delta_time: double, default = 1/100 of the simulation duration
                timestep size between printing solution data [sec]

        """
        if delta_time > 0.0e0:
            self.setkeyword(key="DELT", value=delta_time)
        else:
            msg = [Color.PURPLE, "solution printing timestep size must > 0.", Color.END]
            log_error_message(msg)

    def adaptive_solution_saving(
        self,
        mode: bool,
        value_change: Union[float, None] = None,
        target: Union[str, None] = None,
        steps: Union[int, None] = None,
    ):
        """Set up adaptive solution data saving."""
        """
        Set up adaptive solution data saving.

        Parameters
        ----------
            mode: boolean
                switch adaptive solution saving ON/OFF
            value_change: double, optional
                change in solution variable value between saving additional
                solution data
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
                log_error_message(msg)
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
                    "a reference variable is required",
                    "for value-change adaptive saving.",
                    Color.END,
                ]
                log_error_message(msg)
            elif not isinstance(target, str):
                # not given as string
                msg = [
                    Color.PURPLE,
                    "a reference variable is assigned as a string, e.g., 'OH'.",
                    Color.END,
                ]
                log_error_message(msg)
            elif value_change <= 0.0:
                # non-positive change value
                msg = [
                    Color.PURPLE,
                    "the value change per adaptive solution saving must > 0.",
                    Color.END,
                ]
                log_error_message(msg)
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
            log_error_message(msg)

    def set_ignition_delay(
        self,
        method: str = "T_inflection",
        val: float = 0.0,
        target: str = "",
    ):
        """Set ignition detection criterion."""
        """
        Set ignition detection criterion.

        Parameters
        ----------
            method: string
                ignition definition/detection method
            val: double, optional
                temperature or temperature rise value associated with
                the ignition detection method specified
            target: string, optional
                target species symbol if the 'Species_peak' method is used

        """
        if isinstance(method, str):
            # ignition detection method assigned
            if method == "T_inflection":
                # use inflection points in the temperature profile
                self.setkeyword(key="TIFP", value=True)
                # set flag
                self.report_ignition_delay = True
            elif method == "T_rise":
                # use temperature rise
                if val <= 0.0:
                    msg = [
                        Color.PURPLE,
                        "temperature rise value must > 0.",
                        Color.END,
                    ]
                    log_error_message(msg)
                else:
                    self.setkeyword(key="DTIGN", value=val)
                    # set flag
                    self.report_ignition_delay = True
            elif method == "T_ignition":
                # use temperature value
                if val <= 0.0:
                    msg = [
                        Color.PURPLE,
                        "ignition temperature value must > 0.",
                        Color.END,
                    ]
                    log_error_message(msg)
                else:
                    self.setkeyword(key="TLIM", value=val)
                    # set flag
                    self.report_ignition_delay = True
            elif method == "Species_peak":
                # use species peak location
                if target not in self._specieslist:
                    # no species given
                    msg = [
                        Color.PURPLE,
                        "target species is assigned as a string, e.g., 'OH'.",
                        Color.END,
                    ]
                    log_error_message(msg)
                else:
                    self.setkeyword(key="KLIM", value=target)
                    # set flag
                    self.report_ignition_delay = True
            elif method == "Thermicity_peak":
                # use  peak total thermicity location
                self.setkeyword(key="SLIM", value=True)
                # set flag
                self.report_ignition_delay = True
            else:
                # incorrect ignition detection method given
                msg = [
                    Color.PURPLE,
                    "ignition definition",
                    method,
                    "is not recognized.",
                    Color.END,
                ]
                log_error_message(msg)
                show_ignition_definitions()
        else:
            # use error
            show_ignition_definitions()

    def stop_after_ignition(self):
        """Set the option to stop the simulation after ignition is detected."""
        # stop the simulation after ignition is detected
        self.setkeyword(key="IGN_STOP", value=True)

    def get_ignition_delay(self) -> float:
        """Get the predicted ignition delay time."""
        """Get the predicted ignition delay time from
        the transient reactor simulation.

        Returns
        -------
            ignition_delay_time: double
                ignition delay time [msec] or [CA]

        """
        # initialization
        ignitiondelaytime = c_double(0.0e0)
        # check run status
        status = self.getrunstatus(mode="silent")
        if not self.report_ignition_delay:
            msg = [
                Color.YELLOW,
                "ignition delay has not been set.\n",
                Color.SPACEx6,
                "please set the ignition delay definition first.",
                Color.END,
            ]
            log_info_message(msg)
            return ignitiondelaytime.value
        elif status == -100:
            msg = [
                Color.YELLOW,
                "simulation has yet to be run.\n",
                Color.SPACEx6,
                "please run the reactor simulation first.",
                Color.END,
            ]
            log_info_message(msg)
            return ignitiondelaytime.value
        elif status != 0:
            msg = [
                Color.YELLOW,
                "simulation was failed.\n",
                Color.SPACEx6,
                "please correct the error(s) and rerun the reactor simulation.",
                Color.END,
            ]
            log_info_message(msg)
            return ignitiondelaytime.value

        # get the ignition delay time
        # (batch reactor model [sec], engine model [CA])
        ierr = chemkin_wrapper.chemkin.KINAll0D_GetIgnitionDelay(ignitiondelaytime)
        if ierr != 0:
            msg = [
                Color.MAGENTA,
                "potential bad ignition delay time value,",
                "error code =",
                str(ierr),
                "\n",
                Color.SPACEx6,
                "please check the text output and",
                "revisit the reactor/solver settings.",
                Color.END,
            ]
            log_warning_message(msg)
            return ignitiondelaytime.value
        # check reactor model
        if self._reactortype.value == self.ReactorTypes.get("Batch", 1):
            # check ignition delay time value
            if ignitiondelaytime.value <= 0.0:
                msg = [
                    Color.MAGENTA,
                    "potential bad ignition delay time value.\n",
                    Color.SPACEx6,
                    "please check the text output and",
                    "revisit the reactor/solver settings.",
                    Color.END,
                ]
                log_warning_message(msg)
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
                    "please check the text output and",
                    "revisit the reactor/solver settings.",
                    Color.END,
                ]
                log_warning_message(msg)
                return ignitiondelaytime.value
            else:
                # return ignition distance [cm]
                return ignitiondelaytime.value
        else:
            return ignitiondelaytime.value

    def refresh_starting_mixture(self, new_mixture: Mixture):
        """Reset the initial gas-phase mixture condition."""
        """
        Reset the initial gas-phase mixture condition inside the batch reactor.

        Parameters
        ----------
            new_mixture: Mixture object
                new mixture as the initial reactor gas condition
        """
        if isinstance(new_mixture, Mixture):
            self.pressure = new_mixture.pressure
            self.temperature = new_mixture.temperature
            self.volume = new_mixture.volume
            self.y = new_mixture.y
        else:
            # error condition
            msg = [
                Color.PURPLE,
                "the parameter must be a Mixture object.",
                Color.END,
            ]
            error_and_exit(msg)

    def set_volume_profile(
        self, x: npt.NDArray[np.double], vol: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor volume profile."""
        """
        Specify reactor volume profile.

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
            log_error_message(msg)
            return 10
        else:
            keyword = "VPRO"
            ierr = self.setprofile(key=keyword, x=x, y=vol)
            return ierr

    def set_pressure_profile(
        self, x: npt.NDArray[np.double], pres: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor pressure profile."""
        """
        Specify reactor pressure profile.

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
            log_error_message(msg)
            return 10
        else:
            keyword = "PPRO"
            ierr = self.setprofile(key=keyword, x=x, y=pres)
            return ierr

    def set_surfacearea_profile(
        self, x: npt.NDArray[np.double], area: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor reactive surface area profile."""
        """
        Specify reactor reactive surface area profile.

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
        ierr = self.setprofile(key=keyword, x=x, y=area)
        return ierr

    def set_reactortype_keywords(self):
        """Set reactor type keywords under the Full-Keywords mode."""
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
        """Set reactor condition keywords."""
        self.setkeyword(key="PRES", value=self._pressure.value / P_ATM)
        self.setkeyword(key="TEMP", value=self._temperature.value)
        self.setkeyword(key="TIME", value=self._endtime.value)
        # initial mole fraction
        _, species_lines = self.createspeciesinputlines(
            self._solvertype.value, threshold=1.0e-12, molefrac=self.reactormixture.x
        )
        for line in species_lines:
            self.setkeyword(key=line, value=True)

    def stop_output(self):
        """Stop the reactor model from writing text output to file."""
        self.suppress_output = True

    def validate_inputs(self) -> int:
        """Check the required inputs."""
        """Check the required inputs before running the reactor simulation.

        Returns
        -------
            error code: integer

        """
        ierr = 0
        # required inputs:
        if self._numb_requiredinput <= 0:
            # no required input
            return ierr
        else:
            if len(self._inputcheck) < self._numb_requiredinput:
                msg = [Color.PURPLE, "some required inputs are missing.", Color.END]
                log_error_message(msg)
            # verify required inputs one by one
            for k in self._requiredlist:
                if k not in self._inputcheck:
                    ierr += 1
                    msg = [Color.PURPLE, "missing required input:", k, Color.END]
                    log_error_message(msg)

            return ierr

    def __process_keywords_withfullinputs(self) -> int:
        """Process input keywords."""
        """Process input keywords for the batch reactor model
        under the Full-Keyword mode.

        Returns
        -------
            Error code: integer

        """
        ierr = 0
        set_verbose(True)
        # verify required inputs
        ierr = self.validate_inputs()
        if ierr != 0:
            msg = [Color.PURPLE, "missing required input keywords.", Color.END]
            log_error_message(msg)
            return ierr
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
                ierrc = chemkin_wrapper.chemkin.KINAll0D_SetProfilePoints(ipoints)
                # setup reactor model working arrays
                if ierrc == 0:
                    ierrc = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                        self._mylout, self._chemset_index
                    )
                ierr += ierrc
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "profile data generation error, error code =",
                str(ierr),
                Color.END,
            ]
            log_error_message(msg)
            return ierr
        # prepare initial conditions
        # initial mass fraction
        y_init = self.reactormixture.y
        if self.has_surface_chemistry:
            # set surface species fractions
            site_init, bulk_init = self.set_init_surface_coverage()
        else:
            # surface sites
            site_init = np.zeros(1, dtype=np.double)
            # bulk activities (not applicable)
            bulk_init = np.zeros_like(site_init, dtype=np.double)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("Batch"):
            ierrc = chemkin_wrapper.chemkin.KINAll0D_SetupBatchInputs(
                self._chemset_index,
                self._endtime,
                self._temperature,
                self._pressure,
                self._volume,
                self._heat_loss_rate,
                self._reactivearea,
                y_init,
                site_init,
                bulk_init,
            )
            ierr += ierrc
            if ierrc != 0:
                log_error_message(
                    [Color.PURPLE, "failed to set up basic reactor keywords", Color.END]
                )
                return ierrc

        # set reactor type
        self.set_reactortype_keywords()
        # reactor initial/estimated condition
        self.set_reactorcondition_keywords()
        if ierr == 0 and self.has_surface_chemistry:
            # set surface related keywords
            self.set_surface_chemistry_keywords()
        if ierr == 0 and self._numbprofiles > 0:
            # get keyword lines of all profiles
            err_profile, _, prof_lines = self.createprofileinputlines()
            ierr += err_profile
            if err_profile == 0:
                # set the profile keywords
                for pkey in prof_lines:
                    self.setkeyword(key=pkey, value=True)
            else:
                msg = [
                    Color.PURPLE,
                    "failed to set up profile keywords,",
                    "error code =",
                    str(err_profile),
                    Color.END,
                ]
                log_error_message(msg)
                return err_profile
        # solve integrated heat release rate due to chemical reactions
        if self.EnergyTypes.get("ENERGY") == self._energytype.value:
            self.setkeyword(key="QRGEQ", value=True)
        # add the END keyword
        self.setkeyword(key="END", value=True)
        # create input lines from additional user-specified keywords
        ierr, nlines = self.createkeywordinputlines()
        if ierr == 0:
            if verbose():
                msg = [
                    Color.YELLOW,
                    str(nlines),
                    "input lines are added.",
                    Color.END,
                ]
                log_info_message(msg)
        else:
            msg = [
                Color.PURPLE,
                "failed to create additional keyword lines.",
                Color.END,
            ]
            log_error_message(msg)

        return ierr

    def __run_model_withfullinputs(self) -> int:
        """Run the batch reactor model."""
        """Run the batch reactor model after the keywords are processed
        under the Full-Keyword mode. All keywords must be assigned.

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
        ierr = chemkin_wrapper.chemkin.KINAll0D_CalculateInput(
            self._mylout, self._chemset_index, longline, nlines, linelength
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to set up reactor keywords in Full mode,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            log_error_message(msg)

        return ierr

    def __process_keywords(self) -> int:
        """Process input keywords."""
        """
        Process input keywords for the batch reactor model.

        Returns
        -------
            Error code: integer

        """
        ierr = 0
        ierrc = 0
        err_key = 0
        err_inputs = 0
        # set_verbose(True)
        # verify required inputs
        ierr = self.validate_inputs()
        if ierr != 0:
            msg = [Color.PURPLE, "missing required input keywords.", Color.END]
            log_error_message(msg)
            return ierr
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
                ierrc = chemkin_wrapper.chemkin.KINAll0D_SetProfilePoints(ipoints)
                # setup reactor model working arrays
                if ierrc == 0:
                    ierrc = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                        self._mylout, self._chemset_index
                    )
                ierr += ierrc
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "profile data generation error, error code =",
                str(ierr),
                Color.END,
            ]
            log_error_message(msg)
            return ierr
        # prepare initial conditions
        # initial mass fraction
        y_init = self.reactormixture.y
        if self.has_surface_chemistry:
            # set surface species fractions
            site_init, bulk_init = self.set_init_surface_coverage()
        else:
            # surface sites
            site_init = np.zeros(1, dtype=np.double)
            # bulk activities (not applicable)
            bulk_init = np.zeros_like(site_init, dtype=np.double)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("Batch"):
            ierrc = chemkin_wrapper.chemkin.KINAll0D_SetupBatchInputs(
                self._chemset_index,
                self._endtime,
                self._temperature,
                self._pressure,
                self._volume,
                self._heat_loss_rate,
                self._reactivearea,
                y_init,
                site_init,
                bulk_init,
            )
            ierr += ierrc
            if ierrc != 0:
                msg = [
                    Color.PURPLE,
                    "failed to set up basic reactor keywords,",
                    "error code =",
                    str(ierrc),
                    Color.END,
                ]
                log_error_message(msg)
                return ierrc
            # heat transfer (use additional keywords)
            # solver parameters (use additional keywords)
            # output controls (use additional keywords)
            # ROP (use additional keywords)
            # sensitivity (use additional keywords)
            # ignition delay (use additional keywords)
            # solve integrated heat release rate due to chemical reactions
            if self.EnergyTypes.get("ENERGY") == self._energytype.value:
                ierrc = chemkin_wrapper.chemkin.KINAll0D_IntegrateHeatRelease()
                ierr += ierrc
                if ierrc != 0:
                    msg = [
                        Color.PURPLE,
                        "failed to set up heat release keyword,",
                        "error code =",
                        str(ierrc),
                        Color.END,
                    ]
                    log_error_message(msg)
                    return ierrc

        if ierr == 0 and self.has_surface_chemistry:
            # set surface related keywords
            self.set_surface_chemistry_keywords()

        if ierr == 0 and self._numbprofiles > 0:
            for p in self._profiles_list:
                key = bytes(p.profilekey, "utf-8")
                npoints = c_int(p.size)
                x = p.pos
                y = p.value
                err_profile = chemkin_wrapper.chemkin.KINAll0D_SetProfileParameter(
                    key, npoints, x, y
                )
                ierr += err_profile
            if err_profile != 0:
                msg = [
                    Color.PURPLE,
                    "failed to set up profile keywords,",
                    "error code =",
                    str(err_profile),
                    Color.END,
                ]
                log_error_message(msg)
                return err_profile
        if ierr == 0:
            # set additional keywords
            # create input lines from additional user-specified keywords
            err_inputs, nlines = self.createkeywordinputlines()
            if err_inputs == 0:
                # process additional keywords in _keyword_index and _keyword_lines
                for s in self._keyword_lines:
                    # convert string to byte
                    line = bytes(s, "utf-8")
                    # set additional keyword one by one
                    err_key = chemkin_wrapper.chemkin.KINAll0D_SetUserKeyword(line)
                    err_inputs += err_key
                if err_inputs == 0:
                    if verbose():
                        msg = [
                            Color.YELLOW,
                            str(nlines),
                            "additional input lines are added.",
                            Color.END,
                        ]
                        log_info_message(msg)
                else:
                    msg = [
                        Color.PURPLE,
                        "failed to create additional input lines,",
                        "error code =",
                        str(err_inputs),
                        Color.END,
                    ]
                    log_error_message(msg)
            else:
                msg = [
                    Color.PURPLE,
                    "failed to process additional keywords, error code =",
                    str(err_inputs),
                    Color.END,
                ]
                log_error_message(msg)
        #
        ierr = ierr + err_inputs + err_key

        return ierr

    def __run_model(self) -> int:
        """Run the batch reactor model."""
        """
        Run the batch reactor model after the keywords are processed.

        Returns
        -------
            error code: integer

        """
        # run the simulation without keyword inputs
        ierr = chemkin_wrapper.chemkin.KINAll0D_Calculate(self._chemset_index)
        return ierr

    def run(self) -> int:
        """Perform steps to run the batch reactor model."""
        """Perform simulation steps to run a Chemkin batch reactor model.

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
        log_info_message(msg)
        if not check_chemistryset(self._chemset_index.value):
            # Chemkin-CFD-API is not initialized: reinitialize Chemkin-CFD-API
            msg = [Color.YELLOW, "initializing Chemkin ...", Color.END]
            log_info_message(msg)
            ret_val = chemkin_wrapper.chemkin.KINInitialize(
                self._chemset_index, c_int(0)
            )
            if ret_val != 0:
                msg = [
                    Color.RED,
                    "Chemkin-CFD-API initialization failed;",
                    "code =",
                    str(ret_val),
                    Color.END,
                ]
                log_critical_message(msg)
                exit()
            else:
                chemistryset_initialized(self._chemset_index.value)

        # output initialization
        log_info_message([Color.YELLOW, "clearing output ...", Color.END])

        # keyword processing
        msg = [
            Color.YELLOW,
            "processing and generating keyword inputs ...",
            Color.END,
        ]
        log_info_message(msg)
        if Keyword.no_fullkeyword:
            # use API calls
            ret_val = (
                self.__process_keywords()
            )  # each reactor model subclass to perform its own keyword processing
        else:
            # use full keywords
            ret_val = self.__process_keywords_withfullinputs()
        if ret_val != 0:
            msg = [
                Color.RED,
                "generating the keyword inputs,",
                "error code =",
                str(ret_val),
                Color.END,
            ]
            log_critical_message(msg)
            return ret_val
        log_info_message([Color.YELLOW, "Processing keywords complete", Color.END])

        # run reactor model
        msg = [Color.YELLOW, "running reactor simulation ...", Color.END]
        log_info_message(msg)
        # suppress text output to file
        if self.suppress_output:
            ierr = chemkin_wrapper.chemkin.KINAll0D_SuppressOutput()
            if ierr != 0:
                msg = [Color.YELLOW, "failed to turn off text output.", Color.END]
                log_info_message(msg)
        # reset run status
        self.setrunstatus(code=-100)
        #
        if Keyword.no_fullkeyword:
            # use API calls
            ret_val = self.__run_model()
        else:
            # use full keywords
            ret_val = self.__run_model_withfullinputs()
        # update run status
        self.setrunstatus(code=ret_val)
        msg = ["simulation completed,", "status =", str(ret_val), Color.END]
        if ret_val == 0:
            msg.insert(0, Color.GREEN)
            log_info_message(msg)
        else:
            msg.insert(0, Color.RED)
            log_critical_message(msg)

        return ret_val

    def get_solution_size(self) -> tuple[int, int]:
        """Get number of solution time points."""
        """Get the number of reactors and the number of solution points.

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
            log_warning_message(msg)
            exit()
        elif status != 0:
            msg = [
                Color.PURPLE,
                "simulation was failed.\n",
                Color.SPACEx6,
                "please correct the error(s) and rerun the reactor simulation.",
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # number of reactor
        nreac = c_int(0)
        # number of time points in the solution
        npoints = c_int(0)
        # get solution size of the batch reactor
        ierr = chemkin_wrapper.chemkin.KINAll0D_GetSolnResponseSize(nreac, npoints)
        nreactors = nreac.value
        if ierr == 0 and nreactors == self._nreactors:
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
                str(ierr),
                Color.END,
            ]
            log_error_message(msg)
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
            log_error_message(msg)
            exit()

    def process_solution(self):
        """Post-process solution to extract the raw solution variable data."""
        # check existing raw data
        if self.getrawsolutionstatus():
            msg = [
                Color.YELLOW,
                "the solution has been processed before,",
                "any existing solution data will be deleted from the memory.",
                Color.END,
            ]
            log_info_message(msg)

        msg = [Color.YELLOW, "post-processing raw solution data ...", Color.END]
        log_info_message(msg)
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
            error_and_exit(msg)
        else:
            self._numbsolutionpoints = npoints
        # create arrays to hold the raw solution data
        time = np.zeros(self._numbsolutionpoints, dtype=np.double)
        pres = np.zeros_like(time, dtype=np.double)
        temp = np.zeros_like(time, dtype=np.double)
        vol = np.zeros_like(time, dtype=np.double)
        # create a species mass fraction array to
        # hold the solution species fraction profiles
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
        ierr = chemkin_wrapper.chemkin.KINAll0D_GetGasSolnResponse(
            icreac, icnpts, icnspec, time, temp, pres, vol, frac
        )
        if ierr != 0:
            msg = [
                Color.RED,
                "failed to fetch the raw solution data from memory,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            critical_and_exit(msg)
        # store the raw solution data in a dictionary
        # time
        self._solution_rawarray["time"] = copy.deepcopy(time)
        # for PFR the independent variable is distance [cm]
        if self._reactortype.value == self.ReactorTypes.get("PFR", 3):
            self._solution_rawarray["distance"] = copy.deepcopy(time)
        # temperature
        self._solution_rawarray["temperature"] = copy.deepcopy(temp)
        # pressure
        self._solution_rawarray["pressure"] = copy.deepcopy(pres)
        # volume
        self._solution_rawarray["volume"] = copy.deepcopy(vol)
        # species mass fractions
        self.parsespeciessolutiondata(frac)
        # if the energy equation is solved, get gas-phase heat release rate
        if self._energytype.value == self.EnergyTypes.get("ENERGY", 1):
            gas_heatrelease = np.zeros_like(time, dtype=np.double)
            gas_heatrealease_rate = np.zeros_like(time, dtype=np.double)
            ierr = chemkin_wrapper.chemkin.KINAll0D_GetGasHeatReleaseRate(
                icreac, icnpts, gas_heatrelease, gas_heatrealease_rate
            )
            if ierr != 0:
                msg = [
                    Color.RED,
                    "failed to fetch the raw solution heat release rate from memory,",
                    "error code =",
                    str(ierr),
                    Color.END,
                ]
                critical_and_exit(msg)
            #
            # accumulated gas-phase chemistry heat release [erg]
            self._solution_rawarray["gasheatrelease"] = copy.deepcopy(gas_heatrelease)
            # gas-phase chemistry heat rrelease rate [erg/sec]
            self._solution_rawarray["gashrr"] = copy.deepcopy(gas_heatrealease_rate)
            del gas_heatrelease, gas_heatrealease_rate
        # create solution mixture
        ierr = self.create_solution_mixtures(frac)
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "forming solution mixtures",
                "error code =",
                str(ierr),
                Color.END,
            ]
            error_and_exit(msg)
        # calculate total termicity [1/sec]
        thermicity = np.zeros_like(time, dtype=np.double)
        for i, m in enumerate(self._solution_mixturearray):
            thermicity[i] = m.thermicity()
        # total thermicity
        self._solution_rawarray["thermicity"] = copy.deepcopy(thermicity)
        # get surface coverage
        if self.has_surface_chemistry:
            self.process_surface_solution(nreac)
        # set up single point solution parameters
        if self.report_ignition_delay:
            delay = self.get_ignition_delay()
            if delay > 0.0e0:
                self._solution_parameters["ignition_delay"] = delay
        # clean up
        del time, pres, temp, vol, thermicity, frac

    def process_surface_solution(self, nreac: int):
        """Post-process surface solution."""
        """Post-process surface solution to extract the raw solution variable data
        package the steady-state solution into a mixture object.

        Parameters
        ----------
            nreac: integer
                number of reactors

        """
        #
        icreac = c_int(nreac)
        icnpts = c_int(self._numbsolutionpoints)
        # get material information
        mat_names = self.get_material_names()
        # get the site and the bulk species sizes of the mechanism
        max_sites = self.get_total_site_species()
        max_bulks = self.get_total_bulk_species()
        # set up holding arrays
        # create arrays to hold the raw solution data
        time = np.zeros(self._numbsolutionpoints, dtype=np.double)
        # site fractions
        if max_sites > 0:
            s_frac = np.zeros(
                (max_sites, self._numbsolutionpoints),
                dtype=np.double,
                order="F",
            )
        else:
            s_frac = np.zeros((1, self._numbsolutionpoints), dtype=np.double, order="F")
        # bulk activities
        if max_bulks > 0:
            b_act = np.zeros(
                (max_bulks, self._numbsolutionpoints),
                dtype=np.double,
                order="F",
            )
        else:
            b_act = np.zeros((1, self._numbsolutionpoints), dtype=np.double, order="F")

        ierr_s = chemkin_wrapper.chemkin.KINAll0D_GetSurfaceSolnResponse(
            icreac, icnpts, c_int(max_sites), c_int(max_bulks), time, s_frac, b_act
        )
        if ierr_s != 0:
            msg = [
                Color.MAGENTA,
                "failed to fetch the raw surface solution data from memory,",
                "error code =",
                str(ierr_s),
                Color.END,
            ]
            log_error_message(msg)
        # get surface material temperature [K]
        n_material = self.surface_chemistry.number_materials
        surf_temp = np.zeros(
            (n_material, self._numbsolutionpoints),
            dtype=np.double,
            order="F",
        )
        ierr_t = chemkin_wrapper.chemkin.KINAll0D_GetSurfaceTemperature(
            icreac, icnpts, c_int(n_material), surf_temp
        )
        # store surface coverage in the solution mixtures
        # loop over all time points
        soln_a = np.zeros(max(max_sites, max_bulks), dtype=np.double)
        for i, m in enumerate(self._solution_mixturearray):
            # loop over all materials of the surface mechanism
            s_count = 0
            b_count = 0
            for n, mname in enumerate(mat_names):
                # set material temperature
                if ierr_t != 0:
                    # make material temperature equal to the gas temperature
                    s_temp = m.temperature
                else:
                    s_temp = surf_temp[n][i]
                m.surface_chemistry.set_material_temperature(mname, s_temp)
                # number of site species on this material
                n_sites = self.surface_chemistry.number_site_species(mname)
                # number of bulk species on this material
                n_bulks = self.surface_chemistry.number_bulk_species(mname)
                if ierr_s == 0:
                    if n_sites > 0:
                        soln_a[:] = 0.0
                        for k in range(n_sites):
                            j = k + s_count
                            soln_a[k] = s_frac[j][i]
                        m.surface_chemistry.set_site_frac(mname, soln_a[0:n_sites])
                    s_count += n_sites

                    if n_bulks > 0:
                        soln_a[:] = 0.0
                        for k in range(n_bulks):
                            j = k + b_count
                            soln_a[k] = b_act[j][i]
                        m.surface_chemistry.set_bulk_frac(mname, soln_a[0:n_bulks])
                    b_count += n_bulks
        # fill the raw solution array
        ss = np.zeros_like(time, dtype=np.double)
        if max_sites > 0:
            # get site symbols
            s_symbols = self.get_site_species_names()
            for snames in s_symbols:
                ss[:] = 0.0
                for k, symbol in enumerate(snames):
                    for i in range(self._numbsolutionpoints):
                        ss[i] = s_frac[k][i]
                    self._solution_rawarray[symbol] = copy.deepcopy(ss)
        if max_bulks > 0:
            # get site symbols
            b_symbols = self.get_bulk_species_names()
            for bnames in b_symbols:
                ss[:] = 0.0
                for k, symbol in enumerate(bnames):
                    for i in range(self._numbsolutionpoints):
                        ss[i] = b_act[k][i]
                    self._solution_rawarray[symbol] = copy.deepcopy(ss)
        # if the energy equation is solved, get the total surface heat release rate
        if self._energytype.value == self.EnergyTypes.get("ENERGY", 1):
            surf_heatrelease_rate = np.zeros_like(time, dtype=np.double)
            ierr = chemkin_wrapper.chemkin.KINAll0D_GetSurfaceHeatReleaseRate(
                icreac, icnpts, surf_heatrelease_rate
            )
            if ierr != 0:
                msg = [
                    Color.RED,
                    "failed to fetch the raw solution surface heat release rate",
                    "from memory,",
                    "error code =",
                    str(ierr),
                    Color.END,
                ]
                critical_and_exit(msg)
            # total surface chemistry heat release rate [erg/s]
            self._solution_rawarray["surfhrr"] = copy.deepcopy(surf_heatrelease_rate)

        del time, s_frac, b_act

    def get_solution_variable_profile(self, varname: str) -> npt.NDArray[np.double]:
        """Get the profile of the solution variable specified."""
        """
        Get the profile of the solution variable specified.

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
            log_info_message(msg)
            return np.zeros(0, dtype=np.double)
        # check variable name
        vname = varname.rstrip()
        if vname.lower() in self._solution_tags:
            # is a property variable?
            vname = vname.lower()
        else:
            is_species = False
            if vname in self._specieslist:
                # a gas species
                is_species = True
            if self.has_surface_chemistry:
                for names in self.get_site_species_names():
                    if vname in names:
                        # a site species
                        is_species = True
                        break

                for names in self.get_bulk_species_names():
                    if vname in names:
                        # a bulk species
                        is_species = True
                        break
            if not is_species:
                # not a species?
                msg = [
                    Color.PURPLE,
                    "variable name",
                    vname,
                    "is not part of the known solution variable.\n",
                    Color.SPACEx6,
                    "and has to be derived from other variable(s).",
                    Color.END,
                ]
                log_error_message(msg)
                exit()

        # create variable arrays to hold the solution profile
        var = np.zeros(self._numbsolutionpoints, dtype=np.double)
        # get variable profile from the raw solution data
        var = self._solution_rawarray.get(vname)
        return var

    def create_solution_mixtures(self, specfrac: npt.NDArray[np.double]) -> int:
        """Create a list of Mixtures."""
        """Create a list of Mixtures that represent the gas inside
        the reactor at a solution point.

        Parameters
        ----------
            specfrac: 2D double array, dimensions
            = [number_species, numb_solution_point]
                species fractions of all time points [species_fraction, time_point]

        Returns
        -------
            ierror: integer
                 error code

        """
        if not self.getrawsolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'getsolution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            log_info_message(msg)
            return 1
        # create a temporary Mixture object to hold the mixture properties
        # at current solution point
        smixture = copy.deepcopy(self.reactormixture)
        # create variable arrays to hold the solution profile
        species = []
        # create a species fraction array to
        # hold the solution species fraction profiles
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
                smixture.y = frac
            else:
                # mole fractions
                smixture.x = frac
            # add to the solution mixture list
            self._solution_mixturearray.append(copy.deepcopy(smixture))
        # clean up
        species.clear()
        del pres, temp, vol, frac, species, smixture
        return 0

    def get_solution_mixture(self, time: float) -> Mixture:
        """Get the mixture representing the solution at the given time."""
        """Get the mixture representing the solution state inside
        the reactor at the given time.

        Parameters
        ----------
            time: double
                time point value [sec]

        Returns
        -------
            mixturetarget: Mixture object
                a Mixture representing the gas properties in the reactor
                at the specific time

        """
        # check status
        if not self.getmixturesolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'process_solution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            log_info_message(msg)
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
        """Get the mixture representing the solution at the given grid index."""
        """Get the mixture representing the solution state inside
        the reactor at the given solution point index.

        Parameters
        ----------
            solution_index: integer
                0-base solution time point index

        Returns
        -------
            mixturetarget: Mixture object
                a Mixture representing the gas properties in
                the reactor at the specific time

        """
        # check status
        if not self.getmixturesolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'process_solution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            log_info_message(msg)
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
            log_error_message(msg)
            exit()
        # get the mixture
        mixturetarget = copy.deepcopy(self._solution_mixturearray[solution_index])
        return mixturetarget

    def get_last_solution_mixture(self) -> Mixture:
        """Get the solution mixture at the last time point."""
        """
        Get the mixture representing the solution state inside the reactor
        at the last solution point index.

        Returns
        -------
            mixturetarget: Mixture object
                a Mixture representing the gas properties in the reactor
                at the last time/position point
        """
        # check status
        if not self.getmixturesolutionstatus():
            msg = [
                Color.YELLOW,
                "please use 'process_solution' method",
                "to post-process the raw solution data first.",
                Color.END,
            ]
            log_info_message(msg)
            exit()
        # set index to the last solution point
        solution_index = self._numbsolutionpoints - 1
        # get the mixture
        mixturetarget = copy.deepcopy(self._solution_mixturearray[solution_index])
        return mixturetarget


class GivenPressureBatchReactorFixedTemperature(BatchReactors):
    """Chemkin 0-D transient closed homogeneous reactor model."""

    """Chemkin 0-D transient closed homogeneous reactor model
    with given reactor pressure (CONP) and reactor temperature (TGIV).
    """

    def __init__(self, reactor_condition: Stream, label: str = "CONPT"):
        """Initialize a Given Pressure Batch Reactor object with given temperature."""
        """Initialize a Given Pressure Batch Reactor object
        with given reactor temperature.

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties
                inside the batch reactor
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
        self._ninlets = np.zeros(1, dtype=np.int32)
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
        ierr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if ierr == 0:
            # setup reactor model working arrays
            ierr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._mylout, self._chemset_index
            )
            ierr *= 10
        if ierr != 0:
            msg = [
                Color.RED,
                "failed to initialize the reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            critical_and_exit(msg)
        # if full-keyword mode is turned ON
        if not Keyword.no_fullkeyword:
            # populate the reactor setup keywords
            self.set_reactortype_keywords()

    def set_temperature_profile(
        self, x: npt.NDArray[np.double], temp: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor temperature profile."""
        """
        Specify reactor temperature profile.

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
        ierr = self.setprofile(key=keyword, x=x, y=temp)
        return ierr


class GivenPressureBatchReactorEnergyConservation(BatchReactors):
    """Chemkin 0-D transient closed homogeneous reactor model."""

    """Chemkin 0-D transient closed homogeneous reactor model
    with given reactor pressure (CONP) and solving the energy equation (ENRG).
    """

    def __init__(self, reactor_condition: Stream, label: str = "CONP"):
        """Initialize a Given Pressure Batch Reactor object with Energy Equation."""
        """Initialize a Given Pressure Batch Reactor object that
        solves the Energy Equation.

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties
                inside the batch reactor
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
        self._ninlets = np.zeros(1, dtype=np.int32)
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
        # raw solution data structure
        # GasHRR: heat release rate [erg/sec] due to gas-phase chemsitry
        # GasHeatRelease: accumulated heat release [erg] due to gas-phase chemistry
        self._solution_tags.append("gashrr")
        self._solution_tags.append("gasheatrelease")
        if self.has_surface_chemistry:
            self._solution_tags.append("surfhrr")
        # set up basic batch reactor parameters
        ierr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if ierr == 0:
            # setup reactor model working arrays
            ierr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._mylout, self._chemset_index
            )
            ierr *= 10
        if ierr != 0:
            msg = [
                Color.RED,
                "failed to initialize the batch reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            critical_and_exit(msg)
        # if full-keyword mode is turned ON
        if not Keyword.no_fullkeyword:
            # populate the reactor setup keywords
            self.set_reactortype_keywords()

    @property
    def heat_loss_rate(self) -> float:
        """Get heat loss rate."""
        """
        Get heat loss rate from the reactor to the surroundings.

        Returns
        -------
            heat_loss_rate: double
                heat loss rate [cal/sec]

        """
        return self._heat_loss_rate.value

    @heat_loss_rate.setter
    def heat_loss_rate(self, value: float):
        """Set the heat loss rate."""
        """
        Set the heat loss rate from the reactor to the surroundings.

        Parameters
        ----------
            value: double, default = 0.0
                heat loss rate [cal/sec]

        """
        self._heat_loss_rate = c_double(value)
        if not Keyword.no_fullkeyword:
            self.setkeyword(key="QLOS", value=value)

    @property
    def heat_transfer_coefficient(self) -> float:
        """Get heat transfer coefficient."""
        """
        Get heat transfer coefficient between the reactor and the surroundings.

        Returns
        -------
            heat_transfer_coefficient: double
                heat transfer coefficient [cal/cm2-K-sec]

        """
        return self._heat_transfer_coefficient

    @heat_transfer_coefficient.setter
    def heat_transfer_coefficient(self, value: float = 0.0e0):
        """Set heat transfer coefficient."""
        """
        Set heat transfer coefficient between the reactor and the surroundings.

        Parameters
        ----------
            value: double, optional, default = 0.0
                heat transfer coefficient [cal/cm2-K-sec]

        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer coefficient must > 0.", Color.END]
            error_and_exit(msg)
        else:
            self._heat_transfer_coefficient = value
            # set the corresponding keyword
            self.setkeyword(key="HTC", value=value)

    @property
    def ambient_temperature(self) -> float:
        """Get ambient temperature."""
        """
        Get ambient temperature.

        Returns
        -------
            temperature: double
                ambient temperature [K]

        """
        return self._ambient_temperature

    @ambient_temperature.setter
    def ambient_temperature(self, value: float = 0.0e0):
        """Set ambient temperature."""
        """
        Set ambient temperature.

        Parameters
        ----------
            value: double, optional, default = 300.0
                ambient temperature [K]

        """
        if value <= 0.0e0:
            msg = [Color.PURPLE, "ambient temperature must > 0.", Color.END]
            error_and_exit(msg)
        else:
            self._ambient_temperature = value
            # set the corresponding keyword
            self.setkeyword(key="TAMB", value=value)

    @property
    def heat_transfer_area(self) -> float:
        """Get heat transfer area."""
        """
        Get heat transfer area between the reactor and the surroundings.

        Returns
        -------
            area: double, optional, default = 0.0
                heat transfer area [cm2]

        """
        return self._heat_transfer_area

    @heat_transfer_area.setter
    def heat_transfer_area(self, value: float = 0.0e0):
        """Set heat transfer area."""
        """
        Set heat transfer area between the reactor and the surroundings."

        Parameters
        ----------
            value: double, optional, default = 0.0
                heat transfer area [cm2]

        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer area must >= 0.", Color.END]
            error_and_exit(msg)
        else:
            self._heat_transfer_area = value
            # set the corresponding keyword
            self.setkeyword(key="AREAQ", value=value)

    def set_heat_transfer_area_profile(
        self, x: npt.NDArray[np.double], area: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor heat transfer area profile."""
        """
        Specify reactor heat transfer area profile.

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
            log_error_message(msg)
            return 10
        else:
            keyword = "AEXT"
            ierr = self.setprofile(key=keyword, x=x, y=area)
            return ierr

    def set_heat_loss_profile(
        self, x: npt.NDArray[np.double], qloss: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor heat loss rate profile."""
        """
        Specify reactor heat loss rate profile.

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            qloss: 1D double array
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
            log_error_message(msg)
            return 10
        else:
            keyword = "QPRO"
            ierr = self.setprofile(key=keyword, x=x, y=qloss)
            return ierr


class GivenVolumeBatchReactorFixedTemperature(BatchReactors):
    """Chemkin 0-D transient closed homogeneous reactor model."""

    """Chemkin 0-D transient closed homogeneous reactor model
    with given reactor volume (CONV) and reactor temperature (TGIV).
    """

    def __init__(self, reactor_condition: Stream, label: str = "CONVT"):
        """Initialize a Given Volume Batch Reactor object with fixed temperature."""
        """Initialize a Given Volume Batch Reactor object with
        given reactor temperature.

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties
                inside the batch reactor
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
        self._ninlets = np.zeros(1, dtype=np.int32)
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
        ierr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if ierr == 0:
            # setup reactor model working arrays
            ierr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._mylout, self._chemset_index
            )
            ierr *= 10
        if ierr != 0:
            msg = [
                Color.RED,
                "failed to initialize the reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            critical_and_exit(msg)
        # if full-keyword mode is turned ON
        if not Keyword.no_fullkeyword:
            # populate the reactor setup keywords
            self.set_reactortype_keywords()

    def set_temperature_profile(
        self, x: npt.NDArray[np.double], temp: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor temperature profile."""
        """
        Specify reactor temperature profile.

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
        ierr = self.setprofile(key=keyword, x=x, y=temp)
        return ierr


class GivenVolumeBatchReactorEnergyConservation(BatchReactors):
    """Chemkin 0-D transient closed homogeneous reactor model."""

    """Chemkin 0-D transient closed homogeneous reactor model
    with given reactor volume (CONV) and
    solving the energy equation (ENRG).
    """

    def __init__(self, reactor_condition: Stream, label: str = "CONV"):
        """Initialize a Given Volume Batch Reactor object with energy equation."""
        """Initialize a Given Volume Batch Reactor object that
        solves the Energy Equation.

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties
                inside the batch reactor
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
        self._ninlets = np.zeros(1, dtype=np.int32)
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
        # raw solution data structure
        # GasHRR: heat release rate [erg/sec] due to gas-phase chemsitry
        # GasHeatRelease: accumulated heat release [erg] due to gas-phase chemistry
        self._solution_tags.append("gashrr")
        self._solution_tags.append("gasheatrelease")
        if self.has_surface_chemistry:
            self._solution_tags.append("surfhrr")
        # set up basic batch reactor parameters
        ierr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if ierr == 0:
            # setup reactor model working arrays
            ierr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._mylout, self._chemset_index
            )
            ierr *= 10
        if ierr != 0:
            msg = [
                Color.RED,
                "failed to initialize the reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            critical_and_exit(msg)
        # if full-keyword mode is turned ON
        if not Keyword.no_fullkeyword:
            # populate the reactor setup keywords
            self.set_reactortype_keywords()

    @property
    def heat_loss_rate(self) -> float:
        """Get heat loss rate."""
        """
        Get heat loss rate from the reactor to the surroundings.

        Returns
        -------
            heat_loss_rate: double
                heat loss rate [cal/sec]

        """
        return self._heat_loss_rate.value

    @heat_loss_rate.setter
    def heat_loss_rate(self, value: float):
        """Set the heat loss rate."""
        """
        Set the heat loss rate from the reactor to the surroundings (required)."

        Parameters
        ----------
            value: double, default = 0.0
                heat loss rate [cal/sec]

        """
        self._heat_loss_rate = c_double(value)
        if not Keyword.no_fullkeyword:
            self.setkeyword(key="QLOS", value=value)

    @property
    def heat_transfer_coefficient(self) -> float:
        """Get heat transfer coefficient."""
        """
        Get heat transfer coefficient between the reactor and the surroundings.

        Returns
        -------
            heat_transfer_coefficient: double
                heat transfer coefficient [cal/cm2-K-sec]

        """
        return self._heat_transfer_coefficient

    @heat_transfer_coefficient.setter
    def heat_transfer_coefficient(self, value: float = 0.0e0):
        """Set heat transfer coefficient."""
        """
        Set heat transfer coefficient between the reactor and the surroundings.

        Parameters
        ----------
            value: double, optional, default = 0.0
                heat transfer coefficient [cal/cm2-K-sec]

        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer coefficient must > 0.", Color.END]
            error_and_exit(msg)
        else:
            self._heat_transfer_coefficient = value
            # set the corresponding keyword
            self.setkeyword(key="HTC", value=value)

    @property
    def ambient_temperature(self) -> float:
        """Get ambient temperature."""
        """
        Get ambient temperature.

        Returns
        -------
            temperature: double
                ambient temperature [K]

        """
        return self._ambient_temperature

    @ambient_temperature.setter
    def ambient_temperature(self, value: float = 0.0e0):
        """Set ambient temperature."""
        """
        Set ambient temperature.

        Parameters
        ----------
            value: double, optional, default = 300.0
                ambient temperature [K]

        """
        if value <= 0.0e0:
            msg = [Color.PURPLE, "ambient temperature must > 0.", Color.END]
            error_and_exit(msg)
        else:
            self._ambient_temperature = value
            # set the corresponding keyword
            self.setkeyword(key="TAMB", value=value)

    @property
    def heat_transfer_area(self) -> float:
        """Get heat transfer area."""
        """
        Get heat transfer area between the reactor and the surroundings.

        Returns
        -------
            area: double, optional, default = 0.0
                heat transfer area [cm2]

        """
        return self._heat_transfer_area

    @heat_transfer_area.setter
    def heat_transfer_area(self, value: float = 0.0e0):
        """Set heat transfer area."""
        """
        Set heat transfer area between the reactor and the surroundings.

        Parameters
        ----------
            value: double, optional, default = 0.0
                heat transfer area [cm2]

        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer area must > 0.", Color.END]
            error_and_exit(msg)
        else:
            self._heat_transfer_area = value
            # set the corresponding keyword
            self.setkeyword(key="AREAQ", value=value)

    def set_heat_transfer_area_profile(
        self, x: npt.NDArray[np.double], area: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor heat transfer area profile."""
        """
        Specify reactor heat transfer area profile.

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
            log_error_message(msg)
            return 10
        else:
            keyword = "AEXT"
            ierr = self.setprofile(key=keyword, x=x, y=area)
            return ierr

    def set_heat_loss_profile(
        self, x: npt.NDArray[np.double], qloss: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor heat loss rate profile."""
        """
        Specify reactor heat loss rate profile.

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            qloss: 1D double array
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
            log_error_message(msg)
            return 10
        else:
            keyword = "QPRO"
            ierr = self.setprofile(key=keyword, x=x, y=qloss)
            return ierr


def calculate_effective_activation_energy(
    initial_mixture: Mixture,
    duration: float = 0.01,
    perturb: float = 0.02,
    model: str = "CONV",
) -> tuple[float, float]:
    """Estimate the effective reduced activation energy."""
    """
    Estimate the effective reduced activation energy of
    a gas-phase combustion mechanism at the given mixture condition.

    Parameters
    ----------
        initial_mixture: Mixture object
            the condition of the reference gas mixture
        duration: double, optional, default = 0.01 [sec]
            reactor simulation time [sec]
        perturb: double, optional, default = 0.02 (2%)
            the size (fraction) of the perturbation
        model: string, {"CONP", "CONV"}
            the reactor model to be used to estimate the activation energy

    Returns
    -------
        reduced_act_energy: double
            reduced effective activation energy [-]
        activation_energy: double
            the effective activation energy [cal/mol]
    """
    if not isinstance(initial_mixture, Stream):
        msg = [
            Color.PURPLE,
            "the first parameter must be a Mixture object",
            Color.END,
        ]
        log_error_message(msg)
        exit()
    # check reactor volume
    if not initial_mixture.check_volume():
        # set volume to 1 [cm3] if it is not provided
        initial_mixture.volume = 1.0
    #
    if model.upper() == "CONV":
        rxtor = GivenVolumeBatchReactorEnergyConservation(initial_mixture, label=model)
    elif model.upper() == "CONP":
        rxtor = GivenPressureBatchReactorEnergyConservation(
            initial_mixture, label=model
        )
    else:
        msg = [
            Color.PURPLE,
            "invalid reactor model, use either",
            "'CONP' for constant pressure reactor model or",
            "'CONV' for constant volume reactor model.",
            Color.END,
        ]
        log_error_message(msg)
        exit()
    # initialization
    reduced_act_energy = 0.0e0
    activation_energy = 0.0e0
    # perturbation
    temp_base = initial_mixture.temperature
    dtemp = temp_base * perturb
    # set the run parameters
    temps: list[float] = []
    temps.append(temp_base + dtemp)
    temps.append(temp_base - dtemp)
    # reactor set up
    # initial reactor volume [cm3]
    rxtor.volume = initial_mixture.volume
    # simulation end time [sec]
    rxtor.time = duration
    # turn ON adaptive solution saving
    rxtor.adaptive_solution_saving(mode=True, value_change=100, target="TEMPERATURE")
    # set ignition delay
    rxtor.set_ignition_delay(method="T_inflection")
    # tolerances are given in tuple: (absolute tolerance, relative tolerance)
    rxtor.tolerances = (1.0e-10, 1.0e-8)
    # run the parameter study
    delaytime: list[float] = []
    for t in temps:
        # reset the initial reactor temperature
        rxtor.temperature = t
        # run the reactor model
        runstatus = rxtor.run()
        #
        if runstatus == 0:
            # get ignition delay time [msec]
            delay = rxtor.get_ignition_delay()
            delaytime.append(delay * 1.0e-3)
            print(f"ignition delay time = {delay} [msec]")
        else:
            # run failed
            msg = [
                Color.RED,
                ">>> Run failed. <<<\n",
                "for case: temperature =",
                str(t),
                Color.END,
            ]
            log_error_message(msg)
            exit()
    # verify the ignition delay times
    for case in range(len(delaytime)):
        time = delaytime[case]
        temp = temps[case]
        if time <= 0.0e0:
            # if get this, most likely the END time is too short
            msg = [
                Color.PURPLE,
                "ignition delay time not found, for case temperature =",
                str(temp),
                ", most likely the simulation time is too short.",
                Color.END,
            ]
            log_error_message(msg)
            exit()
    # estimate the effective reduced activation energy
    termn = np.log(delaytime[0]) - np.log(delaytime[1])
    termd = 1.0e0 / temps[0] - 1.0e0 / temps[1]
    reduced_act_energy = termn / termd / temp_base
    reduced_act_energy = max(reduced_act_energy, 0.0e0)
    # effect activation energy [cal/mol]
    activation_energy = reduced_act_energy * R_GAS_CAL * temp_base
    # clean up
    del rxtor
    return reduced_act_energy, activation_energy

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

"""Chemkin steady-state solver controlling parameters."""

from typing import Union

import numpy as np

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.reactormodel import Keyword


class SteadyStateSolver:
    """Common steady-state solver controlling parameters"""

    def __init__(self):
        # steady-state solver control parameter class
        # mostly just keyword processing
        # >>> steady-state search algorithm:
        # absolute tolerance for the steady-state solution
        self.SSabsolute_tolerance = 1.0e-9
        # relative tolerance for the steady-state solution
        self.SSrelative_tolerance = 1.0e-4
        # max number of iterations per steady state search
        self.SSmaxiteration = 100
        # number of steady-state searches before evaluating new Jacobian matrix
        self.SSJacobianage = 20
        # max number of calls to pseudo transient algorithm
        self.maxpseudotransient = 100
        # number of pseudo transient "steps" before calling the steady-state search algorithm
        self.numbinitialpseudosteps = 0
        # upper bound of the temperature value during iteration
        self.maxTbound = 5000.0  # [K]
        # floor value (lower bound) of the gas species mass fraction during iteration
        self.speciesfloor = -1.0e-14
        # reset negative gas species fraction to the given value in intermediate solution
        self.species_positive = 0.0e0
        # use legacy steady-state solver algorithm
        self.use_legacy_technique = False
        # use damping in search: 0 = OFF; 1 = ON
        self.SSdamping = 1
        # absolute perturbation for Jacobian evaluation
        self.absolute_perturbation = 0.0e0
        # relative perturbation for Jacobian evaluation
        self.relative_perturbation = 0.0e0
        # >>> pseudo transient (time stepping) algorithm:
        # absolute tolerance for the time stepping solution
        self.TRabsolute_tolerance = 1.0e-9
        # relative tolerance for the time stepping solution
        self.TRrelative_tolerance = 1.0e-4
        # max number of iterations per pseudo time step before cutting the time step size
        self.TRmaxiteration = 25
        # max number of pseudo time steps before increasing the time step size
        self.timestepsizeage = 25
        # minimum time step size allowed
        self.TRminstepsize = 1.0e-10  # [sec]
        # maximum time step size allowed
        self.TRmaxstepsize = 1.0e-2  # [sec]
        # time step size increasing factor
        self.TRupfactor = 2.0
        # time step size decreasing factor
        self.TRdownfactor = 2.2
        # number of pseudo time steps before evaluating new Jacobian matrix
        self.TRJacobianage = 20
        # initial stride and number of steps per pseudo time stepping call
        # for fixed-temperature solution
        self.TRstride_fixT = 1.0e-6  # [sec]
        self.TRnumbsteps_fixT = 100
        # for energy equation solution
        self.TRstride_ENRG = 1.0e-6  # [sec]
        self.TRnumbsteps_ENRG = 100
        # solver message output level: 0 ~ 2
        self.print_level = 1
        # steady-state solver keywords
        self.SSsolverkeywords: dict[str, Union[int, float, str, bool]] = {}

    @property
    def steady_state_tolerances(self) -> tuple[float, float]:
        """Get tolerance for the steady-state search algorithm

        Returns
        -------
            tuple, [absolute_tolerance, relative_tolerance]
                absolute_tolerance: double
                    absolute tolerance
                relative_tolerance: double
                    relative tolerance

        """
        return (self.SSabsolute_tolerance, self.SSrelative_tolerance)

    @steady_state_tolerances.setter
    def steady_state_tolerances(self, tolerances: tuple[float, float]):
        """Set the absolute and the relative tolerances
        for the steady-state solution search algorithm

        Parameters
        ----------
            tolerances: tuple, [absolute_tolerance, relative_tolerance]
                absolute_tolerance: double
                    absolute tolerance for steady-state search algorithm
                relative_tolerance: double
                    relative tolerance for steady-state search algorithm

        """
        ierr = 0
        if tolerances[0] > 0.0:
            self.SSsolverkeywords["ATOL"] = tolerances[0]
            self.SSabsolute_tolerance = tolerances[0]
        else:
            ierr = 1

        if tolerances[1] > 0.0:
            self.SSsolverkeywords["RTOL"] = tolerances[1]
            self.SSrelative_tolerance = tolerances[1]
        else:
            ierr = 1

        if ierr > 0:
            msg = [Color.PURPLE, "tolerance must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    @property
    def time_stepping_tolerances(self) -> tuple[float, float]:
        """Get tolerance for the pseudo time stepping solution algorithm

        Returns
        -------
            tuple, [absolute_tolerance, relative_tolerance]
                absolute_tolerance: double
                    absolute tolerance for time stepping algorithm
                relative_tolerance: double
                    relative tolerance for time stepping algorithm

        """
        return (self.TRabsolute_tolerance, self.TRrelative_tolerance)

    @time_stepping_tolerances.setter
    def time_stepping_tolerances(self, tolerances: tuple[float, float]):
        """Set the absolute and the relative tolerances
        for the pseudo time stepping solution algorithm

        Parameters
        ----------
            tolerances: tuple, [absolute_tolerance, relative_tolerance]
                absolute_tolerance: double
                    absolutie tolerance for the pseudo time stepping
                relative_tolerance: double
                    relative tolerance for the pseudo time stepping

        """
        ierr = 0
        if tolerances[0] > 0.0:
            self.SSsolverkeywords["ATIM"] = tolerances[0]
            self.TRabsolute_tolerance = tolerances[0]
        else:
            ierr = 1

        if tolerances[1] > 0.0:
            self.SSsolverkeywords["RTIM"] = tolerances[1]
            self.TRrelative_tolerance = tolerances[1]
        else:
            ierr = 1

        if ierr > 0:
            msg = [Color.PURPLE, "tolerance must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_max_pseudo_transient_call(self, maxtime: int):
        """Set the maximum number of call to the pseudo transient algorithm
        in an attempt to find the steady-state solution

        Parameters
        ----------
            maxtime: integer
                max number of pseudo transient calls/attempts

        """
        if maxtime >= 1:
            self.SSsolverkeywords["MAXTIME"] = maxtime
            self.maxpseudotransient = maxtime
        else:
            msg = [Color.PURPLE, "parameter must >= 1.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_max_timestep_iteration(self, maxiteration: int):
        """Set the maximum number of iterations per time step when performing the pseudo transient algorithm

        Parameters
        ----------
            maxtime: integer
                max number of iterations per pseudo time step

        """
        if maxiteration >= 1:
            self.SSsolverkeywords["TRMAXITER"] = maxiteration
            self.TRmaxiteration = maxiteration
        else:
            msg = [Color.PURPLE, "parameter must >= 1.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_max_search_iteration(self, maxiteration: int):
        """Set the maximum number of iterations per search when performing
        the steady-state search algorithm.

        Parameters
        ----------
            maxtime: integer
                max number of iterations per steady-state search

        """
        if maxiteration >= 1:
            self.SSsolverkeywords["SSMAXITER"] = maxiteration
            self.SSmaxiteration = maxiteration
        else:
            msg = [Color.PURPLE, "parameter must >= 1.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_initial_timesteps(self, initsteps: int):
        """Set the number of pseudo time steps to be performed to establish a "better"
        set of guessed solution before start the actual steady-state solution search.

        Parameters
        ----------
            initsteps: integer
                number of initial pseudo time steps

        """
        if initsteps >= 1:
            self.SSsolverkeywords["ISTP"] = initsteps
            self.numbinitialpseudosteps = initsteps
        else:
            msg = [Color.PURPLE, "parameter must >= 1.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_species_floor(self, floor_value: float):
        """Set the minimum species fraction value allowed during steady-state solution search.

        Parameters
        ----------
            floor_value: double
                minimum species fraction value

        """
        if np.abs(floor_value) < 1.0:
            self.SSsolverkeywords["SFLR"] = floor_value
            self.speciesfloor = floor_value
        else:
            msg = [Color.PURPLE, "species floor value must < 1.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_temperature_ceiling(self, ceilingvalue: float):
        """Set the maximum temperature value allowed during steady-state solution search.

        Parameters
        ----------
            ceilingvalue: double
                maximum temperature value

        """
        if ceilingvalue > 300.0:
            self.SSsolverkeywords["TBND"] = ceilingvalue
            self.maxTbound = ceilingvalue
        else:
            msg = [Color.PURPLE, "temperature value must > 300.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_species_reset_value(self, resetvalue: float):
        """Set the positive value to reset any negative species fraction in
        intermediate solutions during iterations.

        Parameters
        ----------
            resetvalue: double
                positive value to reset negative species fraction

        """
        if resetvalue >= 0.0:
            self.SSsolverkeywords["SPOS"] = resetvalue
            self.species_positive = resetvalue
        else:
            msg = [Color.PURPLE, "species fraction value must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_max_pseudo_timestep_size(self, dtmax: float):
        """Set the maximum time step sizes allowed by the pseudo time stepping solution.

        Parameters
        ----------
            dtmax: double
                maximum time step size allowed

        """
        if dtmax > 0.0:
            self.SSsolverkeywords["DTMX"] = dtmax
            self.TRmaxstepsize = dtmax
        else:
            msg = [Color.PURPLE, "time step size must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_min_pseudo_timestep_size(self, dtmin: float):
        """Set the minimum time step size allowed by the pseudo time stepping solution.

        Parameters
        ----------
            dtmin: double
                minimum time step size allowed

        """
        if dtmin > 0.0:
            self.SSsolverkeywords["DTMN"] = dtmin
            self.TRminstepsize = dtmin
        else:
            msg = [Color.PURPLE, "time step size must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_pseudo_timestep_age(self, age: int):
        """Set the minimum number of time steps taken before allowing time step size increase.

        Parameters
        ----------
            age: integer
                min age of the pseudo time step size

        """
        if age > 0:
            self.SSsolverkeywords["IRET"] = age
            self.timestepsizeage = age
        else:
            msg = [Color.PURPLE, "number of time step must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_Jacobian_age(self, age: int):
        """Set the number of steady-state searches before re-evaluate the Jacobian matrix.

        Parameters
        ----------
            age: integer
                age of the steady-state Jacobian matrix

        """
        if age > 0:
            self.SSsolverkeywords["NJAC"] = age
            self.SSJacobianage = age
        else:
            msg = [Color.PURPLE, "number of time step must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_pseudo_Jacobian_age(self, age: int):
        """Set the number of time steps taken before re-evaluate the Jacobian matrix.

        Parameters
        ----------
            age: integer
                age of the pseudo time step Jacobian matrix

        """
        if age > 0:
            self.SSsolverkeywords["TJAC"] = age
            self.TRJacobianage = age
        else:
            msg = [Color.PURPLE, "number of time step must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_damping_option(self, ON: bool):
        """Turn ON (True) or OFF (False) the damping option of the steady-state solver.

        Parameters
        ----------
            ON: boolean
                turn On the damping option

        """
        if isinstance(ON, bool):
            if ON:
                self.SSdamping = 1
            else:
                self.SSdamping = 0
            self.SSsolverkeywords["TWOPNT_DAMPING_OPTIN"] = self.SSdamping
        else:
            msg = [Color.PURPLE, "parameter must be either True or False.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_legacy_option(self, ON: bool):
        """Turn ON (True) or OFF (False) the legacy steady-state solver.

        Parameters
        ----------
            ON: boolean
                turn On the legacy solver

        """
        if isinstance(ON, bool):
            self.use_legacy_technique = ON
            if ON:
                self.SSsolverkeywords["USE_LEGACY_TECHNIQUE"] = "4X"
        else:
            msg = [Color.PURPLE, "parameter must be either True or False.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_print_level(self, level: int):
        """Set the level of information to be provided by the steady-state solver
        to the text output.

        Parameters
        ----------
            level: integer, {0, 1, 2}
                solver message details level (0 ~ 2)

        """
        if level in [0, 1, 2]:
            self.SSsolverkeywords["PRNT"] = level
            self.print_level = level
        else:
            msg = [Color.PURPLE, "print level must be either 0, 1, or 2.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_pseudo_timestepping_parameters(
        self, numb_steps: int = 100, step_size: float = 1.0e-6, stage: int = 1
    ):
        """Set the parameters for the pseudo time stepping process of the steady state solver.

        Parameters
        ----------
            numb_step: integer, default = 100
                the number of pseudo time steps to be taken during each time stepping process
            step_size: double, default = 1.0e-6 [sec]
                the initial time step size for each time stepping process
            stage: integer, {1, 2}
                the stage the time stepping process is in
                1 = fixed temperature stage
                2 = solving energy equation

        """
        if stage in [1, 2]:
            this_key = "TIM" + str(stage)
            this_phrase = this_key + Keyword.fourspaces + str(numb_steps)
            self.SSsolverkeywords[this_phrase] = step_size
        else:
            msg = [Color.PURPLE, "the stage must be either 1 or 2.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

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

"""Plug Flow Reactor (PFR) model."""

import copy
from ctypes import c_double, c_int

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core import chemkin_wrapper
from ansys.chemkin.core.batchreactors.batchreactor import BatchReactors
from ansys.chemkin.core.chemistry import (
    check_chemistryset,
    chemistryset_initialized,
    force_activate_chemistryset,
    verbose,
)
from ansys.chemkin.core.color import Color as Color
from ansys.chemkin.core.inlet import Stream
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.reactormodel import Keyword


class PlugFlowReactor(BatchReactors):
    """Plug Flow Reactor (PFR) model with energy equation."""

    def __init__(self, inlet: Stream, label: str = "PFR"):
        """Initialize a generic PFR object."""
        """
        Initialize a generic PFR object.

        Parameters
        ----------
            inlet: Stream object
                an inlet stream representing the gas properties and
                the flow rate at the PFR entrance
            label: string, optional
                reactor name

        """
        # check Inlet
        if isinstance(inlet, Stream):
            # initialization
            super().__init__(inlet, label)
        else:
            # wrong argument type
            msg = [Color.RED, "the first argument must be an Inlet object.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()

        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("PFR", 3))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("CONP", 1))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # defaults for all plug flow reactor models
        self._nreactors = 1
        self._npsrs = c_int(self._nreactors)
        self._ninlets = np.zeros(1, dtype=np.int32)
        # number of zones
        self._nzones = c_int(0)
        # use API mode for PFR simulations
        Keyword.no_fullkeyword = True
        # FORTRAN file unit of the text output file
        self._mylout = c_int(157)
        #
        # starting position [cm]
        self.startposition = c_double(0.0)
        # reactor length [cm]
        self.reactorlength = c_double(0.0)
        # reactor diameter [cm]
        self.reactordiameter = c_double(0.0)
        # cross-sectional flow area [cm2]
        self.reactorflowarea = 0.0
        # flow area given for the inlet
        if inlet._haveflowarea:
            self.reactorflowarea = inlet._flowarea
            # compute reactor diameter from the flow area
            dia = np.sqrt(4.0 * self.reactorflowarea / np.pi)
            self.reactordiameter = c_double(dia)
            if self.reactorflowarea <= 0.0:
                msg = [Color.YELLOW, "inlet flow area is not set.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
        # inlet mass flow rate [g/sec]
        self._massflowrate = c_double(0.0)
        self._flowrate = 0.0
        if inlet._flowratemode < 0:
            # no given in the inlet
            msg = [
                Color.PURPLE,
                "inlet flow rate is not set.\n",
                Color.SPACEx6,
                "please specify the flow rate of the 'Inlet' object.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self.inletflowratemode = inlet._flowratemode
            self._flowrate = inlet._inletflowrate[inlet._flowratemode]
            self.inletflowrate = copy.deepcopy(inlet._inletflowrate)
            if self._flowrate <= 0.0:
                msg = [
                    Color.PURPLE,
                    "inlet flow rate is not set correctly.\n",
                    Color.SPACEx6,
                    "please specify the flow rate of the 'Inlet' object.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        # solver parameters
        self._absolute_tolerance = 1.0e-12
        self._relative_tolerance = 1.0e-6
        # required inputs: (1) reactor length (2) flow area
        self._numb_requiredinput = 2
        self._requiredlist = ["XEND", "AREAF"]
        # always calculate the residence time
        self.setkeyword(key="RTIME", value="ON")
        # solve the momentum equation in most cases
        # turn the momentum equation OFF
        # when the velocity or the pressure profile along the reactor is given
        self.setkeyword(key="MOMEN", value="ON")
        # profile points
        self._profilesize = int(0)

    @property
    def length(self) -> float:
        """Get reactor length."""
        """
        Get reactor length.

        Returns
        -------
            length: double
                reactor length [cm]

        """
        return self.reactorlength.value

    @length.setter
    def length(self, length: float = 0.0e0):
        """Set reactor length."""
        """
        Set reactor length.

        Parameters
        ----------
            length: double
                reactor length [cm]

        """
        if length <= 0.0e0:
            msg = [Color.PURPLE, "reactor length must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self._inputcheck.append("XEND")
            self.reactorlength = c_double(length)

    def set_start_position(self, x0: float):
        """Set the PFR simulation starting position."""
        """
        Set the PFR simulation starting position, the
        default reactor inlet: x0 = 0.0

        Parameters
        ----------
            x0: double, default = 0.0
                starting position

        """
        if x0 >= self.reactorlength.value:
            msg = [Color.PURPLE, "starting position must < reactor length.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        elif x0 <= 0.0:
            msg = [Color.PURPLE, "reactor diameter must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self.startposition = c_double(x0)

    @property
    def diameter(self) -> float:
        """Get reactor diameter."""
        """
        Get reactor diameter.

        Returns
        -------
            diam: double
                Reactor diameter [cm]

        """
        return self.reactordiameter.value

    @diameter.setter
    def diameter(self, diam: float):
        """Set the PFR diameter."""
        """
        Set the PFR diameter.

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
            self._inputcheck.append("AREAF")
            self.reactordiameter = c_double(diam)
            # set flow area at the inlet
            self.reactormixture._haveflowarea = True
            area = np.pi * diam * diam / 4.0
            self.reactormixture._flowarea = area
            self.reactorflowarea = area

    def set_diameter_profile(
        self, x: npt.NDArray[np.double], diam: npt.NDArray[np.double]
    ) -> int:
        """Specify plug-flow reactor diameter profile."""
        """
        Specify plug-flow reactor diameter profile.

        Parameters
        ----------
            x: 1-D double array
                position value of the profile data [cm]
            diam: 1-D double array
                PFR diameter value of the profile data [cm]

        Returns
        -------
            error code: integer

        """
        keyword = "DPRO"
        ierr = self.setprofile(key=keyword, x=x, y=diam)
        if ierr == 0:
            self._inputcheck.append("AREAF")
            # set flow area at the inlet
            self.reactordiameter = c_double(diam[0])
            self.reactormixture._haveflowarea = True
            area = np.pi * diam * diam / 4.0
            self.reactormixture._flowarea = area
            self.reactorflowarea = area
        return ierr

    @property
    def flowarea(self) -> float:
        """Cross-sectional flow area of the PFR."""
        """
        Cross-sectional flow area of the PFR [cm2].

        Returns
        -------
            area: double
                cross-sectional flow rate [cm2]

        """
        return self.reactorflowarea

    @flowarea.setter
    def flowarea(self, area: float):
        """Set the cross-sectional flow area of the PFR."""
        """
        Set the cross-sectional flow area of the PFR.

        Parameters
        ----------
            area: double
                cross-sectional flow area [cm2]

        """
        if area <= 0.0:
            msg = [Color.PURPLE, "cross-sectional flow area must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            # set the flow area keyword
            self._inputcheck.append("AREAF")
            self.setkeyword(key="AREAF", value=area)
            self.reactorflowarea = area
            # set flow area at the inlet
            diam = np.sqrt(4.0 * area / np.pi)
            self.reactordiameter = c_double(diam)
            self.reactormixture._haveflowarea = True
            self.reactormixture._flowarea = area

    def set_flowarea_profile(
        self, x: npt.NDArray[np.double], area: npt.NDArray[np.double]
    ) -> int:
        """Specify plug-flow reactor cross-sectional flow area profile."""
        """
        Specify plug-flow reactor cross-sectional flow area profile.

        Parameters
        ----------
            x: 1-D double array
                position value of the profile data [cm]
            area: 1-D double array
                PFR cross-sectional flow area value of the profile data [cm2]

        Returns
        -------
            error code: integer

        """
        keyword = "AFLO"
        ierr = self.setprofile(key=keyword, x=x, y=area)
        if ierr == 0:
            self._inputcheck.append("AREAF")
            self.setkeyword(key="AREAF", value=area[0])
            # set flow area at the inlet
            diam = np.sqrt(4.0 * area[0] / np.pi)
            self.reactordiameter = c_double(diam)
            self.reactormixture._haveflowarea = True
            self.reactormixture._flowarea = area[0]
            self.reactorflowarea = area[0]
        return ierr

    def set_inlet_viscosity(self, visc: float):
        """Set the gas mixture viscocity at the PFR inlet."""
        """
        Set the gas mixture viscocity at the PFR inlet.

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
            # set the mixture viscosity keyword
            self.setkeyword(key="VISC", value=visc)

    def set_solver_max_timestep_size(self, size: float):
        """Set the maximum time step size allowed."""
        """Set the maximum time step size allowed by the transient solver.

        Parameters
        ----------
            size: double
                maximum solver step size [cm]

        """
        if size > 0.0e0:
            self.setkeyword(key="DXMX", value=size)
        else:
            msg = [Color.PURPLE, "solver timestep size must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_pseudo_surface_velocity(self, vel: float):
        """Set the pseudo surface velocity at the active surface."""
        """Set the pseudo surface velocity at the reactive surface
        to improve convergence due to surface chemistry stiffness.

        Note: set this parameter only when having convergence issue
        with surface chemistry

        Parameters
        ----------
            vel: double, default = 0.0
                pseudo surface velocity [cm/sec]

        """
        if vel > 0.0e0:
            self.setkeyword(key="PSV", value=vel)
        else:
            msg = [Color.PURPLE, "pseudo velocity must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def mass_flowrate(self) -> float:
        """Get plug flow reactor inlet mass flow rate."""
        """
        Get plug flow reactor inlet mass flow rate [g/sec].

        Returns
        -------
            massflowrate: double
                mass flow rate [g/sec]

        """
        return self.reactormixture.mass_flowrate

    @property
    def velocity(self) -> float:
        """Get plug flow reactor inlet velocity [cm/sec]."""
        """
        Returns
        -------
            vel: double
                inlet velocity [cm/sec]

        """
        return self.reactormixture.velocity

    @property
    def vol_flowrate(self) -> float:
        """Get plug flow reactor inlet volumetric flow rate."""
        """
        Get plug flow reactor inlet volumetric flow rate [cm3/sec].

        Returns
        -------
            volflowrate: double
                volumetric flow rate [cm3/sec]

        """
        return self.reactormixture.vol_flowrate

    @property
    def sccm(self) -> float:
        """Get inlet volumetric flow rate in SCCM."""
        """Get plug flow reactor inlet volumetric flow rate in
        SCCM [standard cm3/min].

        Returns
        -------
            volflowrate: double
                volumetric flow rate in SCCM [standard cm3/min]

        """
        return self.reactormixture.sccm

    def __process_keywords(self) -> int:
        """Process input keywords for the PFR model."""
        """
        Process input keywords for the PFR model.

        Returns
        -------
            error code: integer

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
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
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
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return ierr
        # prepare inlet conditions
        # get inlet mass flow rate
        self._massflowrate = c_double(self.mass_flowrate)
        # inlet mass fraction
        y_init = self.reactormixture.y
        # surface sites (not applicable)
        site_init = np.zeros(1, dtype=np.double)
        # bulk activities (not applicable)
        bulk_init = np.zeros_like(site_init, dtype=np.double)
        # set reactor inlet conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("PFR", 3):
            ierrc = chemkin_wrapper.chemkin.KINAll0D_SetupPFRInputs(
                self._chemset_index,
                self.startposition,
                self.reactorlength,
                self._temperature,
                self._pressure,
                self._heat_loss_rate,
                self.reactordiameter,
                site_init,
                bulk_init,
                self._massflowrate,
                y_init,
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
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return ierrc
            # turn OFF the momentum equation when pressure profile is set
            if self._numbprofiles > 0:
                if "PPRO" in self._profiles_list or "VELPRO" in self._profiles_list:
                    self.setkeyword(key="MOMEN", value="OFF")
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
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    return ierrc

        if ierr == 0 and self._numbprofiles > 0:
            for p in self._profiles_list:
                key = bytes(p.profilekey, "utf-8")
                npoints = c_int(p.size)
                x = p.pos
                y = p.value
                ierr_prof = chemkin_wrapper.chemkin.KINAll0D_SetProfileParameter(
                    key, npoints, x, y
                )
                ierr += ierr_prof
            if ierr_prof != 0:
                msg = [
                    Color.PURPLE,
                    "failed to set up profile keywords,",
                    "error code =",
                    str(ierr_prof),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return ierr_prof
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
                if err_inputs == 0:
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
                        str(err_inputs),
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
            else:
                msg = [
                    Color.PURPLE,
                    "failed to process additional keywords, error code =",
                    str(err_inputs),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
        #
        ierr = ierr + err_inputs + err_key

        return ierr

    def __run_model(self) -> int:
        """Run the PFR model after the keywords are processed."""
        """
        Run the PFR model after the keywords are processed.

        Returns
        -------
            error code: integer

        """
        # run the simulation without keyword inputs
        ierr = chemkin_wrapper.chemkin.KINAll0D_Calculate(self._chemset_index)
        return ierr

    def run(self) -> int:
        """Perform the common steps to run a Chemkin reactor model."""
        """
        Perform the common steps to run a Chemkin reactor model.

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
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            return ret_val
        logger.debug("Processing keywords complete")

        # run reactor model
        msg = [Color.YELLOW, "running reactor simulation ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
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
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        else:
            msg.insert(0, Color.RED)
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)

        return ret_val


class PFREnergyConservation(PlugFlowReactor):
    """Plug Flow Reactor (PFR) model with energy equation."""

    def __init__(self, inlet, label: str = "PFR"):
        """Initialize a PFR object that solves the Energy Equation."""
        """
        Initialize a PFR object that solves the Energy Equation.

        Parameters
        ----------
            inlet: Inlet object
                an inlet stream representing the gas properties and
                the flow rate at the PFR entrance
            label: string, optional
                reactor name

        """
        # check Inlet
        if isinstance(inlet, Stream):
            # initialization
            super().__init__(inlet, label)
        else:
            # wrong argument type
            msg = [Color.RED, "the first argument must be an Inlet object.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()

        # set reactor type
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        self._heat_transfer_coefficient = 0.0e0
        self._ambient_temperature = 3.0e2
        # external heat transfer area per reactor length [cm2/cm]
        self._heat_transfer_area = 0.0e0
        # set up basic PFR parameters
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
                "failed to initialize the plug-flow reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()

    @property
    def heat_loss_rate(self) -> float:
        """Get heat loss rate from the reactor to the surroundings."""
        """
        Get heat loss rate from the reactor to the surroundings.

        Returns
        -------
            qloss: double
                heat loss rate [cal/sec-cm]

        """
        return self._heat_loss_rate.value

    @heat_loss_rate.setter
    def heat_loss_rate(self, value: float):
        """Set the heat loss rate."""
        """Set the heat loss rate per length from the reactor
        to the surroundings (required).

        Parameters
        ----------
            value: double, default = 0.0
                heat loss rate [cal/sec-cm]

        """
        self._heat_loss_rate = c_double(value)
        if not Keyword.no_fullkeyword:
            self.setkeyword(key="QLOS", value=value)

    @property
    def heat_transfer_coefficient(self) -> float:
        """Get heat transfer coefficient between the reactor and the surroundings."""
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
        """Set heat transfer coefficient between the reactor and the surroundings."""
        """
        Set heat transfer coefficient between the reactor and the surroundings.

        Parameters
        ----------
            value: double, default = 0.0
                heat transfer coefficient [cal/cm2-K-sec]

        """
        if value < 0.0e0:
            msg = [Color.PURPLE, "heat transfer coefficient must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
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
            ambient_temperature: double
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
            value: double, default = 300.0
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
        """Get heat transfer area."""
        """Get heat transfer area per length between the reactor
        and the surroundings.

        Returns
        -------
            heat_transfer_area: double
                heat transfer area [cm2/cm]

        """
        return self._heat_transfer_area

    @heat_transfer_area.setter
    def heat_transfer_area(self, value: float = 0.0e0):
        """Set heat transfer area."""
        """Set heat transfer area per length between the reactor
        and the surroundings.

        Parameters
        ----------
            value: double, default = 0.0
                heat transfer area [cm2/cm]

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
        """Specify reactor heat transfer area per reactor length profile."""
        """
        Specify reactor heat transfer area per reactor length profile.

        Parameters
        ----------
            x: 1-D double array
                position value of the profile data [cm]
            area: 1-D double array
                heat transfer area value of the profile data [cm2/cm]

        Returns
        -------
            error code: integer

        """
        keyword = "AEXT"
        ierr = self.setprofile(key=keyword, x=x, y=area)
        return ierr

    def set_heat_loss_profile(
        self, x: npt.NDArray[np.double], qloss: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor heat loss rate per length profile."""
        """
        Specify reactor heat loss rate per length profile.

        Parameters
        ----------
            x: 1-D double array
                position value of the profile data [cm]
            qloss: 1-D double array
                heat loss rate value of the profile data [cal/sec-cm]

        Returns
        -------
            error code: integer

        """
        keyword = "QPRO"
        ierr = self.setprofile(key=keyword, x=x, y=qloss)
        return ierr

    def set_velocity_profile(
        self, x: npt.NDArray[np.double], vel: npt.NDArray[np.double]
    ) -> int:
        """Specify axial velocity profile along the plug-flow reactor."""
        """
        Specify axial velocity profile along the plug-flow reactor.

        Parameters
        ----------
            x: 1-D double array
                position value of the profile data [cm]
            vel: 1-D double array
                axial velocity value of the profile data [cm/sec]

        Returns
        -------
            error code: integer

        """
        keyword = "VELPRO"
        ierr = self.setprofile(key=keyword, x=x, y=vel)
        return ierr


class PFRFixedTemperature(PlugFlowReactor):
    """Plug Flow Reactor (PFR) model with given temperature."""

    def __init__(self, inlet, label: str = "PFR"):
        """Initialize a PFR object with given temperature profile."""
        """Initialize a PFR object with given temperature profile
        along the length of the reactor.

        Parameters
        ----------
            inlet: Stream object
                an inlet stream representing the gas properties and
                the flow rate at the PFR entrance
            label: string, optional
                reactor name

        """
        # check Inlet
        if isinstance(inlet, Stream):
            # initialization
            super().__init__(inlet, label)
        else:
            # wrong argument type
            msg = [
                Color.PURPLE,
                "the first argument must be an Inlet object",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # set reactor type
        self._energytype = c_int(self.EnergyTypes.get("GivenT", 2))
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
                "failed to initialize the plug-flow reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()

    def set_temperature_profile(
        self, x: npt.NDArray[np.double], temp: npt.NDArray[np.double]
    ) -> int:
        """Specify reactor temperature profile."""
        """
        Specify reactor temperature profile.

        Parameters
        ----------
            x: 1-D double array
                position value of the profile data [cm]
            temp: 1-D double array
                temperature value of the profile data [K]

        Returns
        -------
            error code: integer

        """
        keyword = "TPRO"
        ierr = self.setprofile(key=keyword, x=x, y=temp)
        return ierr

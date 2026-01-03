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

"""Spark Ignition (SI) engine model."""

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
    verbose,
)
from ansys.chemkin.core.color import Color as Color
from ansys.chemkin.core.engines.engine import Engine
from ansys.chemkin.core.inlet import Stream
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.reactormodel import Keyword


class SIengine(Engine):
    """Spark Ignition (SI) engine model."""

    def __init__(self, reactor_condition: Stream, label: Union[str, None] = None):
        """Initialize a spark-ignition Engine object."""
        """
        Initialize a spark-ignition Engine object.

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside
                the engine cylinder/zone
            label: string, optional
                engine reactor name

        """
        # set default number of zone(s)
        # 2 zones: the unburned and the burned zones
        nzones = 2
        # set default label
        if label is None:
            label = "SI"

        # use the first zone to initialize the engine model
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("SI", 5))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("ICEN", 3))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # defaults for all closed homogeneous reactor models
        # 2 zones: the unburned and the burned zones
        self._nreactors = nzones
        self._npsrs = c_int(1)
        self._ninlets = np.zeros(1, dtype=np.int32)
        # number of zones
        self._nzones = c_int(nzones)
        # use API mode for SI simulations
        Keyword.no_fullkeyword = True
        # FORTRAN file unit of the text output file
        self._mylout = c_int(156)
        # profile points
        self._profilesize = int(0)
        # burn mass profile mode
        # 0: unset
        # 1: Wiebe function with n, b, SOI, burn duration
        # 2: anchor points 10%, 50%, and 90% mass burned CAs
        # 3: burned mass fraction profile, SOI, burn duration
        self._burnmode: int = 0
        self.sparktiming = -180.0
        self.burnduration = 0.0
        self.wieben = 2.0
        self.Wiebeb = 5.0
        self.massburnedCA10 = -180.0
        self.massburnedCA50 = -180.0
        self.massburnedCA90 = -180.0
        # combustion efficiency
        self.burnefficiency = 1.0
        # numbwer of points in the mass burned fraction profile
        self.MBpoints = 0
        self.MBangles: list[float] = []
        self.MBfractions: list[float] = []
        # set up basic SI engine model parameters
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
            # setup SI engine model working arrays
            ierr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._mylout, self._chemset_index
            )
            ierr *= 10
        if ierr != 0:
            msg = [
                Color.RED,
                "failed to initialize the SI engine model",
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

    def wiebe_parameters(self, n: float, b: float):
        """Set Wiebe function parameters."""
        """
        Set Wiebe function parameters.

        .. math::

            Wiebe = 1 - exp^{-b[(CA-SOC)/duration]^(n+1)]}

        Parameters
        ----------
            n: double
                exponent parameter of the Wiebe function [-]
            b: double
                multiplier parameter of the Wiebe function [-]

        """
        # check input values
        if n <= 0.0 or b <= 0.0:
            msg = [
                Color.PURPLE,
                "Wiebe function parameters n and b must > 0.0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        #
        if self._burnmode > 0:
            msg = [
                Color.YELLOW,
                "previous burned mass profile setup will be overridden.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        #
        self._burnmode = 1
        self.wieben = n
        self.wiebeb = b

    def set_burn_timing(self, soc: float, duration: float = 0.0):
        """Set SI engine SOC timing."""
        """
        Set SI engine start of combustion (SOC) timing.

        Parameters
        ----------
            soc: double
                start of combustion in crank angle [degree]
            duration: double
                burn duration in crank angles [degree]

        """
        if soc <= self.ivc_ca:
            msg = [
                Color.PURPLE,
                "start of combustion CA must > start of simulation CA",
                str(self.ivc_ca),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if duration <= 0.0:
            msg = [Color.PURPLE, "mass burned duration must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        #
        self.sparktiming = soc
        self.burnduration = duration

    def set_burn_anchor_points(self, ca10: float, ca50: float, ca90: float):
        """Set the SI mass burned profile using the anchor points."""
        """
        Set the SI mass burned profile using the anchor points.

        Parameters
        ----------
            ca10: double
                crank angle of 10% mass burned [degree]
            ca50: double
                crank angle of 50% mass burned [degree]
            ca90: double
                crank angle of 90% mass burned [degree]

        """
        if ca10 > ca50:
            msg = [
                Color.PURPLE,
                "the anchor points must be given in ascending order.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if ca90 < ca50:
            msg = [
                Color.PURPLE,
                "the anchor points must be given in ascending order.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if ca10 <= self.ivc_ca:
            msg = [
                Color.PURPLE,
                "anchor point CAs must > start of simulation CA",
                str(self.ivc_ca),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        #
        if self._burnmode > 0:
            msg = [
                Color.YELLOW,
                "previous burned mass profile setup will be overridden.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        #
        self._burnmode = 2
        self.massburnedCA10 = ca10
        self.massburnedCA50 = ca50
        self.massburnedCA90 = ca90

    def set_mass_burned_profile(
        self, crankangles: npt.NDArray[np.double], fractions: npt.NDArray[np.double]
    ) -> int:
        """Specify SI engine mass burned fraction profile."""
        """
        Specify SI engine mass burned fraction profile.

        Parameters
        ----------
            crankangles: 1-D double array
                normalized crank angles of the profile data [degree]
                the crank angles must 0 <= and <= 1
            fractions: 1-D double array
                mass burned fraction of the profile data [-]

        Returns
        -------
            error code: integer

        """
        # set the mass burned profile
        ierror = 0
        self.MBpoints = len(crankangles)
        if len(fractions) != self.MBpoints:
            msg = [Color.PURPLE, "data arrays must have the same size.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierror = 1
        elif self.MBpoints > 1:
            self.MBangles = copy.deepcopy(crankangles)
            self.MBfractions = copy.deepcopy(fractions)
            self._burnmode = 3
        else:
            msg = [Color.PURPLE, "profile must have more than 1 data pair.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierror = 2
        return ierror

    def set_combustion_efficiency(self, efficiency: float):
        """Set the overall combustion efficiency."""
        """
        Set the overall combustion efficiency.

        Parameters
        ----------
            efficiency: double, default = 1.0
                combustion efficiency [-]

        """
        # check value
        if efficiency < 0.0 or efficiency > 1.0:
            msg = [Color.PURPLE, "efficiency must > 0.0 and <= 1.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # set keyword
        self.burnefficiency = efficiency
        self.setkeyword(key="BEFF", value=efficiency)

    def set_burned_products_minimum_mole_fraction(self, bound: float):
        """Set the minimum gas species mole fraction value."""
        """Set the minimum gas species mole fraction value from the flame sheet
        to be injected to the burned zone.

        Parameters
        ----------
            bound: double
                minimum species mole fraction value [-]

        """
        if bound > 0.0:
            # set keyword
            self.setkeyword(key="EQMN", value=bound)
        else:
            msg = [Color.PURPLE, "species fraction value must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_wiebe_keywords(self) -> int:
        """Set the Wiebe function parameters."""
        """Set the Wiebe function parameters keywords
        for the SI engine model.

        Returns
        -------
            error code: integer

        """
        ierror = 0
        if self._burnmode == 1:
            # set start of combustion time
            self.setkeyword(key="BINI", value=self.sparktiming)
            # set burn duration
            self.setkeyword(key="BDUR", value=self.burnduration)
            # set Wiebe parameter b
            self.setkeyword(key="WBFB", value=self.wiebeb)
            # set Wiebe parameter n
            self.setkeyword(key="WBFN", value=self.wieben)
        else:
            msg = [
                Color.PURPLE,
                "incorrect burned mass profile set up,",
                "error code = 10",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierror = 10
        return ierror

    def set_burn_anchor_points_keywords(self) -> int:
        """Set the mass burned profile anchor points."""
        """Set the mass burned profile anchor points keywords
        for the SI engine model.

        Returns
        -------
            error code: integer

        """
        ierror = 0
        if self._burnmode == 2:
            # set 10% mass burned crank angle
            self.setkeyword(key="CASC", value=self.massburnedCA10)
            # set 50% mass burned crank angle
            self.setkeyword(key="CAAC", value=self.massburnedCA50)
            # set 90% mass burned crank angle
            self.setkeyword(key="CAEC", value=self.massburnedCA90)
        else:
            msg = [
                Color.PURPLE,
                "incorrect burned mass profile set up,",
                "error code = 11",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierror = 11
        return ierror

    def set_burn_profile_keywords(self) -> int:
        """Set the mass burned fraction profile."""
        """Set the mass burned fraction profile keywords
        for the SI engine model.

        Returns
        -------
            error code: integer

        """
        ierror = 0
        if self._burnmode == 3 and self.MBpoints > 0:
            # set start of combustion time
            self.setkeyword(key="BINI", value=self.sparktiming)
            # set burn duration
            self.setkeyword(key="BDUR", value=self.burnduration)
            # set number of burned mass fraction profile data points
            self.setkeyword(key="NBFP", value=self.MBpoints)
            # set mass burned fraction profile keywords
            for i in range(self.MBpoints):
                # set mass burned fraction profile
                keyline = (
                    "BFP"
                    + Keyword.fourspaces
                    + str(self.MBangles[i])
                    + Keyword.fourspaces
                    + str(self.MBfractions[i])
                )
                self.setkeyword(key=keyline, value=True)
        else:
            msg = [
                Color.PURPLE,
                "incorrect burned mass profile set up,",
                "error code = 12",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierror = 12

        return ierror

    def __process_keywords(self) -> int:
        """Process input keywords."""
        """Process input keywords for the SI engine model.

        Returns
        -------
            error code: integer

        """
        ierr = 0
        ierrc = 0
        err_key = 0
        err_inputs = 0
        err_profile = 0
        # set_verbose(True)
        # verify required inputs
        ierr = self.validate_inputs()
        # check start of combustion CA and burn duration
        if self._burnmode != 2:
            if self.sparktiming <= self.ivc_ca:
                msg = [
                    Color.PURPLE,
                    "missing 'start of combustion' parameter",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr += 1
            if self.burnduration <= 0.0:
                msg = [
                    Color.PURPLE,
                    "missing 'burn duration' parameter",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr += 1
        if ierr != 0:
            msg = [Color.PURPLE, "missing required input keywords", Color.END]
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
                "failed to set up profile keywords,",
                "error code =",
                str(err_profile),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return ierr
        # prepare initial conditions
        # initial mass fraction
        y_init = self.reactormixture.y
        # connecting rod length to crank radius ratio
        lolr = c_double(self.connectrodlength / self.crankradius)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("SI", 5):
            # insert the ICEN keywords
            self.setkeyword(key="ICEN", value=True)
            #
            ierrc = chemkin_wrapper.chemkin.KINAll0D_SetupHCCIInputs(
                self._chemset_index,
                c_double(self.ivc_ca),
                c_double(self.evo_ca),
                c_double(self.enginespeed),
                c_double(self.compressratio),
                c_double(self.borediam),
                c_double(self.enginestroke),
                lolr,
                self._temperature,
                self._pressure,
                self._heat_loss_rate,
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
            # set SI engine parameter
            if self._burnmode == 1:
                # use Wiebe function to specify the mass burned profile
                ierrc = self.set_wiebe_keywords()
                ierr += ierrc
                if ierrc != 0:
                    msg = [Color.PURPLE, "setting Wiebe function keywords.", Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
            elif self._burnmode == 2:
                # use anchor points to specify the mass burned profile
                ierrc = self.set_burn_anchor_points_keywords()
                ierr += ierrc
                if ierrc != 0:
                    msg = [Color.PURPLE, "setting anchor point keywords.", Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
            elif self._burnmode == 3:
                # use normalized profile to specify the mass burned profile
                ierrc = self.set_burn_profile_keywords()
                ierr += ierrc
                if ierrc != 0:
                    msg = [
                        Color.PURPLE,
                        "setting burned mass profile keywords.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
            else:
                msg = [
                    Color.RED,
                    "burned mass rate is not defined for the SI engine simulation.\n",
                    Color.SPACEx6,
                    "Chemkin SI engine model provides three methods:\n",
                    Color.SPACEx6,
                    Color.SPACE,
                    "1. Wiebe function\n",
                    Color.SPACEx6,
                    Color.SPACE,
                    "2. anchor points\n",
                    Color.SPACEx6,
                    Color.SPACE,
                    "3. piece-wise linear normalized CA-burned fraction profile\n",
                    Color.SPACEx6,
                    "please see Chemkin Theory manual for details.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()

            # heat transfer (use additional keywords)
            # solver parameters (use additional keywords)
            # output controls (use additional keywords)
            # ROP (use additional keywords)
            # sensitivity (use additional keywords)
            # ignition delay (use additional keywords)
            # solve integrated heat release rate due to chemical reactions
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
        # check if the wall heat transfer model is set up
        if ierr == 0 and self._wallheattransfer:
            self.set_heat_transfer_keywords()
        #
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
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
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
        """Run the SI engine model."""
        """Run the SI engine model after the keywords are processed.

        Returns
        -------
            error code: integer

        """
        # run the simulation without keyword inputs
        ierr = chemkin_wrapper.chemkin.KINAll0D_Calculate(self._chemset_index)
        return ierr

    def run(self) -> int:
        """Run Chemkin SI engine model method."""
        """
        Run Chemkin SI engine model method.

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
        logger.debug("Clearing output")

        # keyword processing
        msg = [Color.YELLOW, "processing and generating keyword inputs ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        #
        if Keyword.no_fullkeyword:
            # use API calls
            ret_val = (
                self.__process_keywords()
            )  # each reactor model subclass to perform its own keyword processing
        else:
            # use full keywords
            ret_val = -1
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
        msg = [Color.YELLOW, "running SI engine simulation ...", Color.END]
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

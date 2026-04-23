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

"""Chemkin transient perfectly stirred reactor model."""

from ctypes import c_double, c_int
from typing import Union

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
from ansys.chemkin.core.constants import P_ATM, R_GAS_CAL
from ansys.chemkin.core.inlet import (
    Stream,
    clone_stream,
)
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.mixture import equilibrium
from ansys.chemkin.core.reactormodel import Keyword

class TransientPSR(BatchReactors):
    """Generic transient perfectly-stirred reactor model."""

    def __init__(self, reactor_condition: Stream, label: str = ""):
        """Initialize a generic 0-D transient PSR object."""
        """
        Initialize a generic 0-D transient perfectly-stirred reactor object.
        The transient PSR is a batch reactor that has multiple inlets and
        one outlet. The outlet mass flow rate equals the total inlet mass
        flow rate from all connected external inlets.

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside
                the transient PSR
            label: string, optional
                reactor name

        """
        # initialize the base module
        super().__init__(reactor_condition, label)
        #
        # Perfectly-Stirred Reactor (PSR) model
        self._reactortype = c_int(2)
        # default number of reactors
        self._nreactors = 1
        self._npsrs = c_int(self._nreactors)
        self._ninlets = np.zeros(self._nreactors, dtype=np.int32)
        self._nzones = c_int(0)
        # reactor volume [cm3]
        if reactor_condition._vol > 0.0:
            self._volume = c_double(reactor_condition._vol)
        else:
            self._volume = c_double(0.0)
        # reactor residence time [sec]
        self._residencetime = c_double(0.0)
        # reactive surface area [cm2]
        self._reactivearea = c_double(0.0)
        # simulation time [sec] (not in use)
        self._endtime = c_double(0.0)
        # set reactor number (single reactor)
        self.ireac = c_int(1)
        # single reactor (default) or reactor network
        self.standalone = True
        # FORTRAN file unit of the text output file
        self._mylout = c_int(159)
        # inlet information
        # number of external inlets
        self.numbexternalinlets = 0
        # dict of external inlet objects {inlet label: inlet object}
        self.externalinlets: dict[str, Stream] = {}
        # total mass flow rate into this reactor [g/sec]
        self.totalmassflowrate = 0.0
        # profile points
        self._profilesize = int(0)

    @property
    def residence_time(self) -> float:
        """Get reactor residence time."""
        """
        Returns
        -------
            residence_time: double
                apparent PSR residence time [sec]

        """
        return self._residencetime.value

    @residence_time.setter
    def residence_time(self, value: float):
        """Set reactor residence time (required)."""
        """
        Parameters
        ----------
            value: double
                reactor residence time [sec]

        """
        if value > 0.0e0:
            # set reactor residence time
            self._residencetime = c_double(value)
        else:
            msg = [Color.PURPLE, "residence time must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_inlet(self, extinlet: Stream):
        """Add an external inlet to the reactor."""
        """
        Parameters
        ----------
            extinlet: Stream object
                external inlet to the open reactor

        """
        # check Inlet
        if not isinstance(extinlet, Stream):
            # wrong argument type
            msg = [Color.RED, "the argument must be a Stream object", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # current external inlet count
        count = self.numbexternalinlets + 1
        if extinlet.label is None:
            inletname = self.label + "_inlet_" + str(count)
        else:
            # inlet has label
            inletname = self.label + "_" + extinlet.label
        # check inlet name uniqueness
        if inletname in self.externalinlets:
            # append '_dup' to the given inlet name when
            inletname += "_dup"
            msg = [
                Color.YELLOW,
                "inlet",
                inletname,
                "already connected.\n",
                Color.SPACEx6,
                "will append '_dup' to the original name.\n",
                Color.SPACEx6,
                "the new inlet name is",
                inletname,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        # check inlet flow rate
        if extinlet._flowratemode < 0:
            # no given in the inlet
            msg = [
                Color.PURPLE,
                "inlet flow rate is not set.\n",
                Color.SPACEx6,
                "specify flow rate of the 'Inlet' object",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            if extinlet._use_flow_rate_profile:
                # use the first profile data point value
                flowrate = extinlet.convert_to_mass_flowrate()
                if flowrate <= 0.0:
                    msg = [
                        Color.YELLOW,
                        "initial inlet flow rate of stream",
                        extinlet.label,
                        "<= 0.\n",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
                    if count == 1:
                        # first inlet cannot have zero flow rate
                        # at the start of simulation
                        exit()
            else:
                flowrate = extinlet.mass_flowrate
                if flowrate <= 0.0:
                    msg = [
                        Color.PURPLE,
                        "inlet flow rate of stream",
                        extinlet.label,
                        "<= 0.\n",
                        Color.SPACEx6,
                        "specify flow rate of the 'Stream'.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
        # add the inlet object to the inlet dict of the reactor
        self.externalinlets[inletname] = extinlet
        self.numbexternalinlets = count
        self.totalmassflowrate += flowrate
        #
        msg = [Color.YELLOW, "new inlet", inletname, "is added.", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)

    def reset_inlet(self, new_stream: Stream):
        """Reset the properties of an existing external inlet."""
        """Reset the properties of an existing external inlet from the reactor
        by the inlet name.

        Parameters
        ----------
            new_stream: Stream object
                the updated inlet properties with the same stream label

        """
        # check input
        if not isinstance(new_stream, Stream):
            msg = [Color.PURPLE, "the argument must be a Stream.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # update the named inlet from the externalinlets dict
        missed = True
        # construct the full inlet name
        inletname = self.label + "_" + new_stream.label
        # loop over the external inlet dict
        for iname, inlet in self.externalinlets.items():
            if inletname == iname:
                # found matching inlet
                missed = False
                # take out the original mass flow rate contribution
                self.totalmassflowrate -= inlet.mass_flowrate
                # add the new mass flow rate contribution
                self.totalmassflowrate += new_stream.mass_flowrate
                # update the inlet properties
                clone_stream(new_stream, inlet)

        if missed:
            msg = [Color.PURPLE, "inlet", new_stream.label, "is not found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def remove_inlet(self, name: str):
        """Delete an existing external inlet from the reactor by the inlet name."""
        """
        Parameters
        ----------
            name: string
                external inlet name/label

        """
        # check input
        if not isinstance(name, str):
            msg = [Color.PURPLE, "the argument must be a string.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # delete the named inlet from the externalinlets dict
        missed = False
        inletname = self.label + "_" + name
        if inletname in self.externalinlets:
            # existing inlet
            extinlet = self.externalinlets.pop(inletname, None)
            if extinlet is None:
                missed = True
            elif not isinstance(extinlet, Stream):
                # some internal messed up
                missed = True
                msg = [Color.RED, name, "is not an inlet.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()
            else:
                # decrease the external inlet count by 1
                self.numbexternalinlets -= 1
                # take out the mass flow rate contribution
                self.totalmassflowrate -= extinlet.mass_flowrate
                #
                msg = [
                    Color.YELLOW,
                    "inlet",
                    name,
                    "is removed from reactor",
                    self.label,
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
        else:
            # not in the external inlet dict
            missed = True
        if missed:
            msg = [Color.PURPLE, "inlet", name, "is not found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def net_mass_flowrate(self) -> float:
        """Get the net external inlet mass flow rate."""
        """
        Returns
        -------
            massflowrate: double
                net/total external mass flow rate into the reactor [g/sec]

        """
        return self.totalmassflowrate

    @property
    def net_vol_flowrate(self) -> float:
        """Get the net external volumetric flow rate."""
        """
        Returns
        -------
            volflowrate: double
                net/total external volumetric flow rate into the reactor [cm3/sec]

        """
        vrate = 0.0e0
        inletlist = list(self.externalinlets.keys())
        for inl in inletlist:
            # get inlet volumetric flow rate
            vrate += self.externalinlets[inl].vol_flowrate
        del inletlist
        return vrate

    @property
    def number_external_inlets(self) -> int:
        """Get the number of external inlets to the reactor."""
        """
        Returns
        -------
            ninlet: integer
                total number of external inlets to the reactor

        """
        return self.numbexternalinlets

    def set_inlet_keywords(self) -> int:
        """Set up inlet keywords."""
        """
        Returns
        -------
            error code: integer

        """
        ierr = 0
        ierrc = 0
        # find any inlet flow rate profile data
        for i, inlet in enumerate(self.externalinlets.values()):
            if inlet._use_flow_rate_profile:
                prof_key = inlet.flow_rate_profile_keys.get(inlet._flowratemode, "NA")
                if prof_key == "NA":
                    msg = [
                        Color.PURPLE,
                        "failed to find flow rate profile for inlet",
                        inlet.label,
                        ".\n",
                        Color.SPACEx6,
                        "the profile type =",
                        str(inlet._flowratemode),
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    ierrc = i + 1
                if ierrc == 0:
                    time, flowrate = inlet.flowrate_profile.get(
                        prof_key, ([0.0], [0.0])
                    )
                    tag = prof_key
                    if self.numbexternalinlets > 1:
                        # insert inlet stream name
                        tag += Keyword.fourspaces
                        # l = inlet.label
                        sub_tag = "Inlet" + str(i + 1) + "_Reactor1"
                        tag += sub_tag
                    # profile data point
                    # set up the profile data keyword
                    _ = self.setprofile(key=tag, x=time, y=flowrate, label=True)
                else:
                    exit()

        # loop over all external inlets into the reactor
        i_inlet = 0
        flowrate_sum = 0.0
        check_total_flowrate = True
        #
        for key, inlet in self.externalinlets.items():
            # flow rate profile is given for this inlet
            if inlet._use_flow_rate_profile:
                _, y = inlet.flowrate_profile.get(prof_key, ([0.0], [0.0]))
                flowrate = y[0]
                flowrate_sum += flowrate
                # use full keyword mode for the inlet
                key_mode = Keyword.no_fullkeyword
                Keyword.no_fullkeyword = False
                tag = ""
                if self.numbexternalinlets > 1:
                    # use the assigned inlet stream name
                    # tag = key.rstrip()
                    # the reactor name and the inlet name are
                    # hard coded in chemkin-CFD-API
                    tag = "Inlet" + str(i_inlet + 1) + "_Reactor1"
                # inlet mole fraction
                _, species_lines = self.create_inletspeciesinputlines(
                    inlet_name=tag,
                    threshold=1.0e-10,
                    molefrac=inlet.x,
                )
                for line in species_lines:
                    self.setkeyword(key=line, value=True)
                Keyword.no_fullkeyword = key_mode
                # do not check total mass flow rate
                if inlet._flowratemode != 0:
                    check_total_flowrate = False
            # get inlet mass flow rate
            elif inlet._flowratemode == 0:
                # mass flow rate is specified for this inlet
                flowrate = inlet.mass_flowrate
                flowrate_sum += flowrate
                # inlet temperature
                t_inlet = inlet.temperature
                # inlet mass fraction
                y_inlet = inlet.y
                #
                if np.isclose(0.0, flowrate, atol=1.0e-6):
                    msg = [Color.PURPLE, "inlet", key, "has zero flow rate", Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    ierrc = 100 + i_inlet + 1
                else:
                    # set inlet inputs
                    ierrc = chemkin_wrapper.chemkin.KINAll0D_SetupPSRInletInputs(
                        self._chemset_index,
                        self.ireac,
                        c_int(i_inlet + 1),
                        c_double(flowrate),
                        c_double(t_inlet),
                        y_inlet,
                    )
                ierr += ierrc
            else:
                # volumetric flow rate is specified for this inlet
                # use full keyword mode for the inlet
                key_mode = Keyword.no_fullkeyword
                Keyword.no_fullkeyword = False
                tag = ""
                if self.numbexternalinlets > 1:
                    # use the assigned inlet stream name
                    # tag = key.rstrip()
                    # the reactor name and the inlet name are
                    # hard coded in chemkin-CFD-API
                    tag = "Inlet" + str(i_inlet + 1) + "_Reactor1"
                # inlet mole fraction
                _, species_lines = self.create_inletspeciesinputlines(
                    inlet_name=tag,
                    threshold=1.0e-10,
                    molefrac=inlet.x,
                )
                for line in species_lines:
                    self.setkeyword(key=line, value=True)

                if inlet._flowratemode == 1:
                    this_key = "VDOT" + Keyword.fourspaces + tag
                    self.setkeyword(key=this_key.rstrip(), value=inlet.vol_flowrate)
                elif inlet._flowratemode == 3:
                    this_key = "SCCM" + Keyword.fourspaces + tag
                    self.setkeyword(key=this_key.rstrip(), value=inlet.sccm)
                if inlet._t_set == 1:
                    # inlet temperature is provided
                    this_key = "TINL" + Keyword.fourspaces + tag
                    self.setkeyword(key=this_key.rstrip(), value=inlet.temperature)
                # do not check total mass flow rate
                check_total_flowrate = False
                Keyword.no_fullkeyword = key_mode
            i_inlet += 1
        # check number of external inlet
        if i_inlet == 0:
            msg = [Color.PURPLE, "PSR has no external inlet.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr += 10
        elif i_inlet != self.numbexternalinlets:
            msg = [
                Color.PURPLE,
                "inconsistent number of external inlets.\n",
                Color.SPACEx6,
                "expected number of inlets:",
                str(self.numbexternalinlets),
                "\n",
                Color.SPACEx6,
                "actual number of inlets:",
                str(i_inlet),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr += 11

        # check total mass flow rate
        if ierr == 0 and check_total_flowrate:
            # check total mass flow rate
            if not np.isclose(flowrate_sum, self.totalmassflowrate, atol=1.0e-6):
                msg = [
                    Color.PURPLE,
                    "inconsistent inlet mass flow rate value.\n",
                    Color.SPACEx6,
                    "expected total mass flow rate:",
                    str(self.totalmassflowrate),
                    "\n",
                    Color.SPACEx6,
                    "actual total mass flow rate:",
                    str(flowrate_sum),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr += 12
        return ierr

    def set_initial_conditions(
        self, option: str, guess_temp: Union[float, None] = None
    ):
        """Reset the initial reactor gas mixture."""
        """Reset the initial reactor gas mixture properties by using
        the equilibrium states. This is useful for simulating premixed
        combustion problems.

        Parameters
        ----------
            option: str, {"TP", "HP", "TT"}
                options for additional transformation of the
                guessed mixture composition.
                "HP" indicates the new guessed mixture is the
                equilibrium state with the same enthalpy.
                "TP" indicates the new guessed mixture is the
                equilibrium state at the new given guess_temp
                "TT" indicates the new guessed mixture has the
                composition but at the new given guess_temp
            guess_temp: double, optional
                new mixture temperature [K] used by options "TP" and "TT"

        """
        if option.upper() == "HP":
            # use the constant enthalpy equilibirum mixture as the new guess
            newmixture = equilibrium(self.reactormixture, opt=5)
            # update the guess mixture properties
            self.reactormixture.temperature = newmixture.temperature
            self.temperature = newmixture.temperature
            self.reactormixture.x = newmixture.x
            del newmixture
        else:
            if guess_temp is None:
                # use the current mixture temperature
                msg = [
                    Color.PURPLE,
                    "new gas temperature is not provided,\n",
                    "the original temperature",
                    str(self.reactormixture.temperature),
                    "is applied.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            elif guess_temp < 250.0:
                # use the current mixture temperature
                msg = [
                    Color.PURPLE,
                    "new gas temperature value is invalid,\n",
                    "the original temperature",
                    str(self.reactormixture.temperature),
                    "is applied.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            else:
                # use the given temperature
                self.reactormixture.temperature = guess_temp
                self.temperature = guess_temp

            if option.upper() == "TP":
                # use the constant temperature equilibirum mixture as the new guess
                newmixture = equilibrium(self.reactormixture, opt=1)
                # update the guess mixture composition
                self.reactormixture.x = newmixture.x
                del newmixture

    def reset_initial_temperature(self, temp: float):
        """Reset the initial reactor gas temperature."""
        """
        Parameters
        ----------
            temp: double
                initial reactor gas temperature [K]

        """
        if temp < 250.0:
            # bad value, use the current mixture temperature
            msg = [
                Color.PURPLE,
                "new gas temperature value is invalid,\n",
                "the original temperature",
                str(self.reactormixture.temperature),
                "is applied.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            exit()
        else:
            # use the given temperature
            self.reactormixture.temperature = temp
            self.temperature = temp

    def reset_initial_composition(
        self,
        fraction: Union[npt.NDArray[np.double], list[tuple[str, float]]],
        mode: str = "mole",
    ):
        """Reset the initial reactor gas composition."""
        """
        Parameters
        ----------
            fraction: 1-D double array, dimension = number_species or PyChemkin recipe
                initial reactor gas composition
            mode: string, {"mole", "mass"}
                the given fractions are mole or mass fractions

        """
        if mode.lower() == "mole":
            # set mole fraction
            self.reactormixture.x = fraction
        elif mode.lower() == "mass":
            # set mass fraction
            self.reactormixture.y = fraction
        else:
            # error condition
            msg = [
                Color.PURPLE,
                "the mode of the new composition is invalid,",
                'should be either "mole" or "mass".',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            exit()

    def validate_inputs(self) -> int:
        """Validate the keywords."""
        """Validate the keywords specified by the user before
        running the simulation.

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
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
            # verify required inputs one by one
            for k in self._requiredlist:
                if k not in self._inputcheck:
                    ierr += 1
                    msg = [Color.PURPLE, "missing required input", k, Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
            return ierr

    def __process_keywords(self) -> int:
        """Process input keywords for the reactor model."""
        """
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
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return ierr
        # check external inlet
        if self.numbexternalinlets <= 0 or self.totalmassflowrate <= 0.0:
            # no external inlet for an open reactor
            if self.standalone:
                ierr = 100
                msg = [
                    Color.PURPLE,
                    "missing external inlet for an open reactor.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return ierr
        else:
            # set up inlets
            err_inputs = self.set_inlet_keywords()
            if err_inputs != 0:
                print(f"error code = {err_inputs}")
            ierr += err_inputs
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
            if ierrc != 0:
                msg = [
                    Color.PURPLE,
                    "profile data generation error, error code =",
                    str(ierrc),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return ierrc
        # prepare estimated reactor conditions
        # estimated reactor mass fraction
        y_init = self.reactormixture.y
        if self.has_surface_chemistry:
            # set surface species fractions
            site_init, bulk_init = self.set_init_surface_coverage()
        else:
            # surface sites
            site_init = np.zeros(1, dtype=np.double)
            # bulk activities (not applicable)
            bulk_init = np.zeros_like(site_init, dtype=np.double)
        # set estimated reactor conditions and geometry parameters
        if self._reactortype.value == 2:
            ierrc = chemkin_wrapper.chemkin.KINAll0D_SetupPSRReactorInputs(
                self._chemset_index,
                self.ireac,
                self._endtime,
                self._temperature,
                self._pressure,
                self._volume,
                self._heat_loss_rate,
                self._residencetime,
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
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return ierrc
            # heat transfer (use additional keywords)
            # solver parameters (use additional keywords)
            # output controls (use additional keywords)
            # ROP (use additional keywords)
            # sensitivity (use additional keywords)
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

        if ierr == 0 and self._numbprofiles > 0 and self.standalone:
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
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return err_profile
        if ierr == 0 and self.has_surface_chemistry:
            # set surface related keywords
            self.set_surface_chemistry_keywords()
        if ierr == 0:
            # create input lines from additional user-specified keywords
            if self.standalone:
                # single PSR
                err_inputs, nlines = self.createkeywordinputlines()
            else:
                # PSR is in a cluster
                id_tag = str(self.ireac.value)
                err_inputs, nlines = self.createkeywordinputlines_with_tag(id_tag)
            if err_inputs == 0 and nlines > 0:
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
            elif err_inputs == 0:
                # do nothing
                pass
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
        """Run the reactor model after the keywords are processed."""
        """
        Returns
        -------
            Error code: integer

        """
        # run the simulation without keyword inputs
        ierr = chemkin_wrapper.chemkin.KINAll0D_Calculate(self._chemset_index)
        return ierr

    def run(self) -> int:
        """Perform common steps to run a Chemkin reactor model."""
        """
        Returns
        -------
            Error code: integer

        """
        # initialize the PSR model
        # set up basic PSR parameters
        # number of external inlets
        self._ninlets[0] = self.numbexternalinlets
        #
        # activate the Chemistry set associated with the Reactor instance
        force_activate_chemistryset(self._chemset_index.value)
        #
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
                Color.PURPLE,
                "failed to initialize the transient perfectly-stirred reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
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
            msg = [
                Color.RED,
                "full keyword option not available for PSR models.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            ret_val = 100
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
        logger.debug("processing keywords complete")

        # run reactor model
        msg = [Color.YELLOW, "running reactor simulation ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        # suppress text output to file
        if self.suppress_output:
            ierr = chemkin_wrapper.chemkin.KINAll0D_SuppressOutput()
        if Keyword.no_fullkeyword:
            # use API calls
            ret_val = self.__run_model()
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


class TransientPSRSetVolumeFixedTemperature(TransientPSR):
    """Transient PSR model with reactor volume and fixed reactor temperature."""

    """
    PSR model with given reactor volume (CONV)
    and given reactor temperature (GivenT).
    """

    def __init__(self, reactor_condition: Stream, label: str = ""):
        """Create a transient constant volume perfectly-stirred reactor (PSR)."""
        """
        Create a transient perfectly-stirred reactor (PSR)
        with given reactor volume (CONV)
        and reactor temperature (GivenT).

        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                inlet name/label
        """
        # initialize the base module
        super().__init__(reactor_condition, label)
        # transient
        self._solvertype = c_int(self.SolverTypes.get("Transient", int(1)))
        # specify PSR "reactor volume"
        self._problemtype = c_int(1)
        # given reactor temperature
        self._energytype = c_int(2)
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


class TransientPSRSetVolumeEnergyConservation(TransientPSR):
    """Transient PSR model with reactor volume and solve energy equation."""

    """
    PSR model with given reactor volume (CONV)
    and solve the energy equation (ENERGY).
    """

    def __init__(self, reactor_condition: Stream, label: str = ""):
        """Create a transient constant volume perfectly-stirred reactor (PSR)."""
        """
        Create a transient perfectly-stirred reactor (PSR)
        with given reactor volume (CONV)
        and solve energy equation (ENERGY).

        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                inlet name/label

        """
        # initialize the base module
        super().__init__(reactor_condition, label)
        # transient
        self._solvertype = c_int(self.SolverTypes.get("Transient", int(1)))
        # specify PSR "reactor volume"
        self._problemtype = c_int(1)
        # solve energy equation
        self._energytype = c_int(1)
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

    @property
    def heat_loss_rate(self) -> float:
        """Get heat loss rate from the reactor to the surroundings."""
        """
        Returns
        -------
            qloss: double
                heat loss rate [cal/sec-cm]

        """
        return self._heat_loss_rate.value

    @heat_loss_rate.setter
    def heat_loss_rate(self, value: float):
        """Set the heat loss rate from the reactor to the surroundings."""
        """
        Parameters
        ----------
            value: double, default = 0.0
                heat loss rate [cal/sec-cm]

        """
        self._heat_loss_rate = c_double(value)

    @property
    def heat_transfer_coefficient(self) -> float:
        """Get heat transfer coefficient between the reactor and the surroundings."""
        """
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
        Parameters
        ----------
            value: double, default = 0.0
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
        """Get ambient temperature."""
        """
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
        """Get heat transfer area between the reactor and the surroundings."""
        """
        Returns
        -------
            heat_transfer_area: double
                heat transfer area [cm2]

        """
        return self._heat_transfer_area

    @heat_transfer_area.setter
    def heat_transfer_area(self, value: float = 0.0e0):
        """Set heat transfer area between the reactor and the surroundings."""
        """
        Parameters
        ----------
            value: double, default = 0.0
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


class TransientPSRSetResTimeFixedTemperature(TransientPSR):
    """Transient PSR model with given reactor residence time and reactor temperature."""

    """
    PSR model with given reactor residence time (CONP)
    and reactor temperature (GivenT).

    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor volume and density are varying in this case.
    """

    def __init__(self, guessedmixture: Stream, label: Union[str, None] = None):
        """Create a transient constant pressure perfectly-stirred reactor (PSR)."""
        """
        Create a transient constant pressure perfectly-stirred reactor
        with given reactor residence time (CONP)
        and reactor temperature (GivenT).

        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                inlet name/label

        """
        if label is None:
            label = "Trans_PSR"
        # initialization
        super().__init__(guessedmixture, label)
        # specify PSR "residence time"
        self._problemtype = c_int(2)
        # given temperature
        self._energytype = c_int(2)
        # solver parameters
        self._absolute_tolerance = 1.0e-12
        self._relative_tolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]


class TransientPSRSetResTimeEnergyConservation(TransientPSR):
    """Transient PSR with given reactor residence time and solve energy equation."""

    """
    PSR model with given reactor residence time (CONP)
    and solve energy equation (ENERGY).

    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor volume and density are varying in this case.
    """

    def __init__(self, guessedmixture: Stream, label: Union[str, None] = None):
        """Create a transient constant pressure perfectly-stirred reactor (PSR)."""
        """
        PSR model with given reactor residence time (CONP)
        and solve energy equation (ENERGY).

        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                inlet name/label

        """
        if label is None:
            label = "Trans_PSR"
        # initialization
        super().__init__(guessedmixture, label)
        # specify PSR "residence time"
        self._problemtype = c_int(2)
        # given temperature
        self._energytype = c_int(1)
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        self._heat_transfer_coefficient = 0.0e0
        self._ambient_temperature = 3.0e2
        # external heat transfer [cm2]
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

    @property
    def heat_loss_rate(self) -> float:
        """Get heat loss rate from the reactor to the surroundings."""
        """
        Returns
        -------
            qloss: double
                heat loss rate to the surroundings [cal/sec]

        """
        return self._heat_loss_rate.value

    @heat_loss_rate.setter
    def heat_loss_rate(self, value: float):
        """Set the heat loss rate from the reactor to the surroundings."""
        """
        Parameters
        ----------
            value: double, default = 0.0
                heat loss rate [cal/sec]

        """
        self._heat_loss_rate = c_double(value)

    @property
    def heat_transfer_coefficient(self) -> float:
        """Get heat transfer coefficient between the reactor and the surroundings."""
        """
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
        Parameters
        ----------
            value: double, default = 0.0
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
        """Get ambient temperature."""
        """
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
        """Get heat transfer area between the reactor and the surroundings."""
        """
        Returns
        -------
            heat_transfer_area: double
                heat transfer area [cm2]

        """
        return self._heat_transfer_area

    @heat_transfer_area.setter
    def heat_transfer_area(self, value: float = 0.0e0):
        """Set heat transfer area between the reactor and the surroundings."""
        """
        Parameters
        ----------
            value: double, default = 0.0
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

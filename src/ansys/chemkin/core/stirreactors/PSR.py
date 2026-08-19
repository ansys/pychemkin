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

"""Perfectly stirred reactor (PSR) model."""

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
from ansys.chemkin.core.inlet import Stream
from ansys.chemkin.core.mixture import equilibrium
from ansys.chemkin.core.reactormodel import Keyword
from ansys.chemkin.core.stirreactors.openreactor import OpenReactor
from ansys.chemkin.core.utilities import (
    log_critical_message,
    log_error_message,
    log_info_message,
)
from ansys.chemkin.core.validation import validate_minimum_value


class PerfectlyStirredReactor(OpenReactor):
    """Generic perfectly-stirred reactor model."""

    def __init__(self, guessedmixture: Stream, label: Union[str, None] = None):
        """Instantiate a steady-state PSR object."""
        """Initialize a steady-state constant pressure
        perfectly-stirred reactor (PSR) object.

        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                reactor name

        """
        if label is None:
            label = "PSR"
        # initialization
        super().__init__(guessedmixture, label)
        # reactor pressure [dynes/cm2] (c_double)
        # self._pressure
        # reactor temperature [K] (c_double)
        # self._temperature
        # reactor volume [cm3]
        if guessedmixture._vol > 0.0:
            self._volume = c_double(guessedmixture._vol)
        else:
            self._volume = c_double(0.0)
        # reactor residence time [sec]
        self._residencetime = c_double(0.0)
        # reactive surface area [cm2]
        self._reactivearea = c_double(0.0)
        # simulation time [sec] (not in use)
        self._endtime = c_double(0.0)
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        # single reactor (default) or reactor network
        self.standalone = True
        # check required inputs
        self._numb_requiredinput = 0
        self._requiredlist: list[str] = []
        self._inputcheck: list[str] = []
        # default number of reactors
        self._nreactors = 1
        self._npsrs = c_int(self._nreactors)
        self._ninlets = np.zeros(self._nreactors, dtype=np.int32)
        self._nzones = c_int(0)
        # default reactor type settings
        # Perfectly-Stirred Reactor (PSR) model
        self._reactortype = c_int(2)
        # Steady-State PSR only
        self._solvertype = c_int(self.solver_types.get("SteadyState", 2))
        # default options
        self._problemtype = c_int(self.problem_types.get("SETVOL", 1))
        self._energytype = c_int(self.energy_types.get("ENERGY", 1))
        # set reactor number (single reactor)
        self.ireac = c_int(1)

    @property
    def area(self) -> float:
        """Get reactive surface area."""
        """
        Returns
        -------
            reactive_area: double
                reactive surace area [cm2]

        """
        return self._reactivearea.value

    @area.setter
    def area(self, value: float = 0.0e0):
        """Set reactive surface area (optional)."""
        """
        Parameters
        ----------
            value: double, default = 0.0
                reactive surface area [cm2]

        """
        validate_minimum_value(value, 0.0e0, "reactive surface area must >= 0.")
        self._reactivearea = c_double(value)

    @property
    def volume(self) -> float:
        """Get reactor volume."""
        """
        Returns
        -------
            volume: douoble
                reactor volume [cm3]

        """
        return self._volume.value

    @volume.setter
    def volume(self, value: float):
        """Set reactor volume (required)."""
        """
        Parameters
        ----------
            value: double, default = 1.0
                reactor volume [cm3]

        """
        validate_minimum_value(
            value,
            np.nextafter(0.0, 1.0),
            "reactor volume must > 0.",
        )
        # set reactor volume
        self._volume = c_double(value)
        # set initial mixture volume
        self.reactormixture.volume = value

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
        validate_minimum_value(
            value,
            np.nextafter(0.0, 1.0),
            "residence time must > 0.",
        )
        # set reactor residence time
        self._residencetime = c_double(value)

    def set_inlet_keywords(self) -> int:
        """Set up inlet keywords."""
        """
        Returns
        -------
            error code: integer

        """
        ierr = 0
        # loop over all external inlets into the reactor
        i_inlet = 0
        flowrate_sum = 0.0
        check_total_flowrate = True
        #
        for key, inlet in self.externalinlets.items():
            # get inlet mass flow rate
            if inlet._flowratemode == 0:
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
                    log_error_message(msg)
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
                    # the reactor name and the inlet name are hard coded
                    # in chemkin-CFD-API
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
            log_error_message(msg)
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
            log_error_message(msg)
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
                log_error_message(msg)
                ierr += 12
        return ierr

    def set_reactor_index(self, reactorindex: int):
        """Assign the reactor index/number."""
        """Assign the reactor index/number of the current reactor
        in the reactor network.

        This method should be called by the PSR cluster/network Class/Module
        For single PSR the reactor index is always 1 (default)

        Parameters
        ----------
            reactorindex: integer, default = 1
                reactor index/number in the reactor network

        """
        if reactorindex > 0:
            self.ireac = c_int(reactorindex)
            self.standalone = False

    def set_estimate_conditions(
        self, option: str, guess_temp: Union[float, None] = None
    ):
        """Reset the initial/guessed reactor gas mixture."""
        """Reset the initial/guessed reactor gas mixture properties to
        improve the steady-state solution finding performance.

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
                log_info_message(msg)
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
                log_info_message(msg)
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

    def reset_estimate_temperature(self, temp: float):
        """Reset the estimated reactor gas temperature."""
        """
        Parameters
        ----------
            temp: double
                estimated reactor gas temperature [K]

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
            log_info_message(msg)
            exit()
        else:
            # use the given temperature
            self.reactormixture.temperature = temp
            self.temperature = temp

    def reset_estimate_composition(
        self,
        fraction: Union[npt.NDArray[np.double], list[tuple[str, float]]],
        mode: str = "mole",
    ):
        """Reset the estimated reactor gas composition."""
        """
        Parameters
        ----------
            fraction: 1-D double array, dimension = number_species or PyChemkin recipe
                estimated reactor gas composition
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
            log_info_message(msg)
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
                log_error_message(msg)
            # verify required inputs one by one
            for k in self._requiredlist:
                if k not in self._inputcheck:
                    ierr += 1
                    msg = [Color.PURPLE, "missing required input", k, Color.END]
                    log_error_message(msg)
            return ierr

    def set_ss_solver_keywords(self):
        """Add steady-state solver parameter keywords to the keyword list."""
        # steady-state solver parameter given
        if len(self.ss_solverkeywords) > 0:
            #
            for k, v in self.ss_solverkeywords.items():
                self.setkeyword(k, v)

    def cluster_process_keywords(self) -> int:
        """Process input keywords for the reactor model in cluster network mode."""
        """
        Returns
        -------
            Error code: integer

        """
        ierr = self.__process_keywords()
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
            log_error_message(msg)
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
                log_error_message(msg)
                return ierr
        else:
            # set up inlets
            err_inputs = self.set_inlet_keywords()
            if err_inputs != 0:
                print(f"error code = {err_inputs}")
            ierr += err_inputs
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
                log_error_message(msg)
                return ierrc
            # heat transfer (use additional keywords)
            # solver parameters (use additional keywords)
            # output controls (use additional keywords)
            # ROP (use additional keywords)
            # sensitivity (use additional keywords)

        if ierr == 0 and self._numbprofiles > 0 and self.standalone:
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
        if ierr == 0 and self.has_surface_chemistry:
            # set surface related keywords
            self.set_surface_chemistry_keywords()
        if ierr == 0:
            # set additional keywords
            self.set_ss_solver_keywords()
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
                log_error_message(msg)
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
                "failed to initialize the perfectly-stirred reactor model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            log_error_message(msg)
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
            msg = [
                Color.RED,
                "full keyword option not available for PSR models.",
                Color.END,
            ]
            log_critical_message(msg)
            ret_val = 100
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
        log_info_message([Color.YELLOW, "processing keywords complete", Color.END])

        # run reactor model
        msg = [Color.YELLOW, "running reactor simulation ...", Color.END]
        log_info_message(msg)
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
            log_info_message(msg)
        else:
            msg.insert(0, Color.RED)
            log_critical_message(msg)

        return ret_val

    def process_solution(self) -> Stream:
        """Post-process solution."""
        """Post-process solution to extract the raw solution variable data
        package the steady-state solution into a mixture object.

        Returns
        -------
            smixture: Stream object
                gas stream representing the steady-state solution

        """
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
        # create a Stream object to hold the mixture properties of current solution
        smixture = copy.deepcopy(self.reactormixture)
        # create a species mass fraction array to hold the steady-state solution
        frac = np.empty(self.numbspecies, dtype=np.double)
        # get raw solution data
        temp = c_double(0.0)
        pres = c_double(0.0)
        ierr = chemkin_wrapper.chemkin.KINAll0D_GetSolution(temp, pres, frac)
        if ierr != 0:
            msg = [
                Color.RED,
                "failed to fetch the raw solution data from memory,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            log_critical_message(msg)
            exit()
        # set solution mixture
        self._numbsolutionpoints = self._nreactors
        # steady-state pressure solution [dynes/cm2]
        smixture.pressure = pres.value
        # steady-state temperature solution [K]
        smixture.temperature = temp.value
        # set mixture composition
        if self._speciesmode == "mass":
            # mass fractions
            smixture.y = frac
        else:
            # mole fractions
            smixture.x = frac
        # get surface coverage
        if self.has_surface_chemistry:
            self.process_surface_solution(smixture)

        # get reactor outlet mass flow rate [g/sec]
        exitmassflowrate = c_double(0.0)
        ierr = chemkin_wrapper.chemkin.KINAll0D_GetExitMassFlowRate(exitmassflowrate)
        if ierr == 0:
            smixture.mass_flowrate = max(0.0, exitmassflowrate.value)
        else:
            smixture.mass_flowrate = 0.0
            msg = [
                Color.RED,
                "failed to get the total outlet mass flow rate,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            log_critical_message(msg)
            exit()
        # clean up
        del frac
        # return the solution mixture
        return smixture

    def process_surface_solution(self, smixture: Stream):
        """Post-process surface solution."""
        """Post-process surface solution to extract the raw solution variable data
        package the steady-state solution into a mixture object.

        Parameters
        ----------
            smixture: Stream object
                gas stream representing the steady-state solution

        """
        # get material information
        mat_names = self.get_material_names()
        # get the site and the bulk species sizes of the mechanism
        max_sites = self.get_total_site_species()
        max_bulks = self.get_total_bulk_species()
        # set up holding arrays
        if max_sites > 0:
            s_frac = np.empty(max_sites, dtype=np.double)
            ierr_s = chemkin_wrapper.chemkin.KINAll0D_GetSurfaceSolution(s_frac)
            if ierr_s != 0:
                msg = [
                    Color.MAGENTA,
                    "failed to fetch the raw surface solution data from memory,",
                    "error code =",
                    str(ierr_s),
                    Color.END,
                ]
                log_error_message(msg)
        else:
            ierr_s = -1
            s_frac = np.zeros(1, dtype=np.double)
        if max_bulks > 0:
            b_rate = np.empty(max_bulks, dtype=np.double)
            ierr_b = chemkin_wrapper.chemkin.KINAll0D_GetBulkGrowthRate(b_rate)
            if ierr_b != 0:
                msg = [
                    Color.RED,
                    "failed to fetch the bulk growth rate solution data from memory,",
                    "error code =",
                    str(ierr_b),
                    Color.END,
                ]
                log_error_message(msg)
        else:
            ierr_b = -1
            b_rate = np.zeros(1, dtype=np.double)
        # store surface coverage in the solution mixture
        # loop over all materials of the surface mechanism
        s_count = 0
        b_count = 0
        for mname in mat_names:
            # set material temperature
            if self.surface_chemistry.check_material_temperature(mname) == 0:
                # material temperature equals to the gas temperature
                surf_temp = smixture.temperature
            else:
                surf_temp = self.surface_chemistry.get_material_temperature(mname)
            smixture.surface_chemistry.set_material_temperature(mname, surf_temp)
            # number of site species on this material
            n_sites = self.surface_chemistry.number_site_species(mname)
            # number of bulk species on this material
            n_bulks = self.surface_chemistry.number_bulk_species(mname)
            if ierr_s == 0:
                is_start = s_count
                is_end = is_start + n_sites
                smixture.surface_chemistry.set_site_frac(mname, s_frac[is_start:is_end])
            s_count += n_sites

            if ierr_b == 0:
                ib_start = b_count
                ib_end = ib_start + n_bulks
                smixture.surface_chemistry.set_bulk_growth_rates(
                    mname, b_rate[ib_start:ib_end]
                )
            b_count += n_bulks
        del s_frac, b_rate


class PSRSetResTimeEnergyConservation(PerfectlyStirredReactor):
    """PSR model with given reactor residence time and solve energy equation."""

    """
    PSR model with given reactor residence time (CONP)
    and solve energy equation (ENERGY).

    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor volume and density are varying in this case.
    """

    def __init__(self, guessedmixture: Stream, label: Union[str, None] = None):
        """Create a steady-state constant pressure perfectly-stirred reactor (PSR)."""
        """
        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                inlet name/label

        """
        if label is None:
            label = "PSR"
        # initialization
        super().__init__(guessedmixture, label)
        # specify residence
        self._problemtype = c_int(self.problem_types.get("SETTAU", 2))
        self._energytype = c_int(self.energy_types.get("ENERGY", 1))
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        self._heat_transfer_coefficient = 0.0e0
        self._ambient_temperature = 3.0e2
        # external heat transfer [cm2]
        self._heat_transfer_area = 0.0e0

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
        validate_minimum_value(
            value,
            0.0e0,
            "heat transfer coefficient must > 0.",
        )
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
        validate_minimum_value(
            value,
            np.nextafter(0.0, 1.0),
            "ambient temperature must > 0.",
        )
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
        validate_minimum_value(value, 0.0e0, "heat transfer area must >= 0.")
        self._heat_transfer_area = value
        # set the corresponding keyword
        self.setkeyword(key="AREAQ", value=value)


class PSRSetVolumeEnergyConservation(PerfectlyStirredReactor):
    """PSR model with given reactor volume and solve energy equation."""

    """
    PSR model with given reactor volume (CONV)
    and solve energy equation (ENERGY).

    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor residence time and density are varying in this case.
    """

    def __init__(self, guessedmixture: Stream, label: Union[str, None] = None):
        """Create a steady-state constant pressure perfectly-stirred reactor (PSR)."""
        """
        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                inlet name/label

        """
        if label is None:
            label = "PSR"
        # initialization
        super().__init__(guessedmixture, label)
        # specify volume
        self._problemtype = c_int(self.problem_types.get("SETVOL", 1))
        self._energytype = c_int(self.energy_types.get("ENERGY", 1))
        # heat transfer parameters
        self._heat_loss_rate = c_double(0.0e0)
        self._heat_transfer_coefficient = 0.0e0
        self._ambient_temperature = 3.0e2
        # external heat transfer [cm2]
        self._heat_transfer_area = 0.0e0

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
        validate_minimum_value(
            value,
            0.0e0,
            "heat transfer coefficient must > 0.",
        )
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
        validate_minimum_value(
            value,
            np.nextafter(0.0, 1.0),
            "ambient temperature must > 0.",
        )
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
        validate_minimum_value(value, 0.0e0, "heat transfer area must >= 0.")
        self._heat_transfer_area = value
        # set the corresponding keyword
        self.setkeyword(key="AREAQ", value=value)


class PSRSetResTimeFixedTemperature(PerfectlyStirredReactor):
    """PSR model with given reactor residence time and reactor temperature."""

    """
    PSR model with given reactor residence time (CONP)
    and reactor temperature (GivenT).

    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor volume and density are varying in this case.
    """

    def __init__(self, guessedmixture: Stream, label: Union[str, None] = None):
        """Create a steady-state constant pressure perfectly-stirred reactor (PSR)."""
        """
        Create a steady-state constant pressure perfectly-stirred reactor (PSR).

        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                inlet name/label

        """
        if label is None:
            label = "PSR"
        # initialization
        super().__init__(guessedmixture, label)
        # specify residence time
        self._problemtype = c_int(self.problem_types.get("SETTAU", 2))
        self._energytype = c_int(self.energy_types.get("GivenT", 2))


class PSRSetVolumeFixedTemperature(PerfectlyStirredReactor):
    """PSR model with given reactor volume and reactor temperature."""

    """
    PSR model with given reactor volume (CONV)
    and reactor temperature (GivenT).

    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor residence time and density are varying in this case.
    """

    def __init__(self, guessedmixture: Stream, label: Union[str, None] = None):
        """Create a steady-state constant pressure perfectly-stirred reactor (PSR)."""
        """
        Parameters
        ----------
            guessedmixture: Mixture object
                a mixture representing the estimated/guessed gas properties of the PSR
            label: string, optional
                inlet name/label

        """
        if label is None:
            label = "PSR"
        # initialization
        super().__init__(guessedmixture, label)
        # specify volume
        self._problemtype = c_int(self.problem_types.get("SETVOL", 1))
        self._energytype = c_int(self.energy_types.get("GivenT", 2))

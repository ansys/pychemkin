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

"""Chemkin reactor inlet/stream utilities."""

import copy
from typing import Union

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core.chemistry import Chemistry
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.constants import P_ATM
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.mixture import (
    Mixture,
    cal_mixture_temperature_from_enthalpy,
    compare_mixtures,
)


class Stream(Mixture):
    """Generic inlet stream object for Chemkin open reactor models."""

    # The "Stream" class is an extension of the "Mixture" class

    def __init__(self, chem, label: Union[str, None] = None):
        """Initialize an inlet object with a given chemistry set."""
        """Initialize an inlet object with a given chemistry set
        for open reactor models. The Stream object contains
        pressure, tempersture, flow rate, and gas composition
        of the inlet gas mixture.

        Parameters
        ----------
            chem: Chemistry object
            label: string, optional
                inlet name

        """
        super().__init__(chem)
        # 0=mass flow rate/1=volumetric flow rate/2=velocity/3=SCCM
        # flag
        self._flowratemode = -1  # not given
        self._inletflowrate = [0.0] * 4
        # types of flow rate allowed
        # mass flow rate FLRT [g/sec]
        self._massflowrate = 0.0
        # volumetric flow rate VDOT [cm3/sec]
        self._volflowrate = 0.0
        # gas velocity VEL [cm/sec] for plug flow reactor model
        self._velocity = 0.0
        # standard (198.15K, 1atm) cubic centimeters per minute SCCM [standard cm3/min]
        self._sccm = 0.0
        # inlet velocity gradient [1/sec] (for premixed, oppdif, amd spin)
        self._velgrad = 0.0
        # flow area (for velocity in plug flow reactor model)
        self._haveflowarea = False
        # cross-sectional flow area [cm2]
        self._flowarea = 1.0
        # flag for using flow rate profile
        # -1: no inlet flow rate profile is used
        # 0: mass flow rate profile "FPRO"
        # 1: volumetric flow rate profile "VDOTPRO"
        # 2: velocity profile "VELPRO"
        # 3: SCCM profile "SCCMPRO"
        self._use_flow_rate_profile = False
        self._inlet_flow_rate_profile_mode = -1
        self.flow_rate_profile_keys: dict[int, str] = {
            int(0): "FPRO",
            int(1): "VDOTPRO",
            int(2): "VELPRO",
            int(3): "SCCMPRO",
        }
        self.flowrate_profile: dict[
            str, tuple[npt.NDArray[np.double], npt.NDArray[np.double]]
        ] = {}
        # flag for inlet stream temperature profile "TINPRO"
        self._use_inlet_temperature_profile = False
        # set inlet label
        if label is None:
            label = "inlet"
        self._label = label

    def convert_to_mass_flowrate(self) -> float:
        """Convert different types of flow rate value to mass flow rate."""
        """
        Convert different types of flow rate value to mass flow rate.

        Returns
        -------
            mrate: double
                mass flow rate [g/sec]

        """
        #
        if self._flowratemode == 1:
            # volumetric flow rate
            # get inlet gas mixture density
            mrate = self.rho * self._volflowrate
            return mrate
        elif self._flowratemode == 2:
            # velocity
            if self._haveflowarea:
                mrate = self.rho * self._flowarea * self._velocity
                return mrate
            else:
                # no flow area
                msg = [
                    Color.PURPLE,
                    "flow area is not given for this inlet.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

        elif self._flowratemode == 3:
            # SCCM
            chemid = self._chemset_index.value
            # set standard condition
            p = P_ATM  # [atm]
            t = 298.15  # [K]
            # set mass fractions
            frac = copy.deepcopy(self.y)
            # molecular masses
            wt = self._wt
            # get gas density at the standard condition
            standard_den = Mixture.density(chemid, p, t, frac, wt, mode="mass")
            mrate = standard_den * self._sccm / 60.0
            del frac
            return mrate
        else:
            msg = [Color.PURPLE, "unknown flow rate units.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def convert_to_vol_flowrate(self) -> float:
        """Convert different types of flow rate value to volumetric flow rate."""
        """
        Convert different types of flow rate value to volumetric flow rate.

        Returns
        -------
            vrate: double
                volmetric flow rate [cm3/sec]

        """
        #
        if self._flowratemode == 0:
            # mass flow rate
            # get inlet gas mixture density
            vrate = self._massflowrate / self.rho
            return vrate
        elif self._flowratemode == 2:
            # velocity
            if self._haveflowarea:
                vrate = self._flowarea * self._velocity
                return vrate
            else:
                # no flow area
                msg = [
                    Color.PURPLE,
                    "flow area is not given for this inlet.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

        elif self._flowratemode == 3:
            # SCCM
            chemid = self._chemset_index.value
            # set standard condition
            p = P_ATM  # [atm]
            t = 298.15  # [K]
            # set mass fractions
            frac = copy.deepcopy(self.y)
            # molecular masses
            wt = self._wt
            # get gas density at the standard condition
            standard_den = Mixture.density(chemid, p, t, frac, wt, mode="mass")
            mrate = standard_den * self._sccm / 60.0
            vrate = mrate / self.rho
            del frac
            return vrate
        else:
            msg = [Color.PURPLE, "unknown flow rate units.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def convert_to_sccm(self) -> float:
        """Convert different types of flow rate value to SCCM."""
        """
        Convert different types of flow rate value to SCCM.

        Returns
        -------
            sccm: double
                volumetric flow rate in SCCM [standard cm3/min]

        """
        #
        chemid = self._chemset_index.value
        # set standard condition
        p = P_ATM  # [atm]
        t = 298.15  # [K]
        # set mass fractions
        frac = copy.deepcopy(self.y)
        # molecular masses
        wt = self._wt
        # get gas density at the standard condition
        standard_den = Mixture.density(chemid, p, t, frac, wt, mode="mass")
        del frac
        #
        if self._flowratemode == 0:
            # mass flow rate
            sccm = self._massflowrate / standard_den * 60.0
            return sccm
        elif self._flowratemode == 1:
            # volumetric flow rate
            # get inlet gas mixture density
            mrate = self.rho * self._volflowrate
            sccm = mrate / standard_den * 60.0
            return sccm
        elif self._flowratemode == 2:
            # velocity
            if self._haveflowarea:
                mrate = self.rho * self._flowarea * self._velocity
                sccm = mrate / standard_den * 60.0
                return sccm
            else:
                # no flow area
                msg = [
                    Color.PURPLE,
                    "flow area is not given for this inlet.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        else:
            msg = [Color.PURPLE, "unknown flow rate units.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def check_flow_rate_mode(self) -> int:
        """Get the type of flow rate is provided to this inlet."""
        """
        Get the type of flow rate unit is provided to this inlet.

        Returns
        -------
            flow_rate_mode: integer
                type of flow rate (unit) is assigned
        """
        return self._flowratemode

    def notify_inlet_profile(self, mode: str):
        """Display warning that flow rate profile is specified."""
        """
        Display warning message to instruct the user to use the flow rate
        profile methods when the inlet flow rate profile is specified for
        this inlet.

        Parameters
        ----------
            mode: string
                type of flow rate profile unit
        """
        msg = [
            Color.YELLOW,
            "Flow rate profile is assigned to this inlet",
            self.label,
            ".\n",
            f"Please use the 'get_{mode}_from_profile(time)' method",
            f"to get the {mode} at the given time [sec].",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)

    @property
    def flowarea(self) -> float:
        """Get inlet flow area."""
        """
        Get inlet flow area.

        Returns
        -------
            flowarea: double
                cross-sectional flow area [cm2]

        """
        if self._haveflowarea:
            return self._flowarea
        else:
            msg = [Color.PURPLE, "flow area is not given for this inlet.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @flowarea.setter
    def flowarea(self, farea: float):
        """Set inlet cross-sectional flow area."""
        """
        Set inlet cross-sectional flow area.

        Parameters
        ----------
            farea: double
                cross-sectional flow area [cm2]

        """
        if farea <= 0.0:
            msg = [Color.PURPLE, "invalid flow area value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self._haveflowarea = True
        self._flowarea = farea

    @property
    def mass_flowrate(self) -> float:
        """Get inlet mass flow rate."""
        """
        Get inlet mass flow rate.

        Returns
        -------
            mflowrate: double
                mass flow rate [g/sec]

        """
        if self._use_flow_rate_profile:
            self.notify_inlet_profile("mass_flowrate")

        if self._flowratemode == 0:
            return self._massflowrate
        else:
            return self.convert_to_mass_flowrate()

    @mass_flowrate.setter
    def mass_flowrate(self, mflowrate: float):
        """Set inlet mass flow rate."""
        """
        Set inlet mass flow rate.

        Parameters
        ----------
            mflowrate: double
                mass flow rate [g/sec]

        """
        if mflowrate <= 0.0:
            msg = [Color.PURPLE, "invalid mass flow rate value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # reset the flow rates
        self._volflowrate = 0.0
        self._velocity = 0.0
        self._sccm = 0.0
        # set flow rate mode to mass flow rate
        self._flowratemode = 0
        self._inletflowrate[self._flowratemode] = mflowrate
        self._massflowrate = mflowrate

    @property
    def vol_flowrate(self) -> float:
        """Get inlet volumetric flow rate."""
        """
        Get inlet volumetric flow rate.

        Returns
        -------
            vflowrate: double
                volumetric flow rate [cm3/sec]

        """
        if self._use_flow_rate_profile:
            self.notify_inlet_profile("vol_flowrate")

        if self._flowratemode == 1:
            return self._volflowrate
        else:
            return self.convert_to_vol_flowrate()

    @vol_flowrate.setter
    def vol_flowrate(self, vflowrate: float):
        """Set inlet volumetric flow rate."""
        """
        Set inlet volumetric flow rate.

        Parameters
        ----------
            vflowrate: double
                volumetric flow rate [cm3/sec]

        """
        if vflowrate <= 0.0:
            msg = [Color.PURPLE, "invalid volumetric flow rate value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # reset the flow rates
        self._massflowrate = 0.0
        self._velocity = 0.0
        self._sccm = 0.0
        # set flow rate mode to volumetric flow rate
        self._flowratemode = 1
        self._inletflowrate[self._flowratemode] = vflowrate
        self._volflowrate = vflowrate

    @property
    def sccm(self) -> float:
        """Get inlet SCCM volumetric flow rate."""
        """
        Get inlet SCCM volumetric flow rate.

        Returns
        -------
            vflowrate: double
                SCCM volumetric flow rate [standard cm3/min]

        """
        if self._use_flow_rate_profile:
            self.notify_inlet_profile("sccm_flowrate")
            return self._sccm
        if self._flowratemode == 3:
            return self._sccm
        else:
            return self.convert_to_sccm()

    @sccm.setter
    def sccm(self, vflowrate: float):
        """Set inlet volumetric flow rate in SCCM."""
        """
        Set inlet volumetric flow rate in SCCM.

        Parameters
        ----------
            vflowrate: double
                SCCM volumetric flow rate [standard cm3/min]

        """
        if vflowrate <= 0.0:
            msg = [Color.PURPLE, "invalid SCCM volumetric flow rate value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # reset the flow rates
        self._massflowrate = 0.0
        self._volflowrate = 0.0
        self._velocity = 0.0
        # set flow rate mode to volumetric flow rate
        self._flowratemode = 3
        self._inletflowrate[self._flowratemode] = vflowrate
        self._sccm = vflowrate

    @property
    def velocity(self) -> float:
        """Get inlet gas velocity."""
        """
        Get inlet gas velocity.

        Returns
        -------
            vel: double
                velocity [cm/sec]

        """
        if self._use_flow_rate_profile:
            self.notify_inlet_profile("velocity")

        if self._flowratemode == 2:
            return self._velocity
        else:
            if self._haveflowarea:
                # have flow area
                if self._flowratemode == 1:
                    vrate = self._volflowrate
                else:
                    vrate = self.convert_to_vol_flowrate()
                # convert volumetric flow rate to velocity
                return vrate / self._flowarea
            else:
                # flow area not defined
                msg = [
                    Color.PURPLE,
                    "flow area is not given for this inlet.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    @velocity.setter
    def velocity(self, vel: float):
        """Set inlet velocity."""
        """
        Set inlet velocity.

        Parameters
        ----------
            vel: velocity [cm/sec]

        """
        if vel <= 0.0:
            msg = [Color.PURPLE, "invalid inlet velocity value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # reset the flow rates
        self._massflowrate = 0.0
        self._volflowrate = 0.0
        self._sccm = 0.0
        # set flow rate mode to velocity
        self._flowratemode = 2
        self._inletflowrate[self._flowratemode] = vel
        self._velocity = vel

    @property
    def velocity_gradient(self) -> float:
        """Get inlet gas axial velocity gradient."""
        """Get inlet gas axial velocity gradient (for premixed, oppdif, and spin)
        or radial velocity spreading rate (v_r/r) at the inlet.

        Returns
        -------
            velgrad: double
                velocity gradient [1/sec]

        """
        return self._velgrad

    @velocity_gradient.setter
    def velocity_gradient(self, velgrad: float):
        """Set inlet axial velocity gradient."""
        """
        Set inlet axial velocity gradient.

        Parameters
        ----------
            velgrad: double
                axial velocity gradient [1/sec]
        :return: None

        """
        if velgrad <= 0.0:
            msg = [
                Color.PURPLE,
                "invalid inlet radial velocity spreading rate value.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # set velocity gradient
        self._velgrad = velgrad

    @property
    def label(self) -> str:
        """Get the label of the Stream."""
        """
        Get the label of the Stream.

        Returns
        -------
            label: string
                label of the Stream

        """
        return self._label

    @label.setter
    def label(self, name: str):
        """Set the label of the Stream."""
        """
        Parameters
        ----------
            name: string
                label of the Stream

        """
        self._label = name

    def profile_out_of_range(self, time: float, max: float, min: float):
        """Report value out of the profile data range."""
        """
        Issue an error message when the given time value is out of the
        time profile data range.

        Parameters
        ----------
            time: double
                the given time value [sec]
            max: double
                maximum time value in the profile data
            min: double
                minimum time value in the profile data
        """
        msg = [
            Color.PURPLE,
            "the given time value",
            str(time),
            "is out of the profile data range of [",
            str(min),
            ",",
            str(max),
            "].",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)

    def set_mass_flowrate_profile(
        self, x: npt.NDArray[np.double], flow_rate: npt.NDArray[np.double]
    ):
        """Specify inlet mass flow rate profile."""
        """
        Specify inlet mass flow rate profile for transient reactor models.

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm or sec]
            flow_rate: 1D double array
                mass flow arte value of the profile data [g/sec]
        """
        # reset the flow rates
        self._volflowrate = 0.0
        self._velocity = 0.0
        self._sccm = 0.0
        # set flow rate mode to mass flow rate
        self._flowratemode = 0
        # set inlet profile flag
        self._use_flow_rate_profile = True
        self._inletflowrate[self._flowratemode] = flow_rate[0]
        self._massflowrate = flow_rate[0]
        # profile keywords will be set by the set_profile
        if len(self.flowrate_profile.keys()) > 0:
            msg = [
                Color.YELLOW,
                "existing flow rate profile will be deleted.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.flowrate_profile.clear()
        #
        keyword = self.flow_rate_profile_keys.get(self._flowratemode, "FPRO")
        self.flowrate_profile[keyword] = (x, flow_rate)

    def get_mass_flowrate_profile_data(self, time: float = 0.0) -> float:
        """Get the mass flow rate at the given time from the profile data."""
        """
        Return the inlet mass flow rate at the given time value
        from the inlet mass flow rate profile data by interpolation.

        Parameters
        ----------
            time: double
                time [sec]

        Returns
        -------
            mass_flowrate: double
                the interpolated inlet mass flow rate [g/sec] at the given time
        """
        # profile keywords will be set by the set_profile
        if self._use_flow_rate_profile and self._flowratemode == 0:
            # mass flow rate profile data does not exist
            msg = [
                Color.PURPLE,
                "cannot find the inlet mass flow rate profile data,",
                "please verify the profile data set.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 0.0
        this_x, this_values = self.flowrate_profile.get("FPRO", ([0.0], [0.0]))
        if len(this_x) <= 1:
            msg = [
                Color.PURPLE,
                "cannot find the mass flow rate profile data.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 0.0
        if time < np.min(this_x) or time > np.max(this_x):
            self.profile_out_of_range(time, np.max(this_x), np.min(this_x))
        values = np.interp(time, this_x, this_values)
        mass_flowrate = values[0]
        return mass_flowrate


# stream utilities
def clone_stream(source: Stream, target: Stream):
    """Copy the properties of the source Stream to the target Stream."""
    """
    Duplicate the properties of the source Stream to the target Stream.

    Parameters
    ----------
        source: Stream object
            the "source" Stream to be cloned
        target: Stream object
            the "target" Stream to get new properties

    """
    # check Chemistry set
    if source.chemid == target.chemid:
        # clone the Stream properties
        target.temperature = source.temperature
        target.pressure = source.pressure
        target.mass_flowrate = source.mass_flowrate
        target.x = source.x
    else:
        msg = [
            Color.PURPLE,
            "the streams have different Chemistry Sets.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()


def compare_streams(
    stream_a: Stream,
    stream_b: Stream,
    atol: float = 1.0e-10,
    rtol: float = 1.0e-3,
    mode: str = "mass",
) -> tuple[bool, float, float]:
    """Compare properties of stream B against those of stream A."""
    """Compare properties of stream B against those of stream A.
    The stream properties include mixture properties such as pressure [atm],
    temperature [K], and species mass/mole fractions and the stream mass flow rate.
    When the differences in the property values satisfy both the absolute and
    the relative tolerances, this method will return "True", that is,
    stream B is essentially identical to stream A; otherwise,
    "False" will be returned.

    Parameters
    ----------
        stream_a: Stream object
            mixture A, the target stream
        stream_b: Stream object
            stream B, the sample stream
        atol: double, default = 1.0e-10
            the absolute tolerance for the max property differences
        rtol: double, default = 1.0e-3
            the relative tolerance for the max property differences
        mode: string {"mass", "mole"}, default = "mass"
            compare species "mass" or "mole" fractions

    Returns
    -------
        issame: boolean
            the equivalence of the two mixtures
        atol_max: double
            the max absolute difference value
        rtol_max: double
            the max relative difference value

    """
    # check mixtures first
    issame, diff_max, var_max = compare_mixtures(
        stream_a, stream_b, atol=atol, rtol=rtol, mode=mode
    )
    # compare stream mass flow rate
    mflr_diff = abs(stream_a.mass_flowrate - stream_b.mass_flowrate)
    # find relative difference
    mflr_var = mflr_diff / stream_a.mass_flowrate
    # check tolerances
    flowsame = mflr_diff <= atol
    flowsame = flowsame or mflr_var <= rtol
    if not flowsame:
        print(f"=== mass flow rate difference: {mflr_diff}   {mflr_var}\n")
    #
    diff_max = max(diff_max, mflr_diff)
    var_max = max(var_max, mflr_var)
    issame = issame or flowsame
    #
    return issame, diff_max, var_max


def adiabatic_mixing_streams(stream_a: Stream, stream_b: Stream) -> Stream:
    """Create a new Stream object by mixing two streams adiabatically."""
    """Create a new Stream object by mixing two streams adiabatically.
    The enthalpy of the final stream is the sum of the enthalpies of
    the two streams, and so is the final mass flow rate.

    Parameters
    ----------
        stream_a: Stream object
            mixture A, the target stream
        stream_b: Stream object
            stream B, the sample stream

    Returns
    -------
        final_stream: Stream object
            the final stream from combining stream_a and stream_b

    """
    if isinstance(stream_a, Stream) and isinstance(stream_b, Stream):
        if stream_a.chemid == stream_b.chemid:
            final_stream = copy.deepcopy(stream_a)
        else:
            msg = [
                Color.PURPLE,
                "the streams have different Chemistry Sets.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
    else:
        msg = [
            Color.PURPLE,
            "the streams must be Stream objects.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()

    # number of gas species
    numb_species = stream_a._kk
    # initialization
    massfrac = np.zeros(numb_species, dtype=np.double)
    mix_h = 0.0e0
    speciesfrac_sum = 0.0e0
    # total mass flow rate
    total_mass_flow_rate = stream_a.mass_flowrate + stream_b.mass_flowrate
    final_stream.mass_flowrate = total_mass_flow_rate
    # compute the species mass fractions in the combined stream
    for k in range(numb_species):
        massfrac[k] = stream_a.y[k] * stream_a.mass_flowrate
        massfrac[k] += stream_b.y[k] * stream_b.mass_flowrate
        massfrac[k] /= total_mass_flow_rate
        speciesfrac_sum += massfrac[k]

    if speciesfrac_sum > 0.0e0:
        for k in range(numb_species):
            massfrac[k] /= speciesfrac_sum
        final_stream.y = massfrac

    # compute the final stream's enthalpy flux [ergs/g]
    mix_h = stream_a.hml() / stream_a.wtm * stream_a.mass_flowrate
    mix_h += stream_b.hml() / stream_b.wtm * stream_b.mass_flowrate
    mix_h /= total_mass_flow_rate
    # use mean molecular weight of the gas mixture of the combined stream
    # to convert the enthalpy from [erg/g] to [erg/mol]
    mix_h *= final_stream.wtm
    # compute temperature of the final mixture from the mixture enthalpy
    # set the guessed temperature
    t_guessed = 0.0e0
    ierror = cal_mixture_temperature_from_enthalpy(
        mixture=final_stream, h_mixture=mix_h, guesstemperature=t_guessed
    )
    if ierror != 0:
        msg = [
            Color.PURPLE,
            "failed to compute the final stream temperature,",
            "error code =",
            str(ierror),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    print(f"final stream temperature = {final_stream.temperature} [K]")
    return final_stream


def create_stream_from_mixture(
    chem: Chemistry, mixture: Mixture, flow_rate: float = 0.0, mode: str = "mass"
) -> Stream:
    """Create a new Stream object from the given Mixture object."""
    """
    Create a new Stream object from the given Mixture object.

    Parameters
    ----------
        chem: Chemistry object
            the Chemistry Set used to instantiate the mixture
        mixture: Mixture object
            the Mixture whose properties will be used to set up the new Stream
        flow_rate: double, >= 0.0, default = 0.0
            the flow rate/velocity of the new Stream, [g/sec], [cm3/sec],
            [cm/sec], [standard cm3/minute]
        mode: string, {"mass", "vol", "vel", "sccm"}
            the type of flow rate data is given by the flow_rate parameter

    Returns
    -------
        new_stream: Stream object
            the new Stream based on the given Mixture

    """
    # create the Stream object
    new_stream = Stream(chem)
    # transfer the Mixture properties to the new Stream object
    if isinstance(mixture, Mixture):
        new_stream.pressure = mixture.pressure
        new_stream.temperature = mixture.temperature
        new_stream.x = mixture.x
    else:
        msg = [
            Color.PURPLE,
            "the second parameter must be a Mixture object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # set the flow rate of the new Stream object
    if flow_rate > 0.0:
        if mode.lower() == "vol":
            new_stream.vol_flowrate = flow_rate
        elif mode.lower() == "mass":
            new_stream.mass_flowrate = flow_rate
        elif mode.lower() == "sccm":
            new_stream.sccm = flow_rate
        else:
            new_stream.velocity = flow_rate
            msg = [
                Color.MAGENTA,
                'remember to provide "flow area" to this Stream.',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
    else:
        msg = [Color.PURPLE, "flow rate must > 0.", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    return new_stream

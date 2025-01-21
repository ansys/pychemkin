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
    Chemkin reactor inlet utilities.
"""

import copy

from chemkin.color import Color
from chemkin.constants import Patm
from chemkin.logger import logger
from chemkin.mixture import Mixture


class Inlet(Mixture):
    """
    Generic inlet stream consists of the gas species defined in the given chemistry set
    for Chemkin open reactor models
    """

    # The "Inlet" class is an extension of the "Mixture" class

    def __init__(self, chem, label: str | None = None):
        """
        Initialize an inlet object with a given chemistry set for open reactor models

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
        self._massflowrate = 0.0  # mass flow rate FLRT [g/sec]
        self._volflowrate = 0.0  # volumetric flow rate VDOT [cm3/sec]
        self._velocity = 0.0  # gas velocity VEL [cm/sec] for plug flow reactor model
        self._SCCM = 0.0  # standard (198.15K, 1atm) cubic centimeters per minute SCCM [standard cm3/min]
        # inlet velocity gradient [1/sec] (for premixed, oppdif, amd spin)
        self._velgrad = 0.0
        # flow area (for velocity in plug flow reactor model)
        self._haveflowarea = False
        # cross-sectional flow area [cm2]
        self._flowarea = 1.0
        # set inlet label
        if label is None:
            label = "inlet"
        self.label = label

    def convert_to_mass_flowrate(self) -> float:
        """
        convert different types of flow rate value to mass flow rate

        Returns
        -------
            mrate: double
                mass flow rate [g/sec]
        """
        #
        if self._flowratemode == 1:
            # volumetric flow rate
            # get inlet gas mixture density
            mrate = self.RHO * self._volflowrate
            return mrate
        elif self._flowratemode == 2:
            # velocity
            if self._haveflowarea:
                mrate = self.RHO * self._flowarea * self._velocity
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
            chemID = self._chemset_index.value
            # set standard condition
            p = Patm  # [atm]
            t = 298.15  # [K]
            # set mass fractions
            frac = copy.deepcopy(self.Y)
            # molecular masses
            wt = self._WT
            # get gas density at the standard condition
            standard_den = Mixture.density(chemID, p, t, frac, wt, mode="mass")
            mrate = standard_den * self._SCCM / 60.0
            del frac
            return mrate
        else:
            msg = [Color.PURPLE, "unknown flow rate units.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def convert_to_vol_flowrate(self) -> float:
        """
        convert different types of flow rate value to volumetric flow rate

        Returns
        -------
            vrate: double
                volmetric flow rate [cm3/sec]
        """
        #
        if self._flowratemode == 0:
            # mass flow rate
            # get inlet gas mixture density
            vrate = self._massflowrate / self.RHO
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
            chemID = self._chemset_index.value
            # set standard condition
            p = Patm  # [atm]
            t = 298.15  # [K]
            # set mass fractions
            frac = copy.deepcopy(self.Y)
            # molecular masses
            wt = self._WT
            # get gas density at the standard condition
            standard_den = Mixture.density(chemID, p, t, frac, wt, mode="mass")
            mrate = standard_den * self._SCCM / 60.0
            vrate = mrate / self.RHO
            del frac
            return vrate
        else:
            msg = [Color.PURPLE, "unknown flow rate units.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def convert_to_SCCM(self) -> float:
        """
        convert different types of flow rate value to SCCM

        Returns
        -------
            sccm: double
                volumetric flow rate in SCCM [standard cm3/min]
        """
        #
        chemID = self._chemset_index.value
        # set standard condition
        p = Patm  # [atm]
        t = 298.15  # [K]
        # set mass fractions
        frac = copy.deepcopy(self.Y)
        # molecular masses
        wt = self._WT
        # get gas density at the standard condition
        standard_den = Mixture.density(chemID, p, t, frac, wt, mode="mass")
        del frac
        #
        if self._flowratemode == 0:
            # mass flow rate
            sccm = self._massflowrate / standard_den * 60.0
            return sccm
        elif self._flowratemode == 1:
            # volumetric flow rate
            # get inlet gas mixture density
            mrate = self.RHO * self._volflowrate
            sccm = mrate / standard_den * 60.0
            return sccm
        elif self._flowratemode == 2:
            # velocity
            if self._haveflowarea:
                mrate = self.RHO * self._flowarea * self._velocity
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

    @property
    def flowarea(self) -> float:
        """
        Get inlet flow area

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
        """
        Set inlet cross-sectional flow area

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
        """
        Get inlet mass flow rate

        Returns
        -------
            mflowrate: double
                mass flow rate [g/sec]
        """
        if self._flowratemode == 0:
            return self._massflowrate
        else:
            return self.convert_to_mass_flowrate()

    @mass_flowrate.setter
    def mass_flowrate(self, mflowrate: float):
        """
        Set inlet mass flow rate

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
        self._SCCM = 0.0
        # set flow rate mode to mass flow rate
        self._flowratemode = 0
        self._inletflowrate[self._flowratemode] = mflowrate
        self._massflowrate = mflowrate

    @property
    def vol_flowrate(self) -> float:
        """
        Get inlet volumetric flow rate

        Returns
        -------
            vflowrate: double
                volumetric flow rate [cm3/sec]
        """
        if self._flowratemode == 1:
            return self._volflowrate
        else:
            return self.convert_to_vol_flowrate()

    @vol_flowrate.setter
    def vol_flowrate(self, vflowrate: float):
        """
        Set inlet volumetric flow rate

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
        self._SCCM = 0.0
        # set flow rate mode to volumetric flow rate
        self._flowratemode = 1
        self._inletflowrate[self._flowratemode] = vflowrate
        self._volflowrate = vflowrate

    @property
    def sccm(self) -> float:
        """
        Get inlet SCCM volumetric flow rate

        Returns
        -------
            vflowrate: double
                SCCM volumetric flow rate [standard cm3/min]
        """
        if self._flowratemode == 3:
            return self._SCCM
        else:
            return self.convert_to_SCCM()

    @sccm.setter
    def sccm(self, vflowrate: float):
        """
        Set inlet volumetric flow rate in SCCM

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
        self._SCCM = vflowrate

    @property
    def velocity(self) -> float:
        """
        Get inlet gas velocity

        Returns
        -------
            vel: double
                velocity [cm/sec]
        """
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
        """
        Set inlet velocity

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
        self._SCCM = 0.0
        # set flow rate mode to velocity
        self._flowratemode = 2
        self._inletflowrate[self._flowratemode] = vel
        self._velocity = vel

    @property
    def velocity_gradient(self) -> float:
        """
        Get inlet gas axial velocity gradient (for premixed, oppdif, and spin)
        or radial velocity spreading rate (v_r/r) at the inlet.

        Returns
        -------
            velgrad: double
                velocity gradient [1/sec]
        """
        return self._velgrad

    @velocity_gradient.setter
    def velocity_gradient(self, velgrad: float):
        """
        Set inlet axial velocity gradient

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

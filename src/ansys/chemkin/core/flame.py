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

"""Chemkin utilities for steady state, 1-D flame models."""

from ctypes import c_int

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core.color import Color as Color
from ansys.chemkin.core.grid import Grid
from ansys.chemkin.core.inlet import Stream
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.reactormodel import Keyword, ReactorModel
from ansys.chemkin.core.steadystatesolver import SteadyStateSolver


class Flame(ReactorModel, SteadyStateSolver, Grid):
    """Generic steady state, one dimensional flame model."""

    def __init__(self, fuelstream: Stream, label: str):
        """Create a 1-D flame object."""
        """
        Create a 1-D flame object.

        Parameters
        ----------
            fuelstream: Stream object
                unburned fuel (or fuel + air) inlet stream
            label: string, optional
                reactor name

        """
        # check Stream object
        if isinstance(fuelstream, Stream):
            self.label = label
            # initialization
            ReactorModel.__init__(self, reactor_condition=fuelstream, label=self.label)
        else:
            # wrong argument type
            msg = [Color.RED, "the first argument must be a Stream object.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        if self.reactormixture.transport_data == 0:
            # transport property is required by the flame models
            msg = [Color.RED, "transport properties are required.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # initialize steady-state solver
        SteadyStateSolver.__init__(self)
        # initialize mesh quality control
        Grid.__init__(self)
        # get mass flow rate
        self.mass_flow_rate = fuelstream.mass_flowrate
        # use API mode for steady-state flame simulations
        Keyword.no_fullkeyword = False
        # FORTRAN file unit of the text output file
        self._mylout = c_int(160)
        # temperature profile is set
        self.temp_profile_set = False
        # use the grid points in the temperature profile as the initial mesh
        self.grid_t_profile = False
        # energy equation
        self.energytypes = {"ENERGY": 1, "GivenT": 2}
        self._energytype = c_int(1)
        # solver type: always steady state
        self._solvertype = c_int(1)
        # species mass transport property type
        # not set = 0; mixture_averaged = 1; multi-component = 2
        self.transport_mode = 0
        # number of grid points in the solution
        self._numbsolutionpoints = 0
        #
        self._speciesmode = "mass"
        # number of required input
        self._numb_requiredinput = 0

    def set_temperature_profile(
        self, x: npt.NDArray[np.double], temp: npt.NDArray[np.double]
    ) -> int:
        """Specify temperature profile."""
        """
        Specify temperature profile.

        Parameters
        ----------
            x: 1D double array
                position value of the profile data [cm]
            temp: 1D double array
                temperature value of the profile data [K]

        Returns
        -------
            error code: integer

        """
        keyword = "TPRO"
        ierr = self.setprofile(key=keyword, x=x, y=temp)
        if ierr == 0:
            self.temp_profile_set = True
        return ierr

    def use_temp_profiel_initial_mesh(self, on: bool = False):
        """Use the temperature profile grid as the initial grid."""
        """Use the grid points in the user defined initial/estimated
        temperature profile as the initial/starting grid points.

        Parameters
        ----------
            on: boolean, default = False
                use the grid points of the temperature profile
                as the initial grid points

        """
        self.grid_t_profile = on

    def set_convection_differencing_type(self, mode: str):
        """Specify the finite differencing scheme."""
        """Set the finite differencing scheme for the convective terms
        in the transport equations.

        Parameters
        ----------
            mode: string, {"central", "upwind"}
                finite difference discretizing scheme

        """
        if mode.lower() == "central":
            # central differencing scheme
            self.setkeyword(key="CDIF", value=True)
            if "WDIF" in self._keyword_index:
                self.setkeyword(key="WDIF", value=False)
        else:
            # upwind differencing scheme, default
            self.setkeyword(key="WDIF", value=True)
            if "CDIF" in self._keyword_index:
                self.setkeyword(key="CDIF", value=False)

    def set_mesh_keywords(self) -> int:
        """Set mesh related keywords."""
        """
        Set mesh related keywords.

        Returns
        -------
            error code: integer

        """
        ierr = 0
        # set up initial mesh related keywords
        if self.grid_t_profile:
            # use temperature profile grid as the initial grid
            if self.temp_profile_set:
                # verify TPRO is set
                self.setkeyword(key="USE_TPRO_GRID", value=True)
            else:
                msg = [Color.PURPLE, "temperature profile is NOT set.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr = 1
                return ierr
        elif self.numb_grid_profile > 0:
            # use user provided mesh
            if self.grid_profile[0] != self.starting_x:
                # check startig point
                msg = [
                    Color.PURPLE,
                    "first user grid point does not match the starting position.\n",
                    Color.SPACEx6,
                    "expected starting position = ",
                    str(self.starting_x),
                    "\n",
                    Color.SPACEx6,
                    "first user grid point =",
                    str(self.grid_profile[0]),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr = 2
                return ierr
            if self.grid_profile[self.numb_grid_profile - 1] != self.ending_x:
                # check startig point
                msg = [
                    Color.PURPLE,
                    "last user grid point does not match the ending position.\n",
                    Color.SPACEx6,
                    "expected ending position = ",
                    str(self.ending_x),
                    "\n",
                    Color.SPACEx6,
                    "last user grid point =",
                    str(self.grid_profile[self.numb_grid_profile - 1]),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr = 3
                return ierr
            count = 0
            for x in self.grid_profile:
                this_key = "GRID    " + str(x)
                self.setkeyword(this_key, True)
                count += 1
            if count != self.numb_grid_profile:
                msg = [Color.PURPLE, "grid profile has problem.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierr = abs(count - self.numb_grid_profile) * 10
                return ierr
        else:
            # uniform mesh
            self.setkeyword(key="NPTS", value=self.numb_grid_points)

        # set maximum number of grid points allowed
        self.setkeyword(key="NTOT", value=self.max_numb_grid_points)
        # reaction zone center
        if self.reaction_zone_center_x > 0.0:
            self.setkeyword(key="XCEN", value=self.reaction_zone_center_x)
        # reaction zone width
        if self.reaction_zone_width > 0.0:
            self.setkeyword(key="WMIX", value=self.reaction_zone_width)
        # mesh adaptation
        # maximum solution profile gradient
        self.setkeyword(key="GRAD", value=self.gradient)
        # maximum solution profile curvature
        self.setkeyword(key="CURV", value=self.curvature)
        # maximum number of adaptive points can be added at one time
        self.setkeyword(key="NADP", value=self.max_numb_adapt_points)
        return ierr

    def set_ss_solver_keywords(self):
        """Add steady-state solver parameter keywoprds to the keyword list."""
        # steady-state solver parameter given
        if len(self.ss_solverkeywords) > 0:
            #
            for k, v in self.ss_solverkeywords.items():
                self.setkeyword(k, v)

    # mass transport options

    def use_mixture_averaged_transport(self):
        """Use the mixture-averaged transport properties."""
        self.setkeyword(key="MIX", value=True)
        if self.transport_mode == 2:
            # turn OFF the multi-component transport properties
            self.removekeyword(key="MULT")
        self.transport_mode = 1

    def use_multicomponent_transport(self):
        """Use the multi-component transport properties."""
        """
        Use of the multi-component transport properties is recommended when
        the pressure is low.
        """
        self.setkeyword(key="MULT", value=True)
        if self.transport_mode == 1:
            # turn OFF the mixture-averaged transport properties
            self.removekeyword(key="MIX")
        self.transport_mode = 2

    def use_fixed_lewis_number_transport(self, lewis: float = 1.0):
        """Compute the species diffusion coefficient with fixed Lewis number."""
        """Use a fixed Lewis number to compute the species diffusion coefficient
        from mixture conductivity.

        .. math::

            Le = Sc/Pr = \\frac{(\\kappa/\\rho C_{p})}{D}

        Parameters
        ----------
            Lewis: double
                Lewis number

        """
        if lewis > 0.0:
            self.setkeyword(key="LEWIS", value=lewis)
            # turn OFF the multi-component transport properties
            if self.transport_mode == 2:
                self.use_mixture_averaged_transport()
            # turn OFF thermal diffusion
            if "TDIF" in self._keyword_index:
                self.removekeyword(key="TDIF")
        else:
            msg = [Color.PURPLE, "Lewis number > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def use_thermal_diffusion(self, mode: bool = True):
        """Include the thermal diffusion (Doret) effect."""
        """
        The inclusion of thermal diffucivity is recommended when there are significant
        amount of "light" species in the system. Species with molecular weight
        less than 5 g/mol is considered a light species, for example, hydrogen.

        Parameters
        ----------
            mode: boolean {True, False}
                ON/OFF

        """
        self.setkeyword(key="TDIF", value=mode)

    def set_species_boundary_types(self, mode: str = "comp"):
        """Set the species boundary condition type (at inlet and outlet)."""
        """
        When the species value is fixed at the boundary, the "back" diffusion
        of the species into the inlet stream is not considered. That is, when
        the species profile is positive at the inlet, a fixed species value implies
        the species concentration in the inlet stream is actually lower than its
        boundary concentration. When the "flux" option is employed, species mass flux
        will be balanced at the inlet. That is, the species concentration at the inlet
        will be slightly varied from its specified boundary value so that the "net"
        species mass flux at the inlet is zero.

        Parameters
        ----------
            mode: string, {"flux", "comp"}
                keep the species mole/mass fraction at a fixed value or
                keep the species mass flux conserved

        """
        if mode.lower() == "comp":
            # keep the species boundary value fixed at the inlet
            self.setkeyword(key="COMP", value=True)
            if "FLUX" in self._keyword_index:
                self.setkeyword(key="FLUX", value=False)
        else:
            # keep the species mass flux conserved at the inlet
            self.setkeyword(key="FLUX", value=True)
            if "COMP" in self._keyword_index:
                self.setkeyword(key="COMP", value=False)

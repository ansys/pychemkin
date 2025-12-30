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

"""Chemkin open reactor utilities."""

from ctypes import c_int
from typing import Union

from ansys.chemkin.core.color import Color as Color
from ansys.chemkin.core.inlet import Stream, clone_stream
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.mixture import Mixture
from ansys.chemkin.core.reactormodel import Keyword, ReactorModel
from ansys.chemkin.core.steadystatesolver import SteadyStateSolver


class OpenReactor(ReactorModel, SteadyStateSolver):
    """Generic open reactor model."""

    def __init__(self, guessedmixture: Stream, label: Union[str, None] = None):
        """Create a steady-state flow reactor object."""
        """
        Parameters
        ----------
            guessedmixture: Mixture object
                guessed/estimate reactor condition
            label: string, optional
                reactor name

        """
        # check reactor Mixture object
        if isinstance(guessedmixture, (Stream, Mixture)):
            if label is None:
                self.label = "Reactor"
            else:
                self.label = label
            # initialization
            ReactorModel.__init__(
                self, reactor_condition=guessedmixture, label=self.label
            )
        else:
            # wrong argument type
            msg = [Color.RED, "the first argument must be a Mixture object.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # initialize steady-state solver
        SteadyStateSolver.__init__(self)
        # use API mode for steady-state open reactor/flame simulations
        Keyword.no_fullkeyword = True
        # FORTRAN file unit of the text output file
        self._mylout = c_int(158)
        # inlet information
        # number of external inlets
        self.numbexternalinlets = 0
        # dict of external inlet objects {inlet label: inlet object}
        self.externalinlets: dict[str, Stream] = {}
        # total mass flow rate into this reactor [g/sec]
        self.totalmassflowrate = 0.0
        #
        self.solver_types = {"Transient": 1, "SteadyState": 2}
        self.energy_types = {"ENERGY": 1, "GivenT": 2}
        # specify "reactor residence time" = "SETTAU" or
        # specify "reactor volume" = "SETVOLV"
        self.problem_types = {"SETVOL": 1, "SETTAU": 2}

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
            flowrate = extinlet.mass_flowrate
            if flowrate <= 0.0:
                msg = [
                    Color.PURPLE,
                    "inlet flow rate < 0.\n",
                    Color.SPACEx6,
                    "specify flow rate of the 'Stream' object.",
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

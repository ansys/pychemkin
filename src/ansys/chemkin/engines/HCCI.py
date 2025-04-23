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
    Single or multi- zone homogeneous charge compression ignition (HCCI) engine model.
"""

import copy
from ctypes import c_double, c_int
from typing import Union

from ansys.chemkin import chemkin_wrapper
from ansys.chemkin.chemistry import (
    check_chemistryset,
    chemistryset_initialized,
    force_activate_chemistryset,
    set_verbose,
    verbose,
)
from ansys.chemkin.color import Color as Color
from ansys.chemkin.constants import Patm
from ansys.chemkin.engines.engine import Engine
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
from ansys.chemkin.reactormodel import Keyword
import numpy as np


class HCCIengine(Engine):
    """
    Single or multi- zone homogeneous charge compression ignition (HCCI) engine model
    """

    def __init__(
        self,
        reactor_condition: Stream,
        label: str = "",
        nzones: Union[int, None] = None,
    ):
        """
        Initialize a single- or multi- zone homogeneous charge compression ignition (HCCI) engine object

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside the HCCI engine/zone
            label: string, optional
                HCCI engine name
            nzones: integer, optional, default = 1
                number of zones in the HCCI engine model
        """
        # set default number of zone(s): single-zone
        if nzones is None:
            nzones = 1
        # set default label
        if label == "":
            if nzones == 1:
                label = "HCCI"
            else:
                label = "Multi-Zone HCCI"

        # use the first zone to initialize the engine model
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("HCCI", 4))
        self._solvertype = c_int(self.SolverTypes.get("Transient", 1))
        self._problemtype = c_int(self.ProblemTypes.get("ICEN", 3))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # defaults for all closed homogeneous reactor models
        self._nreactors = nzones
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        # number of zones
        self._nzones = c_int(nzones)
        # must use full keyword mode for multi-zone simulations
        if self._nzones.value > 1:
            Keyword.noFullKeyword = False
        # zonal setup mode for the multi-zone engine simulation
        # 0 = single-zone or multi-zone with uniform zonal properties
        # 1 = multi-zone with raw species mole fractions
        # 2 = multi-zone with equivalence ratio
        self._zonalsetupmode: int = 0
        # zonal temperature values
        self.zonetemperature: list[float] = []
        # zonal volume fractions
        self.zonevolume: list[float] = []
        # zonal mass fractions
        self.usezonemass = False
        self.zonemass: list[float] = []
        # zonal wall heat transfer area fraction
        self.zoneHTarea: list[float] = []
        # zonal gas compositions in mole fraction (for zonalsetupmode =1)
        self.zonemolefrac: list[float] = []  # list of mole fraction arrays
        # zonal equivalence ratios (for zonalsetupmode =2)
        self.zoneequivalenceratio: list[float] = []
        # fuel composition for all zones
        self.zonefueldefined: list[tuple] = []
        # oxidizer composition for all zones
        self.zoneoxiddefined: list[tuple] = []
        # product composition for all zones
        self.zoneproductdefined: list[str] = []
        # zonal additive gas composition
        self.zoneaddmolefrac: list[float] = []
        # zonal EGR ratios
        self.zoneEGRR: list[float] = []
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(155)
        # profile points
        self._profilesize = int(0)
        # set up basic HCCI engine model parameters
        iErr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if iErr == 0:
            # setup HCCI engine model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to initialize the HCCI engine model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()

    def get_number_of_zones(self) -> int:
        """
        Get the number of zones used by the current HCCI simulation

        Returns
        -------
            nzones: integer
                number of zones in the HCCI engine model
        """
        return self._nzones.value

    def set_zonal_temperature(self, zonetemp: list[float]):
        """
        set zonal temperatures for muti-zone HCCI engine simulation

        Parameters
        ----------
            zonetemp: list of doubles, dimension = [nzones]
                zonal temperatures [K]
        """
        nzones = self._nzones.value
        if len(zonetemp) != nzones:
            msg = [
                Color.PURPLE,
                "zonal temperature must be a list of float of size",
                str(nzones),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if len(self.zonetemperature) > 0:
            msg = [Color.YELLOW, "zonal temperatures will be reset.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.zonetemperature = []
        # set zonal temperatures
        for t in zonetemp:
            if t > 0.0:
                self.zonetemperature.append(t)
            else:
                msg = [Color.PURPLE, "zonal temperature must > 0.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        # set zonal definition mode
        if self._zonalsetupmode == 0:
            self._zonalsetupmode = 1

    def set_zonal_volume_fraction(self, zonevol: list[float]):
        """
        set zonal volume fractions for muti-zone HCCI engine simulation

        Parameters
        ----------
            zonevol: list of doubles, dimension = [nzones]
                zonal volume fractions [-]
        """
        nzones = self._nzones.value
        if len(zonevol) != nzones:
            msg = [
                Color.PURPLE,
                "zonal volume fractions must be a list of float of size",
                str(nzones),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if len(self.zonevolume) > 0:
            msg = [
                Color.YELLOW,
                "zonal volume fractions will be reset.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.zonevolume = []
        # set zonal volume fractions (will be normalized)
        for v in zonevol:
            if v > 0.0:
                self.zonevolume.append(v)
            else:
                msg = [Color.PURPLE, "zonal volume must > 0.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    def set_zonal_mass_fraction(self, zonemass: list[float]):
        """
        set zonal mass fractions for muti-zone HCCI engine simulation

        Parameters
        ----------
            zonemass: list of doubles, dimension = [nzones]
                zonal mass fractions [-]
        """
        nzones = self._nzones.value
        if len(zonemass) != nzones:
            msg = [
                Color.PURPLE,
                "zonal mass fractions must be a list of float of size",
                str(nzones),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if len(self.zonemass) > 0:
            msg = [
                Color.YELLOW,
                "zonal mass fractions will be reset.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.zonemass = []
        # set zonal mass fractions (will be normalized)
        for v in zonemass:
            if v > 0.0:
                self.zonemass.append(v)
            else:
                msg = [Color.PURPLE, "zonal mass must > 0.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        # set flag
        self.usezonemass = True

    def set_zonal_heat_transfer_area_fraction(self, zonearea: list[float]):
        """
        set zonal wall heat transfer area fractions for muti-zone HCCI engine simulation

        Parameters
        ----------
            zonearea: list of doubles, dimension = [nzones]
                zonal heat transfer area fractions [-]
        """
        nzones = self._nzones.value
        if len(zonearea) != nzones:
            msg = [
                Color.PURPLE,
                "zonal heat transfer area fractions must be a list of float of size",
                str(nzones),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if len(self.zoneHTarea) > 0:
            msg = [
                Color.YELLOW,
                "zonal heat transfer area fractions will be reset.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.zoneHTarea = []
        # set zonal wall heat transfer area fractions (will be normalized)
        for a in zonearea:
            if a >= 0.0:
                self.zoneHTarea.append(a)
            else:
                msg = [Color.PURPLE, "zonal area must >= 0.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    def set_zonal_gas_mole_fractions(self, zonemolefrac: list[float]):
        """
        set zonal gas mole fractions for muti-zone HCCI engine simulation

        Parameters
        ----------
            zonemolefrac: list of 1-D double arrays (number of gas species), dimension = [nzones]
                zonal gas mole fractions [-]
        """
        nzones = self._nzones.value
        if len(zonemolefrac) != nzones:
            msg = [
                Color.PURPLE,
                "zonal gas mole fraction must be a list of",
                str(nzones),
                "mole fraction arrays of size = number_species",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if self._zonalsetupmode == 2:
            msg = [
                Color.YELLOW,
                "raw gas composition will replace equivalence ratio",
                "to set up the zonal gas compositions",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

        if len(self.zonemolefrac) > 0:
            msg = [Color.YELLOW, "zonal gas mole fractions will be reset.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.zonemolefrac = []
        # set zonal definition mode
        self._zonalsetupmode = 1
        # set zonal gas mole fractions
        for x in zonemolefrac:
            # x must be a double array of size = number of gas species
            self.zonemolefrac.append(x)

    def define_fuel_composition(self, recipe: list[tuple[str, float]]):
        """
        set the fuel composition for setting up zonal gas composition by zonal equivalence ratio

        Parameters
        ----------
            recipe: list of tuples formatted as (species, mole fraction) pairs
        """
        if len(self.zonefueldefined) > 0:
            msg = [
                Color.YELLOW,
                "previous fuel definition will be reset.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.fuel_composition = []
        self.fuel_composition = copy.deepcopy(recipe)

    def define_oxid_composition(self, recipe: list[tuple[str, float]]):
        """
        set the oxidizer composition for setting up zonal gas composition by zonal equivalence ratio

        Parameters
        ----------
            recipe: list of tuples formatted as (species, mole fraction) pairs
        """
        if len(self.zoneoxiddefined) > 0:
            msg = [
                Color.YELLOW,
                "previous oxidizer definition will be reset.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.oxid_composition = []
        self.oxid_composition = copy.deepcopy(recipe)

    def define_product_composition(self, products: list[str]):
        """
        set the complete combustion product species for setting up zonal gas composition by zonal equivalence ratio

        Parameters
        ----------
            products: list of strings
                product species symbols
        """
        if len(self.zoneproductdefined) > 0:
            msg = [
                Color.YELLOW,
                "previous product definition will be reset.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.product_composition = []
        self.product_composition = copy.deepcopy(products)

    def define_additive_fractions(self, addfrac: list[float]):
        """
        set zonal additive gas mole fractions for setting up zonal gas composition by zonal equivalence ratio

        Parameters
        ----------
            addfrac: 1-D double array, dimension = [number of gas species]
                additive gas mole fractions [-]
        """
        nzones = self._nzones.value
        if len(addfrac) != nzones:
            msg = [
                Color.PURPLE,
                "zonal additive gas mole fraction must be a list of",
                str(nzones),
                "mole fraction arrays of size = number_species",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if len(self.zoneaddmolefrac) > 0:
            msg = [
                Color.YELLOW,
                "zonal additive gas mole fractions will be reset.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.zoneaddmolefrac = []
        # set zonal gas mole fractions
        for x in addfrac:
            # x must be a double array of size = number of gas species
            self.zoneaddmolefrac.append(x)

    def set_zonal_equivalence_ratio(self, zonephi: list[float]):
        """
        set zonal wall heat transfer area fractions for setting up zonal gas composition by zonal equivalence ratio

        Parameters
        ----------
            zonephi: 1-D double array, dimension = [nzones]
                zonal equivalence ratios [-]
        """
        nzones = self._nzones.value
        if len(zonephi) != nzones:
            msg = [
                Color.PURPLE,
                "zonal equivalence ratio must be a list of float of size",
                str(nzones),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if self._zonalsetupmode == 1:
            msg = [
                Color.YELLOW,
                "equivalence ratio will replace raw gas mole fractions",
                "to set up the zonal gas compositions",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

        if len(self.zoneequivalenceratio) > 0:
            msg = [
                Color.YELLOW,
                "previous zonal equivalence ratios will be reset.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.zoneequivalenceratio = []
        # set zonal definition mode
        self._zonalsetupmode = 2
        # set zonal equivalence ratio
        for p in zonephi:
            if p >= 0.0:
                self.zoneequivalenceratio.append(p)
            else:
                msg = [Color.PURPLE, "zonal equivalence ratio must >= 0.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    def set_zonal_EGR_ratio(self, zoneegr: list[float]):
        """
        set zonal exhaust gas recirculation (EGR) ratios for setting up zonal gas composition by zonal equivalence ratio

        Parameters
        ----------
            zoneegr: 1-D double array, dimension = [nzones]
                zonal EGR ratios [-]
        """
        nzones = self._nzones.value
        if len(zoneegr) != nzones:
            msg = [
                Color.PURPLE,
                "zonal EGR ratio must be a list of float of size",
                str(nzones),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if len(self.zoneEGRR) > 0:
            msg = [Color.YELLOW, "previous zonal EGR ratios will be reset.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        self.zoneEGRR = []
        # set zonal EGR ratio
        for r in zoneegr:
            if r >= 0.0:
                self.zoneEGRR.append(r)
            else:
                msg = [Color.PURPLE, "zonal EGR ratio must >= 0.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    def set_energy_equation_switch_ON_CA(self, switchCA: float):
        """
        Set the crank angle at which the energy equation will be turn ON for
        the rest of the simulation.
        Before this switch crank angle the given temperature profile(s) or value(s)
        is used in the multi-zone HCCI simulation.

        Parameters
        ----------
            switchCA: double
                energy equation activation crank angle [degree]
        """
        if self._nzones.value == 1:
            msg = [
                Color.PURPLE,
                "energy switch CA",
                "is valid for the multi-zone HCCI engine model only.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if switchCA > self.IVCCA:
            # set keyword
            self.setkeyword(key="ASWH", value=switchCA)
        else:
            msg = [
                Color.PURPLE,
                "energy switch on CA must > starting CA",
                str(self.IVCCA),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_zonal_volume_keyword(self):
        """
        Set zonal volume keyword for the multi-zone HCCI engine simulation
        """
        if self._nzones.value == 1:
            # single zone is not allowed here
            msg = [
                Color.PURPLE,
                "single-zone engine model should not use this method",
                "to set up zonal conditions.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        for izone in range(self._nzones.value):
            # set the zonal number string
            addon = str(izone + 1)
            # set zonal volume fraction
            keyline = (
                "VOL"
                + Keyword.fourspaces
                + str(self.zonevolume[izone])
                + Keyword.fourspaces
                + addon
            )
            self.setkeyword(key=keyline, value=True)

    def set_zonal_mass_keyword(self):
        """
        Set zonal mass keyword for the multi-zone HCCI engine simulation
        """
        for izone in range(self._nzones.value):
            # set the zonal number string
            addon = str(izone + 1)
            # set zonal volume fraction
            keyline = (
                "MZMAS"
                + Keyword.fourspaces
                + str(self.zonemass[izone])
                + Keyword.fourspaces
                + addon
            )
            self.setkeyword(key=keyline, value=True)

    def set_zonal_condition_keywords(self):
        """
        Set zonal initial condition keywords under the Full-Keywords mode
        and use raw species mole fractions to set up zonal gas compositions
        for multi-zone HCCI engine simulation
        """
        if self._nzones.value == 1:
            # single zone is not allowed here
            msg = [
                Color.PURPLE,
                "single-zone engine model should not use this method",
                "to set up zonal conditions.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # set zonal pressure (same for all zones)
        self.setkeyword(key="PRES", value=self._pressure.value / Patm)

        if self._zonalsetupmode == 0:
            # set zonal setup mode to 'zonal species mole fraction'
            self._zonalsetupmode = 1
        for izone in range(self._nzones.value):
            # set the zonal number string
            addon = str(izone + 1)
            # set zonal temperature
            keyline = (
                "TEMP"
                + Keyword.fourspaces
                + str(self.zonetemperature[izone])
                + Keyword.fourspaces
                + addon
            )
            self.setkeyword(key=keyline, value=True)
            if self.usezonemass:
                # set zonal mass fractions
                keyline = (
                    "MZMAS"
                    + Keyword.fourspaces
                    + str(self.zonemass[izone])
                    + Keyword.fourspaces
                    + addon
                )
                self.setkeyword(key=keyline, value=True)
            else:
                # set zonal volume fraction
                if len(self.zonevolume) != self._nzones.value:
                    # zonal volumes are not set
                    # use equal zonal volume fractions
                    volfrac = 1.0 / self._nzones.value
                    self.zonevolume = []
                    for i in range(self._nzones.value):
                        self.zonevolume.append(volfrac)

                keyline = (
                    "VOL"
                    + Keyword.fourspaces
                    + str(self.zonevolume[izone])
                    + Keyword.fourspaces
                    + addon
                )
                self.setkeyword(key=keyline, value=True)

            if len(self.zoneHTarea) > 0:
                # set zonal heat transfer area fraction
                keyline = (
                    "MQAFR"
                    + Keyword.fourspaces
                    + str(self.zoneHTarea[izone])
                    + Keyword.fourspaces
                    + addon
                )
                self.setkeyword(key=keyline, value=True)

            # initial mole fraction
            nspecieslines, species_lines = self.createspeciesinputlineswithaddon(
                "XEST",
                threshold=1.0e-12,
                molefrac=self.zonemolefrac[izone],
                addon=addon,
            )
            for line in species_lines:
                self.setkeyword(key=line, value=True)

    def set_zonal_equivalence_ratio_keywords(self):
        """
        Set zonal initial condition keywords under the Full-Keywords mode
        and use equivalence ratios to set up zonal gas compositions
        for multi-zone HCCI engine simulation
        """
        if self._nzones.value == 1:
            # single zone is not allowed here
            msg = [
                Color.PURPLE,
                "single-zone engine model should not use this method",
                "to set up zonal conditions.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # set zonal pressure (same for all zones)
        self.setkeyword(key="PRES", value=self._pressure.value / Patm)

        if self._zonalsetupmode == 0:
            # set zonal setup mode to 'zonal equivalence ratio'
            self._zonalsetupmode = 2
        for izone in range(self._nzones.value):
            # set the zonal number string
            addon = str(izone + 1)
            # set zonal temperature
            keyline = (
                "TEMP"
                + Keyword.fourspaces
                + str(self.zonetemperature[izone])
                + Keyword.fourspaces
                + addon
            )
            self.setkeyword(key=keyline, value=True)
            if self.usezonemass:
                # set zonal mass fractions
                keyline = (
                    "MZMAS"
                    + Keyword.fourspaces
                    + str(self.zonemass[izone])
                    + Keyword.fourspaces
                    + addon
                )
                self.setkeyword(key=keyline, value=True)
            else:
                # set zonal volume fraction
                if len(self.zonevolume) != self._nzones.value:
                    # zonal volumes are not set
                    # use equal zonal volume fractions
                    volfrac = 1.0 / self._nzones.value
                    self.zonevolume = []
                    for i in range(self._nzones.value):
                        self.zonevolume.append(volfrac)

                keyline = (
                    "VOL"
                    + Keyword.fourspaces
                    + str(self.zonevolume[izone])
                    + Keyword.fourspaces
                    + addon
                )
                self.setkeyword(key=keyline, value=True)

            if len(self.zoneHTarea) > 0:
                # set zonal heat transfer area fraction
                keyline = (
                    "MQAFR"
                    + Keyword.fourspaces
                    + str(self.zoneHTarea[izone])
                    + Keyword.fourspaces
                    + addon
                )
                self.setkeyword(key=keyline, value=True)

            # set zonal equivalence ratio
            keyline = (
                "EQUI"
                + Keyword.fourspaces
                + str(self.zoneequivalenceratio[izone])
                + Keyword.fourspaces
                + addon
            )
            self.setkeyword(key=keyline, value=True)
            if len(self.zoneEGRR) > 0:
                # set zonal EGR ratio
                keyline = (
                    "EGRR"
                    + Keyword.fourspaces
                    + str(self.zoneEGRR[izone])
                    + Keyword.fourspaces
                    + addon
                )
                self.setkeyword(key=keyline, value=True)
            if izone == 0:
                # only need to set the definitions once
                # set fuel composition by using recipe
                for s, x in self.fuel_composition:
                    keyline = (
                        "FUEL" + Keyword.fourspaces + s + Keyword.fourspaces + str(x)
                    )
                    self.setkeyword(key=keyline, value=True)
                # set oxidizer composition by using recipe
                for s, x in self.oxid_composition:
                    keyline = (
                        "OXID" + Keyword.fourspaces + s + Keyword.fourspaces + str(x)
                    )
                    self.setkeyword(key=keyline, value=True)
                # define complete combustion products by using species symbol list
                for s in self.product_composition:
                    keyline = "CPROD" + Keyword.fourspaces + s
                    self.setkeyword(key=keyline, value=True)

            # zonal additive mole fraction
            nspecieslines, species_lines = self.createspeciesinputlineswithaddon(
                "ADD",
                threshold=1.0e-12,
                molefrac=self.zoneaddmolefrac[izone],
                addon=addon,
            )
            for line in species_lines:
                self.setkeyword(key=line, value=True)

    def __process_keywords_withFullInputs(self) -> int:
        """
        Process input keywords for the HCCI engine model under the Full-Keyword mode

        Returns
        -------
            Error code: integer
        """
        iErr = 0
        iErrProf = 0
        set_verbose(True)
        # verify required inputs
        iErr = self.validate_inputs()
        if iErr != 0:
            msg = [Color.PURPLE, "missing required input keywords", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return iErr
        # re-size work arrays if profile is used
        if self._numbprofiles > 0:
            # find total profile data points
            numbprofilepoints = 0
            for p in self._profiles_list:
                numbprofilepoints += p.size
            if numbprofilepoints > 0:
                # re-size work arrays
                ipoints = c_int(numbprofilepoints)
                iErrc = chemkin_wrapper.chemkin.KINAll0D_SetProfilePoints(ipoints)
                # setup reactor model working arrays
                if iErrc == 0:
                    iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                        self._myLOUT, self._chemset_index
                    )
                iErr += iErrc
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to set up profile keywords,",
                "error code =",
                str(iErrProf),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return iErr
        # prepare initial conditions
        # initial mass fraction
        Y_init = self.reactormixture.Y
        # connecting rod length to crank radius ratio
        LOLR = c_double(self.connectrodlength / self.crankradius)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("HCCI"):
            iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupHCCIInputs(
                self._chemset_index,
                c_double(self.IVCCA),
                c_double(self.EVOCA),
                c_double(self.enginespeed),
                c_double(self.compressratio),
                c_double(self.borediam),
                c_double(self.enginestroke),
                LOLR,
                self._temperature,
                self._pressure,
                self._heat_loss_rate,
                Y_init,
            )
            iErr += iErrc
            if iErrc != 0:
                logger.error("failed to set up basic reactor keywords")
                return iErrc

        # set reactor type
        self.set_reactortype_keywords()
        # set number of zones
        self.setkeyword(key="NZONE", value=self._nzones.value)
        # set engine parameter
        self.set_engine_keywords()
        # check if the wall heat transfer model is set up
        if self._wallheattransfer:
            self.set_heat_transfer_keywords()
        #
        if self._nzones.value == 1 or self._zonalsetupmode == 0:
            # single-zone HCCI engine initial condition or
            # multi-zone with uniform zonal properties
            if self._nzones.value > 1:
                if self.usezonemass:
                    # set zonal mass fractions (required for the multi-zone model)
                    self.set_zonal_mass_keyword()
                else:
                    if len(self.zonevolume) != self._nzones.value:
                        # zonal volumes are not set
                        # set uniform zonal volume fractions
                        volfrac = 1.0 / self._nzones.value
                        self.zonevolume = []
                        for i in range(self._nzones.value):
                            self.zonevolume.append(volfrac)
                    # set zonal volume fractions (required for the multi-zone model)
                    self.set_zonal_volume_keyword()
            # set uniform cylinder properties
            self.set_enginecondition_keywords()
        else:
            # multi-zone HCCI engine zonal conditions
            if self._zonalsetupmode == 0:
                if self.usezonemass:
                    # set zonal mass fractions (required for the multi-zone model)
                    self.set_zonal_mass_keyword()
                else:
                    if len(self.zonevolume) != self._nzones.value:
                        # zonal volumes are not set
                        # use equal zonal volume fractions
                        volfrac = 1.0 / self._nzones.value
                        self.zonevolume = []
                        for i in range(self._nzones.value):
                            self.zonevolume.append(volfrac)
                # uniform zonal properties
                self.set_enginecondition_keywords()
            elif self._zonalsetupmode == 1:
                # non-uniform zonal properties with zonal raw gas mole fractions provided
                self.set_zonal_condition_keywords()
            else:
                # non-uniform zonal properties with zonal equivalence ratios provided
                self.set_zonal_equivalence_ratio_keywords()
        #
        if iErr == 0 and self._numbprofiles > 0:
            # get keyword lines of all profiles
            iErrProf, nproflines, prof_lines = self.createprofileinputlines()
            iErr += iErrProf
            if iErrProf == 0:
                # set the profile keywords
                for pkey in prof_lines:
                    for line in pkey:
                        self.setkeyword(key=line, value=True)
            else:
                msg = [
                    Color.PURPLE,
                    "failed to set up profile keywords,",
                    "error code =",
                    str(iErrProf),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return iErrProf
        # solve integrated heat release rate due to chemical reactions
        self.setkeyword(key="QRGEQ", value=True)
        # add the END keyword
        self.setkeyword(key="END", value=True)
        # create input lines from additional user-specified keywords
        iErr, nlines = self.createkeywordinputlines()
        if iErr == 0:
            if verbose():
                msg = [
                    Color.YELLOW,
                    str(nlines),
                    "input lines are added.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
        else:
            msg = [
                Color.PURPLE,
                "failed to create additional keyword lines.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

        return iErr

    def __run_model_withFullInputs(self) -> int:
        """
        Run the HCCI engine model after the keywords are processed under the Full-Keyword mode
        All keywords must be assigned

        Returns
        -------
            Error code: integer
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
        iErr = chemkin_wrapper.chemkin.KINAll0D_CalculateInput(
            self._myLOUT, self._chemset_index, longline, nlines, linelength
        )
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to set up reactor keywords in Full mode,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

        return iErr

    def __process_keywords(self) -> int:
        """
        Process input keywords for the HCCI engine model

        Returns
        -------
            Error code: integer
        """
        iErr = 0
        iErrc = 0
        iErrKey = 0
        iErrInputs = 0
        # set_verbose(True)
        # verify required inputs
        iErr = self.validate_inputs()
        if iErr != 0:
            msg = [Color.PURPLE, "missing required input keywords", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return iErr
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
                iErrc = chemkin_wrapper.chemkin.KINAll0D_SetProfilePoints(ipoints)
                # setup reactor model working arrays
                if iErrc == 0:
                    iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                        self._myLOUT, self._chemset_index
                    )
                iErr += iErrc
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "profile data generation error, error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return iErr
        # prepare initial conditions
        # initial mass fraction
        Y_init = self.reactormixture.Y
        # connecting rod length to crank radius ratio
        LOLR = c_double(self.connectrodlength / self.crankradius)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("HCCI"):
            iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupHCCIInputs(
                self._chemset_index,
                c_double(self.IVCCA),
                c_double(self.EVOCA),
                c_double(self.enginespeed),
                c_double(self.compressratio),
                c_double(self.borediam),
                c_double(self.enginestroke),
                LOLR,
                self._temperature,
                self._pressure,
                self._heat_loss_rate,
                Y_init,
            )
            iErr += iErrc
            if iErrc != 0:
                msg = [
                    Color.PURPLE,
                    "failed to set up basic reactor keywords,",
                    "error code =",
                    str(iErrc),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return iErrc
            # heat transfer (use additional keywords)
            # solver parameters (use additional keywords)
            # output controls (use additional keywords)
            # ROP (use additional keywords)
            # sensitivity (use additional keywords)
            # ignition delay (use additional keywords)
            # solve integrated heat release rate due to chemical reactions
            iErrc = chemkin_wrapper.chemkin.KINAll0D_IntegrateHeatRelease()
            iErr += iErrc
            if iErrc != 0:
                msg = [
                    Color.PURPLE,
                    "failed to set up heat release keyword,",
                    "error code =",
                    str(iErrc),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return iErrc
        else:
            pass
        # check if the wall heat transfer model is set up
        if iErr == 0 and self._wallheattransfer:
            self.set_heat_transfer_keywords()
        #
        if iErr == 0 and self._numbprofiles > 0:
            for p in self._profiles_list:
                key = bytes(p.profilekey, "utf-8")
                npoints = c_int(p.size)
                x = p.pos
                y = p.value
                iErrProf = chemkin_wrapper.chemkin.KINAll0D_SetProfileParameter(
                    key, npoints, x, y
                )
                iErr += iErrProf
            if iErrProf != 0:
                msg = [
                    Color.PURPLE,
                    "failed to set up profile keywords,",
                    "error code =",
                    str(iErrProf),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return iErrProf
        if iErr == 0:
            # set additional keywords
            # create input lines from additional user-specified keywords
            iErrInputs, nlines = self.createkeywordinputlines()
            if iErrInputs == 0:
                # process additional keywords in _keyword_index and _keyword_lines
                for s in self._keyword_lines:
                    # convert string to byte
                    line = bytes(s, "utf-8")
                    # set additional keyword one by one
                    iErrKey = chemkin_wrapper.chemkin.KINAll0D_SetUserKeyword(line)
                if iErrInputs == 0:
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
                        str(iErrInputs),
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
            else:
                msg = [
                    Color.PURPLE,
                    "failed to process additional keywords, error code =",
                    str(iErrInputs),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
        #
        iErr = iErr + iErrInputs + iErrKey

        return iErr

    def __run_model(self) -> int:
        """
        Run the HCCI engine model after the keywords are processed

        Returns
        -------
            Error code: integer
        """
        # run the simulation without keyword inputs
        iErr = chemkin_wrapper.chemkin.KINAll0D_Calculate(self._chemset_index)
        return iErr

    def run(self) -> int:
        """
        Generic Chemkin run HCCI engine model method

        Returns
        -------
            Error code: integer
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
            retVal = chemkin_wrapper.chemkin.KINInitialize(
                self._chemset_index, c_int(0)
            )
            if retVal != 0:
                msg = [
                    Color.RED,
                    "Chemkin-CFD-API initialization failed;",
                    "code =",
                    str(retVal),
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
        msg = [Color.YELLOW, "processing and generating keyword inputs ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        #
        if self._nzones.value == 1 and Keyword.noFullKeyword:
            # use API calls
            retVal = (
                self.__process_keywords()
            )  # each reactor model subclass to perform its own keyword processing
        else:
            # use full keywords
            retVal = self.__process_keywords_withFullInputs()
        if retVal != 0:
            msg = [
                Color.RED,
                "generating the keyword inputs,",
                "error code =",
                str(retVal),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            return retVal
        logger.debug("Processing keywords complete")

        # run reactor model
        msg = [Color.YELLOW, "running HCCI engine simulation ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        if self._nzones.value == 1 and Keyword.noFullKeyword:
            # single-zone HCCI
            # use API calls
            retVal = self.__run_model()
        else:
            # multi-zone HCCI
            # use full keywords
            retVal = self.__run_model_withFullInputs()
        # update run status
        self.setrunstatus(code=retVal)
        msg = ["simulation completed,", "status =", str(retVal), Color.END]
        if retVal == 0:
            msg.insert(0, Color.GREEN)
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        else:
            msg.insert(0, Color.RED)
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)

        return retVal

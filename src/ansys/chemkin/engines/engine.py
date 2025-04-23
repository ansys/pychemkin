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
    Chemkin Engine Utilities.
"""

import copy
from ctypes import c_double, c_int
from typing import Union

from ansys.chemkin import chemkin_wrapper
from ansys.chemkin.batchreactors.batchreactor import BatchReactors
from ansys.chemkin.color import Color as Color
from ansys.chemkin.constants import Patm
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
from ansys.chemkin.reactormodel import Keyword
import numpy as np


class Engine(BatchReactors):
    """
    Generic engine cylinder model
    """

    def __init__(self, reactor_condition: Stream, label: str):
        """
        Initialize a generic Engine object

        Parameters
        ----------
            reactor_condition: Mixture object
                a mixture representing the initial gas properties inside the engine cylinder/zone
            label: string, optional
                engine reactor name
        """
        super().__init__(reactor_condition, label)
        # engine parameters
        # stroke type
        self._numstroke = 4
        # opposed-piston engine
        self._opposedpistonmode = 0
        # bore diameter [cm]
        self.borediam = 0.0
        # bore cross-sectional area [cm2]
        self.borearea = 0.0
        # stroke [cm]
        self.enginestroke = 0.0
        # crank radius [cm]
        self.crankradius = 0.0
        # connecting rod length [cm]
        self.connectrodlength = 0.0
        # piston pin offset [cm]
        self.pistonoffset = 0.0e0
        # cylinder head surface area [cm2] (won't change)
        self.cylinderheadarea = 0.0
        # piston head surface area [cm2] (won't change)
        self.pistonheadarea = 0.0
        # head areas = cylinder head area + piston head area
        self.headareas = 0.0
        # compression ratio
        self.compressratio = 1.0e0
        # engine speed RPM
        self.enginespeed = 1.0e0
        self.degpersec = 0.0
        self.radpersec = 0.0
        # IVC CA (start of engine simulation when the cylinder becomes a closed system because the intake valve is closed)
        self.IVCCA = -180.0
        # EVO CA (end of engine simulation when the cylinder becomes an open system because the exhaust valve is opened)
        self.EVOCA = 180.0
        # default duration = 1 engine revolution
        self.rundurationCA = 360.0
        # engine wall heat transfer models
        # ICHX: dimensionless correlation "ICHX <a> <b> <c> <Twall>"
        # ICHW: dimensional correlation "ICHW <a> <b> <c> <Twall>"
        # ICHH: Hohenburg correlation "ICHH <a> <b> <c> <d> <e> <Twall>"
        self._WallHeatTransferModels = ["ICHX", "ICHW", "ICHH"]
        # heat transfer model parameters
        self.numbHTmodelparameters = [3, 3, 5]
        self.heattransfermodel: int = -1
        self.heattransferparameters: list[float] = []
        self.cylinderwalltemperature = 298.15  # [K]
        # incylinder gas speed correlations
        # velocity correlation parameters
        # Woschni "GVEL <C11> <C12> <C2> <swirling ratio>"
        self.gasvelocity: list[float] = []
        # Woschni+Huber IMEP "HIMP <IMEP>"
        self.HuberIMEP: Union[float, None] = None  # [atm]
        # flag for wall heat transfer, default = adiabatic
        self._wallheattransfer = False
        # check required inputs
        # number of required parameters:
        # starting CA, ending CA, engine speed, compression ratio, bore, stroke, connecting rod length
        self._numb_requiredinput = 7
        self._requiredlist = [
            "DEG0",
            "DEGE",
            "RPM",
            "CMPR",
            "BORE",
            "STRK",
            "CRLEN",
        ]
        # add engine specific keywords
        Keyword._protectedkeywords.extend(self._requiredlist)

    @staticmethod
    def convert_CA_to_Time(CA: float, startCA: float, RPM: float) -> float:
        """
        Convert the current crank angle value to simulation time

        Parameters
        ----------
            CA: double
                engine crank angle [degree]
            startCA: double
                starting crank angle, IVC timing [degree]
            RPM: double
                engine speed RPM [revolutions per minute]

        Returns
        -------
            time: double
                simulation time [sec]
        """
        if RPM <= 0.0:
            msg = [Color.PURPLE, "engine speed RPM must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 0.0
        #
        time = (CA - startCA) / RPM / 6.0e0
        if time < 0.0:
            msg = [
                Color.PURPLE,
                "given CA is less then the starting CA @ IVC.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 0.0
        else:
            return time

    @staticmethod
    def convert_Time_to_CA(time: float, startCA: float, RPM: float) -> float:
        """
        Convert the current time to crank angle

        Parameters
        ----------
            time: double
                current simulation time [sec]
            startCA: double
                starting crank angle, IVC timing [degree]
            RPM: double
                engine speed RPM [revolutions per minute]

        Returns
        -------
            CA: double
                engine crank angle [degree]
        """
        if time < 0.0:
            msg = [Color.PURPLE, "simulation time must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 0.0
        #
        CA = startCA + time * RPM * 6.0e0
        return CA

    def get_Time(self, CA: float) -> float:
        """
        Convert the current crank angle value to simulation time

        Parameters
        ----------
            CA: double
                current engine crank angle [degree]

        Returns
        -------
            time: double
                simulation time [sec]
        """
        return (CA - self.IVCCA) / self.degpersec

    def get_CA(self, time: float) -> float:
        """
        Convert the current time to crank angle

        Parameters
        ----------
            time: double
                current simulation time [sec]

        Returns
        -------
            CA: double
                engine crank angle [degree]
        """
        return self.IVCCA + time * self.degpersec

    @property
    def starting_CA(self) -> float:
        """
        Get the simulation starting crank angle [degree]
        usually the starting CA ~ the intake valve close (IVC) timing

        Returns
        -------
            IVCCA: double
                intake valve closing (IVC) crank angle [degree]
        """
        return self.IVCCA

    @starting_CA.setter
    def starting_CA(self, startCA: float):
        """
        Set the starting crank angle of engine simulation,
        usually this corresponds to the intake valve close (IVC) timing
        a positive starting CA implies the standard top dead center (TDC) is at 360 degrees CA
        a negative starting CA implies the standard TDC is at 0 degree CA

        Parameters
        ----------
            startCA: double
                starting crank angle [degree]
        """
        # set IVC timing in CA
        self.IVCCA = startCA
        # set keyword
        self._inputcheck.append("DEG0")

    @property
    def ending_CA(self) -> float:
        """
        Get the simulation ending crank angle [degree]
        usually the ending CA ~ the exhaust valve open (EVO) timing

        Returns
        -------
            EVOCA: double
                exhaust valve opening (EVO) crank angle [degree]
        """
        return self.EVOCA

    @ending_CA.setter
    def ending_CA(self, endCA: float):
        """
        Set the ending crank angle of engine simulation,
        usually this corresponds to the exhaust valve open (EVO) timing

        Parameters
        ----------
            endCA: double
                ending crank angle [degree]
        """
        # check EVO timing value
        if endCA <= self.starting_CA:
            msg = [
                Color.PURPLE,
                "ending CA must > starting CA =",
                str(self.starting_CA),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # set EVO timing in CA
        self.EVOCA = endCA
        self.rundurationCA = self.ending_CA - self.starting_CA
        # set keyword
        self._inputcheck.append("DEGE")

    @property
    def duration_CA(self) -> float:
        """
        Get the simulation duration in number of crank angles [degree]

        Returns
        -------
            CA: double
                simulation duration in crank angles [degree]
        """
        return self.rundurationCA

    @duration_CA.setter
    def duration_CA(self, CA: float):
        """
        Set the engine simulation duration in CA

        Parameters
        ----------
            CA: double
                crank angle [degree]
        """
        # check EVO timing value
        if CA <= 0.0:
            msg = [Color.PURPLE, "duration CA must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # set EVO timing in CA
        self.rundurationCA = CA
        self.EVOCA = self.IVCCA + CA
        # set keyword
        self._inputcheck.append("DEGE")

    @property
    def bore(self) -> float:
        """
        Get the engine cylinder bore diameter

        Returns
        -------
            diameter: double
                bore diameter [cm]
        """
        return self.borediam

    @bore.setter
    def bore(self, diameter: float):
        """
        Set the engine cylinder bore diameter

        Parameters
        ----------
            diameter: double
                bore diameter [cm]
        """
        if diameter > 0.0:
            self.borediam = diameter
            self.borearea = np.pi * diameter * diameter / 4.0e0
            # set keyword
            self._inputcheck.append("BORE")
        else:
            msg = [Color.PURPLE, "engine bore diameter must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    @property
    def stroke(self) -> float:
        """
        Get the engine stroke

        Returns
        -------
            stroke: double
                engine stroke [cm]
        """
        return self.enginestroke

    @stroke.setter
    def stroke(self, s: float):
        """
        Set the engine stroke

        Parameters
        ----------
            s: double
                engine stroke [cm]
        """
        if s > 0.0:
            self.enginestroke = s
            self.crankradius = s / 2.0e0
            # set keyword
            self._inputcheck.append("STRK")
        else:
            msg = [Color.PURPLE, "piston stroke must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    @property
    def connecting_rod_length(self) -> float:
        """
        Get the connecting rod length

        Returns
        -------
            length: double
                connecting rod length [cm]
        """
        return self.connectrodlength

    @connecting_rod_length.setter
    def connecting_rod_length(self, s: float):
        """
        Set the engine connecting rod length

        Parameters
        ----------
            s: double
                connecting rod length [cm]
        """
        if s > 0.0:
            self.connectrodlength = s
            # set keyword
            self._inputcheck.append("CRLEN")
        else:
            msg = [Color.PURPLE, "piston connecting rod length must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    @property
    def compression_ratio(self) -> float:
        """
        Get the engine compression ratio

        Returns
        -------
            cratio: double
                compression ratio [-]
        """
        return self.compressratio

    @compression_ratio.setter
    def compression_ratio(self, cratio: float):
        """
        Set the engine compression ratio

        Parameters
        ----------
            cratio: double
                compression ratio [-]
        """
        if cratio > 1.0e0:
            self.compressratio = cratio
            # set keyword
            self._inputcheck.append("CMPR")
        else:
            msg = [Color.PURPLE, "engine compression ratio must > 1.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    @property
    def RPM(self) -> float:
        """
        Get the engine speed in RPM

        Returns
        -------
            speed: double
                engine speed [RPM]
        """
        return self.enginespeed

    @RPM.setter
    def RPM(self, speed: float):
        """
        Set the engine speed in RPM

        Parameters
        ----------
            speed: double
                engine speed [RPM]
        """
        if speed > 0.0:
            self.enginespeed = speed
            self.degpersec = speed * 6.0e0
            self.radpersec = self.degpersec * np.pi / 180.0e0
            # set keyword
            self._inputcheck.append("RPM")
        else:
            msg = [Color.PURPLE, "engine speed RPM must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_cylinder_head_area(self, area: float):
        """
        Set the cylinder head clearance surface area

        Parameters
        ----------
            area: double
                area [cm2]
        """
        if area > 0.0:
            self.cylinderheadarea = area
            self.headareas = area + self.pistonheadarea
            # set keyword
            if "BORE" in self._inputcheck:
                self.setkeyword(key="CYBAR", value=area / self.borearea)
            else:
                msg = [
                    Color.PURPLE,
                    "please set cylinder BORE diameter first.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
        else:
            msg = [Color.PURPLE, "cylinder head surface area must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_piston_head_area(self, area: float):
        """
        Set the piston head top surface area

        Parameters
        ----------
            area: double
                piston head top surface area [cm2]
        """
        if area > 0.0:
            self.pistonheadarea = area
            self.headareas = area + self.cylinderheadarea
            # set keyword
            if "BORE" in self._inputcheck:
                self.setkeyword(key="PSBAR", value=area / self.borearea)
            else:
                msg = [
                    Color.PURPLE,
                    "please set cylinder BORE diameter first.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
        else:
            msg = [Color.PURPLE, "piston head surface area must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_piston_pin_offset(self, offset: float):
        """
        Set the piston pin off-set distance

        Parameters
        ----------
            offset: double
                piston pin offset distance [cm]
        """
        if offset < self.crankradius:
            self.pistonoffset = offset
            # set keyword
            self.setkeyword(key="POLEN", value=offset)
        else:
            msg = [
                Color.PURPLE,
                "piston pin offset distance must < crank radius",
                str(self.crankradius),
                "[cm]",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def get_clearance_volume(self) -> float:
        """
        Get the clearance volume

        Returns
        -------
            cvolume: double
                cylinder clearance/minimum volume [cm3]
        """
        if "CMPR" in self._inputcheck:
            dvolume = self.get_displacement_volume()
            cvolume = dvolume / (self.compressratio - 1.0e0)
        else:
            msg = [
                Color.PURPLE,
                "please set engine compression ratio first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            cvolume = 0.0
        return cvolume

    def get_displacement_volume(self) -> float:
        """
        Get the displacement volume

        Returns
        -------
            dvolume: double
                piston displacement/sweeping volume [cm3]
        """
        return self.enginestroke * self.borearea

    def list_engine_parameters(self):
        """
        List engine parameters for verification
        """
        print("      === engine parameters ===")
        print(f"bore diameter         = {self.borediam} [cm]")
        print(f"stroke                = {self.enginestroke} [cm]")
        print(f"connecting rod length = {self.connectrodlength} [cm]")
        print(f"cylinder head area    = {self.cylinderheadarea} [cm2]")
        print(f"piston head area      = {self.pistonheadarea} [cm2]")
        print(f"piston offset         = {self.pistonoffset} [cm]")
        print(f"compression ratio     = {self.compressratio} [-]")
        print(f"engine speed          = {self.enginespeed} [RPM]")
        print(f"IVC crank angle       = {self.IVCCA} [degree]")
        print(f"EVO crank angle       = {self.EVOCA} [degree]")

    @property
    def CAstep_for_saving_solution(self) -> float:
        """
        Get the number of crank angles between saving the solution data

        Returns
        -------
            delta_CA: double
                solution saving interval in crank angles [degree]
        """
        if "DEGSAVE" in self._keyword_index:
            # defined: find index
            i = self._keyword_index.index("DEGSAVE")
            return self._keyword_list[i].value
        else:
            # return default value (100th of the simulation duration)
            if self.rundurationCA > 0.0e0:
                return self.rundurationCA / 1.0e2
            else:
                # not defined yet
                msg = [
                    Color.PURPLE,
                    "solution saving CA interval is not defined",
                    "because the 'ending CA' is not set.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return 0.0

    @CAstep_for_saving_solution.setter
    def CAstep_for_saving_solution(self, delta_CA: float):
        """
        Set the number of crank angles between saving the solution data

        Parameters
        ----------
            delta_CA: double
                number of crank angles between saving solution data [degree]
        """
        if delta_CA > 0.0e0:
            self.setkeyword(key="DEGSAVE", value=delta_CA)
        else:
            msg = [Color.PURPLE, "solution saving CA interval must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    @property
    def CAstep_for_printing_solution(self) -> float:
        """
        Get the number of crank angles between printing the solution data to the text output file

        Returns
        -------
            delta_CA: double
                solution printing interval in crank angles [degree]
        """
        if "DEGPRINT" in self._keyword_index:
            # defined: find index
            i = self._keyword_index.index("DEGPRINT")
            return self._keyword_list[i].value
        else:
            # return default value (100th of the simulation duration in CA)
            if self.rundurationCA > 0.0e0:
                return self.rundurationCA / 1.0e2
            else:
                # not defined yet
                msg = [
                    Color.PURPLE,
                    "solution printing CA interval is not defined",
                    "because the 'ending CA' is not set.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                return 0.0

    @CAstep_for_printing_solution.setter
    def CAstep_for_printing_solution(self, delta_CA: float):
        """
        Set the timestep size between printing the solution data to the text output file

        Parameters
        ----------
            delta_CA: double
                number of crank angles between printing solution data [degree]
        """
        if delta_CA > 0.0e0:
            self.setkeyword(key="DEGPRINT", value=delta_CA)
        else:
            msg = [Color.PURPLE, "solution printing CA interval must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_minimum_zone_mass(self, minmass: float):
        """
        Set the minimum mass in a zone (for Spark Ignition and Direct Injection engine models)

        Parameters
        ----------
            minmass: double, default = 1.0e-6
                minimum zonal mass [g]
        """
        if minmass > 0.0:
            # set keyword
            self.setkeyword(key="MLMT", value=minmass)
        else:
            msg = [Color.PURPLE, "minimum zonal mass must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_zonal_gas_rate_multiplier(
        self, value: float = 1.0e0, zoneID: Union[int, None] = None
    ):
        """
        Set the value of the gas-phase reaction rate multiplier

        Parameters
        ----------
            value: double, default = 1.0
                gas-phase reaction rate multiplier
            zoneID: integer, optional
                zone index to which the multiplier will be applied
                if not provided, the multiplier will be applied to all zones
        """
        if value < 0.0:
            msg = [Color.PURPLE, "reaction rate multiplier must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        else:
            if zoneID is None:
                self._gasratemultiplier = value
                self.setkeyword(key="GFAC", value=value)
            else:
                # zonal GFAC
                self._gasratemultiplier = value
                keyphrase = (
                    "GFAC"
                    + Keyword.fourspaces
                    + str(value)
                    + Keyword.fourspaces
                    + str(zoneID)
                )
                self.setkeyword(key=keyphrase, value=True)

    def set_wall_heat_transfer(
        self, model: str, HTparameters: list[float], walltemperature: float
    ):
        """
        Set cylinder wall heat transfer model and parameters
        engine wall heat transfer models
        ICHX: dimensionless correlation "ICHX <a> <b> <c> <Twall>"
        ICHW: dimensional correlation "ICHW <a> <b> <c> <Twall>"
        ICHH: Hohenburg correlation "ICHH <a> <b> <c> <d> <e> <Twall>"

        Parameters
        ----------
            model: string
                engine wall heat transfer model
            HTparmeters: list of double
                model parameters correspond to the specified heat transfer model
            walltemperature: double
                cylinder wall/cooling oil temperature [K]
        """
        # check existing heat transfer set up
        if self.heattransfermodel >= 0:
            msg = [
                Color.YELLOW,
                "previously defined wall heat transfer model will be overridden.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        #
        mymodel = model.lower()
        # check model
        if mymodel.rstrip() == "dimensionless":
            self.heattransfermodel = 0
        elif mymodel.rstrip() == "dimensioless":
            self.heattransfermodel = 1
        elif mymodel.rstrip() == "hohenburg":
            self.heattransfermodel = 2
        else:
            msg = [
                Color.PURPLE,
                "engine wall heat transfer model",
                model.rstrip(),
                "is not valid.\n",
                Color.SPACEx6,
                "the valid model options are 'dimensional',",
                "'dimensionless', and 'hohenburg'",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check number of parameters
        if len(HTparameters) != self.numbHTmodelparameters[self.heattransfermodel]:
            msg = [
                Color.PURPLE,
                "incorrect number of parameters in the list.\n",
                model,
                "requires",
                str(self.numbHTmodelparameters[self.heattransfermodel]),
                "parameters\n",
                Color.SPACEx6,
                "check Chemkin Input manual for more information.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self.heattransferparameters = []
        # set model parameters
        self.heattransferparameters = copy.copy(HTparameters)
        # set cylinder wall temperature
        self.cylinderwalltemperature = walltemperature
        # set flag for wall heat transfer
        self._wallheattransfer = True

    def set_gas_velocity_correlation(
        self, gasvelparameters: list[float], IMEP: Union[float, None] = None
    ):
        """
        Set the cylinder gas velocity correlation parameters
        Woschni: "GVEL <C11> <C12> <C2> <swirling ratio>"
        Huber: IMEP "HIMP <IMEP>"

        Parameters
        ----------
            gasvelparameters: list of double
                cylinder gas velocity correlation parameters
            IMEP: double, optional, default = 0.0
                indicated mean effective pressure used by the Huber gas velocity correlation [atm]
        """
        # check existing heat transfer set up
        if self.heattransfermodel < 0:
            msg = [
                Color.PURPLE,
                "please specify the wall heat transfer model first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check existing gas velocity correlation parameter
        if len(self.gasvelocity) > 0:
            msg = [
                Color.YELLOW,
                "previously defined gas velocity correlation will be overridden.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        # check number of parameters
        if len(gasvelparameters) != 4:
            msg = [
                Color.PURPLE,
                "incorrect number of parameters.\n",
                Color.SPACEx6,
                "gas velocity correlation requires 4 parameters\n",
                Color.SPACEx6,
                "<C11> <C12> <C2> <swirl ratio>\n",
                Color.SPACEx6,
                "please check Chemkin Input manual for more information.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self.gasvelocity = []
        # set model parameters
        self.gasvelocity = copy.copy(gasvelparameters)
        # IEMP for the Huber correlation
        if IMEP is not None:
            self.HuberIMEP = IMEP  # [atm]

    def set_heat_transfer_keywords(self):
        """
        Set the engine wall heat transfer related keywords
        """
        # check if the wall heat transfer model is set up
        if not self._wallheattransfer:
            return
        # set wall heat transfer model
        line = self._WallHeatTransferModels[self.heattransfermodel]
        for s in self.heattransferparameters:
            line = line + Keyword.fourspaces + str(s).rstrip()
        # add wall temperature to the end
        line = line + Keyword.fourspaces + str(self.cylinderwalltemperature)
        self.setkeyword(key=line, value=True)
        # set incylinder gas speed correlations
        if self.gasvelocity is None:
            # gas velocity correlation is not set
            return
        line = "GVEL"
        for s in self.gasvelocity:
            line = line + Keyword.fourspaces + str(s).rstrip()
        self.setkeyword(key=line, value=True)
        # Huber IMEP
        if self.HuberIMEP is None:
            # IEMP is not set
            return
        self.setkeyword(key="HIMP", value=self.HuberIMEP)

    def set_engine_keywords(self):
        """
        Set engine parameter keywords under the Full-Keywords mode
        """
        self.setkeyword(key="BORE", value=self.borediam)
        self.setkeyword(key="STRK", value=self.enginestroke)
        self.setkeyword(key="CRLEN", value=self.connectrodlength)
        self.setkeyword(key="CMPR", value=self.compressratio)
        self.setkeyword(key="RPM", value=self.enginespeed)
        if np.isclose(abs(self.pistonoffset), 0.0, atol=1.0e-6):
            self.setkeyword(key="POLEN", value=self.pistonoffset)
        self.setkeyword(key="DEG0", value=self.IVCCA)
        self.setkeyword(key="DEGE", value=self.EVOCA)

    def set_enginecondition_keywords(self):
        """
        Set engine initial condition keywords under the Full-Keywords mode
        """
        self.setkeyword(key="PRES", value=self._pressure.value / Patm)
        self.setkeyword(key="TEMP", value=self._temperature.value)
        # initial mole fraction
        nspecieslines, species_lines = self.createspeciesinputlines(
            self._solvertype.value, threshold=1.0e-12, molefrac=self.reactormixture.X
        )
        for line in species_lines:
            self.setkeyword(key=line, value=True)

    def get_engine_heat_release_CAs(self) -> tuple[float, float, float]:
        """
        Get heat release crank angles from the engine solution

        Returns
        -------
            CA10: double
                Crank rotation angle corresponding to 10% of total heat release
            CA50: double
                Crank rotation angle corresponding to 50% of total heat release
            CA90: double
                Crank rotation angle corresponding to 90% of total heat release
        """
        # heat loss rate per CA [erg/degree] sized = 1 + number of surface materials
        QLossRateCA = np.zeros(1, dtype=np.double)
        # apparent heat release rate per CA [erg/degree]
        AHRR = c_double(0.0)
        # apparent heat release rate per CA from PV-ConGamma [erg/degree]
        AHRRP = c_double(0.0)
        # Crank rotation angle corresponding to 10% of total heat release
        HR10 = c_double(self.IVCCA)
        # Crank rotation angle corresponding to 50% of total heat release
        HR50 = c_double(self.IVCCA)
        # Crank rotation angle corresponding to 90% of total heat release
        HR90 = c_double(self.IVCCA)
        # get heat release rate information from the solution
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetEngineHeatRelease(
            QLossRateCA, AHRR, AHRRP, HR10, HR50, HR90
        )
        if iErr != 0:
            # reset the angles to 0 when error is encountered
            HR10 = c_double(0.0)
            HR50 = c_double(0.0)
            HR90 = c_double(0.0)

        return HR10.value, HR50.value, HR90.value

    def get_engine_solution_size(self, expected: int) -> tuple[int, int]:
        """
        Get the number of zones and the number of solution points

        Parameters
        ----------
            expected: integer
                expected number of zonal + mean solution records

        Returns
        -------
            nzones: integer
                number of zones
            npoints: integer
                number of solution points
        """
        # check run completion
        status = self.getrunstatus(mode="silent")
        if status == -100:
            msg = [Color.MAGENTA, "please run the engine simultion first.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        elif status != 0:
            msg = [
                Color.PURPLE,
                "simulation was failed.\n",
                Color.SPACEx6,
                "please correct the error(s) and rerun the engine simulation.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # number of zone
        nzone = c_int(0)
        # number of time points in the solution
        npoints = c_int(0)
        # get solution size of the batch reactor
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetSolnResponseSize(nzone, npoints)
        nzones = nzone.value
        if iErr == 0 and nzones == expected:
            # return the solution sizes
            self._numbsolutionpoints = (
                npoints.value
            )  # number of time points in the solution profile
            return nzones, self._numbsolutionpoints
        elif expected == nzones:
            # fail to get solution sizes
            msg = [
                Color.PURPLE,
                "failed to get the solution size,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            # incorrect number of zones
            msg = [
                Color.PURPLE,
                "incorrect number of zone.\n",
                Color.SPACEx6,
                "the engine model expects",
                str(expected),
                "zone(s)\n",
                Color.SPACEx6,
                str(nzones),
                "found in the solution.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def process_engine_solution(self, zoneID: Union[int, None] = None):
        """
        Post-process solution to extract the raw solution variable data from
        engine simulation results

        Parameters
        ----------
            zoneID: integer
                zone index
        """
        # check existing raw data
        if self.getrawsolutionstatus():
            msg = [
                Color.YELLOW,
                "the solution has been processed before,",
                "any existing solution data will be deleted from the memory.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

        if zoneID is None:
            zoneID = 1

        # reset raw and mixture solution parameters
        self._numbsolutionpoints = 0
        self._solution_rawarray.clear()
        self._solution_mixturearray.clear()
        # check values
        if self._nreactors > 1:
            expectedzones = self._nreactors + 1
        else:
            expectedzones = self._nreactors

        if zoneID > expectedzones:
            msg = [
                Color.PURPLE,
                "zone index must <= number of zones",
                str(self._nreactors + 1),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            if self._nreactors > 1:
                msg = [
                    Color.YELLOW,
                    "for multi-zone engine models,",
                    "zone index",
                    str(self._nreactors + 1),
                    "indicates the cylinder-averaged solution.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            exit()
        elif zoneID > self._nreactors:
            msg = ["Cylinder-averaged Solution"]
            Color.ckprint("info", msg)
        elif self._nreactors > 1:
            msg = ["Zone", str(zoneID), "Solution"]
            Color.ckprint("info", msg)

        # get solution sizes
        nreac, npoints = self.get_engine_solution_size(expectedzones)

        if npoints == 0 or nreac != expectedzones:
            raise ValueError
        else:
            self._numbsolutionpoints = npoints
        # create arrays to hold the raw solution data
        time = np.zeros(self._numbsolutionpoints, dtype=np.double)
        pres = np.zeros_like(time, dtype=np.double)
        temp = np.zeros_like(time, dtype=np.double)
        vol = np.zeros_like(time, dtype=np.double)
        # create a species mass fraction array to hold the solution species fraction profiles
        frac = np.zeros(
            (
                self.numbspecies,
                self._numbsolutionpoints,
            ),
            dtype=np.double,
            order="F",
        )
        # get raw solution data
        icreac = c_int(zoneID)
        icnpts = c_int(npoints)
        icnspec = c_int(self.numbspecies)
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetGasSolnResponse(
            icreac, icnpts, icnspec, time, temp, pres, vol, frac
        )
        if iErr != 0:
            msg = [
                Color.RED,
                "failed to fetch the raw solution data from memory,",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # store the raw solution data in a dictionary
        # time
        self._solution_rawarray["time"] = copy.deepcopy(time)
        # temperature
        self._solution_rawarray["temperature"] = copy.deepcopy(temp)
        # pressure
        self._solution_rawarray["pressure"] = copy.deepcopy(pres)
        # volume
        self._solution_rawarray["volume"] = copy.deepcopy(vol)
        # species mass fractions
        self.parsespeciessolutiondata(frac)
        # create soolution mixture
        iErr = self.create_solution_mixtures(frac)
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "forming solution mixtures",
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # clean up
        del time, pres, temp, vol, frac

    def process_average_engine_solution(self):
        """
        Post-process the ylinder averaged solution profiles from
        multi-zone engine models
        """
        # set the cylinder averge solution record ("zone") index
        meanzoneID = self._nreactors + 1
        # post-process mean solution
        self.process_engine_solution(zoneID=meanzoneID)

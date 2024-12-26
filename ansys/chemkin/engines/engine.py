import copy
from ctypes import c_double, c_int

import numpy as np

from chemkin import chemkin_wrapper
from chemkin.batchreactors.batchreactor import BatchReactors
from chemkin.chemistry import Patm
from chemkin.color import Color as Color
from chemkin.reactormodel import Keyword


class Engine(BatchReactors):
    """
    Common engine cylinder parameters used by Chemkin engine models
    """

    def __init__(self, reactor_condition, label=None):
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
        self.heattransfermodel = None
        self.heattransferparameters = None
        self.cylinderwalltemperature = None  # [K]
        # incylinder gas speed correlations
        # velocity correlation parameters
        # Woschni "GVEL <C11> <C12> <C2> <swirling ratio>"
        self.gasvelocity = None
        # Woschni+Huber IMEP "HIMP <IMEP>"
        self.HuberIMEP = None  # [atm]
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
    def convertCAtoTime(CA, startCA, RPM):
        """
        Convert the current crank angle value to simulation time
        :param CA: current engine crank angle [degree] (double scalar)
        :param startCA: starting crank angle, IVC timing [degree] (double scalar)
        :param RPM: engine speed RPM [revolutions oer minute] (double scalar)
        :return: time [sec] (double scalar)
        """
        if RPM <= 0.0:
            print(Color.PURPLE + "** engine speed RPM must > 0", end=Color.END)
            return 0.0
        #
        time = (CA - startCA) / RPM / 6.0e0
        if time < 0.0:
            print(
                Color.PURPLE + "** given CA is less than the starting CA @ IVC",
                end=Color.END,
            )
            return 0.0
        else:
            return time

    @staticmethod
    def convertTimetoCA(time, startCA, RPM):
        """
        Convert the current time to crank angle
        :param time: current simulation time [sec] (double scalar)
        :param startCA: starting crank angle, IVC timing [degree] (double scalar)
        :param RPM: engine speed RPM [revolutions oer minute] (double scalar)
        :return: engine crank angle [degree] (double scalar)
        """
        if time < 0.0:
            print(Color.PURPLE + "** simulation time must > 0", end=Color.END)
            return 0.0
        #
        CA = startCA + time * RPM * 6.0e0
        return CA

    def toTime(self, CA):
        """
        Convert the current crank angle value to simulation time
        :param CA: current engine crank angle [degree] (double scalar)
        :return: time [sec] (double scalar)
        """
        return (CA - self.IVCCA) / self.degpersec

    def toCA(self, time):
        """
        Convert the current time to crank angle
        :param time: current simulation time [sec] (double scalar)
        :return: engine crank angle [degree] (double scalar)
        """
        return self.IVCCA + time * self.degpersec

    @property
    def startingCA(self):
        """
        Get the simulation starting crank angle [degree]
        usually the starting CA ~ the intake valve close (IVC) timing
        """
        return self.IVCCA

    @startingCA.setter
    def startingCA(self, startCA):
        """
        Set the starting crank angle of engine simulation,
        usually this corresponds to the intake valve close (IVC) timing
        a positive starting CA implies the standard top dead center (TDC) is at 360 degrees CA
        a negative starting CA implies the standard TDC is at 0 degree CA
        :param startCA: starting crank angle [degree] (double scalar)
        """
        # set IVC timing in CA
        self.IVCCA = startCA
        # set keyword
        self._inputcheck.append("DEG0")

    @property
    def endingCA(self):
        """
        Get the simulation ending crank angle [degree]
        usually the ending CA ~ the exhaust valve open (EVO) timing
        """
        return self.EVOCA

    @endingCA.setter
    def endingCA(self, endCA):
        """
        Set the ending crank angle of engine simulation,
        usually this corresponds to the exhaust valve open (EVO) timing
        :param endCA: ending crank angle [degree] (double scalar)
        """
        # check EVO timing value
        if endCA <= self.startingCA:
            print(
                Color.PURPLE + f"** ending CA must > starting CA = {self.startingCA}",
                end=Color.END,
            )
            return
        # set EVO timing in CA
        self.EVOCA = endCA
        self.rundurationCA = self.endingCA - self.startingCA
        # set keyword
        self._inputcheck.append("DEGE")

    @property
    def durationCA(self):
        """
        Get the simulation duration in number of crank angles [degree]
        """
        return self.rundurationCA

    @durationCA.setter
    def durationCA(self, CA):
        """
        Set the engine simulation duration in CA
        :param CA: crank angle [degree] (double scalar)
        """
        # check EVO timing value
        if CA <= 0.0:
            print(Color.PURPLE + "** duration CA must > 0", end=Color.END)
            return
        # set EVO timing in CA
        self.rundurationCA = CA
        self.EVOCA = self.IVCCA + CA
        # set keyword
        self._inputcheck.append("DEGE")

    @property
    def bore(self):
        """
        Get the engine cylinder bore diameter
        :return: bore diameter [cm] (double scalar)
        """
        return self.borediam

    @bore.setter
    def bore(self, diameter):
        """
        Set the engine cylinder bore diameter
        :param diameter: bore diameter [cm] (double scalar)
        """
        if diameter > 0.0:
            self.borediam = diameter
            self.borearea = np.pi * diameter * diameter / 4.0e0
            # set keyword
            self._inputcheck.append("BORE")
        else:
            print(Color.PURPLE + "** bore diameter must > 0", end=Color.END)

    @property
    def stroke(self):
        """
        Get the engine stroke
        :return: stroke [cm] (double scalar)
        """
        return self.enginestroke

    @stroke.setter
    def stroke(self, s):
        """
        Set the engine stroke
        :param s: stroke [cm] (double scalar)
        """
        if s > 0.0:
            self.enginestroke = s
            self.crankradius = s / 2.0e0
            # set keyword
            self._inputcheck.append("STRK")
        else:
            print(Color.PURPLE + "** stroke must > 0", end=Color.END)

    @property
    def connectingrod(self):
        """
        Get the connecting rod length
        :return: length [cm] (double scalar)
        """
        return self.connectrodlength

    @connectingrod.setter
    def connectingrod(self, s):
        """
        Set the engine connecting rod length
        :param s: lemgth [cm] (double scalar)
        """
        if s > 0.0:
            self.connectrodlength = s
            # set keyword
            self._inputcheck.append("CRLEN")
        else:
            print(Color.PURPLE + "** connecting rod length must > 0", end=Color.END)

    @property
    def compressionratio(self):
        """
        Get the engine compression ratio
        :return: compression ratio [-] (double scalar)
        """
        return self.compressratio

    @compressionratio.setter
    def compressionratio(self, cratio):
        """
        Set the engine compression ratio
        :param cratio]: compression ratio [-] (double scalar)
        """
        if cratio > 1.0e0:
            self.compressratio = cratio
            # set keyword
            self._inputcheck.append("CMPR")
        else:
            print(Color.PURPLE + "** compression ratio must > 1", end=Color.END)

    @property
    def RPM(self):
        """
        Get the engine speed in RPM
        :return: engine speed [RPM] (double scalar)
        """
        return self.enginespeed

    @RPM.setter
    def RPM(self, speed):
        """
        Set the engine speed in RPM
        :param speed: engine speed [RPM] (double scalar)
        """
        if speed > 0.0:
            self.enginespeed = speed
            self.degpersec = speed * 6.0e0
            self.radpersec = self.degpersec * np.pi / 180.0e0
            # set keyword
            self._inputcheck.append("RPM")
        else:
            print(Color.PURPLE + "** RPM must > 0", end=Color.END)

    def setcylinderheadarea(self, area):
        """
        Set the cylinder head clearance surface area
        :param area: area [cm2] (double scalar)
        """
        if area > 0.0:
            self.cylinderheadarea = area
            self.headareas = area + self.pistonheadarea
            # set keyword
            if "BORE" in self._inputcheck:
                self.setkeyword(key="CYBAR", value=area / self.borearea)
            else:
                print(Color.PURPLE + "** please set BORE diameter first", end=Color.END)
        else:
            print(Color.PURPLE + "** surface area must > 0", end=Color.END)

    def setpistonheadarea(self, area):
        """
        Set the piston head clearance surface area
        :param area: area [cm2] (double scalar)
        """
        if area > 0.0:
            self.pistonheadarea = area
            self.headareas = area + self.cylinderheadarea
            # set keyword
            if "BORE" in self._inputcheck:
                self.setkeyword(key="PSBAR", value=area / self.borearea)
            else:
                print(Color.PURPLE + "** please set BORE diameter first", end=Color.END)
        else:
            print(Color.PURPLE + "** surface area must > 0", end=Color.END)

    def setpistonpinoffset(self, offset):
        """
        Set the piston head clearance surface area
        :param offset: piston pin offset distance [cm] (double scalar)
        """
        if offset < self.crankradius:
            self.pistonoffset = offset
            # set keyword
            self.setkeyword(key="POLEN", value=offset)
        else:
            print(
                Color.PURPLE
                + f"** piston pin offset distance must < crank radius {self.crankradius}",
                end=Color.END,
            )

    def getclearancevolume(self):
        """
        Get the clearance volume
        :return: clearance volume [cm3] (double scalar)
        """
        if "CMPR" in self._inputcheck:
            dvolume = self.getdisplacementvolume()
            cvolume = dvolume / (self.compressratio - 1.0e0)
        else:
            print(Color.PURPLE + "** please set RPM first", end=Color.END)
            cvolume = 0.0
        return cvolume

    def getdisplacementvolume(self):
        """
        Get the displacement volume
        :return: displacement volume [cm3] (double scalar)
        """
        return self.enginestroke * self.borearea

    def listengineparameters(self):
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
    def CAstepforsavingsolution(self):
        """
        Get the number of crank angles between saving the solution data [degree] (float scalar)
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
                print(
                    Color.YELLOW
                    + "** solution saving CA is not defined because 'ending CA' is not set",
                    end=Color.END,
                )
                return 0.0

    @CAstepforsavingsolution.setter
    def CAstepforsavingsolution(self, delta_CA):
        """
        Set the number of crank angles between saving the solution data
        :param delta_CA: number of crank angles between saving solution data [degree] (float scalar)
        """
        if delta_CA > 0.0e0:
            self.setkeyword(key="DEGSAVE", value=delta_CA)
        else:
            print(
                Color.PURPLE
                + "** solution saving number of crank angles value must > 0",
                end=Color.END,
            )

    @property
    def CAstepforprintingsolution(self):
        """
        Get the number of crank angles between printing the solution data to the text output file [degree] (float scalar)
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
                print(
                    Color.YELLOW
                    + "** solution printing CA is not defined because 'ending CA' is not set",
                    end=Color.END,
                )
                return 0.0

    @CAstepforprintingsolution.setter
    def CAstepforprintingsolution(self, delta_CA):
        """
        Set the timestep size between printing the solution data to the text output file
        :param delta_CA: number of crank angles between printing solution data [degree] (float scalar)
        """
        if delta_CA > 0.0e0:
            self.setkeyword(key="DEGPRINT", value=delta_CA)
        else:
            print(
                Color.PURPLE
                + "** solution printing number of crank angles value must > 0",
                end=Color.END,
            )

    def setminimumzonemass(self, minmass):
        """
        Set the minimum mass in a zone (for Spark Ignition and Direct Injection engine models)
        the default value is 1.0e-6 [g]
        :param minmass: minimum zonal mass [g] (double scalar)
        """
        if minmass > 0.0:
            # set keyword
            self.setkeyword(key="MLMT", value=minmass)
        else:
            print(Color.PURPLE + "** minimum zonal nmass must > 0.0", end=Color.END)
            exit()

    def setzonalgasratemultiplier(self, value=1.0e0, zoneID=None):
        """
        Set the value of the gas-phase reaction rate multiplier (optional)
        default value = 1.0
        :param value: gas-phase reaction rate multiplier (float scalar)
        :return: None
        """
        if value < 0.0:
            print(
                Color.PURPLE + "** reaction rate multiplier must be >= 0",
                end=Color.END,
            )
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

    def setwallheatransfer(self, model, HTparameters, walltemperature):
        """
        Set cylinder wall heat transfer model and parameters
        engine wall heat transfer models
        ICHX: dimensionless correlation "ICHX <a> <b> <c> <Twall>"
        ICHW: dimensional correlation "ICHW <a> <b> <c> <Twall>"
        ICHH: Hohenburg correlation "ICHH <a> <b> <c> <d> <e> <Twall>"
        :param model: engine wall heat transfer model (string)
        :param HTpaarmeters: model parameters correspond to the specified heat transfer model
        (list of double/float)
        :param walltemperature: cylinder wall/cooling oil temperature [K] (double scalar)
        """
        # check existing heat transfer set up
        if self.heattransfermodel is not None:
            print(
                Color.YELLOW + "** previous heat transfer model will be overwritten",
                end=Color.END,
            )
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
            print(
                Color.PURPLE
                + f"** engine wall heat transfer model {model.rstrip()} is not available"
            )
            print(
                "** the valid model options are 'dimensional', 'dimesionless', and 'hohenburg'",
                end=Color.END,
            )
            exit()
        # check number of parameters
        if len(HTparameters) != self.numbHTmodelparameters[self.heattransfermodel]:
            print(Color.PURPLE + "** incorrect number of parameters")
            print(
                f"** {model} requires {self.numbHTmodelparameters[self.heattransfermodel]} parameters"
            )
            print("** check Chemkin Input manual for more information", end=Color.END)
            exit()
        self.heattransferparameters = []
        # set model parameters
        self.heattransferparameters = copy.copy(HTparameters)
        # set cylinder wall temperature
        self.cylinderwalltemperature = walltemperature
        # set flag for wall heat transfer
        self._wallheattransfer = True

    def setgasvelocitycorrelation(self, gasvelparameters, IMEP=None):
        """
        Set the cylinder gas velocity correlation parameters
        Woschni: "GVEL <C11> <C12> <C2> <swirling ratio>"
        Huber: IMEP "HIMP <IMEP>"
        :param gasvelparameters: cylinder gas velocity correlation parameters (list of double)
        :param IMEP: indicated mean effective pressure used by the Huber gas velocity correlation [atm] (double scalar)
        """
        # check existing heat transfer set up
        if self.heattransfermodel is None:
            print(
                Color.YELLOW + "** please specify the wall heat transfer model first",
                end=Color.END,
            )
            return
        # check existing gas velocity correlation parameter
        if self.gasvelocity is not None:
            print(
                Color.YELLOW
                + "** previous gas velocity correlation will be overwritten",
                end=Color.END,
            )
        # check number of parameters
        if len(gasvelparameters) != 4:
            print(Color.PURPLE + "** incorrect number of parameters")
            print("** gas velocity correlation requires 4 parameters")
            print("** <C11> <C12> <C2> <swirl ratio>")
            print("** check Chemkin Input manual for more information", end=Color.END)
            exit()
        self.gasvelocity = []
        # set model parameters
        self.gasvelocity = copy.copy(gasvelparameters)
        # IEMP for the Huber correlation
        self.HuberIMEP = IMEP  # [atm]

    def setheattransferkeywords(self):
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

    def setenginekeywords(self):
        """
        Set engine parameter keywords under the Full-Keywords mode
        :return: None
        """
        self.setkeyword(key="BORE", value=self.borediam)
        self.setkeyword(key="STRK", value=self.enginestroke)
        self.setkeyword(key="CRLEN", value=self.connectrodlength)
        self.setkeyword(key="CMPR", value=self.compressratio)
        self.setkeyword(key="RPM", value=self.enginespeed)
        if self.pistonoffset > 0.0:
            self.setkeyword(key="POLEN", value=self.pistonoffset)
        self.setkeyword(key="DEG0", value=self.IVCCA)
        self.setkeyword(key="DEGE", value=self.EVOCA)

    def setengineconditionkeywords(self):
        """
        Set engine initial condition keywords under the Full-Keywords mode
        :return: None
        """
        self.setkeyword(key="PRES", value=self._pressure.value / Patm)
        self.setkeyword(key="TEMP", value=self._temperature.value)
        # initial mole fraction
        nspecieslines, species_lines = self.createspeciesinputlines(
            self._solvertype.value, threshold=1.0e-12, molefrac=self.reactormixture.X
        )
        for line in species_lines:
            self.setkeyword(key=line, value=True)

    def getengineheatrelease(self):
        """
        Get heat release crank angles from the engine solution
        """
        # heat loss rate per CA [erg/degree] sized = 1 + number of surface materials
        QLossRateCA = np.zeros(1, dtype=np.double)
        # apparent heat release rate per CA [erg/degree]
        AHRR = c_double(0.0)
        # apparent heat release rate per CA from PV-ConGamma [erg/degree]
        AHRRP = c_double(0.0)
        # Crank rotation angle for 10% of total heat release
        HR10 = c_double(self.IVCCA)
        # Crank rotation angle for 50% of total heat release
        HR50 = c_double(self.IVCCA)
        # Crank rotation angle for 90% of total heat release
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

    def getenginesolutionsize(self, expected):
        """
        Get the number of zones and the number of solution points
        :param expected: expected number of zonal + mean solution records (integer scalar)
        :return: nzones = number of zones, npoints = number of solution points (integer scalar, integer scalar)
        """
        # check run completion
        status = self.getrunstatus(mode="silent")
        if status == -100:
            print(
                Color.YELLOW + "** please run the reactor simulation first",
                end=Color.END,
            )
            raise
        elif status != 0:
            print(Color.YELLOW + "** simulation was failed")
            print(
                "** please correct the error and rerun the reactor simulation",
                end=Color.END,
            )
            raise
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
            print(
                Color.PURPLE + f"** failed to get solution size, error code = {iErr}",
                end=Color.END,
            )
            return nzones, 0
        else:
            # incorrect number of zones
            print(Color.PURPLE + f"** incorrect number of zones = {nzones}")
            print(f"   the model expects {expected} zones", end=Color.END)
            return nzones, 0

    def processenginesolution(self, zoneID=None):
        """
        Post-process solution to extract the raw solution variable data from
        engine simulation results
        :param zoneID: zone index (integer scalar)
        """
        # check existing raw data
        if self.getrawsolutionstatus():
            print(
                Color.YELLOW
                + "** solution has been processed before, existing solution data will be deleted",
                end=Color.END,
            )

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
            print(
                Color.PURPLE
                + f"** zone index must <= number of zones {self._nreactors}",
                end=Color.END,
            )
            exit()
        elif zoneID > self._nreactors:
            print(Color.YELLOW + "** cylinder-average solution", end=Color.END)
        elif self._nreactors > 1:
            print(Color.YELLOW + f"** zone {zoneID} solution", end=Color.END)

        # get solution sizes
        nreac, npoints = self.getenginesolutionsize(expectedzones)

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
            print(
                Color.PURPLE
                + f"** error: failed to fetch the raw solution data from memory, error code = {iErr}",
                end=Color.END,
            )
            raise
        # store the ratw solution data in a dictionary
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
        iErr = self.createsolutionmixtures(frac)
        if iErr != 0:
            print(
                Color.PURPLE + "** error: packaging solution mixtures",
                end=Color.END,
            )
            raise
        # clean up
        del time, pres, temp, vol, frac

    def processaverageenginesolution(self):
        """
        Post-process the ylinder averaged solution profiles from
        multi-zone engine models
        """
        # set the cylinder averge solution record ("zone") index
        meanzoneID = self._nreactors + 1
        # post-process mean solution
        self.processenginesolution(zoneID=meanzoneID)

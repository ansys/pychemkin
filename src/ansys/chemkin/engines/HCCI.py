import copy
from ctypes import c_double, c_int
import logging

import numpy as np

from chemkin import chemkin_wrapper
from chemkin.chemistry import (
    Patm,
    checkchemistryset,
    chemistrysetinitialized,
    setverbose,
)
from chemkin.color import Color as Color
from chemkin.engines.engine import Engine
from chemkin.reactormodel import Keyword

logger = logging.getLogger(__name__)


class HCCIengine(Engine):
    """
    Single or multi- zone homogeneous charge compression ignition (HCCI) engine
    """

    def __init__(self, reactor_condition, label=None, nzones=None):
        """
        Single or multi- zone homogeneous charge compression ignition (HCCI) engine
        :param nzones: number of zones in the HCCI engine model (integer scalar)
        """
        # set default number of zone(s): single-zone
        if nzones is None:
            nzones = 1
        # set default label
        if label is None:
            if nzones == 1:
                label = "HCCI"
            else:
                label = "Multi-Zone HCCI"

        # use the first zone to initialize the engine model
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("HCCI"))
        self._solvertype = c_int(self.SolverTypes.get("Transient"))
        self._problemtype = c_int(self.ProblemTypes.get("ICEN"))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY"))
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
        self._zonalsetupmode = 0
        # zonal temperature values
        self.zonetemperature = []
        # zonal volume fractions
        self.zonevolume = []
        # zonal mass fractions
        self.usezonemass = False
        self.zonemass = []
        # zonal wall heat transfer area fraction
        self.zoneHTarea = []
        # zonal gas compositions in mole fraction (for zonalsetupmode =1)
        self.zonemolefrac = []  # list of mole fraction arrays
        # zonal equivalence ratios (for zonalsetupmode =2)
        self.zoneequivalenceratio = []
        # fuel composition for all zones
        self.zonefueldefined = []
        # oxidizer composition for all zones
        self.zoneoxiddefined = []
        # product composition for all zones
        self.zoneproductdefined = []
        # zonal additive gas composition
        self.zoneaddmolefrac = []
        # zonal EGR ratios
        self.zoneEGRR = []
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(155)
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
            print(
                Color.RED
                + f"** error initializing the HCCI engine model: {self.label:s}"
            )
            print(f"   error code = {iErr:d}", end=Color.END)
            exit()

    def getnumberofzones(self):
        """
        Get the number of zones used by the current HCCI simulation
        """
        return self._nzones.value

    def setzonaltemperature(self, zonetemp):
        """
        set zonal temperatures for muti-zone HCCI engine simulation
        :param zonetemp: zonal temperatures (list of doubles)
        """
        nzones = self._nzones.value
        if len(zonetemp) != nzones:
            print(
                Color.PURPLE
                + f"** zonal temperature must be a list double/float of size {nzones}",
                end=Color.END,
            )
            exit()

        if len(self.zonetemperature) > 0:
            print(Color.YELLOW + "** zonal temperatures will be reset", end=Color.END)
        self.zonetemperature = []
        # set zonal temperatures
        for t in zonetemp:
            self.zonetemperature.append(t)
        # set zonal definition mode
        if self._zonalsetupmode == 0:
            self._zonalsetupmode = 1

    def setzonalvolume(self, zonevol):
        """
        set zonal volume fractions for muti-zone HCCI engine simulation
        :param zonevol: zonal volume fractions (list of doubles)
        """
        nzones = self._nzones.value
        if len(zonevol) != nzones:
            print(
                Color.PURPLE
                + f"** zonal volume fraction must be a list of double/float of size {nzones}",
                end=Color.END,
            )
            exit()

        if len(self.zonevolume) > 0:
            print(
                Color.YELLOW + "** zonal volume farctions will be reset", end=Color.END
            )
        self.zonevolume = []
        # set zonal volume fractions (will be normalized)
        for v in zonevol:
            self.zonevolume.append(v)

    def setzonalmassfraction(self, zonemass):
        """
        set zonal mass fractions for muti-zone HCCI engine simulation
        :param zonemass: zonal mass fractions (list of doubles)
        """
        nzones = self._nzones.value
        if len(zonemass) != nzones:
            print(
                Color.PURPLE
                + f"** zonal mass fraction must be a list of double/float of size {nzones}",
                end=Color.END,
            )
            exit()

        if len(self.zonemass) > 0:
            print(Color.YELLOW + "** zonal mass farctions will be reset", end=Color.END)
        self.zonemass = []
        # set zonal mass fractions (will be normalized)
        for v in zonemass:
            self.zonemass.append(v)
        # set flag
        self.usezonemass = True

    def setzonalheattransferarea(self, zonearea):
        """
        set zonal wall heat transfer area fractions for muti-zone HCCI engine simulation
        :param zonearea: zonal heat transfer area fractions (list of doubles)
        """
        nzones = self._nzones.value
        if len(zonearea) != nzones:
            print(
                Color.PURPLE
                + f"** zonal area fraction must be a list of double/float of size {nzones}",
                end=Color.END,
            )
            exit()

        if len(self.zoneHTarea) > 0:
            print(Color.YELLOW + "** zonal area farctions will be reset", end=Color.END)
        self.zoneHTarea = []
        # set zonal wall heat transfer area fractions (will be normalized)
        for a in zonearea:
            self.zoneHTarea.append(a)

    def setzonalmolefractions(self, zonemolefrac):
        """
        set zonal gas mole fractions for muti-zone HCCI engine simulation
        :param zonemolefrac: zonal gas mole fractions (list of double arrays of size = number of gas species)
        """
        nzones = self._nzones.value
        if len(zonemolefrac) != nzones:
            print(
                Color.PURPLE
                + "** zonal gas mole fraction must be a list of "
                + f" {nzones} double/float mole fraction arrays",
                end=Color.END,
            )
            exit()

        if self._zonalsetupmode == 2:
            print(
                Color.YELLOW
                + "** raw gas composition will replace equivalence ratio"
                + " to set up the zonal gas compositions",
                end=Color.END,
            )

        if len(self.zonemolefrac) > 0:
            print(
                Color.YELLOW + "** zonal gas mole farctions will be reset",
                end=Color.END,
            )
        self.zonemolefrac = []
        # set zonal definition mode
        self._zonalsetupmode = 1
        # set zonal gas mole fractions
        for x in zonemolefrac:
            # x must be a double array of size = number of gas species
            self.zonemolefrac.append(x)

    def definefuel(self, recipe):
        """
        set the fuel composition for setting up zonal gas composition by zonal equivalence ratio
        :param recipe: list of tuples for (species, mole fraction) pairs
        """
        if len(self.zonefueldefined) > 0:
            print(
                Color.YELLOW + "** previous fuel definition will be reset",
                end=Color.END,
            )
        self.definefuel = []
        self.definefuel = copy.deepcopy(recipe)

    def defineoxid(self, recipe):
        """
        set the oxidizer composition for setting up zonal gas composition by zonal equivalence ratio
        :param recipe: list of tuples for (species, mole fraction) pairs
        """
        if len(self.zoneoxiddefined) > 0:
            print(
                Color.YELLOW + "** previous oxidizer definition will be reset",
                end=Color.END,
            )
        self.defineoxid = []
        self.defineoxid = copy.deepcopy(recipe)

    def defineproduct(self, products):
        """
        set the complete combustion product species for setting up zonal gas composition by zonal equivalence ratio
        :param products: list of product species symbols
        """
        if len(self.zoneproductdefined) > 0:
            print(
                Color.YELLOW + "** previous product definition will be reset",
                end=Color.END,
            )
        self.defineproduct = []
        self.defineproduct = copy.deepcopy(products)

    def defineaddfractions(self, addfrac):
        """
        set zonal additive gas mole fractions for setting up zonal gas composition by zonal equivalence ratio
        :param addfrac: additive gas mole fractions (list of double arrays of size = number of gas species)
        """
        nzones = self._nzones.value
        if len(addfrac) != nzones:
            print(
                Color.PURPLE
                + "** zonal additive gas mole fraction must be a list of "
                + f"{nzones} double/float mole fraction arrays",
                end=Color.END,
            )
            exit()

        if len(self.zoneaddmolefrac) > 0:
            print(
                Color.YELLOW + "** zonal additive gas mole farctions will be reset",
                end=Color.END,
            )
        self.zoneaddmolefrac = []
        # set zonal gas mole fractions
        for x in addfrac:
            # x must be a double array of size = number of gas species
            self.zoneaddmolefrac.append(x)

    def setzonalequivalenceratio(self, zonephi):
        """
        set zonal wall heat transfer area fractions for setting up zonal gas composition by zonal equivalence ratio
        :param zonephi: zonal equivalence ratios (list of doubles)
        """
        nzones = self._nzones.value
        if len(zonephi) != nzones:
            print(
                Color.PURPLE
                + f"** zonal equivalence ratio must be a list of double/float of size {nzones}",
                end=Color.END,
            )
            exit()

        if self._zonalsetupmode == 1:
            print(
                Color.YELLOW
                + "** equivalence ratio will replace raw gas mole fractions"
                + " to set up the zonal gas compositions",
                end=Color.END,
            )

        if len(self.zoneequivalenceratio) > 0:
            print(
                Color.YELLOW + "** previous zonal equivalence ratio will be reset",
                end=Color.END,
            )
        self.zoneequivalenceratio = []
        # set zonal definition mode
        self._zonalsetupmode = 2
        # set zonal equivalence ratio
        for p in zonephi:
            self.zoneequivalenceratio.append(p)

    def setzonalEGRratio(self, zoneegr):
        """
        set zonal exhaust gas recirculation (EGR) ratios for setting up zonal gas composition by zonal equivalence ratio
        :param zoneegr: zonal EGR ratios (list of doubles)
        """
        nzones = self._nzones.value
        if len(zoneegr) != nzones:
            print(
                Color.PURPLE
                + f"** zonal EGR ratio must be a list of double/float of size {nzones}",
                end=Color.END,
            )
            exit()

        if len(self.zoneEGRR) > 0:
            print(
                Color.YELLOW + "** previous zonal EGR ratio will be reset",
                end=Color.END,
            )
        self.zoneEGRR = []
        # set zonal EGR ratio
        for r in zoneegr:
            self.zoneEGRR.append(r)

    def setenergyequationswitchonCA(self, switchCA):
        """
        Set the crank angle at which the energy equation will be turn ON for
        the rest of the simulation.
        Before this switch crank angle the given temperature profile(s) or value(s)
        is used in the multi-zone HCCI simulation.
        :param switchCA: crank angle [degree] (double scalar)
        """
        if self._nzones.value == 1:
            print(
                Color.PURPLE
                + "**  this parameter is for the multi-zone HCCI engine only",
                end=Color.END,
            )
            return

        if switchCA > self.IVCCA:
            # set keyword
            self.setkeyword(key="ASWH", value=switchCA)
        else:
            print(
                Color.PURPLE + f"** energy switch on CA must > IVC {self.IVCCA}",
                end=Color.END,
            )

    def setzonalvolumekeyword(self):
        """
        Set zonal volume keyword for the multi-zone HCCI engine simulation
        """
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

    def setzonalmasskeyword(self):
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

    def setzonalconditionkeywords(self):
        """
        Set zonal initial condition keywords under the Full-Keywords mode
        and use raw species mole fractions to set up zonal gas compositions
        for multi-zone HCCI engine simulation
        """
        if self._nzones.value == 1:
            # single zone is not allowed here
            print(
                Color.YELLOW
                + "** single-zone should not use this method to set up zonal conditions",
                end=Color.END,
            )
            return

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

    def setzonalequivalenceratiokeywords(self):
        """
        Set zonal initial condition keywords under the Full-Keywords mode
        and use equivalence ratios to set up zonal gas compositions
        for multi-zone HCCI engine simulation
        """
        if self._nzones.value == 1:
            # single zone is not allowed here
            print(
                Color.YELLOW
                + "** single-zone should not use this method to set up zonal conditions",
                end=Color.END,
            )
            return

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
                for s, x in self.definefuel:
                    keyline = (
                        "FUEL" + Keyword.fourspaces + s + Keyword.fourspaces + str(x)
                    )
                    self.setkeyword(key=keyline, value=True)
                # set oxidizer composition by using recipe
                for s, x in self.defineoxid:
                    keyline = (
                        "OXID" + Keyword.fourspaces + s + Keyword.fourspaces + str(x)
                    )
                    self.setkeyword(key=keyline, value=True)
                # define complete combustion products by using species symbol list
                for s in self.defineproduct:
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

    def __process_keywords_withFullInputs(self):
        """
        Process input keywords for the reactor model under the Full-Keyword mode
        :return: Error code (integer scalar)
        """
        iErr = 0
        # setverbose(True)
        # verify required inputs
        iErr = self.inputvalidation()
        if iErr != 0:
            print(
                Color.PURPLE + "** missing required input variable",
                end=Color.END,
            )
            return iErr
        # prepare initial conditions
        # initial mass fraction
        Y_init = self.reactormixture.Y
        # connecting rod length to crank radius ratio
        LOLR = c_double(self.connectingrod / self.crankradius)
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
                self._heatlossrate,
                Y_init,
            )
            iErr += iErrc

        # set reactor type
        self.setreactortypekeywords()
        # set number of zones
        self.setkeyword(key="NZONE", value=self._nzones.value)
        # set engine parameter
        self.setenginekeywords()
        # check if the wall heat transfer model is set up
        if self._wallheattransfer:
            self.setheattransferkeywords()
        #
        if self._nzones.value == 1 or self._zonalsetupmode == 0:
            # single-zone HCCI engine initial condition or
            # multi-zone with uniform zonal properties
            if self._nzones.value > 1:
                if self.usezonemass:
                    # set zonal mass fractions (required for the multi-zone model)
                    self.setzonalmasskeyword()
                else:
                    if len(self.zonevolume) != self._nzones.value:
                        # zonal volumes are not set
                        # set uniform zonal volume fractions
                        volfrac = 1.0 / self._nzones.value
                        self.zonevolume = []
                        for i in range(self._nzones.value):
                            self.zonevolume.append(volfrac)
                    # set zonal volume fractions (required for the multi-zone model)
                    self.setzonalvolumekeyword()
            # set uniform cylinder properties
            self.setengineconditionkeywords()
        else:
            # multi-zone HCCI engine zonal conditions
            if self._zonalsetupmode == 0:
                if self.usezonemass:
                    # set zonal mass fractions (required for the multi-zone model)
                    self.setzonalmasskeyword()
                else:
                    if len(self.zonevolume) != self._nzones.value:
                        # zonal volumes are not set
                        # use equal zonal volume fractions
                        volfrac = 1.0 / self._nzones.value
                        self.zonevolume = []
                        for i in range(self._nzones.value):
                            self.zonevolume.append(volfrac)
                # uniform zonal properties
                self.setengineconditionkeywords()
            elif self._zonalsetupmode == 1:
                # non-uniform zonal properties with zonal raw gas mole fractions provided
                self.setzonalconditionkeywords()
            else:
                # non-uniform zonal properties with zonal equivalence ratios provided
                self.setzonalequivalenceratiokeywords()
        #
        if iErr == 0 and self._numbprofiles > 0:
            # get keyword lines of all profiles
            iErrProf, nproflines, prof_lines = self.createprofileinputlines()
            iErr += iErrProf
            if iErrProf == 0:
                # set the profile keywords
                for p in prof_lines:
                    for line in p:
                        self.setkeyword(key=line, value=True)
        # solve integrated heat release rate due to chemical reactions
        self.setkeyword(key="QRGEQ", value=True)
        # add the END keyword
        self.setkeyword(key="END", value=True)
        # create input lines from additional user-specified keywords
        iErr, nlines = self.createkeywordinputlines()
        print(Color.YELLOW + f"** {nlines} input lines are created", end=Color.END)

        return iErr

    def __run_model_withFullInputs(self, **kwargs):
        """
        Run the reactor model after the keywords are processed under the Full-Keyword mode
        All keywords must be assigned
        :param kwargs: command arguments
        :return: error code (integer scalar)
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

        return iErr

    def __process_keywords(self):
        """
        Process input keywords for the reactor model
        :return: Error code (integer scalar)
        """
        iErr = 0
        iErrc = 0
        iErrKey = 0
        iErrInputs = 0
        setverbose(True)
        # verify required inputs
        iErr = self.inputvalidation()
        if iErr != 0:
            print(
                Color.PURPLE + "** missing required input variable",
                end=Color.END,
            )
            return iErr
        # prepare initial conditions
        # initial mass fraction
        Y_init = self.reactormixture.Y
        # connecting rod length to crank radius ratio
        LOLR = c_double(self.connectingrod / self.crankradius)
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
                self._heatlossrate,
                Y_init,
            )
            iErr += iErrc
            # heat transfer (use additional keywords)
            # solver parameters (use additional keywords)
            # output controls (use additional keywords)
            # ROP (use additional keywords)
            # sensitivity (use additional keywords)
            # ignition delay (use additional keywords)
            # solve integrated heat release rate due to chemical reactions
            iErrc = chemkin_wrapper.chemkin.KINAll0D_IntegrateHeatRelease()
            iErr += iErrc
        else:
            pass
        # check if the wall heat transfer model is set up
        if iErr == 0 and self._wallheattransfer:
            self.setheattransferkeywords()
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
        if iErr == 0:
            # set additional keywords
            # create input lines from additional user-specified keywords
            iErrInputs, nlines = self.createkeywordinputlines()
            print(
                Color.YELLOW + f"** {nlines} additional keywords are created",
                end=Color.END,
            )
            if iErrInputs == 0:
                # process additional keywords in _keyword_index and _keyword_lines
                for s in self._keyword_lines:
                    # convert string to byte
                    line = bytes(s, "utf-8")
                    # set additional keyword one by one
                    iErrKey = chemkin_wrapper.chemkin.KINAll0D_SetUserKeyword(line)
            else:
                print(
                    Color.RED
                    + "** error processing additional keywords. error code = {iErrInputs}",
                    end=Color.END,
                )
        #
        iErr = iErr + iErrInputs + iErrKey

        return iErr

    def __run_model(self, **kwargs):
        """
        Run the reactor model after the keywords are processed
        :param kwargs: command arguments
        :return: error code (integer scalar)
        """
        # run the simulation without keyword inputs
        iErr = chemkin_wrapper.chemkin.KINAll0D_Calculate(self._chemset_index)
        return iErr

    def run(self, **kwargs):
        """
        Generic Chemkin run reactor model method
        :param kwargs: arguments from the run command
        :return: Error code (integer scalar)
        """
        logger.debug("Running " + str(self.__class__.__name__) + " " + self.label)
        print(
            Color.YELLOW + f"** running model {self.__class__.__name__} {self.label}..."
        )
        print(
            f"   initialization = {checkchemistryset(self._chemset_index.value)}",
            end=Color.END,
        )
        if not checkchemistryset(self._chemset_index.value):
            # KINetics is not initialized: reinitialize KINetics
            print(Color.YELLOW + "** initializing chemkin...", end=Color.END)
            retVal = chemkin_wrapper.chemkin.KINInitialize(
                self._chemset_index, c_int(0)
            )
            if retVal != 0:
                print(
                    Color.RED + f"** error processing the keywords (code = {retVal:d})",
                    end=Color.END,
                )
                logger.debug(f"Initializing KINetics failed (code={retVal})")
                return retVal
            else:
                chemistrysetinitialized(self._chemset_index.value)

        for kw in kwargs:
            logger.debug("Reactor model argument " + kw + " = " + str(kwargs[kw]))

        # output initialization
        logger.debug("Clearing output")
        self.output = {}

        # keyword processing
        logger.debug("Processing keywords")
        print(Color.YELLOW + "** processing keywords", end=Color.END)
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
            print(
                Color.RED + f"** error processing the keywords (code = {retVal:d})",
                end=Color.END,
            )
            logger.debug(f"Processing keywords failed (code={retVal})")
            return retVal
        logger.debug("Processing keywords complete")

        # run reactor model
        logger.debug("Running model")
        print(Color.YELLOW + "** running model", end=Color.END)
        if self._nzones.value == 1 and Keyword.noFullKeyword:
            # single-zone HCCI
            # use API calls
            retVal = self.__run_model(**kwargs)
        else:
            # multi-zone HCCI
            # use full keywords
            retVal = self.__run_model_withFullInputs(**kwargs)
        # update run status
        self.setrunstatus(code=retVal)

        logger.debug("Running model complete, status = " + str(retVal))

        return retVal

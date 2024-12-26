import copy
from ctypes import c_double, c_int
import logging

import numpy as np

from chemkin import chemkin_wrapper
from chemkin.batchreactors.batchreactor import BatchReactors
from chemkin.chemistry import checkchemistryset, chemistrysetinitialized, setverbose
from chemkin.color import Color as Color
from chemkin.inlet import Inlet
from chemkin.reactormodel import Keyword

logger = logging.getLogger(__name__)


class PlugFlowReactor(BatchReactors):
    """
    Generic Plug Flow Reactor (PFR) model with energy equation
    """

    def __init__(self, inlet, label=None):
        # set default label
        if label is None:
            label = "PFR"
        # check Inlet
        if isinstance(inlet, Inlet):
            # initialization
            super().__init__(inlet, label)
        else:
            # wrong argument type
            print(
                Color.RED + "** the first argument must be an Inlet object",
                end=Color.END,
            )
            raise TypeError

        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("PFR"))
        self._solvertype = c_int(self.SolverTypes.get("Transient"))
        self._problemtype = c_int(self.ProblemTypes.get("CONP"))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY"))
        # defaults for all plug flow reactor models
        self._nreactors = 1
        self._npsrs = c_int(self._nreactors)
        self._ninlets = c_int(0)
        # number of zones
        self._nzones = c_int(0)
        # use API mode for PFR simulations
        Keyword.noFullKeyword = True
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(157)
        #
        # starting position [cm]
        self.startposition = c_double(0.0)
        # reactor length [cm]
        self.reactorlength = c_double(0.0)
        # reactor diameter [cm]
        self.reactordiameter = c_double(0.0)
        # cross-sectional flow area [cm2]
        self.reactorflowarea = 0.0
        # flow area given for the inlet
        if inlet._haveflowarea:
            self.reactorflowarea = inlet._flowarea
            # compute reactor diameter from the flow area
            self.reactordiameter = c_double(np.sqrt(4.0 * self.reactorflowarea / np.pi))
            if self.reactorflowarea <= 0.0:
                print(Color.YELLOW + "** inlet flow area is not set", end=Color.END)
        # inlet mass flow rate [g/sec]
        self._massflowrate = c_double(0.0)
        self._flowrate = 0.0
        if inlet._flowratemode < 0:
            # no given in the inlet
            print(Color.PURPLE + "** inlet flow rate is not set")
            print("   specify flow rate of the 'Inlet' object", end=Color.END)
        else:
            self.inletflowratemode = inlet._flowratemode
            self._flowrate = inlet._inletflowrate[inlet._flowratemode]
            self.inletflowrate = copy.deepcopy(inlet._inletflowrate)
            if self._flowrate <= 0.0:
                print(Color.PURPLE + "** inlet flow rate is not set correctly")
                print("   specify flow rate of the 'Inlet' object", end=Color.END)
        # solver parameters
        self._absolutetolerance = 1.0e-12
        self._relativetolerance = 1.0e-6
        # required inputs: (1) reactor length (2) flow area
        self._numb_requiredinput = 2
        self._requiredlist = ["XEND", "AREAF"]
        # always calculate the residence time
        self.setkeyword(key="RTIME", value="ON")
        # solve the momentum equation in most cases
        # turn the momentum equation OFF
        # when the velocity or the pressure profile along the reactor is given
        self.setkeyword(key="MOMEN", value="ON")

    @property
    def length(self):
        """
        Get reactor length (required) [cm] (float scalar)
        """
        return self.reactorlength.value

    @length.setter
    def length(self, length=0.0e0):
        """
        Set reactor length (required)
        default value = 0.0 cm
        :param length: reactor length [cm] (float scalar)
        """
        if length <= 0.0e0:
            print(
                Color.PURPLE + "** flow reactor length must be > 0",
                end=Color.END,
            )
        else:
            self._inputcheck.append("XEND")
            self.reactorlength = c_double(length)

    def setstartposition(self, x0):
        """
        Set the PFR simulation starting position
        reactor inlet: x0 = 0.0
        :param x0: starting position [cm] (double scalar)
        """
        if x0 >= self.reactorlength.value:
            print(
                Color.PURPLE + "** starting position must < reactor length",
                end=Color.END,
            )
        else:
            self.startposition = c_double(x0)

    @property
    def diameter(self):
        """
        Reactor diameter [cm]
        """
        return self.reactordiameter.value

    @diameter.setter
    def diameter(self, diam):
        """
        Set the PFR diameter
        :param diam: reactor diameter [cm] (double scalar)
        """
        if diam <= 0.0:
            print(
                Color.PURPLE + "** reactor diameter must > 0",
                end=Color.END,
            )
        else:
            self._inputcheck.append("AREAF")
            self.reactordiameter = c_double(diam)
            # set flow area at the inlet
            self.reactormixture._haveflowarea = True
            area = np.pi * diam * diam / 4.0
            self.reactormixture._flowarea = area
            self.reactorflowarea = area

    def setdiameterprofile(self, x, diam):
        """
        Specify plug-flow reactor diameter profile
        :param x: position value of the profile data [cm] (double array)
        :param diam: PFR diameter value of the profile data [cm] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.ReactorTypes.get(self._reactortype.value) != "PFR":
            print(
                Color.PURPLE
                + "** cannot specify reactor diameter of a non plug-flow reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "DPRO"
            iErr = self.setprofile(key=keyword, x=x, y=diam)
            if iErr == 0:
                self._inputcheck.append("AREAF")
                # set flow area at the inlet
                self.reactordiameter = c_double(diam[0])
                self.reactormixture._haveflowarea = True
                area = np.pi * diam * diam / 4.0
                self.reactormixture._flowarea = area
                self.reactorflowarea = area
            return iErr

    @property
    def flowarea(self):
        """
        Cross-sectional flow area of the PFR [cm2]
        """
        return self.reactorflowarea

    @flowarea.setter
    def flowarea(self, area):
        """
        Set the cross-sectional flow area of the PFR
        :param area: flow area [cm2] (double scalar)
        """
        if area <= 0.0:
            print(
                Color.PURPLE + "** flow area must > 0",
                end=Color.END,
            )
        else:
            # set the flow area keyword
            self._inputcheck.append("AREAF")
            self.setkeyword(key="AREAF", value=area)
            self.reactorflowarea = area
            # set flow area at the inlet
            diam = np.sqrt(4.0 * area / np.pi)
            self.reactordiameter = c_double(diam)
            self.reactormixture._haveflowarea = True
            self.reactormixture._flowarea = area

    def setflowareaprofile(self, x, area):
        """
        Specify plug-flow reactor cross-sectional flow area profile
        :param x: position value of the profile data [cm] (double array)
        :param area: PFR flow area value of the profile data [cm2] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.ReactorTypes.get(self._reactortype.value) != "PFR":
            print(
                Color.PURPLE
                + "** cannot specify reactor flow area of a non plug-flow reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "AFLO"
            iErr = self.setprofile(key=keyword, x=x, y=area)
            if iErr == 0:
                self._inputcheck.append("AREAF")
                self.setkeyword(key="AREAF", value=area[0])
                # set flow area at the inlet
                diam = np.sqrt(4.0 * area[0] / np.pi)
                self.reactordiameter = c_double(diam)
                self.reactormixture._haveflowarea = True
                self.reactormixture._flowarea = area[0]
                self.reactorflowarea = area[0]
            return iErr

    def setinletviscosity(self, visc):
        """
        Set the gas mixture viscocity at the PFR inlet
        :param visc: mixture viscosity [g/cm-sec] or [Poise] (double scalar)
        """
        if visc <= 0.0:
            print(
                Color.PURPLE + "** gas mixture viscosity must > 0",
                end=Color.END,
            )
        else:
            # set the mixture viscosity keyword
            self.setkeyword(key="VISC", value=visc)

    def setsolvermaxtimestepsize(self, size):
        """
        Set the maximum time step size allowed by the solver
        :param size: step size [cm]
        """
        if size > 0.0e0:
            self.setkeyword(key="DXMX", value=size)
        else:
            print(
                Color.PURPLE + "** solver timestep size must > 0",
                end=Color.END,
            )

    def setpseudosurfacevelocity(self, vel):
        """
        Set the pseudo surface velocity at the reactive surface
        to improve convergence due to surface chemistry stiffness
        Note: set this parameter only when having convergence issue with surface chemistry
        :param vel: surface velocity [cm/sec]
        """
        if vel > 0.0e0:
            self.setkeyword(key="PSV", value=vel)
        else:
            print(
                Color.PURPLE + "** pseudo velocity must > 0",
                end=Color.END,
            )

    @property
    def massflowrate(self):
        """
        Get plug flow reactor inlet mass flow rate [g/sec]
        """
        return self.reactormixture.massflowrate

    @property
    def velocity(self):
        """
        Get plug flow reactor inlet velocity [cm/sec]
        """
        return self.reactormixture.velocity

    @property
    def volflowrate(self):
        """
        Get plug flow reactor inlet volumetric flow rate [cm3/sec]
        """
        return self.reactormixture.volflowrate

    @property
    def sccm(self):
        """
        Get plug flow reactor inlet volumetric flow rate in SCCM [standard cm3/min]
        """
        return self.reactormixture.sccm

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
        # prepare inlet conditions
        # get inlet mass flow rate
        self._massflowrate = c_double(self.massflowrate)
        # inlet mass fraction
        Y_init = self.reactormixture.Y
        # surface sites (not applicable)
        Site_init = np.zeros(1, dtype=np.double)
        # bulk activities (not applicable)
        Bulk_init = np.zeros_like(Site_init, dtype=np.double)
        # set reactor inlet conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("PFR"):
            iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupPFRInputs(
                self._chemset_index,
                self.startposition,
                self.reactorlength,
                self._temperature,
                self._pressure,
                self._heatlossrate,
                self.reactordiameter,
                Site_init,
                Bulk_init,
                self._massflowrate,
                Y_init,
            )
            iErr += iErrc
            # turn OFF the momentum equation when pressure profile is set
            if self._numbprofiles > 0:
                if "PPRO" in self._profiles_list or "VELPRO" in self._profiles_list:
                    self.setkeyword(key="MOMEN", value="OFF")
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
        if Keyword.noFullKeyword:
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
        if Keyword.noFullKeyword:
            # use API calls
            retVal = self.__run_model(**kwargs)
        else:
            # use full keywords
            retVal = self.__run_model_withFullInputs(**kwargs)
        # update run status
        self.setrunstatus(code=retVal)

        logger.debug("Running model complete, status = " + str(retVal))

        return retVal


class PlugFlowReactor_EngeryEquation(PlugFlowReactor):
    """
    Plug Flow Reactor (PFR) model with energy equation
    """

    def __init__(self, inlet, label=None):
        # set default label
        if label is None:
            label = "PFR"
        # check Inlet
        if isinstance(inlet, Inlet):
            # initialization
            super().__init__(inlet, label)
        else:
            # wrong argument type
            print(
                Color.RED + "** the first argument must be an Inlet object",
                end=Color.END,
            )
            raise TypeError

        # set reactor type
        self._energytype = c_int(self.EnergyTypes.get("ENERGY"))
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        self._heattransfercoefficient = 0.0e0
        self._ambienttemperature = 3.0e2
        # external heat transfer area per reactor length [cm2/cm]
        self._heattransferarea = 0.0e0
        # set up basic batch reactor parameters
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
            # setup reactor model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            print(
                Color.RED + f"** error initializing the reactor model: {self.label:s}"
            )
            print(f"   error code = {iErr:d}", end=Color.END)
            exit()

    @property
    def heatlossrate(self):
        """
        Get heat loss rate from the reactor to the surroundings [cal/sec] (float scalar)
        default value = 0.0 cal/sec
        """
        return self._heatlossrate.value

    @heatlossrate.setter
    def heatlossrate(self, value):
        """
        Set the heat loss rate per length from the reactor to the surroundings (required)
        default value = 0.0 cal/sec-cm
        :param value: heat loss rate [cal/sec-cm] (float scalar)
        """
        self._heatlossrate = c_double(value)
        if not Keyword.noFullKeyword:
            self.setkeyword(key="QLOS", value=value)

    @property
    def heattransfercoefficient(self):
        """
        Get heat transfer coefficient between the reactor and the surroundings [cal/cm2-K-sec] (float scalar)
        default value = 0.0 cal/cm2-K-sec
        """
        return self._heattransfercoefficient

    @heattransfercoefficient.setter
    def heattransfercoefficient(self, value=0.0e0):
        """
        Set heat transfer coefficient between the reactor and the surroundings (optional)
        default value = 0.0 cal/cm2-K-sec
        :param value: heat transfer coefficient [cal/cm2-K-sec] (float scalar)
        """
        if value < 0.0e0:
            print(
                Color.PURPLE + "** simulation end time must be >= 0",
                end=Color.END,
            )
        else:
            self._heattransfercoefficient = value
            # set the corresponding keyword
            self.setkeyword(key="HTC", value=value)

    @property
    def ambienttemperature(self):
        """
        Get ambient temperature (optional) [K] (float scalar)
        default value = 300 K
        """
        return self._ambienttemperature

    @ambienttemperature.setter
    def ambienttemperature(self, value=0.0e0):
        """
        Set ambient temperature (optional)
        default value = 300 K
        :param value: ambient temperature [K] (float scalar)
        """
        if value <= 0.0e0:
            print(
                Color.PURPLE + "** ambient temperature must be > 0",
                end=Color.END,
            )
        else:
            self._ambienttemperature = value
            # set the corresponding keyword
            self.setkeyword(key="TAMB", value=value)

    @property
    def heattransferarea(self):
        """
        Get heat transfer area per length between the reactor and the surroundings (optional) [cm2/cm] (float scalar)
        default value = 0.0 cm2/cm
        """
        return self._heattransferarea

    @heattransferarea.setter
    def heattransferarea(self, value=0.0e0):
        """
        Set heat transfer area per length between the reactor and the surroundings (optional)
        default value = 0.0 cm2
        :param value: heat transfer area [cm2/cm] (float scalar)
        """
        if value < 0.0e0:
            print(
                Color.PURPLE + "** heat transfer area must be >= 0",
                end=Color.END,
            )
        else:
            self._heattransferarea = value
            # set the corresponding keyword
            self.setkeyword(key="AREAQ", value=value)

    def setheattransferareaprofile(self, x, area):
        """
        Specify reactor heat transfer area per reactor length profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param area: heat transfer area value of the profile data [cm2/cm] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.EnergyTypes.get(self._energytype.value) == "GivenT":
            print(
                Color.PURPLE
                + "** cannot specify heat transfer area of a fixed temperature flow reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "AEXT"
            iErr = self.setprofile(key=keyword, x=x, y=area)
            return iErr

    def setheatlossprofile(self, x, Qloss):
        """
        Specify reactor heat loss rate per length profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param Qloss: heat loss rate value of the profile data [cal/sec-cm] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.EnergyTypes.get(self._energytype.value) == "GivenT":
            print(
                Color.PURPLE
                + "** cannot specify heat loss rate of a fixed temperature flow reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "QPRO"
            iErr = self.setprofile(key=keyword, x=x, y=Qloss)
            return iErr

    def setvelocityprofile(self, x, vel):
        """
        Specify axial velocity profile along the plug-flow reactor
        :param x: position value of the profile data [cm] (double array)
        :param vel: axial velocity value of the profile data [cm/sec] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.ReactorTypes.get(self._reactortype.value) != "PFR":
            print(
                Color.PURPLE
                + "** cannot specify axial velocity of a non plug-flow reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "VELPRO"
            iErr = self.setprofile(key=keyword, x=x, y=vel)
            return iErr


class PlugFlowReactor_GivenTemperature(PlugFlowReactor):
    """
    Plug Flow Reactor (PFR) model with given temperature
    """

    def __init__(self, inlet, label=None):
        # set default label
        if label is None:
            label = "PFR"
        # check Inlet
        if isinstance(inlet, Inlet):
            # initialization
            super().__init__(inlet, label)
        else:
            # wrong argument type
            print(
                Color.RED + "** the first argument must be an Inlet object",
                end=Color.END,
            )
            raise TypeError

        # set reactor type
        self._energytype = c_int(self.EnergyTypes.get("GivenT"))
        # set up basic batch reactor parameters
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
            # setup reactor model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            print(
                Color.RED + f"** error initializing the reactor model: {self.label:s}"
            )
            print(f"   error code = {iErr:d}", end=Color.END)
            exit()

    def settemperatureprofile(self, x, temp):
        """
        Specify reactor temperature profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param temp: temperature value of the profile data [K] (double array)
        :return: error code (integer scalar)
        """
        keyword = "TPRO"
        iErr = self.setprofile(key=keyword, x=x, y=temp)
        return iErr

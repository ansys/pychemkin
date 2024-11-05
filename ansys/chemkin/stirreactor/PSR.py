from ctypes import c_double, c_int
import logging

import numpy as np

from chemkin import chemkin_wrapper
from chemkin.chemistry import checkchemistryset, chemistrysetinitialized, setverbose
from chemkin.color import Color as Color
from chemkin.reactormodel import Keyword
from chemkin.stirreactor.openreactor import openreactor

logger = logging.getLogger(__name__)


class perfectlystirredreactor(openreactor):
    """
    Generic perfectly-stirred reactor
    """

    def __init__(self, guessedmixture, label=None, reactorindex=None):
        """
        Create a steady-state constant pressure perfectly-stirred reactor (PSR)
        :param guessedmixture: guessed/estimate reactor condition (Mixture object)
        :param reactorindex: index of current PSR in the cluster (integer scalar)
        :param label: inlet name (string)
        """
        if label is None:
            label = "PSR"
        if reactorindex is None:
            reactorindex = 1
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
        # heat transfer surface area [cm2]
        self.HTarea = 0.0
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        # check required inputs
        self._numb_requiredinput = 0
        self._requiredlist = []
        self._inputcheck = []
        # default number of reactors
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)  # self.numbexternalinlets
        self._nzones = c_int(0)
        # default reactor type settings
        # Perfectly-Stirred Reactor (PSR) model
        self._reactortype = c_int(2)
        # Steady-State PSR only
        self._solvertype = c_int(self.SolverTypes.get("SteadyState", 2))
        # default options
        self._problemtype = c_int(self.ProblemTypes.get("CONP", 1))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # set reactor number
        self.ireac = c_int(reactorindex)

    @property
    def area(self):
        """
        Get reactive surface area (optional) [cm2] (float scalar)
        """
        return self._reactivearea.value

    @area.setter
    def area(self, value=0.0e0):
        """
        Set reactive surface area (optional)
        default value = 0.0 cm2
        :param value: surface area [cm2] (float scalar)
        """
        if value < 0.0e0:
            print(
                Color.PURPLE + "** reactor reactive area must be >= 0",
                end=Color.END,
            )
        else:
            self._reactivearea = c_double(value)

    def setinletkeywords(self):
        """
        set up inlet keywords
        """
        iErr = 0
        # loop over all external inlets into the reactor
        iInlet = 0
        flowrate_sum = 0.0
        #
        for key, inlet in self.externalinlets.items():
            # get inlet mass flow rate
            flowrate = inlet.massflowrate
            flowrate_sum += flowrate
            # inlet temperature
            T_inlet = inlet.temperature
            # inlet mass fraction
            Y_inlet = inlet.Y
            #
            if np.isclose(0.0, flowrate, atol=1.0e-6):
                print(
                    Color.PURPLE + f"** inlet {key} has zero flow rate", end=Color.END
                )
                iErrc = 100 + iInlet + 1
            else:
                iInlet += 1
                # set inlet inputs
                iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupPSRInletInputs(
                    self._chemset_index,
                    self.ireac,
                    c_double(iInlet),
                    c_double(flowrate),
                    T_inlet,
                    Y_inlet,
                )
            iErr += iErrc
        # check number of external inlet
        if iInlet == 0:
            print(Color.PURPLE + "** no external inlet to the PSR", end=Color.END)
            iErr += 10
        elif iInlet != self.numbexternalinlets:
            print(Color.PURPLE + "** inconsistent number of external inlets")
            print(f"   expected number of inlets: {self.numbexternalinlets}")
            print(f"   actual number of inlets: {iInlet}", end=Color.END)
            iErr += 11
        else:
            pass
        # check total mass flow rate
        if iErr == 0:
            # check total mass flow rate
            if np.isclose(flowrate, self.totalmassflowrate, atol=1.0e-6):
                print(Color.PURPLE + "** inconsistent inlet mass flow rate value")
                print(f"   expected total mass flow rate: {self.totalmassflowrate}")
                print(f"   actual total mass flow rate: {flowrate_sum}", end=Color.END)
                iErr += 12
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
        # set up inlets
        iErrInputs = self.setinletkeywords()
        iErr += iErrInputs
        # prepare estimated reactor conditions
        # estimated reactor mass fraction
        Y_init = self.reactormixture.Y
        # surface sites (not applicable)
        Site_init = np.zeros(1, dtype=np.double)
        # bulk activities (not applicable)
        Bulk_init = np.zeros_like(Site_init, dtype=np.double)
        # set estimated reactor conditions and geometry parameters
        if self._reactortype.value == 2:
            iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupPSRReactorInputs(
                self._chemset_index,
                self.ireac,
                self._endtime,
                self._temperature,
                self._pressure,
                self._volume,
                self._heatlossrate,
                self._residencetime,
                self._reactivearea,
                Y_init,
                Site_init,
                Bulk_init,
            )
            iErr += iErrc
            # heat transfer (use additional keywords)
            # solver parameters (use additional keywords)
            # output controls (use additional keywords)
            # ROP (use additional keywords)
            # sensitivity (use additional keywords)
            # ignition delay (use additional keywords)
            # solve integrated heat release rate due to chemical reactions
            if self.EnergyTypes.get("ENERGY") == self._energytype.value:
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
        # initialize the PSR model
        # set up basic PSR parameters
        # number of external inlets
        self._ninlets = c_double(self.numbexternalinlets)
        #
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
        #
        # get ready to run the reactor model
        # initialize KINetics
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
            print(
                Color.RED + "** full keyword option not available for PSR models",
                end=Color.END,
            )
            retVal = 100
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
        # update run status
        self.setrunstatus(code=retVal)

        logger.debug("Running model complete, status = " + str(retVal))

        return retVal


class PSR_SetResTime_EnergyConservation(perfectlystirredreactor):
    """
    PSR model with given reactor reasidence time (CONP)
    and solve energy equation (ENERGY)
    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor volume and density are varying in this case.
    """

    def __init__(self, guessedmixture, label=None, reactorindex=None):
        """
        Create a steady-state constant pressure perfectly-stirred reactor (PSR)
        :param guessedmixture: guessed/estimate reactor condition (Mixture object)
        :param reactorindex: index of current PSR in the cluster (integer scalar)
        :param label: inlet name (string)
        """
        if label is None:
            label = "PSR"
        if reactorindex is None:
            reactorindex = 1
        # initialization
        super().__init__(guessedmixture, label, reactorindex)
        # specify residence
        self._problemtype = c_int(self.ProblemTypes.get("CONP", 1))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        self._heattransfercoefficient = 0.0e0
        self._ambienttemperature = 3.0e2
        # external heat transfer area per reactor length [cm2/cm]
        self._heattransferarea = 0.0e0

    @property
    def residencetime(self):
        """
        Get reactor residence time (required) [sec] (float scalar)
        """
        return self._residencetime.value

    @residencetime.setter
    def residencetime(self, value):
        """
        Set reactor residence time (required)
        default value = 0.0 sec
        :param value: reactor residence time [sec] (float scalar)
        """
        if value > 0.0e0:
            # set reactor residence time
            self._residencetime = c_double(value)
        else:
            print(Color.PURPLE + "** residence time must be > 0", end=Color.END)

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
        Set the heat loss rate from the reactor to the surroundings (required)
        default value = 0.0 cal/sec
        :param value: heat loss rate [cal/sec] (float scalar)
        """
        self._heatlossrate = c_double(value)

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
        Get heat transfer area between the reactor and the surroundings (optional) [cm2] (float scalar)
        default value = 0.0 cm2
        """
        return self._heattransferarea

    @heattransferarea.setter
    def heattransferarea(self, value=0.0e0):
        """
        Set heat transfer area between the reactor and the surroundings (optional)
        default value = 0.0 cm2
        :param value: heat transfer area [cm2] (float scalar)
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


class PSR_SetVolume_EnergyConservation(perfectlystirredreactor):
    """
    PSR model with given reactor volume (CONV)
    and solve energy equation (ENERGY)
    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor residence time and density are varying in this case.
    """

    def __init__(self, guessedmixture, label=None, reactorindex=None):
        """
        Create a steady-state constant pressure perfectly-stirred reactor (PSR)
        :param guessedmixture: guessed/estimate reactor condition (Mixture object)
        :param reactorindex: index of current PSR in the cluster (integer scalar)
        :param label: inlet name (string)
        """
        if label is None:
            label = "PSR"
        if reactorindex is None:
            reactorindex = 1
        # initialization
        super().__init__(guessedmixture, label, reactorindex)
        # specify volume
        self._problemtype = c_int(self.ProblemTypes.get("CONV", 2))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY", 1))
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        self._heattransfercoefficient = 0.0e0
        self._ambienttemperature = 3.0e2
        # external heat transfer area per reactor length [cm2/cm]
        self._heattransferarea = 0.0e0

    @property
    def volume(self):
        """
        Get reactor volume (required) [cm3] (float scalar)
        """
        return self._volume.value

    @volume.setter
    def volume(self, value):
        """
        Set reactor volume (required)
        default value = 0.0 cm3
        :param value: reactor volume [cm3] (float scalar)
        """
        if value > 0.0e0:
            # set reactor volume
            self._volume = c_double(value)
            # set initial mixture volume
            self.reactormixture.volume = value
        else:
            print(Color.PURPLE + "** reactor volume must be > 0", end=Color.END)

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
        Set the heat loss rate from the reactor to the surroundings (required)
        default value = 0.0 cal/sec
        :param value: heat loss rate [cal/sec] (float scalar)
        """
        self._heatlossrate = c_double(value)

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
        Get heat transfer area between the reactor and the surroundings (optional) [cm2] (float scalar)
        default value = 0.0 cm2
        """
        return self._heattransferarea

    @heattransferarea.setter
    def heattransferarea(self, value=0.0e0):
        """
        Set heat transfer area between the reactor and the surroundings (optional)
        default value = 0.0 cm2
        :param value: heat transfer area [cm2] (float scalar)
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


class PSR_SetResTime_FixedTemperature(perfectlystirredreactor):
    """
    PSR model with given reactor reasidence time (CONP)
    and reactor temperature (GivenT)
    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor volume and density are varying in this case.
    """

    def __init__(self, guessedmixture, label=None, reactorindex=None):
        """
        Create a steady-state constant pressure perfectly-stirred reactor (PSR)
        :param guessedmixture: guessed/estimate reactor condition (Mixture object)
        :param reactorindex: index of current PSR in the cluster (integer scalar)
        :param label: inlet name (string)
        """
        if label is None:
            label = "PSR"
        if reactorindex is None:
            reactorindex = 1
        # initialization
        super().__init__(guessedmixture, label, reactorindex)
        # specify residence time
        self._problemtype = c_int(self.ProblemTypes.get("CONP", 1))
        self._energytype = c_int(self.EnergyTypes.get("GivenT", 2))

    @property
    def residencetime(self):
        """
        Get reactor residence time (required) [sec] (float scalar)
        """
        return self._residencetime.value

    @residencetime.setter
    def residencetime(self, value):
        """
        Set reactor residence time (required)
        default value = 0.0 sec
        :param value: reactor residence time [sec] (float scalar)
        """
        if value > 0.0e0:
            # set reactor residence time
            self._residencetime = c_double(value)
        else:
            print(Color.PURPLE + "** residence time must be > 0", end=Color.END)


class PSR_SetVolume_FixedTemperature(perfectlystirredreactor):
    """
    PSR model with given reactor volume (CONV)
    and reactor temperature (GivenT)
    rho_PSR * Vol_PSR / residence_time = mass_flow_rate
    The reactor pressure and the inlet mass flow rate are always given (fixed)
    so the reactor residence time and density are varying in this case.
    """

    def __init__(self, guessedmixture, label=None, reactorindex=None):
        """
        Create a steady-state constant pressure perfectly-stirred reactor (PSR)
        :param guessedmixture: guessed/estimate reactor condition (Mixture object)
        :param reactorindex: index of current PSR in the cluster (integer scalar)
        :param label: inlet name (string)
        """
        if label is None:
            label = "PSR"
        if reactorindex is None:
            reactorindex = 1
        # initialization
        super().__init__(guessedmixture, label, reactorindex)
        # specify volume
        self._problemtype = c_int(self.ProblemTypes.get("CONV", 2))
        self._energytype = c_int(self.EnergyTypes.get("GivenT", 2))

    @property
    def volume(self):
        """
        Get reactor volume (required) [cm3] (float scalar)
        """
        return self._volume.value

    @volume.setter
    def volume(self, value):
        """
        Set reactor volume (required)
        default value = 0.0 cm3
        :param value: reactor volume [cm3] (float scalar)
        """
        if value > 0.0e0:
            # set reactor volume
            self._volume = c_double(value)
            # set initial mixture volume
            self.reactormixture.volume = value
        else:
            print(Color.PURPLE + "** reactor volume must be > 0", end=Color.END)

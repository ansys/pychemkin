import copy
from ctypes import c_double, c_int
import logging

from chemkin import chemkin_wrapper
from chemkin.chemistry import checkchemistryset, chemistrysetinitialized, setverbose
from chemkin.color import Color as Color
from chemkin.engines.engine import Engine
from chemkin.reactormodel import Keyword

logger = logging.getLogger(__name__)


class SIengine(Engine):
    """
    Spark Ignition (SI) engine model
    """

    def __init__(self, reactor_condition, label=None):
        # set default number of zone(s)
        # 2 zones: the unburned and the burned zones
        nzones = 2
        # set default label
        if label is None:
            label = "SI"

        # use the first zone to initialize the engine model
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("SI"))
        self._solvertype = c_int(self.SolverTypes.get("Transient"))
        self._problemtype = c_int(self.ProblemTypes.get("ICEN"))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY"))
        # defaults for all closed homogeneous reactor models
        # 2 zones: the unburned and the burned zones
        self._nreactors = nzones
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        # number of zones
        self._nzones = c_int(nzones)
        # use API mode for SI simulations
        Keyword.noFullKeyword = True
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(156)
        # burn mass profile mode
        # 0: unset
        # 1: Wiebe function with n, b, SOI, burn duration
        # 2: anchor points 10%, 50%, and 90% mass burned CAs
        # 3: burned mass fraction profile, SOI, burn duration
        self._burnmode = 0
        self.sparktiming = -180
        self.burnduration = 0.0
        self.wieben = 2.0
        self.Wiebeb = 5.0
        self.massburnedCA10 = -180.0
        self.massburnedCA50 = -180.0
        self.massburnedCA90 = -180.0
        # combustion efficiency
        self.burnefficiency = 1.0
        # numbwer of points in the mass burned fraction profile
        self.MBpoints = 0
        self.MBangles = None
        self.MBfractions = None
        # set up basic SI engine model parameters
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
            # setup SI engine model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            print(
                Color.RED + f"** error initializing the SI engine model: {self.label:s}"
            )
            print(f"   error code = {iErr:d}", end=Color.END)
            exit()

    def wiebeparameters(self, n, b):
        """
        Set Wiebe function parameters
        Wiebe = 1 - exp{-b[(CA-SOC)/duration]^(n+1)]}
        :param n: exponent (double scalar)
        :param b: multiplier
        """
        # check input values
        if n <= 0.0 or b <= 0.0:
            print(
                Color.PURPLE + "** Wiebe function parameters n and b must > 0",
                end=Color.END,
            )
            exit()
        #
        if self._burnmode > 0:
            print(
                Color.YELLOW
                + "** previous burned mass profile setup will be overwritten",
                end=Color.END,
            )
        #
        self._burnmode = 1
        self.wieben = n
        self.wiebeb = b

    def setburntiming(self, SOC, duration):
        """
        Set SI engine burn timing
        :param SOC: start of combustion in crank angle [degree] (double scalar)
        :param duration: burn duration in crank angles [degree] (double scalar)
        """
        if SOC <= self.IVCCA:
            print(
                Color.PURPLE
                + f"** start of combustion CA must > start of simulation CA {self.IVCCA}",
                end=Color.END,
            )
            exit()
        if duration <= 0.0:
            print(Color.PURPLE + "** ass burn duration must > 0", end=Color.END)
            exit()
        #
        self.sparktiming = SOC
        self.burnduration = duration

    def setanchorpoints(self, CA10, CA50, CA90):
        """
        Set the SI mass burned profile using the anchor points
        :param CA10: crank angle of 10% mass burned [degree] (double scalar)
        :param CA50: crank angle of 50% mass burned [degree] (double scalar)
        :param CA90: crank angle of 90% mass burned [degree] (double scalar)
        """
        if CA10 > CA50:
            print(
                Color.PURPLE + "** the anchor points must in ascending order",
                end=Color.END,
            )
            exit()
        if CA90 < CA50:
            print(
                Color.PURPLE + "** the anchor points must in ascending order",
                end=Color.END,
            )
            exit()
        if CA10 <= self.IVCCA:
            print(
                Color.PURPLE
                + f"** anchor point CA must > start of simulation CA {self.IVCCA}",
                end=Color.END,
            )
            exit()
        #
        if self._burnmode > 0:
            print(
                Color.YELLOW
                + "** previous burned mass profile setup will be overwritten",
                end=Color.END,
            )
        #
        self._burnmode = 2
        self.massburnedCA10 = CA10
        self.massburnedCA50 = CA50
        self.massburnedCA90 = CA90

    def setmassburnedprofile(self, crankangles, fractions):
        """
        Specify SI engine mass burned fraction profile
        :param crankangles: normalized crank angles of the profile data [degree]
        the crank angles must 0 <= and <= 1 (double array)
        :param fractions: mass burned fraction of the profile data [-] (double array)
        :return: error code (integer scalar)
        """
        # set the mass burned profile
        iError = 0
        self.MBpoints = len(crankangles)
        if len(fractions) != self.MBpoints:
            print(
                Color.PURPLE + "** data arrays must have the same size", end=Color.END
            )
            iError = 1
        elif self.MBpoints > 1:
            self.MBangles = copy.deepcopy(crankangles)
            self.MBfractions = copy.deepcopy(fractions)
            self._burnmode = 3
        else:
            print(
                Color.PURPLE + "** profile must have more than 1 data pair",
                end=Color.END,
            )
            iError = 2
        return iError

    def setcombustionefficiency(self, efficiency):
        """
        Set the overall combustion efficiency
        :param efficiency: combustion efficiency (double scalar)
        """
        # check value
        if efficiency < 0.0 or efficiency > 1.0:
            print(Color.PURPLE + "** efficiency must 0.0 <= and <= 1.0", end=Color.END)
            exit()
        # set keyword
        self.burnefficiency = efficiency
        self.setkeyword(key="BEFF", value=efficiency)

    def setburnedproductminfraction(self, bound):
        """
        Set the minimum gas species mole fraction value from the flame sheet
        to be injected to the burned zone
        :param bound: minimum species mole fraction value (double scalar)
        """
        if bound > 0.0:
            # set keyword
            self.setkeyword(key="EQMN", value=bound)
        else:
            print(Color.PURPLE + "** species fraction value must > 0.0", end=Color.END)
            exit()

    def setwiebekeywords(self):
        """
        Set the Wiebe function parameters keywords for the SI engine model
        :return: error code (integer scalar)
        """
        iError = 0
        if self._burnmode == 1:
            # set start of combustion time
            self.setkeyword(key="BINI", value=self.sparktiming)
            # set burn duration
            self.setkeyword(key="BDUR", value=self.burnduration)
            # set Wiebe parameter b
            self.setkeyword(key="WBFB", value=self.wiebeb)
            # set Wiebe parameter n
            self.setkeyword(key="WBFN", value=self.wieben)
        else:
            print(
                Color.PURPLE + "** incorrect burned mass profile setup option",
                end=Color.END,
            )
            iError = 10
        return iError

    def setanchorpointskeywords(self):
        """
        Set the mass burned porfile anchor points keywords for the SI engine model
        :return: error code (integer scalar)
        """
        iError = 0
        if self._burnmode == 2:
            # set 10% mass burned crank angle
            self.setkeyword(key="CASC", value=self.massburnedCA10)
            # set 50% mass burned crank angle
            self.setkeyword(key="CAAC", value=self.massburnedCA50)
            # set 90% mass burned crank angle
            self.setkeyword(key="CAEC", value=self.massburnedCA90)
        else:
            print(
                Color.PURPLE + "** incorrect burned mass profile setup option",
                end=Color.END,
            )
            iError = 11
        return iError

    def setburnprofilekeywords(self):
        """
        Set the mass burned fraction profile keywords for the SI engine model
        :return: error code (integer scalar)
        """
        iError = 0
        if self._burnmode == 3 and self.MBpoints > 0:
            # set start of combustion time
            self.setkeyword(key="BINI", value=self.sparktiming)
            # set burn duration
            self.setkeyword(key="BDUR", value=self.burnduration)
            # set number of burned mass fraction profile data points
            self.setkeyword(key="NBFP", value=self.MBpoints)
            # set mass burned fraction profile keywords
            for i in range(self.MBpoints):
                # set mass burned fraction profile
                keyline = (
                    "BFP"
                    + Keyword.fourspaces
                    + str(self.MBangles[i])
                    + Keyword.fourspaces
                    + str(self.MBfractions[i])
                )
                self.setkeyword(key=keyline, value=True)
        else:
            print(
                Color.PURPLE + "** incorrect burned mass profile setup option",
                end=Color.END,
            )
            iError = 12

        return iError

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
        if self._reactortype.value == self.ReactorTypes.get("SI"):
            # insert the ICEN keywords
            self.setkeyword(key="ICEN", value=True)
            #
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
            # set SI engine parameter
            if self._burnmode == 1:
                # use Wiebe function to specify the mass burned profile
                iErrc = self.setwiebekeywords()
                iErr += iErrc
            elif self._burnmode == 2:
                # use anchor points to specify the mass burned profile
                iErrc = self.setanchorpointskeywords()
                iErr += iErrc
            elif self._burnmode == 3:
                # use normalized profile to specify the mass burned profile
                iErrc = self.setburnprofilekeywords()
                iErr += iErrc
            else:
                print(
                    Color.RED
                    + "** burned mass rate is not set up for the SI engine simulation"
                )
                print("  Chemkin SI engine model provides three methods:")
                print("  1. Wiebe function")
                print("  2. anchor points")
                print("  3. piece-wise linear normalized CA-burned fraction profile")
                print("  please see Chemkin Theory manual for details", end=Color.END)
                exit()

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
        if Keyword.noFullKeyword:
            # use API calls
            retVal = (
                self.__process_keywords()
            )  # each reactor model subclass to perform its own keyword processing
        else:
            # use full keywords
            retVal = -1
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

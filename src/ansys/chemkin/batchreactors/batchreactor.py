import copy
from ctypes import c_double, c_int
import logging

import numpy as np

from chemkin import chemkin_wrapper
from chemkin.chemistry import (
    Patm,
    checkchemistryset,
    chemistrysetinitialized,
    findinterpolateparameters,
    setverbose,
    showignitiondefinition,
)
from chemkin.color import Color as Color
from chemkin.mixture import interpolatemixtures
from chemkin.reactormodel import Keyword
from chemkin.reactormodel import ReactorModel as reactor

logger = logging.getLogger(__name__)


class BatchReactors(reactor):
    """
    Generic model of Chemkin 0-D transient closed homogeneous reactor models
    """

    # set possible types in batch reactors
    ReactorTypes = {"Batch": 1, "PSR": 2, "PFR": 3, "HCCI": 4, "SI": 5, "DI": 6}
    SolverTypes = {"Transient": 1, "SteadyState": 2}
    EnergyTypes = {"ENERGY": 1, "GivenT": 2}
    ProblemTypes = {"CONP": 1, "CONV": 2, "ICEN": 3}

    def __init__(self, reactor_condition, label):
        # initialize the base module
        super().__init__(reactor_condition, label)
        #
        # reactor parameters (required)
        self._volume = c_double(0.0e0)
        self._endtime = c_double(0.0e0)
        self._reactivearea = c_double(0.0e0)
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        # solver parameters
        self._absolutetolerance = 1.0e-12
        self._relativetolerance = 1.0e-6
        # check required inputs
        self._numb_requiredinput = 0
        self._requiredlist = []
        self._inputcheck = []
        # default number of reactors
        self._nreactors = 1

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
            # set volume keyword (not set by the setup calls)
            self.setkeyword(key="VOL", value=value)
        else:
            print(Color.PURPLE + "** reactor volume must be > 0", end=Color.END)

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
            # set internal surface area keyword for PFR
            if (
                self._reactortype.value == self.ReactorTypes.get("PFR")
                or Keyword.noFullKeyword
            ):
                self.setkeyword(key="AREA", value=value)

    @property
    def tolerances(self):
        """
        Get solver tolerances, absolute tolerance, relative tolerance (float scalar, float scalar)
        """
        return self._absolutetolerance, self._relativetolerance

    def settolerances(self, absolute_tolerance=None, relative_tolerance=None):
        """
        Set solver tolerances
        :param absolute_tolerance: absolute tolerance (float scalar)
        :param relative_tolerance: relative tolerance (float scalar)
        """
        # set absolute tolerance
        if absolute_tolerance is not None:
            self._absolutetolerance = max(absolute_tolerance, 1.0e-20)
            # set keywords
            self.setkeyword(key="ATOL", value=self._absolutetolerance)
        # set relative tolerance
        if relative_tolerance is not None:
            self._relativetolerance = max(relative_tolerance, 1.0e-20)
            # set keywords
            self.setkeyword(key="RTOL", value=self._relativetolerance)

    @property
    def forcenonnegative(self):
        """
        Get the status of the forcing non-negative option of the transient solver (boolean scalar)
        """
        if "NNEG" in self._keyword_index:
            # defined: find index
            i = self._keyword_index.index("NNEG")
            return self._keyword_list[i].value
        else:
            # not defined: return default value
            return False

    @forcenonnegative.setter
    def forcenonnegative(self, mode=False):
        """
        Set the forcing non-negative solution option
        :param mode: True/False (turn the option ON/OFF) (boolean scalar)
        """
        # set keyword
        self.setkeyword(key="NNEG", value=mode)

    def setsolverinitialtimestepsize(self, size):
        """
        Set the initial time step size to be used by the solver
        :param size: step size [sec] or [cm]
        """
        if size > 0.0e0:
            self.setkeyword(key="HO", value=size)
        else:
            print(
                Color.PURPLE + "** solver timestep size must > 0",
                end=Color.END,
            )

    def setsolvermaxtimestepsize(self, size):
        """
        Set the maximum time step size allowed by the solver
        :param size: step size [sec] or [cm]
        """
        if size > 0.0e0:
            self.setkeyword(key="STPT", value=size)
        else:
            print(
                Color.PURPLE + "** solver timestep size must > 0",
                end=Color.END,
            )

    @property
    def timestepforsavingsolution(self):
        """
        Get the timestep size between saving the solution data [sec] (float scalar)
        """
        if "DTSV" in self._keyword_index:
            # defined: find index
            i = self._keyword_index.index("DTSV")
            return self._keyword_list[i].value
        else:
            # return default value (100th of the end time)
            if self._endtime.value > 0.0e0:
                return self._endtime.value / 1.0e2
            else:
                # not defined yet
                print(
                    Color.YELLOW
                    + "** solution saving timestep is not defined because 'end time' is not set",
                    end=Color.END,
                )
                return 0.0

    @timestepforsavingsolution.setter
    def timestepforsavingsolution(self, delta_time):
        """
        Set the timestep size between saving the solution data
        :param delta_time: timestep size between saving solution data [sec] (float scalar)
        """
        if delta_time > 0.0e0:
            self.setkeyword(key="DTSV", value=delta_time)
        else:
            print(
                Color.PURPLE + "** solution saving timestep value must > 0",
                end=Color.END,
            )

    @property
    def timestepforprintingsolution(self):
        """
        Get the timestep size between printing the solution data to the text output file [sec] (float scalar)
        """
        if "DELT" in self._keyword_index:
            # defined: find index
            i = self._keyword_index.index("DELT")
            return self._keyword_list[i].value
        else:
            # return default value (100th of the end time)
            if self._endtime.value > 0.0e0:
                return self._endtime.value / 1.0e2
            else:
                # not defined yet
                print(
                    Color.YELLOW
                    + "** solution printing timestep is not defined because 'end time' is not set",
                    end=Color.END,
                )
                return 0.0

    @timestepforprintingsolution.setter
    def timestepforprintingsolution(self, delta_time):
        """
        Set the timestep size between printing the solution data to the text output file
        :param delta_time: timestep size between printing solution data [sec] (float scalar)
        """
        if delta_time > 0.0e0:
            self.setkeyword(key="DELT", value=delta_time)
        else:
            print(
                Color.PURPLE + "** solution printing timestep value must > 0",
                end=Color.END,
            )

    def adaptivesolutionsaving(self, mode, value_change=None, target=None, steps=None):
        """
        Set up adaptive solution data saving
        :param mode: True/False (switch adaptive solution saving ON/OFF) (boolean scalr)
        :param value_change: change in solution variable value between saving additional solution data (float scalar)
        :param target: the target variable that is used by the value_change option (string scalar)
        :param steps: number of solver time steps between saving additional solution data (integer scalar)
        :return: None
        """
        # turn ON/OFF the adaptive solution saving option
        self.setkeyword(key="ADAP", value=mode)
        self.setkeyword(key="NADAP", value=not mode)
        # set options
        if steps is not None:
            # use number of solver time steps option:
            if steps <= 0.0:
                # non-positive number of steps
                print(
                    Color.PURPLE
                    + "** the number of steps per adaptive solution saving must be > 0 ",
                    end=Color.END,
                )
            else:
                # set parameters
                self.setkeyword(key="ADAP", value=True)
                self.setkeyword(key="ASTEPS", value=int(steps))
        elif value_change is not None:
            # use change in the solution variable value option:
            if target is None:
                # target variable is not given
                print(
                    Color.PURPLE
                    + "** a reference variable is required for value-change adaptive saving",
                    end=Color.END,
                )
            elif not isinstance(target, str):
                # not given as string
                print(
                    Color.PURPLE
                    + "** a reference variable is assigned as a string, e.g., 'OH' ",
                    end=Color.END,
                )
            elif value_change <= 0.0:
                # non-positive change value
                print(
                    Color.PURPLE
                    + "** the value change per adaptive solution saving must be > 0 ",
                    end=Color.END,
                )
            else:
                # set parameters
                self.setkeyword(key="ADAP", value=True)
                self.setkeyword(key="AVAR", value=target)
                self.setkeyword(key="AVALUE", value=value_change)
        elif mode:
            # use error
            print(
                Color.PURPLE
                + "** need to specify either the number of steps or the <change value>-<target variable> pair",
                end=Color.END,
            )

    def setignitiondelay(self, method=None, val=None, target=None):
        """
        Set ignition detection criterion
        :param method: ignition definition/detection method (string scalar)
        :param val: temperature or temperature rise value associated with the ignition detection method specified (float scalar)
        :param target: target species symbol if the 'Species_peak' method is used (string scalar)
        :return: None
        """
        if method is None:
            # method is not given: use the default detection method: inflection points in the temperature profile
            method = "T_inflection"
        if isinstance(method, str):
            # ignition detection method assigned
            if method == "T_inflection":
                # use inflection points in the temperature profile
                self.setkeyword(key="TIFP", value=True)
            elif method == "T_rise":
                # use temperature rise
                if val <= 0.0:
                    print(
                        Color.PURPLE + "** temperature rise must be > 0 ",
                        end=Color.END,
                    )
                else:
                    self.setkeyword(key="DTIGN", value=val)
            elif method == "T_ignition":
                # use temperature value
                if val <= 0.0:
                    print(
                        Color.PURPLE + "** ignition temperature must be > 0 ",
                        end=Color.END,
                    )
                else:
                    self.setkeyword(key="TLIM", value=val)
            elif method == "Species_peak":
                # use species peak location
                if not isinstance(target, str):
                    # no species given
                    print(
                        Color.PURPLE
                        + "** target species is assigned as a string, e.g., 'OH' ",
                        end=Color.END,
                    )
                else:
                    self.setkeyword(key="KLIM", value=target)
            else:
                # incorrect ignition detection method given
                print(
                    Color.PURPLE
                    + f"** ignition definition specified {method:s} is not recognized",
                    end=Color.END,
                )
                showignitiondefinition()
        else:
            # use error
            showignitiondefinition()

    def stopafterignition(self):
        """
        Set the option to stop the simulation after ignition is detected
        :return: None
        """
        # stop the simulation after ignition is detected
        self.setkeyword(key="IGN_STOP", value=True)

    def getignitiondelay(self):
        """
        Get the predicted ignition delay time from the transient reactor simulation
        :return: ignition delay time [msec] or [CA] (double scalar)
        """
        # initialization
        ignitiondelaytime = c_double(0.0e0)
        # check run status
        status = self.getrunstatus(mode="silent")
        if status == -100:
            print(
                Color.YELLOW + "** please run the reactor simulation first",
                end=Color.END,
            )
            return ignitiondelaytime.value
        elif status != 0:
            print(Color.YELLOW + "** simulation was failed")
            print(
                "** please correct the error and rerun the reactor simulation",
                end=Color.END,
            )
            return ignitiondelaytime.value

        # get the ignition delay time (batch reactor model [sec], engine model [CA])
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetIgnitionDelay(ignitiondelaytime)
        if iErr != 0:
            print(Color.PURPLE + "** potential bad ignition delay time value")
            print(
                "** please check the run status and revisit the reactor settings",
                end=Color.END,
            )
        # check reactor model
        if self._reactortype.value == self.ReactorTypes.get("Batch"):
            # check ignition delay time value
            if ignitiondelaytime.value <= 0.0:
                print(Color.PURPLE + "** potential bad ignition delay time value")
                print(
                    "** please check the run status and revisit the reactor settings",
                    end=Color.END,
                )
                return ignitiondelaytime.value
            else:
                # convert ignition  delay time from [sec] to [msec]
                return ignitiondelaytime.value * 1.0e3
        elif self._reactortype.value in [
            self.ReactorTypes.get("HCCI"),
            self.ReactorTypes.get("SI"),
            self.ReactorTypes.get("DI"),
        ]:
            # engine models
            print(Color.YELLOW + "** mean ignition delay time in CA", end=Color.END)
            return ignitiondelaytime.value
        elif self._reactortype.value == self.ReactorTypes.get("PFR"):
            print(Color.YELLOW + "** ignition delay in [cm]", end=Color.END)
            # check ignition distance value
            if ignitiondelaytime.value <= 0.0:
                print(Color.PURPLE + "** potential bad ignition distance value")
                print(
                    "** please check the run status and revisit the reactor settings",
                    end=Color.END,
                )
                return ignitiondelaytime.value
            else:
                # return ignition distance [cm]
                return ignitiondelaytime.value

    def setvolumeprofile(self, x, vol):
        """
        Specify reactor volume profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param vol: volume value of the profile data [cm3] (double array)
        :return: error code (integer scalar)
        """
        if (
            BatchReactors.ProblemTypes.get(self._problemtype.value) == "CONP"
            and BatchReactors.EnergyTypes.get(self._energytype.value) == "GivenT"
        ):
            print(
                Color.PURPLE
                + "** cannot constrain volume of a given pressure fixed temperature batch reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "VPRO"
            iErr = self.setprofile(key=keyword, x=x, y=vol)
            return iErr

    def setpressureprofile(self, x, pres):
        """
        Specify reactor pressure profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param pres: pressure value of the profile data [dynes/cm2] (double array)
        :return: error code (integer scalar)
        """
        if (
            BatchReactors.ProblemTypes.get(self._problemtype.value) == "CONV"
            and BatchReactors.EnergyTypes.get(self._energytype.value) == "GivenT"
        ):
            print(
                Color.PURPLE
                + "** cannot constrain pressure of a given volume fixed temperature batch reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "PPRO"
            iErr = self.setprofile(key=keyword, x=x, y=pres)
            return iErr

    def setsurfaceareaprofile(self, x, area):
        """
        Specify reactor reactive surface area profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param area: reactive surface area value of the profile data [cm2] (double array)
        :return: error code (integer scalar)
        """
        keyword = "AINT"
        iErr = self.setprofile(key=keyword, x=x, y=area)
        return iErr

    def setreactortypekeywords(self):
        """
        Set reactor type keywords under the Full-Keywords mode
        :return: None
        """
        # keyword headers
        # set solver types
        if self._solvertype.value == self.SolverTypes.get("Transient"):
            self.setkeyword(key="TRAN", value=True)
        else:
            self.setkeyword(key="STST", value=True)
        # set reactor related keywords
        if self._reactortype.value == self.ReactorTypes.get("Batch"):
            # batch reactors
            # set problem type
            if self._problemtype.value == self.ProblemTypes.get("CONP"):
                self.setkeyword(key="CONP", value=True)
            else:
                self.setkeyword(key="CONV", value=True)
            # set energy equation
            if self._energytype.value == self.EnergyTypes.get("ENERGY"):
                self.setkeyword(key="ENRG", value=True)
            else:
                self.setkeyword(key="TGIV", value=True)
        elif self._reactortype.value == self.ReactorTypes.get("PSR"):
            # PSR
            # set energy equation
            if self._energytype.value == self.EnergyTypes.get("ENERGY"):
                self.setkeyword(key="ENRG", value=True)
            else:
                self.setkeyword(key="TGIV", value=True)
        elif self._reactortype.value == self.ReactorTypes.get("PFR"):
            # PFR
            self.setkeyword(key="PLUG", value=True)
            # set energy equation
            if self._energytype.value == self.EnergyTypes.get("ENERGY"):
                self.setkeyword(key="ENRG", value=True)
            else:
                self.setkeyword(key="TGIV", value=True)
        else:
            # IC engine models
            self.setkeyword(key="ICEN", value=True)
            self.setkeyword(key="TRAN", value=True)
            # set energy equation
            self.setkeyword(key="ENRG", value=True)

    def setreactorconditionkeywords(self):
        """
        Set reactor initial/estimated condition keywords under the Full-Keywords mode
        :return: None
        """
        self.setkeyword(key="PRES", value=self._pressure.value / Patm)
        self.setkeyword(key="TEMP", value=self._temperature.value)
        self.setkeyword(key="TIME", value=self._endtime.value)
        # initial mole fraction
        nspecieslines, species_lines = self.createspeciesinputlines(
            self._solvertype.value, threshold=1.0e-12, molefrac=self.reactormixture.X
        )
        for line in species_lines:
            self.setkeyword(key=line, value=True)

    def inputvalidation(self):
        iErr = 0
        # required inputs:
        if self._numb_requiredinput <= 0:
            # no required input
            return iErr
        else:
            if len(self._inputcheck) < self._numb_requiredinput:
                print(
                    Color.PURPLE + "** some required inputs are missing",
                    end=Color.END,
                )
            # verify required inputs one by one
            for k in self._requiredlist:
                if k not in self._inputcheck:
                    iErr += 1
                    print(
                        Color.RED + f"** missing required input '{k:s}'",
                        end=Color.END,
                    )
            return iErr

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
        # surface sites (not applicable)
        Site_init = np.zeros(1, dtype=np.double)
        # bulk activities (not applicable)
        Bulk_init = np.zeros_like(Site_init, dtype=np.double)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("Batch"):
            iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupBatchInputs(
                self._chemset_index,
                self._endtime,
                self._temperature,
                self._pressure,
                self._volume,
                self._heatlossrate,
                self._reactivearea,
                Y_init,
                Site_init,
                Bulk_init,
            )
            iErr += iErrc

        # set reactor type
        self.setreactortypekeywords()
        # reactor initial/estimated condition
        self.setreactorconditionkeywords()
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
        # surface sites (not applicable)
        Site_init = np.zeros(1, dtype=np.double)
        # bulk activities (not applicable)
        Bulk_init = np.zeros_like(Site_init, dtype=np.double)
        # set reactor initial conditions and geometry parameters
        if self._reactortype.value == self.ReactorTypes.get("Batch"):
            iErrc = chemkin_wrapper.chemkin.KINAll0D_SetupBatchInputs(
                self._chemset_index,
                self._endtime,
                self._temperature,
                self._pressure,
                self._volume,
                self._heatlossrate,
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

    def getsolutionsize(self):
        """
        Get the number of reactors and the number of solution points
        :return: nreactor = number of reactors, npoints = number of solution points (integer scalar, integer scalar)
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
        # number of reactor
        nreac = c_int(0)
        # number of time points in the solution
        npoints = c_int(0)
        # get solution size of the batch reactor
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetSolnResponseSize(nreac, npoints)
        nreactors = nreac.value
        if iErr == 0 and nreactors == self._nreactors:
            # return the solution sizes
            self._numbsolutionpoints = (
                npoints.value
            )  # number of time points in the solution profile
            return self._nreactors, self._numbsolutionpoints
        elif self._nreactors == nreactors:
            # fail to get solution sizes
            print(
                Color.PURPLE + f"** failed to get solution size, error code = {iErr}",
                end=Color.END,
            )
            return nreactors, 0
        else:
            # incorrect number of reactor (batch reactor is single reactor)
            print(Color.PURPLE + f"** incorrect number of reactor = {nreactors}")
            print(f"   the model expects {self._nreactors} reactor/zone", end=Color.END)
            return nreactors, 0

    def processsolution(self):
        """
        Post-process solution to extract the raw solution variable data
        """
        # check existing raw data
        if self.getrawsolutionstatus():
            print(
                Color.YELLOW
                + "** solution has been processed before, existing solution data will be deleted",
                end=Color.END,
            )

        # reset raw and mixture solution parameters
        self._numbsolutionpoints = 0
        self._solution_rawarray.clear()
        self._solution_mixturearray.clear()
        # get solution sizes
        nreac, npoints = self.getsolutionsize()
        # check values
        if npoints == 0 or nreac != self._nreactors:
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
        icreac = c_int(nreac)
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

    def getsolutionvariableprofile(self, varname):
        """
        Get the profile of the solution variable specified
        :param varname: name of the solution variable (string)
        :return: solution value profile (1d double array)
        """
        if not self.getrawsolutionstatus():
            print(
                Color.YELLOW
                + "** please use getsolution method to process the solution first"
            )
            return 1
        # check variable name
        vname = varname.rstrip()
        if vname.lower() in self._solution_tags:
            # is a property variable?
            vname = vname.lower()
        else:
            if vname not in self._specieslist:
                # is not a species?
                print(Color.PURPLE + f"** Error: variable name {vname}")

        # create variable arrays to hold the solution profile
        var = np.zeros(self._numbsolutionpoints, dtype=np.double)
        # get variable profile from the raw solution data
        var = self._solution_rawarray.get(vname)
        return var

    def createsolutionmixtures(self, specfrac):
        """
        Create a list of Mixtures that represent the gas inside the reactor at a solution point
        :param specfrac: species fractions of all time points numb_species x numb_solution_point (2D double array)
        :return: iError error flag (integer scalar)
        """
        if not self.getrawsolutionstatus():
            print(
                Color.YELLOW
                + "** please use getsolution method to process the raw solution first"
            )
            return 1
        # create a temporary Mixture object to hold the mixture properties at current solution point
        smixture = copy.deepcopy(self.reactormixture)
        # create variable arrays to hold the solution profile
        species = []
        # create a species fraction array to hold the solution species fraction profiles
        frac = np.zeros(self.numbspecies, dtype=np.double)
        # get solution variable profile from the raw solution arrays
        pres = self.getsolutionvariableprofile("pressure")
        temp = self.getsolutionvariableprofile("temperature")
        vol = self.getsolutionvariableprofile("volume")
        # loop over all species
        for sp in self._specieslist:
            species.append(self.getsolutionvariableprofile(sp))
        # loop over all solution points
        for i in range(self._numbsolutionpoints):
            # get mixture properties at the current solution point
            # pressure [dynes/cm2]
            smixture.pressure = pres[i]
            # temperature [K]
            smixture.temperature = temp[i]
            # mixture voluke [cm3]
            smixture.volume = vol[i]
            # species composition
            for k in range(self.numbspecies):
                frac[k] = specfrac[k, i]
            # set mixture composition
            if self._speciesmode == "mass":
                # mass fractions
                smixture.Y = frac
            else:
                # mole fractions
                smixture.X = frac
            # add to the solution mixture list
            self._solution_mixturearray.append(copy.deepcopy(smixture))
        # clean up
        species.clear()
        del pres, temp, vol, frac, species, smixture
        return 0

    def getsolutionmixture(self, time):
        """
        Get the mixture representing the solution state inside the reactor at the given time
        :param time: time point value [sec] (double scalar)
        :return: state a Mixture object representing the mixture properties of the reactor at time
        """
        # check status
        if not self.getmixturesolutionstatus():
            print(
                Color.YELLOW
                + "** error: use processsolution first to process the raw solution data",
                end=Color.END,
            )
            raise
        # get the time point array
        timearray = self.getsolutionvariableprofile("time")
        # find the interpolation parameters
        ileft, ratio = findinterpolateparameters(time, timearray)
        # find the mixture
        if ratio == 0.0e0:
            # get the mixtures
            mixtureleft = copy.deepcopy(self._solution_mixturearray[ileft])
            return mixtureleft
        elif ratio == 1.0e0:
            # get the mixtures
            mixtureright = copy.deepcopy(self._solution_mixturearray[ileft + 1])
            return mixtureright
        else:
            # get the mixtures
            mixtureleft = copy.deepcopy(self._solution_mixturearray[ileft])
            mixtureright = copy.deepcopy(self._solution_mixturearray[ileft + 1])
            # interpolate the mixture properties
            mixturetarget = interpolatemixtures(mixtureleft, mixtureright, ratio)
            # clean up
            del mixtureleft, mixtureright
            #
            return mixturetarget

    def getsolutionmixtureatindex(self, solution_index):
        """
        Get the mixture representing the solution state inside the reactor at the given solution point index
        :param solution_index: 0-base time point index (integer scalar)
        :return: state a Mixture object representing the mixture properties of the reactor at the solution time point
        """
        # check status
        if not self.getmixturesolutionstatus():
            print(
                Color.YELLOW
                + "** error: use processsolution first to process the raw solution data",
                end=Color.END,
            )
            raise
        # check index
        if solution_index > self._numbsolutionpoints - 1:
            print(
                Color.PURPLE
                + f"** error: the given index {solution_index} > the max time points index {self._numbsolutionpoints-1}"
            )
            print("   the solution index is 0-based", end=Color.END)
            raise
        # get the mixture
        mixturetarget = copy.deepcopy(self._solution_mixturearray[solution_index])
        return mixturetarget


class GivenPressureBatchReactor_FixedTemperature(BatchReactors):
    """
    Chemkin 0-D transient closed homogeneous reactor model
    with given reactor pressure (CONP) and reactor temperature (TGIV)
    """

    def __init__(self, reactor_condition, label=None):
        # set default label
        if label is None:
            label = "CONPT"
        # initialize the base module
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("Batch"))
        self._solvertype = c_int(self.SolverTypes.get("Transient"))
        self._problemtype = c_int(self.ProblemTypes.get("CONP"))
        self._energytype = c_int(self.EnergyTypes.get("GivenT"))
        # defaults for all closed homogeneous reactor models
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        # solver parameters
        self._absolutetolerance = 1.0e-12
        self._relativetolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]
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
        # if full-keyword mode is turned ON
        if not Keyword.noFullKeyword:
            # populate the reactor setup keywords
            self.setreactortypekeywords()

    @property
    def time(self):
        """
        Get simulation end time (required) [sec] (float scalar)
        """
        return self._endtime.value

    @time.setter
    def time(self, value=0.0e0):
        """
        Set simulation end time (required)
        default value = 0.0 sec
        :param value: simulation end time [sec] (float scalar)
        """
        if value <= 0.0e0:
            print(
                Color.PURPLE + "** simulation end time must be > 0",
                end=Color.END,
            )
        else:
            self._inputcheck.append("TIME")
            self._endtime = c_double(value)

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


class GivenPressureBatchReactor_EnergyConservation(BatchReactors):
    """
    Chemkin 0-D transient closed homogeneous reactor model
    with given reactor pressure (CONP) and
    solving the energy equation (ENRG)
    """

    def __init__(self, reactor_condition, label=None):
        # set default label
        if label is None:
            label = "CONP"
        # initialize the base module
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("Batch"))
        self._solvertype = c_int(self.SolverTypes.get("Transient"))
        self._problemtype = c_int(self.ProblemTypes.get("CONP"))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY"))
        # defaults for all closed homogeneous reactor models
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        self._heattransfercoefficient = 0.0e0
        self._ambienttemperature = 3.0e2
        self._heattransferarea = 0.0e0
        # solver parameters
        self._absolutetolerance = 1.0e-12
        self._relativetolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]
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
        # if full-keyword mode is turned ON
        if not Keyword.noFullKeyword:
            # populate the reactor setup keywords
            self.setreactortypekeywords()

    @property
    def time(self):
        """
        Get simulation end time (required) [sec] (float scalar)
        """
        return self._endtime.value

    @time.setter
    def time(self, value=0.0e0):
        """
        Set simulation end time (required)
        default value = 0.0 sec
        :param value: simulation end time [sec] (float scalar)
        """
        if value <= 0.0e0:
            print(
                Color.PURPLE + "** simulation end time must be > 0",
                end=Color.END,
            )
        else:
            self._inputcheck.append("TIME")
            self._endtime = c_double(value)

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

    def setheattransferareaprofile(self, x, area):
        """
        Specify reactor heat transfer area profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param area: heat transfer area value of the profile data [cm2] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.EnergyTypes.get(self._energytype.value) == "GivenT":
            print(
                Color.PURPLE
                + "** cannot specify heat transfer area of a fixed temperature batch reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "AEXT"
            iErr = self.setprofile(key=keyword, x=x, y=area)
            return iErr

    def setheatlossprofile(self, x, Qloss):
        """
        Specify reactor heat loss rate profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param Qloss: heat loss rate value of the profile data [cal/sec] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.EnergyTypes.get(self._energytype.value) == "GivenT":
            print(
                Color.PURPLE
                + "** cannot specify heat loss rate of a fixed temperature batch reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "QPRO"
            iErr = self.setprofile(key=keyword, x=x, y=Qloss)
            return iErr


class GivenVolumeBatchReactor_FixedTemperature(BatchReactors):
    """
    Chemkin 0-D transient closed homogeneous reactor model
    with given reactor volume (CONV) and reactor temperature (TGIV)
    """

    def __init__(self, reactor_condition, label=None):
        # set default label
        if label is None:
            label = "CONVT"
        # initialize the base module
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("Batch"))
        self._solvertype = c_int(self.SolverTypes.get("Transient"))
        self._problemtype = c_int(self.ProblemTypes.get("CONV"))
        self._energytype = c_int(self.EnergyTypes.get("GivenT"))
        # defaults for all closed homogeneous reactor models
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        # solver parameters
        self._absolutetolerance = 1.0e-12
        self._relativetolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]
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
        # if full-keyword mode is turned ON
        if not Keyword.noFullKeyword:
            # populate the reactor setup keywords
            self.setreactortypekeywords()

    @property
    def time(self):
        """
        Get simulation end time (required) [sec] (float scalar)
        """
        return self._endtime.value

    @time.setter
    def time(self, value=0.0e0):
        """
        Set simulation end time (required)
        default value = 0.0 sec
        :param value: simulation end time [sec] (float scalar)
        """
        if value <= 0.0e0:
            print(
                Color.PURPLE + "** simulation end time must be > 0",
                end=Color.END,
            )
        else:
            self._inputcheck.append("TIME")
            self._endtime = c_double(value)

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


class GivenVolumeBatchReactor_EnergyConservation(BatchReactors):
    """
    Chemkin 0-D transient closed homogeneous reactor model
    with given reactor volume (CONV) and
    solving the energy equation (ENRG)
    """

    def __init__(self, reactor_condition, label=None):
        # set default label
        if label is None:
            label = "CONV"
        # initialize the base module
        super().__init__(reactor_condition, label)
        # set reactor type
        self._reactortype = c_int(self.ReactorTypes.get("Batch"))
        self._solvertype = c_int(self.SolverTypes.get("Transient"))
        self._problemtype = c_int(self.ProblemTypes.get("CONV"))
        self._energytype = c_int(self.EnergyTypes.get("ENERGY"))
        # defaults for all closed homogeneous reactor models
        self._nreactors = 1
        self._npsrs = c_int(1)
        self._ninlets = c_int(0)
        self._nzones = c_int(0)
        # heat transfer parameters
        self._heatlossrate = c_double(0.0e0)
        self._heattransfercoefficient = 0.0e0
        self._ambienttemperature = 3.0e2
        self._heattransferarea = 0.0e0
        # solver parameters
        self._absolutetolerance = 1.0e-12
        self._relativetolerance = 1.0e-6
        # required inputs: (1) end time
        self._numb_requiredinput = 1
        self._requiredlist = ["TIME"]
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
        # if full-keyword mode is turned ON
        if not Keyword.noFullKeyword:
            # populate the reactor setup keywords
            self.setreactortypekeywords()

    @property
    def time(self):
        """
        Get simulation end time (required) [sec] (float scalar)
        """
        return self._endtime.value

    @time.setter
    def time(self, value=0.0e0):
        """
        Set simulation end time (required)
        default value = 0.0 sec
        :param value: simulation end time [sec] (float scalar)
        """
        if value <= 0.0e0:
            print(
                Color.PURPLE + "** simulation end time must be > 0",
                end=Color.END,
            )
        else:
            self._inputcheck.append("TIME")
            self._endtime = c_double(value)

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

    def setheattransferareaprofile(self, x, area):
        """
        Specify reactor heat transfer area profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param area: heat transfer area value of the profile data [cm2] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.EnergyTypes.get(self._energytype.value) == "GivenT":
            print(
                Color.PURPLE
                + "** cannot specify heat transfer area of a fixed temperature batch reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "AEXT"
            iErr = self.setprofile(key=keyword, x=x, y=area)
            return iErr

    def setheatlossprofile(self, x, Qloss):
        """
        Specify reactor heat loss rate profile
        :param x: position value of the profile data [cm or sec] (double array)
        :param Qloss: heat loss rate value of the profile data [cal/sec] (double array)
        :return: error code (integer scalar)
        """
        if BatchReactors.EnergyTypes.get(self._energytype.value) == "GivenT":
            print(
                Color.PURPLE
                + "** cannot specify heat loss rate of a fixed temperature batch reactor",
                end=Color.END,
            )
            return 10
        else:
            keyword = "QPRO"
            iErr = self.setprofile(key=keyword, x=x, y=Qloss)
            return iErr

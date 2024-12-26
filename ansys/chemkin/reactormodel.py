import copy
import ctypes
from ctypes import c_double, c_int
import logging

import numpy as np

from . import chemkin_wrapper
from .chemistry import Patm, checkchemistryset, chemistrysetinitialized, verbose
from .color import Color
from .mixture import Mixture

logger = logging.getLogger(__name__)


#
# Base class for keyword data
#
class Keyword:
    """
    A Chemkin keyword
    """

    # supported Chemkin keyword data types
    _keyworddatatypes = ["bool", "int", "float", "str"]
    _valuetypes = (bool, int, float, str)
    # required keywords that are given as reactor properties or as mixture properties
    # and will be set by using the KINAll0D_SetupBatchInputs call
    _protectedkeywords = [
        "CONP",
        "CONV",
        "TRAN",
        "STST",
        "TGIV",
        "ENRG",
        "PRES",
        "TEMP",
        "TAU",
        "TIME",
        "XEND",
        "FLRT",
        "VDOT",
        "SCCM",
        "VDOT",
        "DIAM",
        "AREA",
        "REAC",
        "GAS",
        "INIT",
        "XEST",
        "SURF",
        "ACT",
        "TINL",
        "FUEL",
        "OXID",
        "PROD",
        "ASEN",
        "ATLS",
        "RTLS",
        "EPST",
        "EPSS",
    ]
    gasspecieskeywords = ["REAC", "XEST", "FUEL", "OXID"]
    flowratekeywords = ["FLRT", "VDOT", "VEL", "SCCM"]
    profilekeywords = [
        "TPRO",
        "PPRO",
        "VPRO",
        "QPRO",
        "AINT",
        "AEXT",
        "DPRO",
        "FPRO",
        "SCCMPRO",
        "VDOTPRO",
        "VELPRO",
        "TINPRO",
        "AFLO",
    ]
    fourspaces = "    "
    # Under the default API-call mode, important keywords (the _protectedkeywords) are set by direct API calls,
    # the rest of the keywords can be set by keyword input lines (i.e., using the setkeyword method).
    # Under the full-keyword mode, all keywords and their parameters are set by keyword input lines, and specifying
    # those _protectedkeywords via the setkeyword method are required.
    noFullKeyword = True  # default: API-call mode

    def __init__(self, phrase, value, data_type):
        """
        Initialize the Chemkin keyword
        :param phrase: Chemkin keyword phrase (string)
        :param value: value assigned to the Chemkin keyword (type indicated by the data_type scalar)
        :param data_type: data type of the value (string: 'int', 'float', 'string', or 'bool')
        """
        self._set = False
        iErr = 0
        # check value data type
        if data_type not in Keyword._keyworddatatypes:
            # the declared data type is not supported
            print(
                Color.PURPLE + f"** unsupported data type specified {data_type}",
                end=Color.END,
            )
            if not isinstance(value, (bool, int, float, str)):
                # value does not match the declared data type
                print(
                    Color.PURPLE + f"** variable has different data type {type(value)}",
                    end=Color.END,
                )
            iErr = 1
        # block the protected keywords
        if Keyword.noFullKeyword:
            if phrase.upper() in Keyword._protectedkeywords:
                print(
                    Color.PURPLE
                    + f"** use reactor property setter to assign '{phrase}' value"
                )
                print(
                    "   for example, to set the reactor volume use: 'MyBatchReactor.volume = 100'",
                    end=Color.END,
                )
                iErr = 2
        if iErr > 0:
            return
        self._key = phrase.upper()  # Chemkin keyword phrase
        self._value = value  # value assigned to the keyword
        self._data_type = data_type  # data type of the values
        self._prefix = ""  # a prefix to the keyword that can be used
        # to comment out/disable the keyword by setting it to '!'
        self._set = True

    @staticmethod
    def setfullkeywords(mode):
        """
        All keywords and their parameters must be specified by using the setkeyword method
        and will be passed to the reactor model for further processing
        :param mode: True/False turn the full keyword mode ON/OFF (boolean scalar)
        :return: None
        """
        if mode:
            # turn ON the full keyword mode (no checking on protected keywords)
            Keyword.noFullKeyword = False
        else:
            Keyword.noFullKeyword = True

    def show(self):
        """
        return the Chemkin keyword and its parameter value
        :return: None
        """
        if self._set:
            if isinstance(self._value, (int, float)):
                print(f'** keyword "{self._key:s}": value = {self._value}')
            elif isinstance(self._value, bool):
                if self._value:
                    print(f'** keyword "{self._key:s}": value = True')
                else:
                    print(
                        f'** keyword "{self._prefix:1s}{self._key:s}": value = Disabled'
                    )
            else:
                print(f'** keyword "{self._key:s}": value = {self._value}')
        else:
            print(f'** keyword "{self._key:s}": value not set')

    def resetvalue(self, value):
        """
        Reset the parameter value of an existing keyword
        :param value: keyword parameter (int, float, string, or bool scalar depending on the keyword)
        :return: None
        """
        if isinstance(value, Keyword._valuetypes):
            if isinstance(value, bool):
                if value:
                    # true: keep the keyword active
                    self._prefix = ""
                    self._value = True
                else:
                    # false: disable the keyword
                    self._prefix = "!"
                    self._value = False
            else:
                # integer, float, or string parameter
                self._value = value
        else:
            print(
                Color.PURPLE
                + f"** value has a wrong data type {type(value)}, value will not be reset"
            )
            print(
                f"   data type expected by keyword {self._key:s} is {self._data_type:s}",
                end=Color.END,
            )

    def parametertype(self):
        """
        get parameter type of the keyword
        :return: parameter data type (string: 'int', 'float', 'string', or 'bool')
        """
        return type(self._value)

    @property
    def value(self):
        """
        Get parameter value of the keyword
        :return: parameter value (integer, floating, string, or boolean scalar depending on the keyword)
        """
        # extract the keyword value
        if self._data_type == "bool":
            mode = self._prefix != "!"
            return mode
        else:
            return self._value

    def getvalue_as_string(self):
        """
        Create the keyword input line for Chemkin applications
        :return: line length, line (integer scalar, string)
        """
        # initialization
        line = ""
        linelength = 0
        # assembly the keyword line
        if self._data_type == "bool":
            # boolean keyword (active or disabled by '!')
            line = self._prefix + self._key
        else:
            # integer, double, or string parameter
            line = self._prefix + self._key + Keyword.fourspaces + str(self._value)

        linelength = len(line)
        return linelength, line


#
# This keyword type is used to distinguish keywords that act as on/off switches by their presence
#
class BooleanKeyword(Keyword):
    """
    Chemkin boolean keyword
    """

    def __init__(self, phrase):
        """
        set up a Chemkin keyword with a boolean parameter or with no parameter
        :param phrase: Chemkin keyword phrase (string)
        """
        value = True
        super().__init__(phrase, value, "bool")


#
# This keyword type is used to hold integer keyword types (not sure if there actually are any of these)
#
class IntegerKeyword(Keyword):
    """
    A Chemkin integer keyword
    """

    def __init__(self, phrase, value=0):
        """
        set up a Chemkin keyword with an integer parameter
        :param phrase: Chemkin keyword phrase (string)
        :param value: parameter value (integer scalar)
        """
        super().__init__(phrase, value, "int")


#
# This keyword type is used to hold real keyword types
#
class RealKeyword(Keyword):
    """
    A Chemkin real keyword
    """

    def __init__(self, phrase, value=0.0e0):
        """
        set up a Chemkin keyword with a real number (floating number) parameter
        :param phrase: Chemkin keyword phrase (string)
        :param value: parameter value (floating number scalar)
        """
        super().__init__(phrase, value, "float")


#
# This keyword type is used to hold string keyword types
#
class StringKeyword(Keyword):
    """
    A Chemkin string keyword
    """

    def __init__(self, phrase, value=""):
        """
        set up a Chemkin keyword with a string parameter
        :param phrase: Chemkin keyword phrase (string)
        :param value: parameter value (string)
        """
        if len(value) <= 0:
            print(Color.PURPLE + "** no string parameter given", end=Color.END)
            return
        super().__init__(phrase, value, "str")


class Profile:
    """
    Generic Chemkin profile keyword class
    """

    def __init__(self, key, x, y):
        """
        Create a profile object
        :param key: profile keyword (string scalar)
        :param x: position of the profile data points (double array)
        :param y: variable value of the profile data (double array)
        """
        # initialization
        self._profilekeyword = ""
        self._status = 0
        # check
        if key.upper() in Keyword.profilekeywords:
            self._profilekeyword = key.upper()
        else:
            print(
                Color.PURPLE + "** profile is not available under the reactor model",
                end=Color.END,
            )
            self._status = -1
            return
        # profile data sizes
        xsize = len(x)
        ysize = len(y)
        if xsize == ysize:
            self._size = xsize
            # independent variable (time, location, grid, ...)
            if isinstance(x, np.double):
                self._pos = copy.deepcopy(x)
            else:
                self._pos = np.array(x, dtype=np.double)
            # dependent variable value at the corresponding position
            if isinstance(y, np.double):
                self._val = copy.deepcopy(y)
            else:
                self._val = np.array(y, dtype=np.double)
        else:
            print(
                Color.PURPLE
                + "** the number of positions does not match the number of values"
            )
            print(f"   number of positions = {xsize:d}")
            print(f"   number of values    = {ysize:d}", end=Color.END)
            self._status = -2
            return

    @property
    def size(self):
        """
        Get number of data points in the profile
        :return: number of points (integer scalar)
        """
        return self._size

    @property
    def status(self):
        """
        Get the validity of the profile object
        :return: status (integer scalar)
        """
        return self._status

    @property
    def pos(self):
        """
        Get position values of profiles data
        :return: position [sec, cm] (double array)
        """
        return self._pos

    @property
    def value(self):
        """
        Get variable values of profile data
        :return: variable value (double array)
        """
        return self._val

    @property
    def profilekey(self):
        """
        Get profile keyword
        :return: keyword (string)
        """
        return self._profilekeyword

    def show(self):
        """
        Show the profile data
        :return: None
        """
        print(f"profile size: {self._size:d}")
        print(f" position           {self._profilekeyword:s}  ")
        for i in range(self._size):
            print(f"{self._pos[i]:f}         {self._val[i]}")

    def resetprofile(self, size, x, y):
        # check array size
        if size == self._size:
            # new profile has the same size
            self._pos[:] = x[:]
            self._val[:] = y[:]
        else:
            # new profile has different size
            self._size = size
            # resize the arrays
            self._pos.resize(size, refcheck=False)
            self._val.resize(size, refcheck=False)
        # fill the arrays with new values
        self._pos[:] = x[:]
        self._val[:] = y[:]

    def getprofile_as_string_list(self):
        """
        Create the keyword input lines as a list for Chemkin applications
        :return: number of profile lines, line (integer scalar, string list)
        """
        # initialization
        line = []
        # special treatment for pressure profile
        factor = 1.0e0
        if self._profilekeyword == "PPRO":
            if Keyword.noFullKeyword:
                # use API calls: pressure profile units = dynes/cm2
                pass
            else:
                # use Full Keywords: pressure units = atm
                factor = Patm
        # assembly the profile keyword lines
        for i in range(self._size):
            thisline = ""
            thisline = (
                self._profilekeyword
                + Keyword.fourspaces
                + str(self._pos[i])
                + Keyword.fourspaces
                + str(self._val[i] / factor)
            )
            line.append(thisline)
        return self._size, line


#
# Framework and generic base classes for running Chemkin reactor models,
# defining methods to set chemistry, process keywords, and run
#
class ReactorModel:
    """
    A generic Chemkin reactor model framework
    """

    def __init__(self, reactor_condition, label):
        """
        Initialize the basic parameters of Chemkin reactor model
        :param reactor_condition: mixture containing the initial/estimate reactor pressure, temperature, and gas composition (Mixture object)
        :param label: reactor label/name (string scalar)
        """
        # check mixture
        if not isinstance(reactor_condition, Mixture):
            print(
                Color.RED + "** the first argument must be a Mixture object",
                end=Color.END,
            )
            raise TypeError
        iErr = reactor_condition.validate()
        if iErr != 0:
            print(Color.YELLOW + "** mixture is not fully defined", end=Color.END)
            raise TypeError
        # initialization
        self.label = label
        # chemistry set index
        self._chemset_index = ctypes.c_int(reactor_condition.chemID)
        # mixture
        self.reactormixture = copy.deepcopy(reactor_condition)
        # mixture gas species symbols
        self._specieslist = reactor_condition._specieslist  # gas species symbols
        # mixture temperature [K]
        self._temperature = ctypes.c_double(reactor_condition.temperature)
        # mixture pressure [dynes/cm2]
        self._pressure = ctypes.c_double(reactor_condition.pressure)
        # number of required input
        self._numb_requiredinput = 0
        self._inputcheck = []
        # gas reaction rate multiplier
        self._gasratemultiplier = 1.0e0
        # write text output file
        self._TextOut = True
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(154)
        # write XML solution file
        self._XMLOut = True
        # number of keywords used
        self._numbkeywords = 0
        # list of keyword phrases used for easy searching
        self._keyword_index = []
        # list of keyword objects defined
        self._keyword_list = []
        # list of keyword lines
        # (each line is a string consists of: '<keyword> <parameter>', i.e., _keyword_index + _keyword_parameters)
        self._keyword_lines = []
        # number of keyword lines
        self._numblines = 0
        # length of each keyword line
        self._linelength = []
        # number of profile assigned
        self._numbprofiles = 0
        # list of profile keywords used for easy searching
        self._profiles_index = []
        # list of profile objects defined
        self._profiles_list = []
        # simulation run status
        #  -100 = not yet run
        #     0 = run success
        # other = run failed
        self.runstatus = -100
        # raw solution data structure
        self._solution_tags = [
            "time",
            "distance",
            "temperature",
            "pressure",
            "volume",
            "velocity",
            "flowrate",
        ]
        self.numbspecies = self.reactormixture._KK
        self._speciesmode = "mass"
        self._numbsolutionpoints = 0
        self._solution_rawarray = {}
        self._numbsolutionmixtures = 0
        self._solution_mixturearray = []
        # initialize output buffer
        self.output = {}
        # initialize KINetics
        if not checkchemistryset(self._chemset_index.value):
            # need to initialize KINetics
            print(Color.YELLOW + "** initializing chemkin...", end=Color.END)
            iErr = chemkin_wrapper.chemkin.KINInitialize(self._chemset_index, c_int(0))
            if iErr == 0:
                chemistrysetinitialized(self._chemset_index.value)
            else:
                print(Color.RED + "** fail to initialize KINetics", end=Color.END)

    def usefullkeywords(self, mode):
        """
        Specify all necessary keywords explicitly
        :param mode: True/False turn full keyword mode ON/OFF (boolean scalar)
        :return: None
        """
        Keyword.setfullkeywords(mode)
        if mode:
            print(
                Color.YELLOW
                + f"** reactor {self.label} will be run with full keyword input mode",
                end=Color.END,
            )

    def __findkeywordslot(self, key):
        """
        Find the proper index in the global keyword list to add a new keyword or to modify the keyword parameter
        :param key: Chemkin keyword (string scalar)
        :return: index in the global keyword list for the keyword, whether this is a new keyword (integer scalar, bool scalar)
        """
        # check existing keyword
        if self._numbkeywords == 0:
            return 0, True
        else:
            if key in self._keyword_index:
                return self._keyword_index.index(key), False
            else:
                # new keyword
                return self._numbkeywords, True

    def setkeyword(self, key, value):
        """
        Set a Chemkin keyword and its parameter
        :param key: Chemkin keyword phrase (string scalar)
        :param value: value of the parameter (int, float, string, or bool scalar depending on the keyword)
        :return: None
        """
        # find the keyword
        i, newkey = self.__findkeywordslot(key.upper())
        # add the keyword to the keywords list
        if newkey:
            # a new keyword
            if isinstance(value, str):
                # value is a string
                self._keyword_list.append(StringKeyword(key.upper(), value))
                self._keyword_index.append(key.upper())
            elif isinstance(value, bool):
                # value is a boolean value
                if value:
                    # set the keyword only if the value is True
                    self._keyword_list.append(BooleanKeyword(key.upper()))
                    self._keyword_index.append(key.upper())
                else:
                    # remove the count
                    self._numbkeywords -= 1
            elif isinstance(value, int):
                # value is an integer
                self._keyword_list.append(IntegerKeyword(key.upper(), value))
                self._keyword_index.append(key.upper())
            elif isinstance(value, float):
                # value is a real number
                self._keyword_list.append(RealKeyword(key.upper(), value))
                self._keyword_index.append(key.upper())
            else:
                print(
                    Color.PURPLE + "** invalid keyword value data type",
                    end=Color.END,
                )
                return
            self._numbkeywords += 1
        else:
            # an existing keyword, just update its value
            if isinstance(value, (str, bool)):
                self._keyword_list[i].resetvalue(value)
            elif isinstance(value, (int, float)):
                self._keyword_list[i].resetvalue(value)
            else:
                print(
                    Color.PURPLE + "** invalid keyword value data type",
                    end=Color.END,
                )
                return

    def showkeywordinputlines(self):
        """
        list all currently-defined keywords and their parameters line by line
        :return: None
        """
        # header
        print("** INPUT KEYWORDS: \n")
        print("=" * 40)
        # display the keyword and the parameters line by line
        for k in self._keyword_list:
            n, line = k.getvalue_as_string()
            print(f"{line[:n]:s}")
        print("=" * 40)

    def createkeywordinputlines(self):
        """
        Create keyword input lines for Chemkin applications
        one keyword per line: <keyword>     <parameter>
        :return: Error code, number of lines (integer scalar, integer scalar)
        """
        # initialization
        self._numblines = 0
        self._linelength.clear()
        self._keyword_lines.clear()
        # create the keyword lines from the keyword objects in the keyword list
        for k in self._keyword_list:
            n, line = k.getvalue_as_string()
            self._linelength.append(n)
            self._keyword_lines.append(line)
            self._numblines += 1
        # print the entire keyword input block
        if verbose():
            print("** INPUT KEYWORDS:")
            # print(f'number of keyword input lines: {self._numblines:d} == {self._numbkeywords:d} \n')
            print("=" * 40)
            for line in self._keyword_lines:
                print(line)
            print("=" * 40)
        iErr = self._numbkeywords - self._numblines
        return iErr, self._numblines

    def __findprofileslot(self, key):
        """
        Find the proper index in the global profile list either to add a new profile or to modify the existing profile parameter
        :param key: Chemkin profile keyword (string scalar)
        :return: index in the global profile list for the keyword, whether this is a new profile (integer scalar, bool scalar)
        """
        # check existing keyword
        if self._numbprofiles == 0:
            return 0, True
        else:
            if key in self._profiles_index:
                return self._profiles_index.index(key), False
            else:
                # new keyword
                return self._numbprofiles, True

    def setprofile(self, key, x, y):
        """
        Set a Chemkin profile and its parameter
        :param key: Chemkin profile keyword phrase (string scalar)
        :param x: position values of the profile data (double array)
        :param y: variable values of the profile data (double array)
        :return: Error code (integer scalar)
        """
        #
        iErr = 0
        # find the keyword
        i, newprofile = self.__findprofileslot(key.upper())
        # add the profile to the profiles index list
        if newprofile:
            # a new profile
            self._profiles_list.append(Profile(key.upper(), x, y))
            status = self._profiles_list[i].status
            if status == 0:
                self._profiles_index.append(key.upper())
                self._numbprofiles += 1
            else:
                print(Color.PURPLE + "** fail to create the profile '{key}'")
                print(f"   error code = {status:d}", end=Color.END)
                iErr = status
        else:
            # an existing keyword, just update its value
            xsize = len(x)
            ysize = len(y)
            if xsize == ysize:
                self._profiles_list[i].resetprofile(xsize, x, y)
            else:
                print(
                    Color.PURPLE
                    + "** the number of positions does not match the number of values"
                )
                print(f"   number of positions = {xsize:d}")
                print(f"   number of values    = {ysize:d}", end=Color.END)
                iErr = 1
        return iErr

    def createprofileinputlines(self):
        """
        Create profile keyword input lines for Chemkin applications
        one keyword per line: <profile keyword>     <position>  <value>
        :return: Error code, number of lines, list of string lists (integer scalar, integer scalar, list of string lists)
        """
        # initialization
        numblines = 0
        numbprofiles = 0
        keyword_lines = []
        # create the keyword lines from the keyword objects in the profile list
        for p in self._profiles_list:
            n, lines = p.getprofile_as_string_list()
            keyword_lines.append(lines)
            numblines += n
            numbprofiles += 1
            # print the entire keyword input block per profile
            if verbose():
                print("** PROFILE KEYWORDS:")
                print(f"{n:d} keyword input lines in {p._profilekeyword} profile\n")
                print("=" * 40)
                for line in lines:
                    print(line)
                print("=" * 40)
        # lines: list of strings of a profile ['VPRO x1 v1', 'VPRO x2 v2', ...]
        # keyword_lines: list of lines:  [['VPRO x1 v1', 'VPRO x2 v2', ...], ['PPRO x1 p1', 'PPRO x2 p2', ..] , ... ]
        iErr = numbprofiles - self._numbprofiles
        return iErr, numblines, keyword_lines

    def createspeciesinputlines(self, solvertype, threshold=1.0e-12, molefrac=None):
        """
        Create keyword input lines for initial/estimated species mole fraction inside the batch reactor
        :param solvertype: solver type of the reactor model (integer scalar)
        :param threshold: minimum species mole fraction value to be included in the species keyword (double scalar)
        :param molefrac: species composition in mole fractions (1-D double array)
        :return: Number of keyword lines, list of keyword line strings (integer scalar, list of strings)
        """
        # initial(transient)/estimate(steady-state) composition keyword depends on the solver type
        key = Keyword.gasspecieskeywords[solvertype - 1]
        KSYM = self._specieslist
        lines = []
        numb_lines = 0
        for i in range(len(molefrac)):
            if molefrac[i] > threshold:
                thisline = (
                    key
                    + Keyword.fourspaces
                    + KSYM[i].rstrip()
                    + Keyword.fourspaces
                    + str(molefrac[i])
                )
                lines.append(thisline)
                numb_lines += 1
        return numb_lines, lines

    def createspeciesinputlineswithaddon(
        self, key="XEST", threshold=1.0e-12, molefrac=None, addon=""
    ):
        """
        Create keyword input lines for initial/estimated species mole fraction inside the batch reactor
        :param key: keyword for species value (string)
        :param threshold: minimum species mole fraction value to be included in the species keyword (double scalar)
        :param molefrac: species composition in mole fractions (1-D double array)
        :param addon: add-on string to the species input, usually the reactor/zone number (string)
        :return: Number of keyword lines, list of keyword line strings (integer scalar, list of strings)
        """
        # must use estimate composition keyword 'XEST'
        # (the 'REAC' keyword does not accept reactor/zone number)
        KSYM = self._specieslist
        ksize = len(molefrac)
        if ksize != len(KSYM):
            print(
                Color.PURPLE
                + f"** species mole fraction array has size {ksize}"
                + f" but {len(KSYM)} is expected",
                end=Color.END,
            )
            exit()

        lines = []
        numb_lines = 0
        for i in range(ksize):
            if molefrac[i] > threshold:
                thisline = (
                    key.rstrip()
                    + Keyword.fourspaces
                    + KSYM[i].rstrip()
                    + Keyword.fourspaces
                    + str(molefrac[i])
                    + Keyword.fourspaces
                    + addon.rstrip()
                )
                lines.append(thisline)
                numb_lines += 1
        return numb_lines, lines

    def chemID(self):
        """
        Get chemistry set index
        :return: chemistry set index (integer scalar)
        """
        return self._chemset_index.value

    @property
    def temperature(self):
        """
        Get reactor initial temperature
        :return: temperature [K] (double scalar)
        """
        return self.reactormixture.temperature

    @temperature.setter
    def temperature(self, t):
        """
        (Re)set reactor temperature
        :param t: temperature [K] (double scalar)
        :return: None
        """
        if t <= 1.0e1:
            print(Color.PURPLE + "** invalid temperature value", end=Color.END)
            pass
        self._temperature = c_double(t)
        self.reactormixture.temperature = t

    @property
    def pressure(self):
        """
        Get reactor pressure
        :return: pressure [dynes/cm2] (double scalar)
        """
        return self.reactormixture.pressure

    @pressure.setter
    def pressure(self, p):
        """
        (Re)set reactor pressure
        :param p: pressure [dynes/cm2] (double scalar)
        :return: None
        """
        if p <= 0.0e0:
            print(Color.PURPLE + "** invalid pressure value", end=Color.END)
            pass
        self._pressure = c_double(p)
        self.reactormixture.pressure = p

    @property
    def massfraction(self):
        """
        Get the initial/guessed/estimate gas species mass fractions inside the reactor
        :return: mixture mass fraction (double array)
        """
        return self.reactormixture.Y

    @massfraction.setter
    def massfraction(self, recipe):
        """
        (Re)set the initial/guessed/estimate gas species mass fractions inside the reactor
        :param recipe: mixture composition a list of [('species_symbol', mass_fraction), ('species_symbol', mass_fraction), ...]
        :return: None
        """
        self.reactormixture.Y(recipe)

    @property
    def molefraction(self):
        """
        Get the initial/guessed/estimate gas species mole fractions inside the reactor
        :return: mixture mole fraction (double array)
        """
        return self.reactormixture.X

    @molefraction.setter
    def molefraction(self, recipe):
        """
        (Re)set the initial/guessed/estimate gas species mole fractions inside the reactor
        :param recipe: mixture composition a list of [('species_symbol', mole_fraction), ('species_symbol', mole_fraction), ...]
        :return: None
        """
        self.reactormixture.X(recipe)

    @property
    def concentration(self):
        """
        Get the initial/guessed/estimate gas species molar concentrations inside the reactor
        :return: mixture molar concentration [mole/cm3] (double array)
        """
        return self.reactormixture.concentration

    def listcomposition(self, mode, option=" ", bound=0.0e0):
        """
        List the gas mixture composition inside the reactor
        :param mode: flag indicates the fractions returned are 'mass' or 'mole' fractions
        :param option: flag indicates to list 'all' species or just the species with non-zero fraction (default)
        :param bound: minimum fraction value for the species to be printed (double scalar)
        :return: None
        """
        self.reactormixture.listcomposition(mode=mode, option=option, bound=bound)

    @property
    def gasratemultiplier(self):
        """
        Get the value of the gas-phase reaction rate multiplier (optional)
        :return: gas-phase reaction rate multiplier (float scalar)
        """
        return self._gasratemultiplier

    @gasratemultiplier.setter
    def gasratemultiplier(self, value=1.0e0):
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
            self._gasratemultiplier = value
            self.setkeyword(key="GFAC", value=value)

    @property
    def STD_Output(self):
        """
        Get text output status (optional)
        :return: text output ON=True/OFF=False (boolean scala)
        """
        return self._TextOut

    @STD_Output.setter
    def STD_Output(self, mode):
        """
        Set text output status (optional)
        default value = True: always write to the text output file
        :param mode: True/False (turn ON/turn OFF) (boolean scalar)
        :return: None
        """
        off = not mode
        self.setkeyword(key="NO_SDOUTPUT_WRITE", value=off)
        self._TextOut = mode

    @property
    def XML_Output(self):
        """
        Get XML solution output status (optional)
        :return: XML solution output ON=True/OFF=False (boolean scala)
        """
        return self._XMLOut

    @XML_Output.setter
    def XML_Output(self, mode):
        """
        Set XML solution output status (optional)
        default value = True: always create the XML solution file
        :param mode: True/False (turn ON/turn OFF) (boolean scalar)
        :return: None
        """
        off = not mode
        self.setkeyword(key="NO_XMLOUTPUT_WRITE", value=off)
        self._XMLOut = mode

    def setsensitivityanalysis(
        self,
        mode=True,
        absolute_tolerance=None,
        relative_tolerance=None,
        temperature_threshold=None,
        species_threshold=None,
    ):
        """
        Switch ON/OFF A-factor sensitivity analysis
        :param mode: True/False (turn A-factor sensitivity ON/OFF) (boolean scalar)
        :param absolute_tolerance: absolute tolerance of the sensitivity parameters (float scalar)
        :param relative_tolerance: relative tolerance of the sensitivity parameters (float scalar)
        :param temperature_threshold: threshold normalized temperature sensitivity parameter value to print out to the text output file (float scalar)
        :param species_threshold: threshold normalized species sensitivity parameter value to print out to the text output file (float scalar)
        :return: None
        """
        if "ASEN" in self._keyword_index:
            # already defined
            i = self._keyword_index.index("ASEN")
            if mode:
                # reactivate the keyword if it is disabled
                if self._keyword_list[i]._prefix == "!":
                    self._keyword_list[i]._prefix = ""
                # set tolerances if given
                if absolute_tolerance is not None:
                    self.setkeyword(key="ATLS", value=absolute_tolerance)
                if relative_tolerance is not None:
                    self.setkeyword(key="RTLS", value=relative_tolerance)
                # reset the thresholds
                if temperature_threshold is not None:
                    self.setkeyword(key="EPST", value=temperature_threshold)
                if species_threshold is not None:
                    self.setkeyword(key="EPSS", value=species_threshold)
            else:
                # disable the keyword
                if self._keyword_list[i]._prefix != "!":
                    self._keyword_list[i]._prefix = "!"
        else:
            # not defined
            if mode:
                # enable sensitivity analysis
                self.setkeyword(key="ASEN", value=mode)
                # set sensitivity analysis related parameters
                if absolute_tolerance is not None:
                    self.setkeyword(key="ATLS", value=absolute_tolerance)
                if relative_tolerance is not None:
                    self.setkeyword(key="RTLS", value=relative_tolerance)
                if temperature_threshold is not None:
                    self.setkeyword(key="EPST", value=temperature_threshold)
                if species_threshold is not None:
                    self.setkeyword(key="EPSS", value=species_threshold)
            else:
                # do nothing
                pass

    def setROPanalysis(self, mode=True, threshold=None):
        """
        Switch ON/OFF the ROP (Rate Of Production) analysis
        :param mode: True/False (turn ROP ON/OFF) (boolean scalar)
        :param threshold: threshold ROP value to print out to the text output file (float scalar)
        :return: None
        """
        if "AROP" in self._keyword_index:
            # already defined
            i = self._keyword_index.index("AROP")
            if mode:
                # reactivate the keyword if it is disabled
                if self._keyword_list[i]._prefix == "!":
                    self._keyword_list[i]._prefix = ""
                # reset the threshold
                if threshold is not None:
                    self.setkeyword(key="EPSR", value=threshold)
            else:
                # disable the keyword
                if self._keyword_list[i]._prefix != "!":
                    self._keyword_list[i]._prefix = "!"
        else:
            # not defined
            if mode:
                # enable ROP analysis
                self.setkeyword(key="AROP", value=mode)
                if threshold is not None:
                    self.setkeyword(key="EPSR", value=threshold)
            else:
                # do nothing
                pass

    @property
    def realgas(self):
        """
        Get the real gas EOS status
        True: real gas EOS is turned ON
        :return: True/False (boolean scalar)
        """
        if "RLGAS" in self._keyword_index:
            # already defined
            i = self._keyword_index.index("RLGAS")
            if self._keyword_list[i]._prefix == "!":
                # commented out
                return False
            else:
                # is turned ON
                return True
        else:
            # has not been turned ON
            return False

    def userealgasEOS(self, mode):
        """
        Set the option to turn ON/OFF the real gas model
        :param mode: True/False (turn ROP ON/OFF) (boolean scalar)
        :return: None
        """
        # turn ON/OFF the real gas EOS
        self.setkeyword(key="RLGAS", value=mode)

    def setrealgasmixingmodel(self, model):
        """
        Set the real gas mixing rule/model
        :param model: 0/1 (integer scalar)
        :return: None
        """
        # set the real gas mixing model
        _mixingmodels = ["Van der Waals", "pseudocritical"]
        if model in [0, 1]:
            print(
                Color.YELLOW + f"** the {_mixingmodels[model]:s} mixing model is used",
                end=Color.END,
            )
            self.setkeyword(key="RLMIX", value=model)
        else:
            print(Color.PURPLE + f"** model index {model:d} is not valid")
            print(f"    set model = 0 to use the {_mixingmodels[0]} mixing model")
            print(
                f"    set model = 1 to use the {_mixingmodels[1]} mixing model",
                end=Color.END,
            )

    def setrunstatus(self, code):
        """
        Set the simulation run status
        :param code: status/error code (integer scalar)
        """
        self.runstatus = code

    def getrunstatus(self, mode="silent"):
        """
        Get the reactor model simuation status
        :param mode: 'verbose' or 'silent' (default) option for additional print information (string)
        :return: run status 0=success; -100=not run; other=failed (integer scalar)
        """
        if mode.lower() == "verbose":
            if self.runstatus == -100:
                print("** Simulation not run")
            elif self.runstatus == 0:
                print("** simulation run successfully")
            else:
                print(f"** simulation failed with code {self.runstatus}")

        return self.runstatus

    def __process_keywords(self):
        """
        Generic Chemkin reactor keyword processing method
        :return: Error code (integer scalar)
        """
        iErr = 0
        return iErr

    def __run_model(self, **kwargs):
        """
        Simulation execution procedures specific to a particular Chemkin reactor model
        :param kwargs: arguments from the run command
        :return: Error code (integer scalar)
        """
        iErr = 0
        return iErr

    def run(self, **kwargs):
        """
        Generic Chemkin run reactor model method
        :param kwargs: arguments from the run command
        :return: Error code (integer scalar)
        """
        logger.debug("Running " + str(self.__class__.__name__) + " " + self.label)
        for kw in kwargs:
            logger.debug("Reactor model argument " + kw + " = " + str(kwargs[kw]))
        # output initialization
        logger.debug("Clearing output")
        self.output.clear()
        # keyword processing
        logger.debug("Processing keywords")
        retVal = (
            self.__process_keywords()
        )  # each reactor model subclass to perform its own keyword processing
        logger.debug("Processing keywords complete")
        # run reactor model
        logger.debug("Running model")
        retVal = self.__run_model(**kwargs)
        logger.debug("Running model complete, status = " + str(retVal))

        return retVal

    def setsolutionspeciesfracmode(self, mode="mass"):
        """
        Set the type of species fractions in the solution
        :param mode: species fraction type = 'mass' (default) or 'mole' (string)
        :return: None
        """
        if mode.lower() in ["mole", "mass"]:
            self._speciesmode = mode.lower()
        else:
            # wrong mode value
            print(
                Color.PURPLE
                + "** species fraction mode not found, use 'mass' or 'mole'",
                end=Color.END,
            )

    def getrawsolutionstatus(self):
        """
        Get the status of the post-process
        :return: True = raw solution is ready, False = raw solution is yet to be processed (boolean scalar)
        """
        status = False
        if self._numbsolutionpoints > 0:
            status = True
        return status

    def getmixturesolutionstatus(self):
        """
        Get the status of the post-process
        :return: True = solution mixtures is ready, False = solution mixtures are yet to be processed (boolean scalar)
        """
        status = False
        if len(self._solution_mixturearray) > 0:
            status = True
        return status

    def getsolutionsize(self):
        """
        Get the number of reactors and the number of solution points
        :return: nreactor = number of reactors, npoints = number of solution points (integer scalar, integer scalar)
        """
        pass
        return 1, self._numbsolutionpoints

    def getnumbersolutionpoints(self):
        """
        Get  the number of solution points per reactor
        :return: npoints = number of solution points (integer scalar)
        """
        return self._numbsolutionpoints

    def parsespeciessolutiondata(self, frac):
        """
        Parse the species fraction solution data that are stored in a 2D array (numb_species x numb_solution)
        :param frac: species fraction array (2D double array)
        """
        # create a temporary array to hold the solution data of one species
        y = np.zeros(self._numbsolutionpoints, dtype=np.double)
        #
        for k in range(self.numbspecies):
            y[:] = frac[k, :]
            # add to the raw solution data
            self._solution_rawarray[self._specieslist[k].rstrip()] = copy.deepcopy(y)
            y[:] = 0.0e0
        # clean up
        del y

    def processsolution(self):
        """
        Post-process solution to extract the raw solution variable data
        """
        pass

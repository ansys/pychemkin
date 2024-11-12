import ctypes
from ctypes import POINTER, c_char_p, c_double, c_int
import os
from typing import Dict, List
import webbrowser

import numpy as np

from . import chemkin_wrapper as ck_wrapper
from .color import Color

# == Chemkin module global parameters -- DO NOT MODIFY without asking Chemkin development team members
# Chemkin constants
boltzmann = 1.3806504e-16  # Boltzmann constant [ergs/K] (double scalar)
avogadro = 6.02214179e23  # Avogadro number [1/mole] (double scalar)
PI = 3.1415926535897932384626433832795e0  # PI (double scalar)
Patm = 1.01325e06  # atmospheric pressure [dynes/cm2] (double scalar)
ergsperjoule = 1.0e7  # ergs per joule [ergs/J] (double scalar)
joulespercalorie = 4.184e0  # joules per calorie [J/cal] (double scalar)
ergspercalorie = (
    joulespercalorie * ergsperjoule
)  # ergs per calorie [erg/cal] (double scalar)
ergspereV = 1.602176487e-12  # ergs per eV [erg/volt] (double scalar)
eVperK = ergspereV / boltzmann  # eV per K [volt/K] (double scalar)
RGas = boltzmann * avogadro  # universal gas constant R [ergs/mol-K] (double scalar)
RGasCal = (
    RGas * 1.0e-7 / joulespercalorie
)  # universal gas constant R [cal/mol-K] (double scalar)
# == end of global constants
_symbollength = 16  # Chemkin element/species symbol length
MAX_SPECIES_LENGTH = _symbollength + 1  # Chemkin element/species symbol length + 1
LP_c_char = ctypes.POINTER(ctypes.c_char)  # pointer to C type character array
COMPLETE = 0

_chemsetidentifiers: List = (
    []
)  # string used to identify different chemistry sets in the same project

_chemkinverbose = True  # verbose mode to turn ON/OFF the print statements that do not have the leading '**' characters
_CKdict: Dict = {}  # chemkin hints
_CKInitialized: Dict = {}  # KINetics initialization flag for every Chemistry Set
# == end of global parameters


#
# Chemkin module level methods
#
def verbose() -> bool:
    """
    return the global verbose mode indicating the status (ON/OFF) of printing statements that do not have the leading '**' characters
    :return: status of the verbose mode (logical scalar)
    """
    global _chemkinverbose
    return _chemkinverbose


def setverbose(OnOff):
    """
    set the global verbose mode to turn ON(True) or OFF(False) of printing statements that do not have the leading '**' characters
    :param OnOff: the verbose status (True or False) (logical scalar)
    :return: None
    """
    global _chemkinverbose
    _chemkinverbose = OnOff


# utilities


def whereelementinarray1D(arr, target):
    """
    Find the number of occurrence and the element index in the 1D arr array that matches the target value.
    Using numpy.argwhere might be more efficient. However, the numpy method returns a list of lists of occurrence indices
    while this might be necessary for general applications, it is an overkill for simple 1D array cases.
    :param arr: the reference 1D integer or double array (1D integer or double array)
    :param target: the target value to be matched (integer or double scalar)
    :return: number_of_occurrences, occurrence_index (integer scalar, 1D integer array)
    """
    count = 0
    # check arr array size
    arr_size = len(arr)
    if arr_size == 0:
        # nothing in arr
        return count, []
    temp_index = np.zeros(arr_size, dtype=np.int32)
    value = type(arr[0])(target)
    # find all the matching occurrences
    for m in range(arr_size):
        if arr[m] == value:
            temp_index[count] = m
            count += 1
    if count == 0:
        # target is not in arr
        where_index = []
    else:
        where_index = temp_index[:count]
    return count, where_index


def bisect(ileft, iright, x, xarray):
    """
    Use bisectional method to find the largest index in the xarray of which its value is small or equal to the target x value
    :param ileft: index of xarray that represents the current lower bound
    :param iright: index of xarray that represents the current upper bound
    :param x: target value (double scalar)
    :param xarray: a sorted array containing all x values in strictly ascending order x[i] < x[i+1] (1D double array)
    :return: itarget = the largest index in the xarray of which its value is small or equal to the target x value (integer scalar)
    """
    if (iright - ileft) > 1:
        ihalf = int((ileft + iright) / 2)
        if xarray[ihalf] > x:
            iright = ihalf
        else:
            ileft = ihalf
        itarget = bisect(ileft, iright, x, xarray)
        # print(f"lower bound = {ileft}, upper bound = {iright}, target = {itarget}")
    else:
        itarget = ileft
    return itarget


def findinterpolateparameters(x, xarray):
    """
    Find the index ileft that
       xarray[ileft] <= x <= xarray[iright] where iright = ileft + 1
    :param x: target value (double scalar)
    :param xarray: a sorted array containing all x values in strictly ascending order x[i] < x[i+1] (1D double array)
    :return: itarget = the largest index in the xarray of which its value is small or equal to the target x value (integer scalar) and
    ratio = the distance ratio  (x - xarray[ileft])/(xarray[ileft+1] - xarray[ileft])
    """
    iarraysize = len(xarray)
    if x == xarray[0]:
        # x = xarray[0]
        itarget = 0
        ratio = 0.0e0
        return itarget, ratio
    if x == xarray[iarraysize - 1]:
        # x = xarray[max]
        itarget = iarraysize - 2
        ratio = 1.0e0
        return itarget, ratio
    if (x - xarray[0]) * (x - xarray[iarraysize - 1]) > 0.0e0:
        # x value is out of bound
        print(
            Color.PURPLE
            + (
                f"** Error: the target value x={x} does not fall between {xarray[0]} and {xarray[iarraysize-1]}"
            ),
            end=Color.END,
        )
        raise ValueError
    # bisect method
    ileft = 0
    iright = iarraysize - 1
    itarget = bisect(ileft, iright, x, xarray)
    ratio = (x - xarray[itarget]) / (xarray[itarget + 1] - xarray[itarget])
    return itarget, ratio


def interpolatearray(x, xarray, yarray):
    """
    Find the value in the yarray from the interpolation parameters ileft and ratio
        y = (1-ratio)*yarray[ileft] + ratio*yarray[ileft+1]
        where ileft and ratio are determined from the target x value and the xarray
    :param x: target value (double scalar)
    :param xarray: a sorted array containing all x values in strictly ascending order x[i] < x[i+1] (1D double array)
    :param yarray: variable array (1D double array)
    :return: y the interpolated variable value (double scalar)
    """
    # find the interpolation parameters
    ileft, ratio = findinterpolateparameters(x, xarray)
    # perform the interpolation to find the y value from the yarray
    y = (1.0e0 - ratio) * yarray[ileft]
    y += ratio * yarray[ileft + 1]
    return y


def checkchemistryset(chemID):
    """
    check whether the Chemistry Set is initialized in KINetics
    :param chemID: chemistry set index associated with the Chemistry Set
    """
    global _CKInitialized
    status = _CKInitialized.get(chemID, False)
    return status


def chemistrysetnew(chemID):
    """
    Create a new Chemistry Set initialization flag and set the value to False
    :param chemID: chemistry set index associated with the Chemistry Set
    """
    global _CKInitialized
    _CKInitialized[chemID] = False


def chemistrysetinitialized(chemID):
    """
    Set the Chemistry Set Initialization flag to True
    :param chemID: chemistry set index associated with the Chemistry Set
    """
    global _CKInitialized
    _CKInitialized[chemID] = True


def checkrealgasstatus(chemID):
    """
    Check whether the real-gas cubic EOS is active
    :param chemID: chemistry set index associated with the mixture (integer scalar)
    :return: status (boolean scalar)
    """
    # initialization assuming the real-gas EOS is not ON
    status = False
    #
    chemset_index = c_int(chemID)
    mode = c_int(0)
    iErr = ck_wrapper.chemkin.KINRealGas_CheckRealGasStatus(chemset_index, mode)
    if iErr == 0:
        status = mode.value == 1
    return status


def setcurrentpressure(chemID, pressure):
    """
    Set gas mixture pressure for real-gas EOS calculations
    :param chemID: Chemistry set index (integer scalar)
    :param pressure: gas pressure [dynes/cm2] (double scalar)
    """
    # convert variables
    chemset_index = c_int(chemID)
    p = c_double(pressure)
    iErr = ck_wrapper.chemkin.KINRealGas_SetCurrentPressure(chemset_index, p)
    return iErr


def help(topic=None):
    """
    Provide assistance on finding information about Chemkin keywords
    :return: None
    """
    #
    if topic is None:
        # general information about getting help
        print(
            Color.YELLOW
            + "** For detailed information about all Chemkin keywords and reactor models, "
            + "use 'chemkin.help('manual')'."
        )
        print(
            "   For usage of the real-gas cubic EOS in mixture thermodynamic property calculation, "
            + "use 'chemkin.help('real-gas')'."
        )
        print(
            "   For mixture equilibrium calculation options, use 'chemkin.help('equilibrium')'."
        )
        print(
            Color.YELLOW
            + "   For information about a Chemkin reactor model keyword, "
            + "use 'chemkin.help('keyword')'."
        )
        print(
            "   For batch reactors ignition delay time definitions, use 'chemkin.help('ignition')'.",
            end=Color.END,
        )
    elif topic.lower() in "manual manuals":
        # information about chemkin manuals
        print(
            Color.YELLOW
            + "** chemkin.manuals will open the Chemkin manuals page of the Ansys Help portal "
            + "in a new tab of the default browser."
        )
        print("      * provide Ansys login credentials to access the manuals.")
        print("      * for Chemkin keywords: check out the Input manual.")
        print(
            "      * for Chemkin reactor models: check out the Theory manual.",
            end=Color.END,
        )
        manuals()
    elif topic.lower() in "keyword keywords":
        # information about getting information about specific reactor model keywords
        print(Color.YELLOW + "** For information about a Chemkin reactor model keyword")
        print("      use 'chemkin.keyhints('<keyword>')'.")
        print("      ex: chemkin.keyhints('HTC')")
        print("   For information about keywords related to a phrase")
        print("      use 'chemkin.phrasehints('<phrase>')'.")
        print("      ex: chemkin.phrasehints('tolerance')", end=Color.END)
    elif topic.lower() in "ignition":
        # ignition definition options
        showignitiondefinition()
    elif topic.lower() in "real gas  real-gas":
        # real-gas model usage
        showrealgasusage()
    elif topic.lower() in "equilibrium":
        # equilibrium calculation options
        showequilibriumoption()
    else:
        print(
            Color.PURPLE
            + f"** cannot find help topic: '{topic:s}', the valid topics are"
        )
        print(
            "'manual', 'keyword', 'ignition', 'real gas', and 'equilibrium'",
            end=Color.END,
        )


def showrealgasusage():
    print(
        Color.YELLOW
        + "** the real-gas cubic equation of state is available only when the mechanism contains the 'EOS_' data"
    )
    print(
        "   after the gas-phase mechanism is pre-processed, the pre-processor will indicate "
        + "if the mechanism contains the necessary real-gas data"
    )
    print("   * for real-gas eligible mechanisms,")
    print("     > using the real-gas EOS with mixtures:")
    print("       to check the current activation status of the real-gas EOS, use")
    print("           status = chemkin.checkrealgasstatus()")
    print("              status = True means the real-gas EOS is active")
    print("                     = False means the ideal gas law is active")
    print(
        "       to activate the real-gas cubic EOS in mixture property calculation, use"
    )
    print("              <mixture_object>.userealgascubicEOS()")
    print("       to select the mixing rule after the real-gas EOS is activated, use")
    print("              <mixture_object>.setrealgasmixingrule(rule)")
    print("              rule is the mixing rule option:")
    print("                   0: the Van der Waals mixing method (default)")
    print("                   1: the pseudo-critical property method")
    print(
        "       to deactivate the real-gas cubic EOS in mixture property calculation, use"
    )
    print("              <mixture_object>.useidealgaslaw()")
    print("     > using the real-gas EOS with reactor models:")
    print("       see reactor model keywords: 'RLGAS' and 'RLMIX'")
    print("              ex: chemkin.keywordhints('RLGA'", end=Color.END)


def showequilibriumoption():
    """
    Show the equilibrium calculation usage and options
    :return: None
    """
    print(Color.YELLOW + "** equilibrium calculation usage: ")
    print("      EQ_mixture = chemkin.equilibrium(INIT_mixture, opt)")
    print("      INIT_mixture is the initial mixture (object)")
    print("      EQ_mixture is the final/equilibrium mixture (object)")
    print("      opt is the equilibrium calculation option: ")
    print("           1: SPECIFIED T AND P (default)")
    print("           2: SPECIFIED T AND V")
    print("           4: SPECIFIED P AND V")
    print("           5: SPECIFIED P AND H")
    print("           7: SPECIFIED V AND U")
    print("           8: SPECIFIED V AND H")
    print("** Chapman-Jouguet detonation calculation usage:")
    print("      speed_list, CJ_mixture = chemkin.detonation(INIT_mixture)")
    print("      INIT_mixture is the initial mixture (object)")
    print("      CJ_mixture is the C-J state mixture (object)")
    print("      speed_list is a list consists of two speed values at the C-J state: ")
    print(
        "           [sound_speed, detonation_wave_speed] in cm/sec",
        end=Color.END,
    )


def showignitiondefinition():
    """
    Show the ignition definitions
    :return: None
    """
    # show ignition definition usage
    print(Color.YELLOW + "** ignition definition is assigned as a string, e.g., 'OH' ")
    print("   valid options are: ")
    print("    1: 'T_inflection'")
    print("    2: 'T_rise', <val>")
    print("    3: 'T_ignition', <val>")
    print("    4: 'Species_peak', '<target species>'", end=Color.END)


def manuals():
    """
    Open the Chemkin manuals page of the Ansys Help portal
    :return: None
    """
    # Chemkin manual page
    chemkin_manual_url = (
        "https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/prod_page.html?"
        + "pn=Chemkin&pid=ChemkinPro&lang=en"
    )
    # open the web page on a new tab of the default web browser
    webbrowser.open_new_tab(chemkin_manual_url)


def createmixturerecipefromfractions(chemistryset, frac):
    """
    Build a PyChemkin mixture recipe/formula from a species fraction array (i.e., mixture mole/mass composition).
    This mixture recipe can then be used to create the corresponding Mixture object.
    :param chemistryset: the Chemistry object will be used to create the mixture (Chemistry object)
    :param frac: mole or mass fractions of the mixture (1D double array)
    :return: recipe_size, recipe (integer scalar, list of tuples [(species_symbol, fraction), ... ]
    """
    # initialization
    count = 0
    recipe = []
    # check Chemistry object
    if not isinstance(chemistryset, Chemistry):
        print(
            Color.RED + "** the first argument must be a Chemistry object",
            end=Color.END,
        )
        # print('\033[39m')
        return count, recipe
    # check array size
    numb_species = chemistryset.KK
    if len(frac) != numb_species:
        print(
            Color.RED
            + "** the size of the fractions array does not match the number of species in the chemistry set."
        )
        print(f"the fraction array size should be {numb_species}", end=Color.END)
        return count, recipe
    # build the recipe from frac array
    for k in range(numb_species):
        if frac[k] > 0.0e0:
            speciessymbol = chemistryset.KSymbol[k]
            recipe.append((speciessymbol, frac[k]))
            count += 1
    return count, recipe


# stoichiometric
#
def _nonzeroelementinarray1D(arr, threshold=0):
    """
    Find the number of occurrence and the indices of the non-zero (> 0) element in the array arr.
    Using numpy.nonzero might be more efficient. However, the numpy method returns a list of lists of occurrence indices
    while this might be necessary for general applications, it is an overkill for simple 1D array cases.
    :param arr: the reference array with non-negative integer or double (1D integer or double array)
    :param threshold: the reference value (integer or double scalar) by default the value = 0
    :return: number_of_occurrences, occurrence_index (integer scalar, 1D integer array)
    """
    # find the number of non-zero counts
    nonzero_count = np.count_nonzero(arr)
    if nonzero_count == 0:
        return nonzero_count, []
    nonzero_index = np.zeros(nonzero_count, dtype=np.int32)
    thrd = type(arr[0])(threshold)
    j = 0
    # find all non-zero occurrences
    for m in range(len(arr)):
        if arr[m] > thrd:
            nonzero_index[j] = m
            j += 1
    return nonzero_count, nonzero_index


def calculatestoichiometrics(chemistryset, fuel_molefrac, oxid_molefrac, prod_index):
    """
    calculate the stoichiometric coefficients of the complete combustion reaction of the given fuel and oxidizer mixtures.
    Consider the complete combustion of the fuel + oxidizer mixture:
    (fuel species) + alpha*(oxidizer species) <=> nu(1)*prod(1) + ... + nu(numb_prod)*prod(numb_prod)
    The number of unknowns is equal to the number of elements that make of all the fuel and oxidizer species. And the
    number of product species must be one less than the number of unknowns.
    The unknowns are:
    alpha is the stoichiometric coefficient multiplier of the oxidizer species
    nu(1), ... nu(numb_prod) are the stoichiometric coefficients of the product species
    The conservation of elements yields a set of linear algebraic equations:
    A x = b
    in which x = [ -alpha | nu(1), ...., nu(numb_prod) ]  (a vector of size numb_elem ) can be obtained.
    :param chemistryset: the Chemistry object used to create the fuel and the oxidizer mixtures (Chemistry object)
    :param fuel_molefrac: mole fractions of the fuel mixture (1D double array)
    :param oxid_molefrac: mole fractions of the oxidizer mixture (1D double array)
    :param prod_index: the species indices of the complete combustion products (1D integer array)
    :return: oxidizer_coefficient_multiplier, stoichiometric_coefficients_of_products (double scalar, 1D double array)
    """
    # check the Chemistry object
    iErr = 0
    if not isinstance(chemistryset, Chemistry):
        print(
            Color.RED + "** the first argument must be a Chemistry object",
            end=Color.END,
        )
        return 0.0, np.zeros_like(prod_index)
    # get the number of elements and the number of gas species from the chemistry set
    numb_elem = chemistryset.MM
    numb_species = chemistryset.KK
    # find fuel species array size
    kfuel = len(fuel_molefrac)
    # find oxidizer array size
    koxid = len(oxid_molefrac)
    # check fuel composition array
    if numb_species != kfuel:
        print(
            Color.PURPLE + f"** the fuel species array size must be {numb_species:d}",
            end=Color.END,
        )
        return 0.0e0, np.zeros_like(prod_index)
    # check oxidizer composition array
    if numb_species != koxid:
        print(
            Color.PURPLE
            + f"** the oxidizer species array size must be {numb_species:d}",
            end=Color.END,
        )
        return 0.0, np.zeros_like(prod_index)
    # find number of product species
    numb_prod = len(prod_index)
    # find fuel species index and count
    numb_fuel, fuel_index = _nonzeroelementinarray1D(fuel_molefrac)
    # find oxidizer species index and count
    numb_oxid, oxid_index = _nonzeroelementinarray1D(oxid_molefrac)
    # the same species cannot be fuel and oxidizer at the same time
    for i in oxid_index:
        j, j_index = whereelementinarray1D(fuel_index, i)
        if j != 0:
            print(str(j) + "   " + str(j_index))
            print(
                Color.PURPLE
                + f"** species {chemistryset.KSymbol[i]} cannot be in both the fuel and the oxidizer mixtures",
                end=Color.END,
            )
            iErr += 1
    if iErr > 0:
        return 0.0, np.zeros_like(prod_index)
    # find the actual number of elements in fuel and oxidizer
    elem_tally = np.zeros(numb_elem, dtype=np.int32)
    # elements in the fuel species
    for k in fuel_index:
        for m in range(numb_elem):
            elem_count = chemistryset.SpeciesComposition(m, k)
            if elem_count > 0:
                elem_tally[m] += elem_count
    # elements in the oxidizer species
    for k in oxid_index:
        for m in range(numb_elem):
            elem_count = chemistryset.SpeciesComposition(m, k)
            if elem_count > 0:
                elem_tally[m] += elem_count
    numb_coreelem, coreelem_index = _nonzeroelementinarray1D(elem_tally)
    # check the number of product species
    if numb_prod != (numb_coreelem - 1):
        print(
            Color.PURPLE
            + f"** the number of product species must be {numb_coreelem - 1:3}",
            end=Color.END,
        )
        return 0.0, np.zeros_like(prod_index)
    else:
        # check product elements
        # find elements in product species
        elem_prod = np.zeros(numb_elem, dtype=np.int32)
        for k in prod_index:
            for m in range(numb_elem):
                elem_count = chemistryset.SpeciesComposition(m, k)
                if elem_count > 0:
                    elem_prod[m] += elem_count
        numb_prodelem, prodelem_index = _nonzeroelementinarray1D(elem_prod)
        # check elements in the products and in the fuel and oxidzer mixtures
        elname = ""
        if numb_prodelem == numb_coreelem:
            for m in prodelem_index:
                if m in coreelem_index:
                    pass
                else:
                    elname = chemistryset.elementsymbols[m]
                    print(
                        Color.PURPLE
                        + f"** element {elname:s} in products is not in fuel or oxidizer mixtures",
                        end=Color.END,
                    )
                    return 0.0, np.zeros_like(prod_index)
        else:
            print(
                Color.PURPLE
                + "** the number of product elements must be the same as the number of elements in fuel and oxidizer"
            )
            print(f"   the number of elements in products: {numb_prodelem}")
            print(
                f"   the number of elements in the fuel and the oxidizer: {numb_coreelem}",
                end=Color.END,
            )
            return 0.0, np.zeros_like(prod_index)
    # create arrays of the linear algebraic system
    A = np.zeros((numb_coreelem, numb_coreelem), dtype=np.double)
    b = np.zeros(numb_coreelem, dtype=np.double)
    # construct the (numb_coreelem x 1) b array on the right-hand side
    #   b = [SUM_k(NCF(1,k)*fuel_molefrac(k)), ... SUM_k(NCF(numb_elem,k)*fuel_molefrac(k))]
    for m in range(numb_coreelem):
        b[m] = 0.0e0
        this_elem = coreelem_index[m]
        for k in range(numb_species):
            elem_count = chemistryset.SpeciesComposition(this_elem, k)
            b[m] += elem_count.astype(np.double) * fuel_molefrac[k]
            # first column of A[1:numb_coreelem, 1]
            A[m][0] += elem_count.astype(np.double) * oxid_molefrac[k]
    # construct the sub-matrix on the right of A[1:numb_coreelem, 2:numb_prod]
    for m in range(numb_coreelem):
        this_elem = coreelem_index[m]
        for k in range(numb_prod):
            k_prod = prod_index[k]
            A[m][k + 1] = chemistryset.SpeciesComposition(this_elem, k_prod)
    # solve the linear system: A x = b
    x = np.linalg.solve(A, b)
    alpha = -x[0]
    nu = x[1:numb_coreelem]
    return alpha, nu


class Chemistry:
    """
    define and preprocess Chemkin chemistry set
    """

    realgas_CuEOS = [
        "ideal gas",
        "Van der Waals",
        "Redlich-Kwong",
        "Soave",
        "Aungier",
        "Peng-Robinson",
    ]
    realgas_mixingrules = ["Van der Waals", "pseudocritical"]

    def __init__(self, chem=None, surf=None, therm=None, tran=None, label=None):
        # set flags
        self._index_surf = c_int(0)
        self._index_tran = c_int(0)
        self._error_code = c_int(-1)
        # initialization
        self._chemset_index = c_int(-1)  # chemistry set index
        self._num_elements = c_int(0)  # number of elements
        self._num_gas_species = c_int(0)  # number of gas species
        self._num_gas_reactions = c_int(0)  # number of gas-phase reactions
        self._num_materials = c_int(0)  # number of materials
        self._num_max_site_species = c_int(0)  # total number of surface site species
        self._num_max_bulk_species = c_int(0)  # total number of bulk species
        self._num_max_phases = c_int(0)  # number of phases
        self._num_max_surf_reactions = c_int(0)  # total number of surface reactions
        self._gas_species = {}  # gas species symbols dictionary
        self._elements = {}  # element symbols dictionary
        self._EOS = c_int(0)  # equation of state (default 0 = ideal gas)
        self.userealgas = False  # use ideal gas law by default
        self._AWTdone = 0
        self._WTdone = 0
        self._NCFdone = 0
        # fake initialization
        self._AWT = np.zeros(1, dtype=np.double)
        self._WT = np.zeros(1, dtype=np.double)
        self._KSYMdone = 0
        self.KSymbol = []
        self.ESymbol = []
        self.label = " "
        if label is not None:
            self.label = label
        # constants
        # default input file names
        self._gas_file = ctypes.c_char_p(b"chem.inp")
        self._surf_file = ctypes.c_char_p(b"surf.inp")
        self._therm_file = ctypes.c_char_p(b"therm.dat")
        self._tran_file = ctypes.c_char_p(b"tran.dat")
        # default linking file names
        self._gas_link = ctypes.c_char_p(b"chem.asc")
        self._surf_link = ctypes.c_char_p(b"surf.asc")
        self._tran_link = ctypes.c_char_p(b"tran.asc")
        # summary file for the preprocessing step
        self._summary_out = ctypes.c_char_p(b"Summary.out")
        # set the mechanism input files names if given
        if chem or surf:
            self.set_file_names(chem, surf, therm, tran)
        # check surface mechanism
        if os.path.isfile(self._surf_file.value):
            self._index_surf = ctypes.c_int(1)

        # check transport data file
        if os.path.isfile(self._tran_file.value):
            self._index_tran = ctypes.c_int(1)

    @property
    def chemfile(self):
        """
        Get gas-phase mechanism file name of this chemistry set
        """
        return self._gas_file.value.decode()

    @chemfile.setter
    def chemfile(self, filename):
        """
        Assign the gas-phase mechanism filename
        :param filename: name of the gas-phase mechanism file with the full path
        """
        self._gas_file = ctypes.c_char_p(filename.encode())

    @property
    def thermfile(self):
        """
        Get thermodynamic data filename of this chemistry set
        """
        return self._therm_file.value.decode()

    @thermfile.setter
    def thermfile(self, filename):
        """
        Assign the thermodynamic data filename
        :param filename: name of the thermodynamic data file with the full path
        """
        self._therm_file = ctypes.c_char_p(filename.encode())

    @property
    def tranfile(self):
        """
        Get transport data filename of this chemistry set
        """
        return self._tran_file.value.decode()

    @tranfile.setter
    def tranfile(self, filename):
        """
        Assign the transport data filename
        :param filename: name of the transport data file with the full path
        """
        self._tran_file = ctypes.c_char_p(filename.encode())
        if os.path.isfile(self._tran_file.value):
            self._index_tran = ctypes.c_int(1)
        else:
            self._index_tran = c_int(0)
            raise OSError(
                f"transport data file {self._tran_file.value.decode():s} not found"
            )

    @property
    def summaryfile(self):
        """
        Get the name of the summary file from the preprocessor
        """
        return self._summary_out.value.decode()

    def preprocesstransportdata(self):
        if self._index_tran.value == 0:
            # send a warning message
            print(
                Color.PURPLE
                + '** make sure the gas mechanism contains the "TRANSPORT ALL" block',
                end=Color.END,
            )
            self._index_tran = ctypes.c_int(1)
        else:
            # send the confirmation message
            if self._index_tran.value == 1:
                print(
                    Color.YELLOW
                    + f'** transport data in file: "{self._tran_file.value.decode()}" will be processed',
                    end=Color.END,
                )
            else:
                print(
                    Color.YELLOW
                    + f'** transport data in file: "{self._gas_file.value.decode()}" will be processed',
                    end=Color.END,
                )

    @property
    def surffile(self):
        """
        Get surface mechanism filename of this chemistry set
        """
        return self._surf_file.value.decode()

    @surffile.setter
    def surffile(self, filename):
        """
        Assign the surface mechanism filename
        :param filename: name of the surface mechanism file with the full path
        """
        self._surf_file = ctypes.c_char_p(filename.encode())
        if os.path.isfile(self._surf_file.value):
            self._index_surf = ctypes.c_int(1)
        else:
            self._index_surf = c_int(0)
            raise OSError(
                f"surface mechanism file {self._surf_file.value.decode():s} not found"
            )

    def set_file_names(self, chem=None, surf=None, therm=None, tran=None):
        """
        Assign all input files of the chemistry set
        :param chem: name of the gas mechanism file with the full path
        :param surf: name of the surface mechanism file with the full path
        :param therm: name of the thermodynamic data file with the full path
        :param tran: name of the transport data file with the full path
        :return: None
        """
        self._chemset_index = c_int(-1)
        if chem:
            self._gas_file = ctypes.c_char_p(chem.encode())
        if surf:
            self._surf_file = ctypes.c_char_p(surf.encode())
            if os.path.isfile(self._surf_file.value):
                self._index_surf = ctypes.c_int(1)
            else:
                self._index_surf = ctypes.c_int(0)
                raise OSError(
                    f"surface mechanism file {self._surf_file.value.decode():s} not found"
                )
        else:
            self._index_surf = c_int(0)
        if therm:
            self._therm_file = ctypes.c_char_p(therm.encode())
        if tran:
            self._tran_file = ctypes.c_char_p(tran.encode())
            if os.path.isfile(self._tran_file.value):
                self._index_tran = ctypes.c_int(1)
            else:
                self._index_tran = c_int(0)
                raise OSError(
                    f"transport data file {self._tran_file.value.decode():s} not found"
                )
        else:
            self._index_tran = c_int(0)

    def preprocess(self):
        """
        Run Chemkin preprocessor
        :return: Error code (integer scalar)
        """
        # check minimum set of required files
        if not os.path.isfile(self._gas_file.value):
            print(
                Color.RED
                + f"gas mechanism file {self._gas_file.value.decode():s} not found",
                end=Color.END,
            )
            raise OSError
        if not os.path.isfile(self._therm_file.value):
            print(Color.YELLOW + "** thermodynamic data file not found/specified")
            print(
                f"   make sure the mechanism file {self._gas_file.value.decode()} contains the 'THERM ALL' keyword",
                end=Color.END,
            )

        # verify chemistry set
        # create a new identifier for this chemistry set
        name = self._gas_file.value + self._therm_file.value
        if self._index_tran.value == 1:
            name = name + self._tran_file.value
        if self._index_surf.value == 1:
            name = name + self._surf_file.value
        identifier = name.decode()
        # check if this chemistry set is already processed by this project
        if identifier in _chemsetidentifiers:
            # existing chemistry set
            print(Color.YELLOW + "** chemistry set is already processed")
            myindex = _chemsetidentifiers.index(identifier)
            print(f"** the chemistry set index is: {myindex:d}", end=Color.END)
        else:
            # new chemistry set
            # add the identifier to the chemistry identifiers list
            _chemsetidentifiers.append(identifier)
            # modify linking file names
            myindex = len(_chemsetidentifiers)

        if myindex > 1:
            myfilename = "chem_" + str(myindex - 1) + ".asc"
            self._gas_link = ctypes.c_char_p(myfilename.encode())
            myfilename = "surf_" + str(myindex - 1) + ".asc"
            self._surf_link = ctypes.c_char_p(myfilename.encode())
            myfilename = "tran_" + str(myindex - 1) + ".asc"
            self._tran_link = ctypes.c_char_p(myfilename.encode())
            # modify the summary file for the preprocessing step
            myfilename = "Summary_" + str(myindex - 1) + ".out"
            self._summary_out = ctypes.c_char_p(myfilename.encode())

        # run preprocessor
        self._error_code = ck_wrapper.chemkin.KINPreProcess(
            self._index_surf,
            self._index_tran,
            self._gas_file,
            self._surf_file,
            self._therm_file,
            self._tran_file,
            self._gas_link,
            self._surf_link,
            self._tran_link,
            self._summary_out,
            self._chemset_index,
        )
        if self._error_code == 0:
            iErr = ck_wrapper.chemkin.KINGetChemistrySizes(
                self._chemset_index,
                self._num_elements,
                self._num_gas_species,
                self._num_gas_reactions,
                self._num_materials,
                self._num_max_site_species,
                self._num_max_bulk_species,
                self._num_max_phases,
                self._num_max_surf_reactions,
            )

            if iErr != 0:
                # failed to find mechanism sizes
                print(
                    Color.RED + "** failed to find chemistry parameters",
                    end=Color.END,
                )
                return iErr

            # get species symbols in a dictionary
            # LP_c_char = ctypes.POINTER(ctypes.c_char)
            buff = (LP_c_char * self._num_gas_species.value)()

            for i in range(0, self._num_gas_species.value):
                buff[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp = ctypes.cast(buff, POINTER(LP_c_char))
            iErr = ck_wrapper.chemkin.KINGetGasSpeciesNames(self._chemset_index, pp)
            if iErr != 0:
                # failed to get species symbols
                print(Color.RED + "** failed to get species symbols", end=Color.END)
                return iErr

            for index in range(0, len(buff)):
                strVal = ctypes.cast(buff[index], c_char_p).value
                self._gas_species[strVal] = index

            self._KSYMdone = 1
            del buff
            # get species molecular masses
            self._WT = self.WT
            # get atomic masses
            self._AWT = self.AWT
            # check real-gas model
            self.verifyrealgasmodel()
            # create a new KINetics initialization flag for this Chemistry Set
            chemistrysetnew(self._chemset_index.value)
            # save the chemkin work spaces for later use (by using active())
            self.save()
        else:
            # fail to preprocess the chemistry files
            print(
                Color.RED
                + f"** failed to preprocess the chemistry set, error code = {self._error_code}"
            )
            print(
                f"** the chemistry set index is: {self._chemset_index.value:d}",
                end=Color.END,
            )

        return self._error_code

    def verifyrealgasmodel(self):
        """
        Verify the availability of real-gas data in the mechanism
        """
        EOSModel = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
        try:
            # check if the mechanism contains the real-gas EOS data
            iErr = ck_wrapper.chemkin.KINRealGas_GetEOSMode(
                self._chemset_index, self._EOS, EOSModel
            )
            if iErr == 0:
                if self._EOS.value > 0 and self._EOS.value <= 5:
                    print(
                        Color.YELLOW
                        + f"** real-gas cubic EOS '{EOSModel.value.decode()}' is available",
                        end=Color.END,
                    )
                    del EOSModel
                    return

            self._EOS = c_int(0)
        except OSError:
            # mechanism contains no real-gas data
            self._EOS = c_int(0)
            if verbose():
                print(
                    Color.PURPLE + "** error accessing the real gas information",
                    end=Color.END,
                )

        del EOSModel
        if verbose():
            print(
                Color.YELLOW + "** mechanism is for ideal gas law only",
                end=Color.END,
            )

    @property
    def speciessymbols(self):
        """
        Get list of gas species symbols (string list)
        """
        if self._KSYMdone == 0:
            # recycle existing data
            buff = (LP_c_char * self._num_gas_species.value)()
            for i in range(0, self._num_gas_species.value):
                buff[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp = ctypes.cast(buff, POINTER(LP_c_char))
            iErr = ck_wrapper.chemkin.KINGetGasSpeciesNames(self._chemset_index, pp)
            if iErr == 0:
                self._gas_species.clear()
                for index in range(0, len(buff)):
                    strVal = ctypes.cast(buff[index], c_char_p).value
                    self._gas_species[strVal] = index
                    self._KSYMdone == 1
            else:
                # failed to get species symbols
                print(Color.RED + "** failed to get species symbols", end=Color.END)
                return [""]
            del buff

        # convert string type
        mylist = list(self._gas_species.keys())
        self.KSymbol.clear()
        for s in mylist:
            self.KSymbol.append(s.decode())
        del mylist
        return self.KSymbol

    @property
    def elementsymbols(self):
        """
        Get the list of element symbols (string list)
        """
        buff = (LP_c_char * self._num_elements.value)()
        for i in range(0, self._num_elements.value):
            buff[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
        pp = ctypes.cast(buff, POINTER(LP_c_char))
        iErr = ck_wrapper.chemkin.KINGetElementNames(self._chemset_index, pp)
        if iErr == 0:
            self._elements.clear()
            for index in range(0, len(buff)):
                strVal = ctypes.cast(buff[index], c_char_p).value
                self._elements[strVal] = index
        else:
            # failed to get element symbols
            print(Color.RED + "** failed to get element symbols", end=Color.END)
            return [""]

        del buff
        # convert string type
        mylist = list(self._elements.keys())
        self.ESymbol.clear()
        for s in mylist:
            self.ESymbol.append(s.decode())
        del mylist
        return self.ESymbol

    def getspecindex(self, specname=None):
        """
        Get index of the gas species (integer scalar)
        """
        nn = ctypes.c_char_p(specname.encode())
        specindex = self._gas_species.get(nn.value, -1)
        if specindex <= 0:
            print(
                Color.RED + f"** species symbol not found: {specname:s}",
                end=Color.END,
            )
            exit()
        return specindex

    @property
    def chemID(self):
        """
        Get chemistry set index (integer scalar)
        """
        if self._chemset_index.value >= 0:
            return self._chemset_index.value
        else:
            return -1

    @property
    def surfchem(self):
        """
        Get surface chemistry status
        0 = this chemistry set does NOT include a surface chemistry
        1 = this chemistry set includes a  surface chemistry
        """
        return self._index_surf.value

    @property
    def KK(self):
        """
        Get number of gas species (integer scalar)
        """
        return self._num_gas_species.value

    @property
    def MM(self):
        """
        Get number of elements in the chemistry set (integer scalar)
        """
        return self._num_elements.value

    @property
    def IIGas(self):
        """
        Get number of gas-phase reactions (integer scalar)
        """
        return self._num_gas_reactions.value

    @property
    def AWT(self):
        """
        compute atomic masses [gm/mol] (1D double array)
        """
        if self._AWTdone == 1:
            return self._AWT
        if self._chemset_index.value < 0:
            raise "Please preprocess the chemistry set"
        del self._AWT  # clear the "original" definition in __init__
        self._AWT = np.zeros(self._num_elements.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetAtomicWeights(self._chemset_index, self._AWT)
        if iErr == 0:
            self._AWTdone = 1
        else:
            # failed to find atomic masses
            print(Color.RED + "** failed to get atomic masses", end=Color.END)
        return self._AWT

    @property
    def WT(self):
        """
        compute gas species molecular masses [gm/mol] (1D double array)
        """
        if self._WTdone == 1:
            return self._WT
        if self._chemset_index.value < 0:
            raise "Please preprocess the chemistry set"
        del self._WT  # clear the "original" definition in __init__
        self._WT = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasMolecularWeights(
            self._chemset_index, self._WT
        )
        if iErr == 0:
            self._WTdone = 1
        else:
            # failed to find molecular masses
            print(Color.RED + "** failed to get molecular masses", end=Color.END)
        return self._WT

    def SpeciesCp(self, temp=0.0, pres=None):
        """
        Get species specific heat capacity at constant pressure
        :param temp: Temperature [K]
        :param pres: Pressure [dynes/cm2]
        :return: species specific heat capacities at constant pressure [ergs/mol-K] (1D double array)
        """
        if self._chemset_index.value < 0:
            raise "Please preprocess the chemistry set"
        if temp <= 1.0e0:
            raise ValueError("Temperature value too low")
        # check real-gas
        if checkrealgasstatus(self.chemID):
            if pres is None:
                # pressure is not assigned
                print(
                    Color.PURPLE
                    + "** must provide pressure to evaluate real-gas species properties"
                )
                print(
                    "   Usage: self.SpeciesCp(temperature, pressure)",
                    end=Color.END,
                )
                return [0.0e0]
            else:
                # set current pressure for the real-gas
                setcurrentpressure(self.chemID, pres)
        #
        TT = c_double(temp)
        Cp = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpecificHeat(self._chemset_index, TT, Cp)
        if iErr == 0:
            # convert [ergs/g-K] to [ergs/mol-K]
            # for k in range(len(Cp)):
            #    Cp[k] = Cp[k] * self._WT[k]
            Cp *= self._WT
        else:
            # failed to compute specific heats
            print(Color.RED + "** failed to compute specific heats", end=Color.END)

        return Cp

    def SpeciesCv(self, temp=0.0):
        """
        Get species specific heat capacity at constant volume (ideal gas only)
        :param temp: Temperature [K]
        :return: species specific heat capacities at constant volume [ergs/mol-K] (1D double array)
        """
        Cv = self.SpeciesCp(temp)
        R = RGas
        # for k in range(len(Cp)):
        #    Cv[k] = Cp[k] - R
        Cv -= R

        return Cv

    def SpeciesH(self, temp=0.0):
        """
        Get species enthalpy
        :param temp: Temperature [K]
        :return: species enthalpy [ergs/mol] (1D double array)
        """
        if self._chemset_index.value < 0:
            raise "Please preprocess the chemistry set"
        if temp <= 1.0e0:
            raise ValueError("Temperature value too low")
        TT = c_double(temp)
        H = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpeciesEnthalpy(self._chemset_index, TT, H)
        if iErr == 0:
            # convert [ergs/gm] to [ergs/mol]
            # for k in range(len(H)):
            #    H[k] = H[k], * self._WT[k]
            H *= self._WT
        else:
            # failed to compute enthalpies
            print(
                Color.RED + "** failed to compute species enthalpies",
                end=Color.END,
            )

        return H

    def SpeciesU(self, temp=0.0):
        """
        Get species internal energy
        :param temp: Temperature [K]
        :return: species internal energy [ergs/mol] (1D double array)
        """
        if self._chemset_index.value < 0:
            raise "Please preprocess the chemistry set"
        if temp <= 1.0e0:
            raise ValueError("Temperature value too low")
        TT = c_double(temp)
        U = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpeciesInternalEnergy(
            self._chemset_index, TT, U
        )
        if iErr == 0:
            # convert [ergs/gm] to [ergs/mol]
            # for k in range(len(U)):
            #    U[k] = U[k], * self._WT[k]
            U *= self._WT
        else:
            # failed to compute internal energies
            print(
                Color.RED + "** failed to compute species internal energies",
                end=Color.END,
            )

        return U

    def SpeciesVisc(self, temp=0.0):
        """
        Get species viscosity
        :param temp: Temperature [K]
        :return: species viscosity [gm/cm-sec] (1D double array)
        """
        if self._chemset_index.value < 0:
            raise "Please preprocess the chemistry set"
        if temp <= 1.0e0:
            raise ValueError("Temperature value too low")
        TT = c_double(temp)
        visc = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetViscosity(self._chemset_index, TT, visc)
        if iErr != 0:
            # failed to compute viscosity
            print(
                Color.RED + "** failed to compute species viscosity",
                end=Color.END,
            )

        return visc

    def SpeciesCond(self, temp=0.0):
        """
        Get species conductivity
        :param temp: Temperature [K]
        :return: species conductivity [ergs/cm-K-sec] (1D double array)
        """
        if self._chemset_index.value < 0:
            raise "Please preprocess the chemistry set"
        if temp <= 1.0e0:
            raise ValueError("Temperature value too low")
        TT = c_double(temp)
        cond = np.zeros(self._num_gas_species.value, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetConductivity(self._chemset_index, TT, cond)
        if iErr != 0:
            # failed to compute conductivities
            print(
                Color.RED + "** failed to compute species conductivity",
                end=Color.END,
            )

        return cond

    def SpeciesDiffusionCoeffs(self, press=0.0, temp=0.0):
        """
        Get species diffusion coefficients
        :param press: Pressure [dynes/cm2]
        :param temp: Temperature [K]
        :return: species diffusion coefficients [cm2/sec] (2D double array)
        """
        if self._chemset_index.value < 0:
            raise "Please preprocess the chemistry set"
        if temp <= 1.0e0:
            raise ValueError("Temperature value too low")
        if press <= 1.0e0:
            raise ValueError("Pressure value too low")
        PP = c_double(press)
        TT = c_double(temp)
        dim = (self._num_gas_species.value, self._num_gas_species.value)
        diffusioncoeffs = np.zeros(dim, dtype=np.double, order="F")
        iErr = ck_wrapper.chemkin.KINGetDiffusionCoeffs(
            self._chemset_index, PP, TT, diffusioncoeffs
        )
        if iErr != 0:
            # failed to compute diffusion coefficients
            print(
                Color.RED + "** failed to compute species diffusion coefficients",
                end=Color.END,
            )

        return diffusioncoeffs

    def SpeciesComposition(self, elemindex=-1, specindex=-1):
        """
        Get elemental composition of a species
        :param elemindex: index of the element
        :param specindex: index of the gas species
        :return: number of the element in the given gas species (integer scalar)
        """
        if self._NCFdone == 0:
            # initialize the NCF matrix
            dim = (self._num_elements.value, self._num_gas_species.value)
            self.elementalcomp = np.zeros(dim, dtype=np.int32, order="F")
            # load the NCF matrix
            iErr = ck_wrapper.chemkin.KINGetGasSpeciesComposition(
                self._chemset_index, self.elementalcomp
            )
            if iErr != 0:
                print(
                    Color.RED + "** failed to compute elemental compositions",
                    end=Color.END,
                )
                return 0
            else:
                self._NCFdone = 1

        # check element index
        if elemindex < 0:
            print(Color.RED + "** element not found", end=Color.END)
            raise ValueError
        elif elemindex >= self._num_elements.value:
            print(Color.PURPLE + "** element index out of bound", end=Color.END)
            return 0
        # check species index
        if specindex < 0:
            print(Color.RED + "** species not found", end=Color.END)
            raise ValueError
        elif elemindex >= self._num_gas_species.value:
            print(Color.PURPLE + "** species index out of bound", end=Color.END)
            return 0

        return self.elementalcomp[elemindex][specindex]

    @property
    def EOS(self):
        """
        Get the available real-gas EOS model (integer scalar)
        """
        return self._EOS.value

    def userealgascubicEOS(self):
        """
        Turn ON the real-gas cubic EOS to compute mixture properties if the mechanism contains necessary data
        :return: None
        """
        if self._EOS.value < 1:
            # no real gas EOS data in the mechanism
            print(
                Color.YELLOW + "** mechanism is for ideal gas law only",
                end=Color.END,
            )
            return
        # check real-gas EOS status
        iFlag = c_int(0)
        iErr = ck_wrapper.chemkin.KINRealGas_UseCubicEOS(self._chemset_index, iFlag)
        if iErr != 0:
            print(
                Color.RED + f"** Warning: method returned '{iErr}' error code",
                end=Color.END,
            )
        if iFlag.value == 0:
            print(
                Color.YELLOW
                + f"** real-gas cubic EOS model {Chemistry.realgas_CuEOS[self._EOS.value]} is turned ON",
                end=Color.END,
            )
            self.userealgas = True
        else:
            self.userealgas = False

    def useidealgaslaw(self):
        """
        Turn on the ideal gas law to compute mixture properties
        :return: None
        """
        if self._EOS.value < 1:
            # no real gas EOS data in the mechanism
            print(
                Color.YELLOW + "** mechanism is for ideal gas law only",
                end=Color.END,
            )
            self.userealgas = False
            return
        # check real-gas EOS status
        iFlag = c_int(0)
        iErr = ck_wrapper.chemkin.KINRealGas_UseIdealGasLaw(self._chemset_index, iFlag)
        if iErr != 0:
            print(
                Color.RED + f"** Warning: method returned '{iErr}' error code",
                end=Color.END,
            )
        if iFlag.value == 0:
            print(Color.YELLOW + "** ideal gas law is turned ON", end=Color.END)
            self.userealgas = False

    def getreactionparameters(self):
        """
        Get the reaction rate parameters of all gas-phase reactions
        :return: A-Factor, exponent, activation energy/temperature [mole-cm3-sec-K], [-], [K] (1D double array, 1D double array, 1D double array)
        """
        reactionsize = self.IIGas
        # pre-exponent A factor of all gas-phase reactions in the mechanism in cgs units [mole-cm3-sec-K]
        AFactor = np.zeros(shape=reactionsize, dtype=np.double)
        # temperature exponent of all reactions [-]
        TBeta = np.zeros_like(AFactor, dtype=np.double)
        # activation energy/temperature of all reactions [K]
        AEnergy = np.zeros_like(AFactor, dtype=np.double)
        # get the reaction parameters
        iErr = ck_wrapper.chemkin.KINGetReactionRateParameters(
            self._chemset_index, AFactor, TBeta, AEnergy
        )
        if iErr != 0:
            AFactor[:] = 0.0e0
            TBeta[:] = 0.0e0
            AEnergy[:] = 0.0e0
        return AFactor, TBeta, AEnergy

    def setreactionAFactor(self, reaction_index, AFactor):
        # check inputs
        if reaction_index > self.IIGas or reaction_index < 1:
            print(
                Color.PURPLE
                + f"** reaction index is out of bound, range = [1 ~ {self.IIGas}]",
                end=Color.END,
            )
            return 1
        if AFactor < 0.0e0:
            print(Color.PURPLE + "** A-factor must >= 0", end=Color.END)
            return 2
        # convert the reaction parameters
        ireac = c_int(-reaction_index)  # negative index to "put" A-factor value
        iErr = ck_wrapper.chemkin.KINSetAFactorForAReaction(
            self._chemset_index, ireac, AFactor
        )
        if iErr != 0:
            print(
                Color.PURPLE + f"** error encountered, code = {iErr}",
                end=Color.END,
            )
            return iErr
        return 0

    def getreactionAFactor(self, reaction_index):
        # initialization
        AFactor = 0.0e0
        # check inputs
        if reaction_index > self.IIGas or reaction_index < 1:
            print(
                Color.PURPLE
                + f"** reaction index is out of bound, range = [1 ~ {self.IIGas}]",
                end=Color.END,
            )
            return AFactor
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        # get the A-factor value
        iErr = ck_wrapper.chemkin.KINSetAFactorForAReaction(
            self._chemset_index, ireac, AFactor
        )
        if iErr != 0:
            print(
                Color.PURPLE + f"** error encountered, code = {iErr}",
                end=Color.END,
            )
            return 0.0e0
        return AFactor

    def getgasreactionstring(self, reaction_index):
        """
        Get the reaction string of the gas-phase reaction specified by the reaction index.
        :param reaction_index: (base-1) gas-phase reaction index (integer scalar)
        :return: reaction string of the reaction_index-th reaction (string)
        """
        # initialization
        reactionstring = ""
        if reaction_index > self._num_gas_reactions.value:
            print(
                Color.PURPLE
                + f"** reaction index must be < {self._num_gas_reactions.value}",
                end=Color.END,
            )
            return reactionstring
        elif reaction_index <= 0:
            print(Color.PURPLE + "** reaction index must > 0", end=Color.END)
            return reactionstring
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        iStringSize = c_int(0)
        # get reaction string
        rstring = bytes(" " * 1024, "utf-8")
        iErr = ck_wrapper.chemkin.KINGetGasReactionString(
            self._chemset_index, ireac, iStringSize, rstring
        )
        if iErr != 0:
            print(
                Color.RED + f"** Warning: method returned '{iErr}' error code",
                end=Color.END,
            )
        # convert C string back to string
        # print(rstring.decode()[0:iStringSize.value])              # check
        reactionstring = rstring.decode()[0 : iStringSize.value]
        del rstring
        return reactionstring

    def save(self):
        """
        Store the work spaces of the current Chemistry Set
        if new Chemistry Set will be created later in the same project
        """
        iErr = ck_wrapper.chemkin.KINUpdateChemistrySet(self._chemset_index)
        if iErr == 0:
            print(
                Color.YELLOW + f"** Work spaces saved for Chemistry Set {self.label}",
                end=Color.END,
            )
        else:
            print(
                Color.PURPLE + "** error saving the Chemistry Set work spaces",
                end=Color.END,
            )

    def activate(self):
        """
        Switch to (re-activate) the work spaces of the current Chemistry Set
        when there are multiple Chemistry Sets in the same project
        """
        iErr = ck_wrapper.chemkin.KINSwitchChemistrySet(self._chemset_index)
        if iErr == 0:
            print(
                Color.YELLOW
                + f"** Work spaces for Chemistry Set {self.label} activated",
                end=Color.END,
            )
        else:
            print(
                Color.PURPLE + "** error reactivating the Chemistry Set work spaces",
                end=Color.END,
            )

import copy
import ctypes
from ctypes import c_double, c_int

import numpy as np

from . import chemkin_wrapper as ck_wrapper
from .chemistry import (
    Chemistry,
    calculatestoichiometrics,
    checkchemistryset,
    checkrealgasstatus,
    chemistrysetinitialized,
    setcurrentpressure,
    verbose,
    whereelementinarray1D,
)
from .color import Color


# mixture mixing
def isothermalmixing(recipe, mode, finaltemperature):
    """
    Find the resulting gas mixture properties from mixing a number of gas mixtures at the given mixture temperature
    :param recipe: mixture and its mixing ratio as a list of [(Mixture_object, fraction), ('Mixture_object', fraction), ...]
    :param mode: indicting the fractions given in the recipe are in 'mole' or 'mass' ratios
    :param finaltemperature: temperature of the resulting gas mixture after mixing [K] (double scalar)
    :return: the resulting gas mixture after mixing (a Mixture object)
    """
    # check number of mixtures
    numb_mixture = len(recipe)
    numb_species = 0
    chem_index_check = -1
    # create the final mixture object
    finalmixture = copy.deepcopy(recipe[0][0])
    if numb_mixture == 0:
        # noting there
        print(Color.RED + "** the mixing recipe is empty", end=Color.END)
        return
    # check types
    if isinstance(recipe[0][0], Mixture):
        numb_species = recipe[0][0]._KK
        # reset the compositions
        finalmixture._Yset = 0
        finalmixture._massfrac[:] = 0.0e0
        finalmixture._Xset = 0
        finalmixture._molefrac[:] = 0.0e0
        # set the chemistry set index
        chem_index_check = finalmixture.chemID
        if chem_index_check < 0:
            raise RuntimeError(
                "** Mixture object #0 is not associated to any chemistry set"
            )
    else:
        # incorrect object type
        # delete the finalmixture object
        del finalmixture
        raise RuntimeError("** first component must be a Chemkin Mixture object")
    # check given final mixture temperature
    if finaltemperature <= 10.0:
        # final mixture temperature is not provided
        # delete the finalmixture object
        del finalmixture
        raise RuntimeError("** temperature of the final mixture must be provided")

    # initialization
    mixfrac = np.zeros(numb_mixture, dtype=np.double)
    mixfrac_sum = 0.0e0
    count = 0
    for s, v in recipe:
        # check object type
        if not isinstance(s, Mixture):
            raise RuntimeError("** first component must be a Chemkin Mixture object")
        # check ratio value
        if v <= 0.0e0:
            raise ValueError("** second component must be a positive floating number")
        # check chemistry set consistency
        if s.chemID != chem_index_check:
            raise RuntimeError(
                f"** Mixture object #{count} has incorrect chemistry set setup"
            )
        # get the mixture's mean molar mass g/mol
        mwt = s.WTM
        # find the composition of the final mixture
        speciesfrac_sum = 0.0e0
        if mode.lower() == "mole":
            # molar ratios are given
            # mole ratio
            mixfrac[count] = v
            # compute the composition of the final mixture
            for k in range(numb_species):
                finalmixture._molefrac[k] += s.X[k] * v
                speciesfrac_sum += finalmixture._molefrac[k]
            if speciesfrac_sum > 0.0e0:
                finalmixture._Xset = 1
        else:
            # assume mass ratios are given by default
            # mass ratio
            mixfrac[count] = v / mwt
            # compute the composition of the final mixture
            for k in range(numb_species):
                finalmixture._molefrac[k] += s.X[k] * mixfrac[count]
                speciesfrac_sum += finalmixture._molefrac[k]
            if speciesfrac_sum > 0.0e0:
                finalmixture._Xset = 1
        count += 1

    # normalize the mixing mole ratios
    mixfrac_sum = np.sum(mixfrac)
    mixfrac /= mixfrac_sum
    # normalize the mole fractions of the final mixture
    finalmixture._molefrac /= mixfrac_sum
    # set the temperature of the final mixture (given as input)
    finalmixture.temperature = finaltemperature
    # print(f'final mixture temperature = {finalmixture.temperature:f} [K]')
    return finalmixture


def adiabaticmixing(recipe, mode):
    """
    Find the resulting gas mixture properties from mixing a number of gas mixtures with constant total enthalpy
    :param recipe: mixture and its mixing ratio as a list of [(Mixture_object, fraction), ('Mixture_object', fraction), ...]
    :param mode: indicting the fractions given in the recipe are in 'mole' or 'mass' ratios
    :return: the resulting gas mixture after mixing (a Mixture object)
    """
    # check number of mixtures
    numb_mixture = len(recipe)
    numb_species = 0
    chem_index_check = -1
    # create the final mixture object
    finalmixture = copy.deepcopy(recipe[0][0])
    if numb_mixture == 0:
        # noting there
        print(Color.RED + "** the mixing recipe is empty", end=Color.END)
        return
    # check types
    if isinstance(recipe[0][0], Mixture):
        numb_species = recipe[0][0]._KK
        # reset the compositions
        finalmixture._Yset = 0
        finalmixture._massfrac[:] = 0.0e0
        finalmixture._Xset = 0
        finalmixture._molefrac[:] = 0.0e0
        # set the chemistry set index
        chem_index_check = finalmixture.chemID
        if chem_index_check < 0:
            raise RuntimeError(
                "** Mixture object #0 is not associated to any chemistry set"
            )
    else:
        # incorrect object type
        # delete the finalmixture object
        del finalmixture
        raise RuntimeError("** first component must be a Chemkin Mixture object")

    # initialization
    mixfrac = np.zeros(numb_mixture, dtype=np.double)
    mixfrac_sum = 0.0e0
    mix_h = 0.0e0
    count = 0
    for s, v in recipe:
        # check object type
        if not isinstance(s, Mixture):
            raise RuntimeError("** first component must be a Chemkin Mixture object")
        # check ratio value
        if v <= 0.0e0:
            raise ValueError("** second component must be a positive floating number")
        # check chemistry set consistency
        if s.chemID != chem_index_check:
            raise RuntimeError(
                f"** Mixture object #{count} has incorrect chemistry set setup"
            )
        # get the mixture's mean molar mass g/mol
        mwt = s.WTM
        # find the composition of the final mixture
        speciesfrac_sum = 0.0e0
        if mode.lower() == "mole":
            # molar ratios are given
            # mole ratio
            mixfrac[count] = v
            # compute the composition of the final mixture
            for k in range(numb_species):
                finalmixture._molefrac[k] += s.X[k] * v
                speciesfrac_sum += finalmixture._molefrac[k]
            if speciesfrac_sum > 0.0e0:
                finalmixture._Xset = 1
        else:
            # assume mass ratios are given by default
            # mass ratio
            mixfrac[count] = v / mwt
            # compute the composition of the final mixture
            for k in range(numb_species):
                finalmixture._molefrac[k] += s.X[k] * mixfrac[count]
                speciesfrac_sum += finalmixture._molefrac[k]
            if speciesfrac_sum > 0.0e0:
                finalmixture._Xset = 1

        # compute the final mixture's enthalpy ergs/mol
        mix_h += s.HML() * mixfrac[count]
        count += 1

    # normalize the mixing mole ratios
    mixfrac_sum = np.sum(mixfrac)
    mixfrac /= mixfrac_sum
    # normalize the mole fractions of the final mixture
    finalmixture._molefrac /= mixfrac_sum
    # normalize the total mixture enthalpy ergs/mol (= the enthalpy of the final mixture)
    mix_h /= mixfrac_sum
    # compute temperature of the final mixture from the mixture enthalpy
    # set the guessed temperature
    t_guessed = 0.0e0
    iErr = calculatemixturetemperaturefromenthalpy(
        mixture=finalmixture, mixtureH=mix_h, guesstemperature=t_guessed
    )
    if iErr != 0:
        print(
            Color.RED + f"** Warning: method returned error code {iErr}",
            end=Color.END,
        )
    if verbose():
        print(f"final mixture temperature = {finalmixture.temperature:f} [K]")
    return finalmixture


def calculatemixturetemperaturefromenthalpy(mixture, mixtureH, guesstemperature=None):
    """
    Compute the mixture temperature from the given mixture enthalpy,
    the solved mixture temperature is stored as the temperature attribute of the given gas mixture (i.e., as mixture.temperature)
    :param mixture: a gas mixture (Mixture object)
    :param mixtureH: mixture enthalpy of the given gas mixture (double scalar)
    :param guesstemperature: a guessed value for the mixture temperature at the start of the iteration process (optional) (double scalar)
    :return: error code (integer scalar)
    """
    # check argument
    if not isinstance(mixture, Mixture):
        raise RuntimeError("** the argument must be a Chemkin Mixture object")
    # make a copy of the mixture object
    localmixture = copy.deepcopy(mixture)
    # set converge tolerance
    tolerance = 0.1  # accurate to 0.1 K
    # iteration count limit
    maxcount = 200
    count = 0
    iErr = 0
    dt = 1.0e3
    # set guessed temperature value if given
    if guesstemperature > 0.0e0:
        localmixture.temperature = guesstemperature
    # solve for the temperature by using the Newton's method
    while True:
        # function: H(T) = mixtureH
        # compute value at T = localmixture.temperature
        f = localmixture.HML() - mixtureH
        # compute slope at T = localmixture.temperature
        df = localmixture.CPBL()
        try:
            # compute correction
            dt = f / df
        except ZeroDivisionError:
            # diverge
            print(
                Color.PURPLE + "** failed to find mixture temperature: search diverges",
                end=Color.END,
            )
            iErr = 1
            break

        if abs(dt) <= tolerance:
            # search converges
            break
        # update temperature T
        localmixture.temperature -= dt
        count += 1

    if count >= maxcount:
        # not converging within count limit
        print(
            Color.PURPLE
            + f"** cannot find mixture temperature within {maxcount:d} iterations"
        )
        print(f"** final tolerance = {abs(dt):f} [K]", end=Color.END)
        iErr = 2
    # update temperature
    if iErr != 1:
        mixture.temperature = localmixture.temperature
    # print(f'** iteration count = {count:d}')
    del localmixture
    return iErr


def interpolatemixtures(mixtureleft, mixtureright, ratio):
    """
    Create a new mixture object by interpolating the two mixture objects with a specific weight ratio
    mixture_new = (1 - ratio) * mixtureleft + ratio * mixtureright
    :param mixtureleft: a Mixture object
    :param mixtureright: a Mixture object
    :param ratio: the weight parameters for interpolation, 0 <= ratio <= 1 (double scalar)
    :return: the resulting Mixture object
    """
    # check mixtures
    if not isinstance(mixtureleft, Mixture):
        print(
            Color.RED + "** the first argument must be a Mixture object",
            end=Color.END,
        )
        raise TypeError
    iErr = mixtureleft.validate()
    if iErr != 0:
        print(
            Color.YELLOW + "** the first mixture is not fully defined",
            end=Color.END,
        )
        raise
    #
    if not isinstance(mixtureright, Mixture):
        print(
            Color.RED + "** the second argument must be a Mixture object",
            end=Color.END,
        )
        raise TypeError
    iErr = mixtureright.validate()
    if iErr != 0:
        print(
            Color.YELLOW + "** the second mixture is not fully defined",
            end=Color.END,
        )
        raise
    # check ratio
    if ratio < 0.0e0 or ratio > 1.0e0:
        print(
            Color.RED + "** error: the weight parameter must be 0 <= and <= 1",
            end=Color.END,
        )
        raise ValueError
    ratiom = 1.0e0 - ratio
    # interpolate the mixture properties
    mixturenew = copy.deepcopy(mixtureleft)
    # temperature
    mixturenew.temperature = (
        ratiom * mixtureleft.temperature + ratio * mixtureright.temperature
    )
    # pressure
    mixturenew.pressure = ratiom * mixtureleft.pressure + ratio * mixtureright.pressure
    # volume
    mixturenew.volume = ratiom * mixtureleft.volume + ratio * mixtureright.volume
    # species composition
    fracleft = mixtureleft.Y
    fracright = mixtureright.Y
    frac = np.zeros(len(fracleft), dtype=np.double)
    frac = ratiom * fracleft + ratio * fracright
    mixturenew.Y = frac
    # clean up
    del frac
    #
    return mixturenew


# equilibrium
#
def calculateequilibrium(
    chemID, p, t, frac, wt, mode_in, mode_out, EQOption=None, useRealGas=None
):
    """
    Get the equilibrium mixture composition corresponding to the given initial mixture composition at the given pressure and temperature
    :param chemID: chemistry set index associated with the mixture (integer scalar)
    :param p: mixture pressure in [dynes/cm2] (double scalar)
    :param t: mixture temperature in [K] (double scalar)
    :param frac: mixture composition given by either mass or mole fractions as specified by mode (double array)
    :param wt: molar masses of the species in the mixture in [gm/mol] (double array)
    :param mode_in: flag indicates the frac array is 'mass' or 'mole' fractions
    :param mode_out: flag to indicate the returning composition is in 'mole' or 'mass' fraction
    :param EQOption: equilibrium type (integer scalar) see below
    1.  SPECIFIED T AND P
    2.  SPECIFIED T AND V
    3.  SPECIFIED T AND S
    4.  SPECIFIED P AND V
    5.  SPECIFIED P AND H
    6.  SPECIFIED P AND S
    7.  SPECIFIED V AND U
    8.  SPECIFIED V AND H
    9.  SPECIFIED V AND S
    10. CHAPMAN-JOUGUET DETONATION (H, S, V, T CONTAIN THE UNBURNED STATE, AND TE IS THE BURNED ESTIMATE)
    :param useRealGas: option to turned ON/OFF (1/0) the real-gas cubic EOS if available (integer scalar)
    :return: list of state variables at equilibrium , equilibrium composition given in fractions assigned by mode_out
    (list of double scalars, 1D double array)
    the state variable list consists of
    equilibrium pressure [dynes/cm2], equilibrium temperature [K], and, if Chapmen-Jouguet option is used,
    speed of sound at equilibrium [cm/sec], detonation wave speed [cm/sec]
    """
    # find the equilibrium composition at the mixture pressure and temperature
    # check inputs
    if chemID < 0:
        print(Color.RED + "** invalid chemistry", end=Color.END)
        return [p, t, 0.0, 0.0], [0.0e0]

    if p <= 0.0 or (p * t) <= 0.0:
        print(
            Color.RED + "** invalid pressure and/or temperature value(s)",
            end=Color.END,
        )
        return [p, t], [0.0e0]

    # number species
    kgas = len(frac)
    if kgas != len(wt):
        print(
            Color.RED
            + f"** fraction and molar mass arrays must have the same size: {kgas:d}",
            end=Color.END,
        )
        return [p, t, 0.0, 0.0], [0.0e0]

    # initialization
    # set to constant T-P option by default
    eq_option = c_int(1)
    # use ideal gas law by default
    iRealGas = c_int(0)
    xx_eq = np.zeros(kgas, dtype=np.double)
    #
    if mode_in.lower() == "mole":
        # normalize mass fractions
        iErr, x = Mixture.normalize(frac=frac)
    elif mode_in.lower() == "mass":
        # convert mass fraction to mole fraction and normalize
        x = Mixture.massfractiontomolefraction(massfrac=frac, wt=wt)
    else:
        # fraction type not given or incorrect
        print(
            Color.PURPLE + '** must specify "mole" or "mass" fractions given',
            end=Color.END,
        )
        return [p, t, 0.0, 0.0], xx_eq
    # check equilibrium calculation option
    if EQOption is None:
        pass
    else:
        eq_option = c_int(EQOption)

    # check real gas option
    if useRealGas is None:
        pass
    else:
        if useRealGas:
            iRealGas = c_int(1)
            setcurrentpressure(chemID, pressure=p)

    # convert parameters to c pointers
    _chemset_index = ctypes.c_int(chemID)
    pp = c_double(p)  # pressure scalar
    tt = c_double(t)  # temperature scalar
    xx = np.ctypeslib.as_array(x)  # mole fraction array
    #
    pp_eq = c_double(p)
    tt_eq = c_double(t)
    detonationwavespeed = c_double(0.0e0)
    soundspeed_eq = c_double(0.0e0)
    # perform gas-phase equilibrium calculationk
    if not checkchemistryset(_chemset_index.value):
        # need to initialize KINetics
        print(Color.YELLOW + "** initializing chemkin...", end=Color.END)
        iErr = ck_wrapper.chemkin.KINInitialize(_chemset_index, c_int(0))
        if iErr != 0:
            print(Color.RED + "** fail to initialize KINetics", end=Color.END)
            statevars = [p, t, 0.0, 0.0]
            return statevars, frac
        else:
            chemistrysetinitialized(_chemset_index.value)
    else:
        iErr = 0

    iErr = ck_wrapper.chemkin.KINCalculateEqGasWithOption(
        _chemset_index,
        eq_option,
        iRealGas,
        pp,
        tt,
        xx,
        pp_eq,
        tt_eq,
        soundspeed_eq,
        detonationwavespeed,
        xx_eq,
    )

    if iErr == 0:
        # process solution
        if eq_option.value == 10 and verbose():
            # CHAPMAN-JOUGUET DETONATION
            print(
                f"** detonation wave speed = {detonationwavespeed.value / 1.0e2} [m/sec]"
            )
            print(
                f"** speed of sound at final state = {soundspeed_eq.value / 1.0e2} [m/sec]"
            )

        if mode_out.lower() == "mass":
            # convert mass fraction to mole fraction and normalize
            y_eq = Mixture.molefractiontomassfraction(molefrac=xx_eq, wt=wt)
            statevars = [
                pp_eq.value,
                tt_eq.value,
                soundspeed_eq.value,
                detonationwavespeed.value,
            ]
            return statevars, y_eq
        else:
            # by default, return mass fraction
            # normalize mass fractions
            iErr, x_eq = Mixture.normalize(frac=xx_eq)
            statevars = [
                pp_eq.value,
                tt_eq.value,
                soundspeed_eq.value,
                detonationwavespeed.value,
            ]
            return statevars, x_eq

    else:
        print(Color.RED + "** fail to find the equilibrium state", end=Color.END)
        statevars = [p, t, 0.0, 0.0]
        return statevars, xx


def equilibrium(mixture, opt=None):
    """
    Find the equilibrium state mixture corresponding to the given mixture
    :param mixture: initial gas mixture (Mixture object)
    :param opt: equilibrium type (integer scalar) see below
    1.  SPECIFIED T AND P
    2.  SPECIFIED T AND V
    *3.  SPECIFIED T AND S
    4.  SPECIFIED P AND V
    5.  SPECIFIED P AND H
    *6.  SPECIFIED P AND S
    7.  SPECIFIED V AND U
    8.  SPECIFIED V AND H
    *9.  SPECIFIED V AND S
    *10. CHAPMAN-JOUGUET DETONATION (H, S, V, T CONTAIN THE UNBURNED STATE, AND TE IS THE BURNED ESTIMATE)
    * indicates the options are not available
    :return: gas mixture at the equilibrium state (Mixture object)
    """
    # check argument
    if not isinstance(mixture, Mixture):
        raise RuntimeError("** the argument must be a Chemkin Mixture object")
    # initialization a Mixture object by duplication
    EQState = copy.deepcopy(mixture)
    # reset mass/mole fractions
    EQState._Xset = 0
    EQState._molefrac[:] = 0.0e0
    EQState._Yset = 0
    EQState._massfrac[:] = 0.0e0
    if opt is None:
        # use default option
        option = 1
    else:
        # check option
        if opt in [3, 6, 9, 10]:
            print(
                Color.PURPLE + f"** equilibrium option {opt:d} is not available",
                end=Color.END,
            )
            return mixture
        else:
            option = opt

    if mixture.EOS == 0:
        userealgas = 0
    else:
        userealgas = mixture.userealgas
    # compute the equilibrium state (mass fraction for now)
    eqvars, EQState._massfrac = calculateequilibrium(
        mixture.chemID,
        p=EQState.pressure,
        t=EQState.temperature,
        frac=mixture.Y,
        wt=mixture._WT,
        mode_in="mass",
        mode_out="mass",
        EQOption=option,
        useRealGas=userealgas,
    )
    if np.sum(EQState._massfrac, dtype=np.double) > 0.0e0:
        EQState._Yset = 1
    EQState.pressure = eqvars[0]
    EQState.temperature = eqvars[1]
    return EQState


def detonation(mixture):
    """
    Find the Chapman-Jouguet state mixture and detonation wave speed corresponding to the given mixture
    :param mixture: initial gas mixture (Mixture object)
    :return: list of speeds (list of double scalars), gas mixture at the equilibrium state (Mixture object)
    the speed list consists of
    speed of sound [cm/sec], detonation wave speed [cm/sec]
    """
    # check argument
    if not isinstance(mixture, Mixture):
        raise RuntimeError("** the argument must be a Chemkin Mixture object")
    # initialization a Mixture object by duplication
    EQState = copy.deepcopy(mixture)
    # reset mass/mole fractions
    EQState._Xset = 0
    EQState._molefrac[:] = 0.0e0
    EQState._Yset = 0
    EQState._massfrac[:] = 0.0e0
    # use the C-J option
    option = 10
    if mixture.EOS == 0:
        userealgas = 0
    else:
        userealgas = mixture.userealgas
    # compute the equilibrium state (mass fraction for now)
    eqvars, EQState._massfrac = calculateequilibrium(
        mixture.chemID,
        p=EQState.pressure,
        t=EQState.temperature,
        frac=mixture.Y,
        wt=mixture._WT,
        mode_in="mass",
        mode_out="mass",
        EQOption=option,
        useRealGas=userealgas,
    )
    if np.sum(EQState._massfrac, dtype=np.double) > 0.0e0:
        EQState._Yset = 1
    EQState.pressure = eqvars[0]
    EQState.temperature = eqvars[1]
    return [eqvars[2], eqvars[3]], EQState


class Mixture:
    """
    define a mixture based on the gas species in the given chemistry set
    """

    def __init__(self, chem):
        self._temp = 0.0e0  # mixture temperature [K]
        self._press = 0.0e0  # mixture pressure [dynes/cm2]
        self._vol = 0.0e0  # mixture volume [cm3]
        # flags
        self._Tset = 0
        self._Pset = 0
        self._Xset = 0
        self._Yset = 0
        # chemistry set validation
        if not isinstance(chem, Chemistry):
            raise TypeError("arg must be a chemkin.Chemistry object")
        if chem.chemID < 0:
            raise "invalid chemistry"
        # shorthand for frequently used variables
        self._chemset_index = ctypes.c_int(chem.chemID)  # chemistry set index
        self._KK = chem.KK  # number of gas species
        self._IIgas = chem.IIGas  # number of gas-phase reactions
        self._specieslist = []
        self._specieslist = chem.speciessymbols  # gas species symbols
        self._WT = chem.WT  # gas species molar masses
        # create internal arrays: array size = number of gas species
        self._molefrac = np.zeros(
            self._KK, dtype=np.double
        )  # mixture composition given in mole fractions
        self._massfrac = np.zeros_like(
            self._molefrac
        )  # mixture composition given in mole fractions
        self._concentration = np.zeros_like(self._molefrac)  # concentrations (not used)
        self._SurfaceChem = c_int(
            chem.surfchem
        )  # flag indicating there is surface chemistry (type c_int: 0 = no, 1 = yes)
        self._EOS = c_int(chem.EOS)  # real-gas EOS model in the mechanism
        self.userealgas = chem.userealgas  # status of the real-gas EOS usage

    @property
    def chemID(self):
        """
        Get chemistry set index (integer scalar)
        """
        return self._chemset_index.value

    @property
    def pressure(self):
        """
        Get gas mixture pressure [dynes/cm2] (double scalar)
        """
        if self._Pset == 1:
            return self._press
        else:
            print(
                Color.PURPLE + "** mixture pressure is not provided",
                end=Color.END,
            )
            return 0.0e0

    @pressure.setter
    def pressure(self, p):
        """
        Set gas mixture pressure
        :param p: pressure [dynes/cm2] (double scalar)
        """
        if p <= 0.0:
            print(Color.RED + "** invalid pressure value", end=Color.END)
            return

        self._press = p
        self._Pset = 1

    @property
    def temperature(self):
        """
        Get gas mixture temperature [K] (double scalar)
        """
        if self._Tset == 1:
            return self._temp
        else:
            print(
                Color.PURPLE + "** mixture temperature is not provided",
                end=Color.END,
            )
            return 0.0e0

    @temperature.setter
    def temperature(self, t):
        """
        Set gas mixture temperature
        :param t: temperature [K] (double scalar)
        """
        if t <= 10.0:
            print(Color.RED + "** invalid temperature value", end=Color.END)
            return
        self._temp = t
        self._Tset = 1

    @property
    def volume(self):
        if self._vol > 0.0e0:
            return self._vol
        else:
            print(Color.PURPLE + "** mixture volume is not provided", end=Color.END)
            return 0.0e0

    @volume.setter
    def volume(self, vol):
        if vol <= 0.0e0:
            print(Color.RED + "** invalid volume value", end=Color.END)
            return
        self._vol = vol

    @property
    def X(self):
        """
        Get mixture mole fraction (double array)
        """
        if self._Xset == 1:
            iErr, x = Mixture.normalize(self._molefrac)
            return x
        elif self._Yset == 1:
            iErr, x = self.__YtoX()
            if iErr != 0:
                print(Color.RED + "** error encountered", end=Color.END)
            return x
        else:
            print(Color.RED + "** mixture composition not defined", end=Color.END)
            iErr, x = Mixture.normalize(self._molefrac)
            return x

    @X.setter
    def X(self, recipe):
        """
        Set mixture molar composition
        :param recipe: mixture composition a list of [('species_symbol', mole_fraction), ('species_symbol', mole_fraction), ...]
        """
        if self._Xset == 1:
            # reset the mole fraction array
            self._molefrac[:] = 0.0e0
        if isinstance(recipe[0], tuple):
            for sp, x in recipe:
                if sp in self._specieslist:
                    index = self._specieslist.index(sp)
                else:
                    print(
                        Color.RED + f"** not a gas species: {sp:s}",
                        end=Color.END,
                    )
                    exit()
                if x < 0.0:
                    print(Color.RED + "** negative mole fraction", end=Color.END)
                    exit()
                # set mole fraction
                self._molefrac[index] = x
        elif isinstance(recipe[0], (float, np.double)):
            kgas = len(recipe)
            if kgas == self._KK:
                for k in range(kgas):
                    self._molefrac[k] = max(recipe[k], 0.0e0)
            else:
                print(
                    Color.RED
                    + f"** size of the mole fraction array must equal to the number of species: {self._KK:d}",
                    end=Color.END,
                )
                return
        else:
            print(Color.RED + "** the argument must be")
            print("   (1) a list of tuple : for example, ('O2', 0.21) or")
            print(
                "   (2) a mole fraction array of size = number_of_species",
                end=Color.END,
            )
            return
        # reset mass fraction
        self._Yset = 0
        self._massfrac[:] = 0.0e0
        self._Xset = 1

    @property
    def Y(self):
        """
        Get mixture mass fraction (double array)
        """
        if self._Yset == 1:
            iErr, y = Mixture.normalize(self._massfrac)
            return y
        elif self._Xset == 1:
            iErr, y = self.__XtoY()
            if iErr != 0:
                print(Color.RED + "** error encountered", end=Color.END)
            return y
        else:
            print(Color.RED + "** mixture composition not defined", end=Color.END)
            iErr, y = Mixture.normalize(self._massfrac)
            return y

    @Y.setter
    def Y(self, recipe):
        """
        Set mixture mass composition
        :param recipe: mixture composition a list of [('species_symbol', mass_fraction), ('species_symbol', mass_fraction), ...]
        """
        if self._Yset == 1:
            # reset the mass fraction array
            self._massfrac[:] = 0.0e0
        if isinstance(recipe[0], tuple):
            for sp, y in recipe:
                if sp in self._specieslist:
                    index = self._specieslist.index(sp)
                else:
                    print(
                        Color.RED + f"** not a gas species: {sp:s}",
                        end=Color.END,
                    )
                    exit()
                if y < 0.0:
                    print(Color.RED + "** negative mass fraction", end=Color.END)
                    exit()
                # set mass fraction
                self._massfrac[index] = y
        elif isinstance(recipe[0], (float, np.double)):
            kgas = len(recipe)
            if kgas == self._KK:
                for k in range(kgas):
                    self._massfrac[k] = max(recipe[k], 0.0e0)
            else:
                print(
                    Color.RED
                    + f"** size of the mass fraction array must equal to the number of species: {self._KK:d}",
                    end=Color.END,
                )
                return
        else:
            print(Color.RED + "** the argument must be")
            print("   (1) a list of tuple : for example, ('O2', 0.23) or")
            print(
                "   (2) a mass fraction array of size = number_of_species",
                end=Color.END,
            )
            return
        # reset mole fraction
        self._Xset = 0
        self._molefrac[:] = 0.0e0
        self._Yset = 1

    @property
    def concentration(self):
        """
        Get mixture molar concentrations [mole/cm3] (1D double array)
        """
        if self._Xset == 1:
            # mole fractions are given
            # remove negative values and normalize fractions
            iErr, c = Mixture.normalize(frac=self._molefrac)
            if iErr == 0:
                # compute mean molar mass
                mwt = self.WTM
                # compute density
                den = self.RHO
                fac = den / mwt
                for k in range(self._KK):
                    c[k] *= fac
                self._concentration[:] = c[:]
            return c
        elif self._Yset == 1:
            # mass fractions are given
            # remove negative values and normalize fractions
            iErr, c = Mixture.normalize(frac=self._massfrac)
            if iErr == 0:
                # compute density
                den = self.RHO
                for k in range(self._KK):
                    c[k] = c[k] * den / self._WT[k]
                self._concentration[:] = c[:]
            return c
        else:
            print(
                Color.PURPLE + "** mixture composition is not defined",
                end=Color.END,
            )
            return [0.0]

    @property
    def EOS(self):
        """
        Get the available real-gas EOS model that is provided in the mechanism (integer scalar)
        """
        return self._EOS.value

    @staticmethod
    def normalize(frac):
        """
        Normalize the mixture composition
        :param frac: mixture composition to be normalized
        :return: (error code, normalized fraction array)  (integer scalar, double array)
        """
        # initialization
        sumx = 0.0e0
        KK = len(frac)  # number of entries
        localfrac = copy.deepcopy(frac)  # make a local copy of the frac array
        # remove negative fraction and calculate sum
        for k in range(KK):
            if localfrac[k] > 0.0:
                sumx += localfrac[k]
            else:
                localfrac[k] = 0.0e0
        # check none zero sum
        if sumx > 0.0:
            # normalization
            for k in range(KK):
                localfrac[k] = localfrac[k] / sumx
            return 0, localfrac
        else:
            # fractions summed to zero
            print(Color.PURPLE + "** fractions summed to zero", end=Color.END)
            return 1, frac

    @property
    def WT(self):
        """
        Get species molecular masses [gm/mole] (1D double array)
        """
        return self._WT

    @property
    def WTM(self):
        """
        Get mean molar mass of the gas mixture [gm/mol] (double scalar)
        """
        mwt = 0.0e0
        if self._Xset == 1:
            # mole fractions are given
            # remove negative values and normalize fractions
            iErr, x = Mixture.normalize(frac=self._molefrac)
            if iErr == 0:
                # compute mean molar mass
                for k in range(self._KK):
                    mwt += x[k] * self._WT[k]

            return mwt

        elif self._Yset == 1:
            # mass fractions are given
            # remove negative values and normalize fractions
            iErr, y = Mixture.normalize(frac=self._massfrac)
            if iErr == 0:
                # compute mean molar mass
                for k in range(self._KK):
                    mwt += y[k] / self._WT[k]

                if mwt > 0.0:
                    return 1.0e0 / mwt
                else:
                    # zero mean molar mass
                    return mwt
            else:
                # zero mean molar mass
                return mwt
        else:
            # no fractions given
            print(
                Color.PURPLE + "** mixture composition is not defined",
                end=Color.END,
            )
            return mwt

    def __XtoY(self):
        """
        Convert mole fraction to mass fraction
        :return: (error_code, mass fraction) (integer scalar, double array)
        """
        # compute mean molar mass
        mwt = self.WTM
        if mwt > 0.0e0:
            # remove negative values and normalize fractions
            iErr, y = Mixture.normalize(frac=self._molefrac)
            if iErr == 0:
                # convert mole fractions to mass fractions
                for k in range(self._KK):
                    y[k] = y[k] * self._WT[k] / mwt
                return 0, y
            else:
                return iErr, self._molefrac
        else:
            # zero mean molar mass
            return 2, self._molefrac

    def __YtoX(self):
        """
        Convert mass fraction to mole fraction
        :return: (error_code, mole fraction) (integer scalar, double array)
        """
        # compute mean molar mass
        mwt = self.WTM
        if mwt > 0.0e0:
            # remove negative values and normalize fractions
            iErr, x = Mixture.normalize(frac=self._massfrac)
            if iErr == 0:
                # convert mass fractions to mole fractions
                for k in range(self._KK):
                    x[k] = x[k] * mwt / self._WT[k]
                return 0, x
            else:
                return iErr, self._massfrac
        else:
            # zero mean molar mass
            return 2, self._massfrac

    @staticmethod
    def meanmolarmass(frac, wt, mode):
        """
        Get mean molar mass of the gas mixture
        :param frac: mixture composition in 'mass' or mole fraction as indicated by mode (double array)
        :param wt: species molar mass [gm/mol] (double array)
        :param mode: flag indicates the frac array is 'mass' or 'mole' fractions
        :return: mean molar mass [gm/mol] (double scalar)
        """
        # initialization
        mwt = 0.0e0
        # check sizes
        kgas = len(frac)
        k = len(wt)
        if k != kgas:
            # mismatch input arrays
            print(
                Color.RED
                + f"** {mode} fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return mwt
        if mode.lower() == "mole":
            # mole fractions are given
            # remove negative values and normalize fractions
            iErr, x = Mixture.normalize(frac=frac)
            if iErr == 0:
                # compute mean molar mass
                for k in range(kgas):
                    mwt += x[k] * wt[k]

            return mwt

        elif mode.lower() == "mass":
            # mass fractions are given
            # remove negative values and normalize fractions
            iErr, y = Mixture.normalize(frac=frac)
            # compute mean molar mass
            if iErr == 0:
                for k in range(kgas):
                    mwt += y[k] / wt[k]

            if mwt > 0.0:
                return 1.0e0 / mwt
            else:
                # zero mean molar mass
                return mwt
        else:
            # fractions summed to zero
            print(
                Color.PURPLE + "** mixture composition is not defined",
                end=Color.END,
            )
            return mwt

    @staticmethod
    def molefractiontomassfraction(molefrac, wt):
        """
        Convert mole fraction to mass fraction
        :param molefrac: mixture composition in mole fractions (double array)
        :param wt: species molar mass [gm/mol] (double array)
        :return: mass fraction (double array)
        """
        # check size
        kgas = len(molefrac)
        k = len(wt)
        if k != kgas:
            # mismatch input arrays
            print(
                Color.RED
                + f"** mole fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return molefrac
        # compute mean molar mass
        mwt = Mixture.meanmolarmass(frac=molefrac, wt=wt, mode="mole")
        if mwt > 0.0e0:
            # remove negative values and normalize fractions
            iErr, massfrac = Mixture.normalize(frac=molefrac)
            if iErr == 0:
                # convert mole fractions to mass fractions
                for k in range(kgas):
                    massfrac[k] = massfrac[k] * wt[k] / mwt

            return massfrac
        else:
            # zero mean molar mass
            return molefrac

    @staticmethod
    def massfractiontomolefraction(massfrac, wt):
        """
        Convert mass fraction to mole fraction
        :param massfrac: mixture composition in mass fractions (double array)
        :param wt: species molar mass [gm/mol] (double array)
        :return: mole fraction (double array)
        """
        # check size
        kgas = len(massfrac)
        k = len(wt)
        if k != kgas:
            print(
                Color.RED
                + f"** mass fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return massfrac
        # compute mean molar mass
        mwt = Mixture.meanmolarmass(frac=massfrac, wt=wt, mode="mass")
        if mwt > 0.0e0:
            # remove negative values and normalize fractions
            iErr, molefrac = Mixture.normalize(frac=massfrac)
            if iErr == 0:
                # convert mass fractions to mole fractions
                for k in range(kgas):
                    molefrac[k] = molefrac[k] * mwt / wt[k]

            return molefrac
        else:
            # zero mean molar mass
            return massfrac

    @staticmethod
    def massfractiontoconcentration(chemID, p, t, massfrac, wt):
        """
        Convert mass fractions to molar concentrations
        :param chemID: chemistry set index associated with the mixture (integer scalar)
        :param p: pressure [dynes/cm2] (double scalar)
        :param t: temperature [K] (double scalar)
        :param massfrac: mixture mass fractions (1D double array)
        :param wt: species molecular masses [gm/mole] (1D double array)
        :return: molar concentration -mole/cm3] (1D double array)
        """
        # check size
        kgas = len(massfrac)
        k = len(wt)
        if k != kgas:
            print(
                Color.RED
                + f"** mass fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return massfrac
        # compute density
        den = Mixture.density(chemID, p, t, frac=massfrac, wt=wt, mode="mass")
        if den > 0.0e0:
            # remove negative values and normalize fractions
            iErr, c = Mixture.normalize(frac=massfrac)
            if iErr == 0:
                # convert mass fractions to mole fractions
                for k in range(kgas):
                    c[k] = c[k] * den / wt[k]
            return c
        else:
            # zero mean molar mass
            return massfrac

    @staticmethod
    def molefractiontoconcentration(chemID, p, t, molefrac, wt):
        """
        Convert mole fractions to molar concentrations
        :param chemID: chemistry set index associated with the mixture (integer scalar)
        :param p: pressure [dynes/cm2] (double scalar)
        :param t: temperature [K] (double scalar)
        :param molefrac: mixture mole fractions (1D double array)
        :param wt: species molecular masses [gm/mole] (1D double array)
        :return: molar concentration -mole/cm3] (1D double array)
        """
        # check size
        kgas = len(molefrac)
        k = len(wt)
        if k != kgas:
            print(
                Color.RED
                + f"** mole fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return molefrac
        # compute mean molar mass
        mwt = Mixture.meanmolarmass(frac=molefrac, wt=wt, mode="mole")
        # compute density
        den = Mixture.density(chemID, p, t, frac=molefrac, wt=wt, mode="mole")
        if mwt * den > 0.0e0:
            # remove negative values and normalize fractions
            iErr, c = Mixture.normalize(frac=molefrac)
            if iErr == 0:
                # convert mass fractions to mole fractions
                fac = den / mwt
                for k in range(kgas):
                    c[k] *= fac
            return c
        else:
            # zero mean molar mass
            return molefrac

    def listcomposition(self, mode, option=" ", bound=0.0e0):
        """
        list the mixture composition
        :param mode: flag indicates the fractions returned are 'mass' or 'mole' fractions
        :param option: flag indicates to list 'all' species or just the species with non-zero fraction (default)
        :param bound: minimum fraction value for the species to be printed (double scalar)
        :return: None
        """
        #
        if option.lower() == "all":
            # list all species
            if mode.lower() == "mass":
                print(f"listing mixture composition in {mode} fractions")
                for k in range(self._KK):
                    print(f"{self._specieslist[k]:18} :  {self.Y[k]:e}")
            elif mode.lower() == "mole":
                print(f"listing mixture composition in {mode} fractions")
                for k in range(self._KK):
                    print(f"{self._specieslist[k]:18} :  {self.X[k]:e}")
            else:
                print(
                    Color.PURPLE + f"** unknown fraction mode: {mode}",
                    end=Color.END,
                )
        else:
            # list no-zero components
            if mode.lower() == "mass":
                print(f"listing mixture composition in {mode} fractions")
                for k in range(self._KK):
                    if self.Y[k] > np.max([bound, 0.0e0]):
                        print(f"{self._specieslist[k]:18} :  {self.Y[k]:e}")
            elif mode.lower() == "mole":
                print(f"listing mixture composition in {mode} fractions")
                for k in range(self._KK):
                    if self.X[k] > np.max([bound, 0.0e0]):
                        print(f"{self._specieslist[k]:18} :  {self.X[k]:e}")
            else:
                print(
                    Color.PURPLE + f"** unknown fraction mode: {mode}",
                    end=Color.END,
                )

    @staticmethod
    def density(chemID, p, t, frac, wt, mode):
        """
        Get mass density from the given mixture condition: pressure, temperature, and species composition
        :param chemID: chemistry set index associated with the mixture (integer scalar)
        :param p: mixture pressure in [dynes/cm2] (double scalar)
        :param t: mixture temperature in [K] (double scalar)
        :param frac: mixture composition given by either mass or mole fractions as specified by mode (double array)
        :param wt: molar masses of the species in the mixture in [gm/mol] (double array)
        :param mode: flag indicates the frac array is 'mass' or 'mole' fractions
        :return: mass density in [gm/cm3] (double scalar)
        """
        # check inputs
        if chemID < 0:
            print(Color.RED + "** invalid chemistry", end=Color.END)
            return 0.0e0

        if p <= 0.0 or (p * t) <= 0.0:
            print(
                Color.RED + "** invalid pressure and/or temperature value(s)",
                end=Color.END,
            )
            return 0.0e0

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            print(
                Color.RED
                + f"** fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return 0.0e0

        # initialization
        den_C = c_double(0.0)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.molefractiontomassfraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            iErr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            print(
                Color.PURPLE + '** must specify "mole" or "mass" fractions given',
                end=Color.END,
            )
            return 0.0e0

        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemID)
        pp = c_double(p)  # pressure scalar
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # compute mass density from mass fraction
        iErr = ck_wrapper.chemkin.KINGetMassDensity(chemset_index, tt, pp, yy, den_C)
        if iErr == 0:
            return den_C.value
        else:
            # failed to compute mixture density
            print(Color.RED + "** failed to compute mixture density", end=Color.END)
            return 0.0e0

    @property
    def RHO(self):
        """
        Get mixture mass density [gm/cm3] (double scalar)
        """
        # initialization
        den = 0.0e0
        # check pressure
        if self._Pset == 0:
            print(
                Color.PURPLE + "** mixture pressure [dynes/cm2] is not provided",
                end=Color.END,
            )
            return den
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return den
        #
        if self._Xset == 1:
            # mixture mole fraction given
            den = Mixture.density(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._WT,
                mode="mole",
            )
            return den
        elif self._Yset == 1:
            # mixture mass fraction given
            den = Mixture.density(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._WT,
                mode="mass",
            )
            return den
        else:
            # mixture composition is not provided
            print(
                Color.PURPLE + "** mixture composition is not provided",
                end=Color.END,
            )
            return den

    @staticmethod
    def mixturespecificheat(chemID, p, t, frac, wt, mode):
        """
        Get mixture specific heat capacity from the given mixture condition: pressure, temperature, and species composition
        :param chemID: chemistry set index associated with the mixture (integer scalar)
        :param p: mixture pressure in [dynes/cm2] (double scalar)
        :param t: mixture temperature in [K] (double scalar)
        :param frac: mixture composition given by either mass or mole fractions as specified by mode (double array)
        :param wt: molar masses of the species in the mixture in [gm/mol] (double array)
        :param mode: flag indicates the frac array is 'mass' or 'mole' fractions
        :return: mixture specific heat capacity [erg/mol-K] (double scalar)
        """
        # check inputs
        if chemID < 0:
            print(Color.RED + "** invalid chemistry", end=Color.END)
            return 0.0e0

        if t <= 10.0:
            print(Color.RED + "** invalid temperature value", end=Color.END)
            return 0.0e0

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            print(
                Color.RED
                + f"** fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return 0.0e0

        # initialization
        CpB_C = c_double(0.0)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.molefractiontomassfraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            iErr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            print(
                Color.PURPLE + '** must specify "mole" or "mass" fractions given',
                end=Color.END,
            )
            return 0.0e0
        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemID)
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # real-gas
        if checkrealgasstatus(chemID):
            # real-gas cubic EOS is active, set current pressure that is required by the chemkin real-gas module
            setcurrentpressure(chemID, p)
        # compute mass density from mass fraction
        iErr = ck_wrapper.chemkin.KINGetGasMixtureSpecificHeat(
            chemset_index, tt, yy, CpB_C
        )
        # compute mean molar mass
        mwt = Mixture.meanmolarmass(frac=y, wt=wt, mode="mass")
        if iErr == 0:
            return CpB_C.value * mwt
        else:
            # failed to compute mixture specific heat
            print(
                Color.RED + "** failed to compute mixture specific heat",
                end=Color.END,
            )
            return 0.0e0

    @staticmethod
    def mixtureenthalpy(chemID, p, t, frac, wt, mode):
        """
        Get mixture enthalpy from the given mixture condition: pressure, temperature, and species composition
        :param chemID: chemistry set index associated with the mixture (integer scalar)
        :param p: mixture pressure in [dynes/cm2] (double scalar)
        :param t: mixture temperature in [K] (double scalar)
        :param frac: mixture composition given by either mass or mole fractions as specified by mode (double array)
        :param wt: molar masses of the species in the mixture in [gm/mol] (double array)
        :param mode: flag indicates the frac array is 'mass' or 'mole' fractions
        :return: mixture enthalpy [erg/mol] (double scalar)
        """
        # check inputs
        if chemID < 0:
            print(Color.RED + "** invalid chemistry", end=Color.END)
            return 0.0e0

        if t <= 10.0:
            print(Color.RED + "** invalid temperature value", end=Color.END)
            return 0.0e0

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            print(
                Color.RED
                + f"** fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return 0.0e0

        # initialization
        H_C = c_double(0.0)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.molefractiontomassfraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            iErr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            print(
                Color.PURPLE + '** must specify "mole" or "mass" fractions given',
                end=Color.END,
            )
            return 0.0e0
        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemID)
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # real-gas
        if checkrealgasstatus(chemID):
            # real-gas cubic EOS is active, set current pressure that is required by the chemkin real-gas module
            setcurrentpressure(chemID, p)
        # compute enthalpy from mass fraction
        iErr = ck_wrapper.chemkin.KINGetGasMixtureEnthalpy(chemset_index, tt, yy, H_C)
        # compute mean molar mass
        mwt = Mixture.meanmolarmass(frac=y, wt=wt, mode="mass")
        if iErr == 0:
            return H_C.value * mwt
        else:
            # failed to compute mixture enthalpy
            print(
                Color.RED + "** failed to compute mixture enthalpy",
                end=Color.END,
            )
            return 0.0e0

    @staticmethod
    def rateofproduction(chemID, p, t, frac, wt, mode):
        """
        Get species molar rate of production from the given mixture condition: pressure, temperature, and species composition
        :param chemID: chemistry set index associated with the mixture (integer scalar)
        :param p: mixture pressure in [dynes/cm2] (double scalar)
        :param t: mixture temperature in [K] (double scalar)
        :param frac: mixture composition given by either mass or mole fractions as specified by mode (double array)
        :param wt: molar masses of the species in the mixture in [gm/mol] (double array)
        :param mode: flag indicates the frac array is 'mass' or 'mole' fractions
        :return: species molar rate of production in [mol/cm3-sec] (double array)
        """
        # check inputs
        if chemID < 0:
            print(Color.RED + "** invalid chemistry", end=Color.END)
            return [0.0e0]

        if p <= 0.0 or (p * t) <= 0.0:
            print(
                Color.RED + "** invalid pressure and/or temperature value(s)",
                end=Color.END,
            )
            return [0.0e0]

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            print(
                Color.RED
                + f"** fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return [0.0e0]

        # initialization
        ROP = np.zeros(kgas, dtype=np.double)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.molefractiontomassfraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            iErr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            print(
                Color.PURPLE + '** must specify "mole" or "mass" fractions given',
                end=Color.END,
            )
            return ROP

        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemID)
        pp = c_double(p)  # pressure scalar
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # compute mass density from mass fraction
        iErr = ck_wrapper.chemkin.KINGetGasROP(chemset_index, tt, pp, yy, ROP)
        if iErr == 0:
            return ROP
        else:
            # failed to compute species molar rates of production
            print(
                Color.RED + "** failed to compute species molar rates of production",
                end=Color.END,
            )
            return [0.0e0]

    @staticmethod
    def reactionrates(chemID, numbreaction, p, t, frac, wt, mode):
        """
        Get molar rates of the gas reactions from the given mixture condition: pressure, temperature, and species composition
        :param chemID: chemistry set index associated with the mixture (integer scalar)
        :param numbreaction: number of gas reactions associated with the chemistry set (integer scalar)
        :param p: mixture pressure in [dynes/cm2] (double scalar)
        :param t: mixture temperature in [K] (double scalar)
        :param frac: mixture composition given by either mass or mole fractions as specified by mode (double array)
        :param wt: molar masses of the species in the mixture in [gm/mol] (double array)
        :param mode: flag indicates the frac array is 'mass' or 'mole' fractions
        :return: forward and reverse molar rates of the reactions in [mol/cm3-sec] (double array, double array)
        """
        # check inputs
        if chemID < 0:
            print(Color.RED + "** invalid chemistry", end=Color.END)
            return [0.0e0], [0.0e0]

        if p <= 0.0 or (p * t) <= 0.0:
            print(
                Color.RED + "** invalid pressure and/or temperature value(s)",
                end=Color.END,
            )
            return [0.0e0], [0.0e0]

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            print(
                Color.RED
                + f"** fraction and molar mass arrays must have the same size: {kgas:d}",
                end=Color.END,
            )
            return [0.0e0], [0.0e0]

        # initialization
        Kforward = np.zeros(numbreaction, dtype=np.double)
        Kreverse = np.zeros_like(Kforward, dtype=np.double)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.molefractiontomassfraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            iErr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            print(
                Color.PURPLE + '** must specify "mole" or "mass" fractions given',
                end=Color.END,
            )
            return Kforward, Kreverse

        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemID)
        pp = c_double(p)  # pressure scalar
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # compute mass density from mass fraction
        iErr = ck_wrapper.chemkin.KINGetGasReactionRates(
            chemset_index, tt, pp, yy, Kforward, Kreverse
        )
        if iErr == 0:
            return Kforward, Kreverse
        else:
            # failed to compute reaction rates
            print(Color.RED + "** failed to compute reaction rates", end=Color.END)
            return [0.0e0], [0.0e0]

    def FindEquilibrium(self):
        """
        Create the equilibrium state mixture corresponding to mixture itself
        :return: gas mixture at the equilibrium state (Mixture object)
        """
        # initialization a Mixture object by duplication
        EQState = copy.deepcopy(self)
        # reset mass/mole fractions
        EQState._Xset = 0
        EQState._molefrac[:] = 0.0e0
        EQState._Yset = 0
        EQState._massfrac[:] = 0.0e0
        # compute the equilibrium state (mass fraction for now)
        eqvars, EQState._massfrac = calculateequilibrium(
            self._chemset_index.value,
            p=EQState.pressure,
            t=EQState.temperature,
            frac=self.Y,
            wt=self._WT,
            mode_in="mass",
            mode_out="mass",
        )
        if np.sum(EQState._massfrac, dtype=np.double) > 0.0e0:
            EQState._Yset = 1
        return EQState

    def HML(self):
        """
        Get enthalpy of the mixture
        :return: mixture enthalpy [erg/mol] (double scalar)
        """
        # initialization
        hml = 0.0e0
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return hml

        if self._Xset == 1:
            # mixture mole fraction given
            hml = Mixture.mixtureenthalpy(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._WT,
                mode="mole",
            )
            return hml
        elif self._Yset == 1:
            # mixture mass fraction given
            hml = Mixture.mixtureenthalpy(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._WT,
                mode="mass",
            )
            return hml
        else:
            # mixture composition is not provided
            print(
                Color.PURPLE + "** mixture composition is not provided",
                end=Color.END,
            )
            return hml

    def CPBL(self):
        """
        Get specific heat capacity of the mixture
        :return: mixture specific heat capacity [erg/mol-K] (double scalar)
        """
        # initialization
        cpbl = 0.0e0
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return cpbl
        #
        if self._Xset == 1:
            # mixture mole fraction given
            cpbl = Mixture.mixturespecificheat(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._WT,
                mode="mole",
            )
            return cpbl
        elif self._Yset == 1:
            # mixture mass fraction given
            cpbl = Mixture.mixturespecificheat(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._WT,
                mode="mass",
            )
            return cpbl
        else:
            # mixture composition is not provided
            print(
                Color.PURPLE + "** mixture composition is not provided",
                end=Color.END,
            )
            return cpbl

    def ROP(self):
        """
        Get species molar rate of production from the given mixture condition: pressure, temperature, and species composition
        :return: species molar rate of production in [mol/cm3-sec] (double array)
        """
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return [0.0e0]
        # check pressure
        if self._Pset == 0:
            print(
                Color.PURPLE + "** mixture pressure [dynes/cm2] is not provided",
                end=Color.END,
            )
            return [0.0e0]
        #
        if self._Xset == 1:
            # mixture mole fraction given
            ROP = Mixture.rateofproduction(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._WT,
                mode="mole",
            )
            return ROP
        elif self._Yset == 1:
            # mixture mass fraction given
            ROP = Mixture.rateofproduction(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._WT,
                mode="mass",
            )
            return ROP
        else:
            # mixture composition is not provided
            print(
                Color.PURPLE + "** mixture composition is not provided",
                end=Color.END,
            )
            return [0.0e0]

    def RxnRates(self):
        """
        Get molar rates of the gas reactions from the given mixture condition: pressure, temperature, and species composition
        :return: forward and reverse molar rates of the reactions in [mol/cm3-sec] (double array, double array)
        """
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return [0.0e0]
        # check pressure
        if self._Pset == 0:
            print(
                Color.PURPLE + "** mixture pressure [dynes/cm2] is not provided",
                end=Color.END,
            )
            return [0.0e0], [0.0e0]
        # initialization
        Kforward = np.zeros(self._IIgas, dtype=np.double)
        Kreverse = np.zeros_like(Kforward, dtype=np.double)
        #
        if self._Xset == 1:
            # mixture mole fraction given
            Kforward, Kreverse = Mixture.reactionrates(
                self._chemset_index.value,
                numbreaction=self._IIgas,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._WT,
                mode="mole",
            )
            return Kforward, Kreverse
        elif self._Yset == 1:
            # mixture mass fraction given
            Kforward, Kreverse = Mixture.reactionrates(
                self._chemset_index.value,
                numbreaction=self._IIgas,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._WT,
                mode="mass",
            )
            return Kforward, Kreverse
        else:
            # mixture composition is not provided
            print(
                Color.PURPLE + "** mixture composition is not provided",
                end=Color.END,
            )
            return [0.0e0], [0.0e0]

    def speciesCp(self):
        """
        Get species specific heat capacity at constant pressure
        :return: species specific heat capacities at constant pressure [ergs/mol-K] (1D double array)
        """
        TT = c_double(self.temperature)
        Cp = np.zeros(self._KK, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpecificHeat(self._chemset_index, TT, Cp)
        if iErr == 0:
            # convert [ergs/g-K] to [ergs/mol-K]
            Cp *= self._WT
        else:
            # failed to compute specific heats
            print(Color.RED + "** failed to compute specific heats", end=Color.END)
        return Cp

    def speciesH(self):
        """
        Get species enthalpy
        :return: species enthalpy [ergs/mol] (1D double array)
        """
        TT = c_double(self.temperature)
        H = np.zeros(self._KK, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpeciesEnthalpy(self._chemset_index, TT, H)
        if iErr == 0:
            # convert [ergs/gm] to [ergs/mol]
            H *= self._WT
        else:
            # failed to compute enthalpies
            print(
                Color.RED + "** failed to compute species enthalpies",
                end=Color.END,
            )
        return H

    def speciesVisc(self):
        """
        Get species viscosity
        :return: species viscosity [gm/cm-sec] (1D double array)
        """
        TT = c_double(self.temperature)
        visc = np.zeros(self._KK, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetViscosity(self._chemset_index, TT, visc)
        if iErr != 0:
            # failed to compute viscosities
            print(
                Color.RED + "** failed to compute species viscosity",
                end=Color.END,
            )
        return visc

    def speciesCond(self):
        """
        Get species conductivity
        :return: species conductivity [ergs/cm-K-sec] (1D double array)
        """
        TT = c_double(self.temperature)
        cond = np.zeros(self._KK, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetConductivity(self._chemset_index, TT, cond)
        if iErr != 0:
            # failed to compute conductivities
            print(
                Color.RED + "** failed to compute species conductivity",
                end=Color.END,
            )
        return cond

    def speciesDiffusionCoeffs(self):
        """
        Get species diffusion coefficients
        :return: species diffusion coefficients [cm2/sec] (2D double array)
        """
        PP = c_double(self.pressure)
        TT = c_double(self.temperature)
        dim = (self._KK, self._KK)
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

    def mixtureviscosity(self):
        """
        Get viscosity of the gas mixture
        :return: mixture viscosity [gm/cm-sec] (double scalar)
        """
        # initialization
        visc = c_double(0.0e0)
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return visc.value
        # get mixture viscosity
        TT = c_double(self.temperature)
        iErr = ck_wrapper.chemkin.KINGetMixtureViscosity(
            self._chemset_index, TT, self.Y, visc
        )
        if iErr != 0:
            # error message
            print(
                Color.PURPLE + "** failed to compute mixture viscosity",
                end=Color.END,
            )
        # mixture viscosity in gm/cm-sec
        return visc.value

    def mixtureconductivity(self):
        """
        Get conductivity of the gas mixture
        :return: mixture conductivity [erg/cm-K-sec] (double scalar)
        """
        # initialization
        cond = c_double(0.0e0)
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return cond.value
        # get mixture viscosity
        TT = c_double(self.temperature)
        iErr = ck_wrapper.chemkin.KINGetMixtureConductivity(
            self._chemset_index, TT, self.Y, cond
        )
        if iErr != 0:
            # error message
            print(
                Color.PURPLE + "** failed to compute mixture conductivity",
                end=Color.END,
            )
        # mixture conductivity in ergs/cm-K-sec
        return cond.value

    def mixturediffusioncoeffs(self):
        """
        Get mixture-averaged species diffusion coefficients of the gas mixture
        :return: mixture-averaged diffusion coefficients [cm2/sec] (1D double array)
        """
        # initialization
        diffusioncoeffs = np.zeros(self._KK, dtype=np.double)
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return diffusioncoeffs
        # check pressure
        if self._Pset == 0:
            print(
                Color.PURPLE + "** mixture pressure [dynes/cm2] is not provided",
                end=Color.END,
            )
            return diffusioncoeffs
        # get mixture viscosity
        PP = c_double(self.pressure)
        TT = c_double(self.temperature)
        iErr = ck_wrapper.chemkin.KINGetMixtureDiffusionCoeffs(
            self._chemset_index, PP, TT, self.Y, diffusioncoeffs
        )
        if iErr != 0:
            # error message
            print(
                Color.PURPLE
                + "** failed to compute mixture-averaged diffusion coefficients",
                end=Color.END,
            )
        # mixture-averaged diffusion coefficients in cm2/sec
        return diffusioncoeffs

    def mixturebinarydiffusioncoeffs(self):
        """
        Get multi-component species binary diffusion coefficients of the gas mixture
        :return: binary diffusion coefficients [cm2/sec] (2D double array)
        """
        # initialization
        dim = (self._KK, self._KK)
        binarydiffusioncoeffs = np.zeros(dim, dtype=np.double, order="F")
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return binarydiffusioncoeffs
        # check pressure
        if self._Pset == 0:
            print(
                Color.PURPLE + "** mixture pressure [dynes/cm2] is not provided",
                end=Color.END,
            )
            return binarydiffusioncoeffs
        # get mixture viscosity
        PP = c_double(self.pressure)
        TT = c_double(self.temperature)
        iErr = ck_wrapper.chemkin.KINGetOrdinaryDiffusionCoeffs(
            self._chemset_index, PP, TT, self.Y, binarydiffusioncoeffs
        )
        if iErr != 0:
            # error message
            print(
                Color.PURPLE
                + "** failed to compute multi-component binary diffusion coefficients",
                end=Color.END,
            )
        # mixture multi-component binary diffusion coefficients in cm2/sec
        return binarydiffusioncoeffs

    def mixturethermaldiffusioncoeffs(self):
        """
        Get thermal diffusivity of the gas mixture
        :return: thermal diffusivity [gm/cm-sec] (1D double array)
        """
        # initialization
        thermaldiffusioncoeffs = np.zeros(self._KK, dtype=np.double)
        cond = c_double(0.0e0)  # mixture thermal conductivity
        # check temperature
        if self._Tset == 0:
            print(
                Color.PURPLE + "** mixture temperature [K] is not provided",
                end=Color.END,
            )
            return thermaldiffusioncoeffs
        # check pressure
        if self._Pset == 0:
            print(
                Color.PURPLE + "** mixture pressure [dynes/cm2] is not provided",
                end=Color.END,
            )
            return thermaldiffusioncoeffs
        # get mixture viscosity
        PP = c_double(self.pressure)
        TT = c_double(self.temperature)
        iErr = ck_wrapper.chemkin.KINGetThermalDiffusionCoeffs(
            self._chemset_index, PP, TT, self.Y, thermaldiffusioncoeffs, cond
        )
        if iErr != 0:
            # error message
            print(
                Color.PURPLE
                + "** failed to compute mixture thermal diffusion coefficients",
                end=Color.END,
            )
        # mixture thermal diffusion coefficients in gm/cm-sec
        return thermaldiffusioncoeffs

    def volHRR(self):
        """
        Get volumetric heat release rate
        :return: volumetric heat release rate [ergs/cm3-sec] (double scalar)
        """
        volHRR = 0.0e0
        # get species enthalpy
        TT = c_double(self.temperature)
        H = np.zeros(self._KK, dtype=np.double)
        iErr = ck_wrapper.chemkin.KINGetGasSpeciesEnthalpy(self._chemset_index, TT, H)
        if iErr == 0:
            # convert H from ergs/gm to ergs/mol
            H *= self._WT
        else:
            return volHRR
        # get species molar rate of production mol/cm3-sec
        ROP = self.ROP()
        # volumetric heat release rate = SUM(H_k * ROP_k)  ergs/cm3-sec
        volHRR = np.dot(H, ROP)
        return volHRR

    def massROP(self):
        """
        Get species mass rates of production
        :return: mass rates of production [gm/cm3-sec] (double array)
        """
        # get species molar rate of production mol/cm3-sec
        ROP = self.ROP()
        # species mass rate of production = ROP_k * WT_k
        massROP = ROP * self._WT
        return massROP

    def listROP(self, threshold=0.0):
        """
        list information about species molar production rate in descending order
        :return: sorted species index, sorted ROP values [-, mol/cm3-sec] (1D integer array, 1D double array)
        """
        # get species molar rate of production mol/cm3-sec
        ROP = self.ROP()
        # create a copy of non-zero rates
        temp_ROP = np.zeros_like(ROP, dtype=np.double)
        temp_order = np.zeros_like(ROP, dtype=np.int32)
        # include reactions with non-zero rate only
        count = 0
        for i in range(len(ROP)):
            if abs(ROP[i]) > threshold:
                # non-zero entries
                temp_ROP[count] = ROP[i]
                temp_order[count] = i
                count += 1

        # sort on the temporary array in descending order
        sorted_ROP = np.flip(np.sort(temp_ROP[:count]))
        # find species index
        new_order = np.zeros_like(sorted_ROP, dtype=np.int32)
        for i in range(len(sorted_ROP)):
            count, speciesID = whereelementinarray1D(temp_ROP, sorted_ROP[i])
            # get the species index corresponding to the first occurrence
            # in case there are multiple species having the same rate
            new_order[i] = temp_order[speciesID[0]]
        # print out the list of species with its ROP value in descending order
        if verbose():
            print("non-zero species molar rate of production ")
            print("=" * 50)
            print(" order    species symbol     rate of production")
            print("                             [mol/cm3-sec]")
            for i in range(len(new_order)):
                print(
                    f" {i+1:-2d} {self._specieslist[new_order[i]]:>16}              {sorted_ROP[i]: e}"
                )
        return new_order, sorted_ROP

    def listmassROP(self, threshold=0.0):
        """
        list information about species mass rate of production in descending order
        :return: sorted species index, sorted massROP values [-, gm/cm3-sec] (1D integer array, 1D double array)
        """
        # get species mass rate of production gm/cm3-sec
        massROP = self.massROP()
        # create a copy of non-zero rates
        temp_ROP = np.zeros_like(massROP, dtype=np.double)
        temp_order = np.zeros_like(massROP, dtype=np.int32)
        # include reactions with non-zero rate only
        count = 0
        for i in range(len(massROP)):
            if abs(massROP[i]) > threshold:
                # non-zero entries
                temp_ROP[count] = massROP[i]
                temp_order[count] = i
                count += 1

        # sort on the temporary array in descending order
        sorted_ROP = np.flip(np.sort(temp_ROP[:count]))
        # find species index
        new_order = np.zeros_like(sorted_ROP, dtype=np.int32)
        for i in range(len(sorted_ROP)):
            count, speciesID = whereelementinarray1D(temp_ROP, sorted_ROP[i])
            # get the species index corresponding to the first occurrence
            # in case there are multiple species having the same rate
            new_order[i] = temp_order[speciesID[0]]
        # print out the list of species with its ROP value in descending order
        if verbose():
            print("non-zero species mass rate of production ")
            print("=" * 50)
            print(" order    species symbol     rate of production")
            print("                             [gm/cm3-sec]")
            for i in range(len(new_order)):
                print(
                    f" {i+1:-2d} {self._specieslist[new_order[i]]:>16}              {sorted_ROP[i]: e}"
                )
        return new_order, sorted_ROP

    def listreactionrates(self, threshold=0.0):
        """
        list information about reaction rate in descending order
        :return: sorted reaction index, sorted reaction rate values [-, mol/cm3-sec] (1D integer array, 1D double array)
        """
        # molar rates of reactions
        RF, RR = self.RxnRates()
        # create a copy of non-zero rates
        temp_netRR = np.zeros_like(RF, dtype=np.double)
        temp_order = np.zeros_like(RF, dtype=np.int32)
        # include reactions with non-zero rate only
        count = 0
        for i in range(len(RF)):
            netRR = RF[i] - RR[i]
            if abs(netRR) > threshold:
                # non-zero entries
                temp_netRR[count] = netRR
                temp_order[count] = i
                count += 1

        # sort on the temporary array in descending order
        sorted_RR = np.flip(np.sort(temp_netRR[:count]))
        # find species index
        new_order = np.zeros_like(sorted_RR, dtype=np.int32)
        new_RR = copy.deepcopy(sorted_RR)
        for i in range(len(sorted_RR)):
            # find the instances of this reaction rate
            count, rxnID = whereelementinarray1D(temp_netRR, new_RR[i])
            # get the reaction number corresponding to the first occurrence
            # in case there are multiple reactions having the same rate
            new_order[i] = temp_order[rxnID[0]]
            # remove this instance from the reaction rate array
            new_RR[i] = 0.0e0
        # print out the list of reaction with its net reaction rate value in descending order
        if verbose():
            print("non-zero molar rates of reaction ")
            print("=" * 50)
            print(" order    reaction number    molar rate of reaction")
            print("                             [mol/cm3-sec]")
            for i in range(len(new_order)):
                print(
                    f" {i+1:-2d}          {new_order[i]+1:-4d}              {sorted_RR[i]: e}"
                )
        return new_order, sorted_RR

    def XbyEquivalenceRatio(
        self,
        chemistryset,
        fuel_molefrac,
        oxid_molefrac,
        add_molefrac,
        products,
        equivalenceratio,
        threshold=1.0e-10,
    ):
        """
        Specify the mixture molar composition by providing the equivalence ratio, the mole fractions of the fuel mixture,
        the oxidizer mixture, and the additives mixture, and the list of the complete combustion product species.
        :param chemistryset: the chemistry set used to create the mixtures (Chemistry object)
        :param fuel_molefrac: mole fractions of the fuel mixture (1D double array)
        :param oxid_molefrac: mole fractions of the oxidizer mixture (1D double array)
        :param add_molefrac: mole fractions of the additives mixture (1D double array)
        :param products: list of the complete combustion species symbols (list of strings)
        :param equivalenceratio: equivalence ratio of the final mixture (double scalar)
        :return: Error status (integer scalar)
        """
        # check chemistry set
        if not isinstance(chemistryset, Chemistry):
            print(
                Color.RED + "** the first argument must be a Chemistry object",
                end=Color.END,
            )
            return 1
        # number of gas species in the mechanism
        kspecies = chemistryset.KK
        # find fuel mole array size
        kfuel = len(fuel_molefrac)
        # find oxidizer mole array size
        koxid = len(oxid_molefrac)
        # find additives mole array size
        kadd = len(add_molefrac)
        # check species number consistency
        iErr = 0
        if kspecies != kfuel:
            print(
                Color.RED
                + f"** the fuel mole fraction array must have size {kspecies:d}",
                end=Color.END,
            )
            iErr += 1
        if kspecies != koxid:
            print(
                Color.RED
                + f"** the oxidizer mole fraction array must have size {kspecies:d}",
                end=Color.END,
            )
            iErr += 1
        if kspecies != kadd:
            print(
                Color.RED
                + f"** the additives mole fraction array must have size {kspecies:d}",
                end=Color.END,
            )
            iErr += 1
        if iErr > 0:
            return 2
        # check equivalence ratio value
        if equivalenceratio <= 0.0e0:
            print(
                Color.RED + "** the equivalence ratio must be > 0",
                end=Color.END,
            )
            return 3
        # check product species
        kprod = len(products)
        if kprod == 0:
            print(
                Color.RED + "** no complete combustion product species given",
                end=Color.END,
            )
            return 4
        # find sum of additives fraction
        suma = 0.0e0
        if kadd > 0:
            # remove tiny mole fraction values
            for i in range(len(add_molefrac)):
                if add_molefrac[i] < threshold:
                    add_molefrac[i] = 0.0e0

            suma = np.sum(add_molefrac)
        # find product species index
        prod_index = np.zeros(kprod, dtype=np.int32)
        j = 0
        for s in products:
            prod_index[j] = chemistryset.getspecindex(s)
            j += 1
        # find the stoichiometric coefficients assuming complete combustion
        alpha, nu = calculatestoichiometrics(
            chemistryset, fuel_molefrac, oxid_molefrac, prod_index
        )
        if alpha <= 0.0e0 or nu[0] == 0:
            print(
                Color.RED + "** failed to find the stoichiometric coefficients",
                end=Color.END,
            )
            return 5
        # find the fuel-oxidizer mixture molar composition
        self._molefrac[:] = 0.0e0
        self._molefrac = equivalenceratio * fuel_molefrac + alpha * oxid_molefrac
        # normalize the mole fractions
        sumx = np.sum(self._molefrac)
        if sumx > 0.0e0:
            ratio = (1.0e0 - suma) / sumx
            self._molefrac *= ratio
            # include additives fractions
            if kadd > 0:
                self._molefrac += add_molefrac
            # set the composition flags of the final mixture
            self._Xset = 1
            self._massfrac[:] = 0.0e0
            self._Yset = 0
            return 0
        else:
            print(
                Color.RED + "** failed to find the stoichiometric coefficients",
                end=Color.END,
            )
            return 6

    def YbyEquivalenceRatio(
        self,
        chemistryset,
        fuel_massfrac,
        oxid_massfrac,
        add_massfrac,
        products,
        equivalenceratio,
        threshold=1.0e-10,
    ):
        """
        Specify the mixture molar composition by providing the equivalence ratio, the mole fractions of the fuel mixture,
        the oxidizer mixture, and the additives mixture, and the list of the complete combustion product species.
        :param chemistryset: the chemistry set used to create the mixtures (Chemistry object)
        :param fuel_massfrac: mass fractions of the fuel mixture (1D double array)
        :param oxid_massfrac: mass fractions of the oxidizer mixture (1D double array)
        :param add_massfrac: mass fractions of the additives mixture (1D double array)
        :param products: list of the complete combustion species symbols (list of strings)
        :param equivalenceratio: equivalence ratio of the final mixture (double scalar)
        :return: Error status (integer scalar)
        """
        # check chemistry set
        if not isinstance(chemistryset, Chemistry):
            print(
                Color.RED + "** the first argument must be a Chemistry object",
                end=Color.END,
            )
            return 1
        # convert mass fractions to mole fractions
        fuel_molefrac = Mixture.massfractiontomolefraction(
            massfrac=fuel_massfrac, wt=chemistryset.WT
        )
        oxid_molefrac = Mixture.massfractiontomolefraction(
            massfrac=oxid_massfrac, wt=chemistryset.WT
        )
        add_molefrac = Mixture.massfractiontomolefraction(
            massfrac=add_massfrac, wt=chemistryset.WT
        )
        # find the final mixture mole fractions and set the flags
        iErr = self.XbyEquivalenceRatio(
            chemistryset,
            fuel_molefrac,
            oxid_molefrac,
            add_molefrac,
            products,
            equivalenceratio,
            threshold,
        )
        return iErr

    def getEGRmolefraction(self, EGRratio, threshold=1.0e-8):
        """
        Compute the EGR composition in mole fraction corresponding to this mixture
        """
        # create burned mixture
        burned = self.FindEquilibrium()
        # compute the EGR stream mole fractions
        klength = len(burned.X)
        EGR_molefrac = np.zeros(klength, dtype=np.double)
        for i in range(klength):
            if burned.X[i] > threshold:
                EGR_molefrac[i] = EGRratio * burned.X[i]
        del burned
        return EGR_molefrac

    def validate(self):
        """
        Check whether the mixture is fully defined before being used by other methods
        :return: Error status (integer scalar)
        """
        iErr = 0
        # check mixture temperature
        if self._Tset == 0:
            print(
                Color.YELLOW + "** mixture temperature is not provided",
                end=Color.END,
            )
            iErr = 1
        if self._Pset == 0:
            print(
                Color.YELLOW + "** mixture pressure is not provided",
                end=Color.END,
            )
            iErr = 1
        if self._Xset == 0 and self._Yset == 0:
            print(
                Color.YELLOW + "** mixture composition is not provided",
                end=Color.END,
            )
            iErr = 1
        return iErr

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
                Color.RED + f"** Warning: method returned error code {iErr}",
                end=Color.END,
            )
        if iFlag.value == 0:
            print(
                Color.YELLOW
                + f"** real-gas cubic EOS model {Chemistry.realgas_CuEOS[self._EOS.value]} is turned ON",
                end=Color.END,
            )
            # set default mixing rule to Van der Waals
            mixingrule = 0
            # set default mixing rule to Van der Waals
            self.setrealgasmixingrule(mixingrule)
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
                Color.RED + f"** Warning: method returned error code {iErr}",
                end=Color.END,
            )
        if iFlag.value == 0:
            print(Color.YELLOW + "** ideal gas law is turned ON", end=Color.END)
            self.userealgas = False

    def setrealgasmixingrule(self, rule=0):
        """
        Set the mixing rule to be used for calculating the real-gas mixture properties
        :param rule: mixing rule: 0 for the Van der Waals mixing rule; 1 for the critical properties mixing rule (integer scalar)
        :return: None
        """
        if self._EOS.value < 1:
            # no real gas EOS data in the mechanism
            print(
                Color.YELLOW + "** mechanism is for ideal gas law only",
                end=Color.END,
            )
            return
        # set default mixing rule to Van der Waals
        mixingrule = c_int(rule + 1)
        iFlag = c_int(0)
        iErr = ck_wrapper.chemkin.KINRealGas_SetMixingRule(
            self._chemset_index, mixingrule, iFlag
        )
        if iErr != 0:
            print(
                Color.RED + f"** Warning: method returned error code {iErr}",
                end=Color.END,
            )
        if iFlag.value == 2:
            # real-gas cubic EOS is turned OFF
            print(Color.YELLOW + "** the ideal gas law is used", end=Color.END)
            self.userealgas = False
        elif iFlag.value != 0:
            print(
                Color.PURPLE
                + f"** error setting up real-gas cubic EOS (code = {iFlag.value})",
                end=Color.END,
            )
        else:
            print(
                Color.YELLOW
                + f"** the cubic EOS is used, set mixing rule to '{Chemistry.realgas_mixingrules[rule]}'",
                end=Color.END,
            )
            self.userealgas = True

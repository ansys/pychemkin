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

"""Chemkin Mixture utilities."""

import copy
import ctypes
from ctypes import c_double, c_int

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core import chemkin_wrapper as ck_wrapper
from ansys.chemkin.core.chemistry import (
    Chemistry,
    check_active_chemistryset,
    check_chemistryset,
    check_realgas_status,
    chemistryset_initialized,
    set_current_pressure,
    verbose,
)
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.constants import P_ATM
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.utilities import (
    calculate_stoichiometrics,
    where_element_in_array_1d,
)


class Mixture:
    """define a mixture based on the gas species in the given chemistry set."""

    def __init__(self, chem: Chemistry):
        """Initialize a Mixture object based on the given Chemistry set."""
        """
        Parameters
        ----------
            chem: Chemistry object

        """
        self._temp = 0.0e0  # mixture temperature [K]
        self._press = 0.0e0  # mixture pressure [dynes/cm2]
        self._vol = 0.0e0  # mixture volume [cm3]
        # flags
        self._t_set = 0
        self._p_set = 0
        self._x_set = 0
        self._y_set = 0
        # chemistry set validation
        if not isinstance(chem, Chemistry):
            msg = [
                Color.RED,
                "the argument must be a Chemkin.Chemistry object.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        if chem.chemid < 0:
            msg = [
                Color.RED,
                "invalid chemistry,",
                "please preprocess the chemistry first.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # shorthand for frequently used variables
        self._chemset_index = ctypes.c_int(chem.chemid)  # chemistry set index
        self._kk = chem.kk  # number of gas species
        self._ii_gas = chem.ii_gas  # number of gas-phase reactions
        self._specieslist: list[str] = []
        self._specieslist = chem.species_symbols  # gas species symbols
        self._wt = chem.wt  # gas species molar masses
        # create internal arrays: array size = number of gas species
        self._molefrac = np.zeros(
            self._kk, dtype=np.double
        )  # mixture composition given in mole fractions
        self._massfrac = np.zeros_like(
            self._molefrac
        )  # mixture composition given in mole fractions
        self._concentration = np.zeros_like(self._molefrac)  # concentrations (not used)
        self._surfacechem = c_int(
            chem.surfchem
        )  # flag indicating there is surface chemistry (type c_int: 0 = no, 1 = yes)
        self.transport_data = chem._index_tran.value
        self._eos = c_int(chem.eos)  # real-gas EOS model in the mechanism
        self.userealgas = chem.userealgas  # status of the real-gas EOS usage

    @property
    def chemid(self) -> int:
        """Get chemistry set index."""
        """
        Returns
        -------
            chemid: integer
                chemistry set index associated with this Mixture

        """
        return self._chemset_index.value

    @property
    def kk(self) -> int:
        """Get the number of gas species."""
        """
        Returns
        -------
            num_spec: integer
                number of gas species in the mixture

        """
        return self._kk

    @property
    def pressure(self) -> float:
        """Get gas mixture pressure [dynes/cm2]."""
        """
        Returns
        -------
            pressure: double
                mixture pressure [dynes/cm2]

        """
        if self._p_set == 1:
            return self._press
        else:
            msg = [Color.PURPLE, "mixture pressure is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @pressure.setter
    def pressure(self, p: float):
        """Set gas mixture pressure."""
        """
        Parameters
        ----------
            p: double
                pressure [dynes/cm2]

        """
        if p <= 0.0:
            msg = [Color.PURPLE, "invalid pressure value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        self._press = p
        self._p_set = 1

    @property
    def temperature(self) -> float:
        """Get gas mixture temperature."""
        """
        Returns
        -------
            temperature: double
                temperature [K]

        """
        if self._t_set == 1:
            return self._temp
        else:
            msg = [Color.PURPLE, "mixture temperature is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @temperature.setter
    def temperature(self, t: float):
        """Set gas mixture temperature."""
        """
        Parameters
        ----------
            t: double
                mixture temperature [K]

        """
        if t <= 10.0:
            msg = [Color.PURPLE, "invalid temperature value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self._temp = t
        self._t_set = 1

    @property
    def volume(self) -> float:
        """Get mixture volume."""
        """
        Returns
        -------
            volume: double
                mixture volume [cm3]

        """
        if self._vol > 0.0e0:
            return self._vol
        else:
            msg = [Color.PURPLE, "mixture volume is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @volume.setter
    def volume(self, vol: float):
        """Set mixture volume."""
        """
        Parameters
        ----------
            vol: double
                mixture volume [cm3]

        """
        if vol <= 0.0e0:
            msg = [Color.PURPLE, "invalid volume value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self._vol = vol

    @property
    def x(self) -> npt.NDArray[np.double]:
        """Get mixture mole fraction."""
        """
        Returns
        -------
            x: 1-D double array, dimensdion = number_species
                mixture composition in mole fractions

        """
        if self._x_set == 1:
            ierr, x = Mixture.normalize(self._molefrac)
            return x
        elif self._y_set == 1:
            ierr, x = self.__ytox()
            if ierr != 0:
                msg = [Color.PURPLE, "fraction conversion failed.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            return x
        else:
            msg = [Color.PURPLE, "mixture composition not yet defined.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @x.setter
    def x(self, recipe: list[tuple[str, float]]):
        """Set mixture molar composition."""
        """
        Parameters
        ----------
        recipe: list of tuples, [(species_symbol, fraction), ... ]
            non-zero mixture composition corresponding to
            the given mole fraction array

        """
        if self._x_set == 1:
            # reset the mole fraction array
            self._molefrac[:] = 0.0e0
        if isinstance(recipe[0], tuple):
            for sp, x in recipe:
                if sp in self._specieslist:
                    index = self._specieslist.index(sp)
                else:
                    msg = [Color.PURPLE, sp, "is not a valid gas species.", Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
                if x < 0.0:
                    msg = [Color.PURPLE, "negative mole fraction.", Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
                # set mole fraction
                self._molefrac[index] = x
        elif isinstance(recipe[0], (float, np.double)):
            kgas = len(recipe)
            if kgas == self._kk:
                for k in range(kgas):
                    self._molefrac[k] = max(recipe[k], 0.0e0)
            else:
                msg = [
                    Color.PURPLE,
                    "size of the mole fraction array must equal to",
                    "the number of gas species:",
                    str(self._kk),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        else:
            msg = [
                Color.PURPLE,
                "the argument must be:\n",
                Color.SPACEx6,
                "(1) a list of tuples: [('O2', 0.21), ('N2', 0.79)]\n",
                "or\n",
                Color.SPACEx6,
                "(2) a mole fraction array of size = <number of gas species>",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # reset mass fraction
        self._y_set = 0
        self._massfrac[:] = 0.0e0
        self._x_set = 1

    # alias
    mole_fractions = x

    @property
    def y(self) -> npt.NDArray[np.double]:
        """Get mixture mass fraction."""
        """
        Returns
        -------
            y: 1-D double array, dimensdion = number_species
                mixture composition in mass fractions

        """
        if self._y_set == 1:
            ierr, y = Mixture.normalize(self._massfrac)
            return y
        elif self._x_set == 1:
            ierr, y = self.__xtoy()
            if ierr != 0:
                msg = [Color.PURPLE, "fraction conversion failed.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            return y
        else:
            msg = [Color.PURPLE, "mixture composition not yet defined.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @y.setter
    def y(self, recipe: list[tuple[str, float]]):
        """Set mixture mass composition."""
        """
        Parameters
        ----------
        recipe: list of tuples, [(species_symbol, fraction), ... ]
            non-zero mixture composition corresponding to the given mass fraction array

        """
        if self._y_set == 1:
            # reset the mass fraction array
            self._massfrac[:] = 0.0e0
        if isinstance(recipe[0], tuple):
            for sp, y in recipe:
                if sp in self._specieslist:
                    index = self._specieslist.index(sp)
                else:
                    msg = [Color.PURPLE, sp, "is not a valid gas species.", Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
                if y < 0.0:
                    msg = [Color.PURPLE, "negative mass fraction value.", Color.END]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
                # set mass fraction
                self._massfrac[index] = y
        elif isinstance(recipe[0], (float, np.double)):
            kgas = len(recipe)
            if kgas == self._kk:
                for k in range(kgas):
                    self._massfrac[k] = max(recipe[k], 0.0e0)
            else:
                msg = [
                    Color.PURPLE,
                    "size of the mass fraction array must equal to",
                    "the number of gas species:",
                    str(self._kk),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        else:
            msg = [
                Color.PURPLE,
                "the argument must be:\n",
                Color.SPACEx6,
                "(1) a list of tuples: [('O2', 0.21), ('N2', 0.79)]\n",
                "or\n",
                Color.SPACEx6,
                "(2) a mole fraction array of size = <number of gas species>",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # reset mole fraction
        self._x_set = 0
        self._molefrac[:] = 0.0e0
        self._y_set = 1

    # alias
    mass_fractions = y

    @property
    def concentration(self) -> npt.NDArray[np.double]:
        """Get mixture molar concentrations."""
        """
        Returns
        -------
            c: 1-D double array, dimensdion = number_species
                mixture compisition in molar concentrations [mole/cm3]

        """
        if self._x_set == 1:
            # mole fractions are given
            # remove negative values and normalize fractions
            ierr, c = Mixture.normalize(frac=self._molefrac)
            if ierr == 0:
                # compute mean molar mass
                mwt = self.wtm
                # compute density
                den = self.rho
                fac = den / mwt
                for k in range(self._kk):
                    c[k] *= fac
                self._concentration[:] = c[:]
            return c
        elif self._y_set == 1:
            # mass fractions are given
            # remove negative values and normalize fractions
            ierr, c = Mixture.normalize(frac=self._massfrac)
            if ierr == 0:
                # compute density
                den = self.rho
                for k in range(self._kk):
                    c[k] = c[k] * den / self._wt[k]
                self._concentration[:] = c[:]
            return c
        else:
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def eos(self) -> int:
        """Get the available real-gas EOS model that is provided in the mechanism."""
        """
        Returns
        -------
            EOS: integer
                index of the realgas EOS model defined in the gas-phase
                mechanism input file

        """
        return self._eos.value

    @staticmethod
    def normalize(frac: npt.ArrayLike) -> tuple[int, npt.NDArray[np.double]]:
        """Normalize the mixture composition."""
        """
        Parameters
        ----------
            frac: 1-D double array
                mixture composition to be normalized

        Returns
        -------
            error code: integer
                error code
            localfrac: 1-D double array
                normalized fraction array

        """
        # initialization
        sumx = 0.0e0
        kk = len(frac)  # number of entries
        localfrac = copy.deepcopy(frac)  # make a local copy of the frac array
        # remove negative fraction and calculate sum
        for k in range(kk):
            if localfrac[k] > 0.0:
                sumx += localfrac[k]
            else:
                localfrac[k] = 0.0e0
        # check none zero sum
        if sumx > 0.0:
            # normalization
            for k in range(kk):
                localfrac[k] = localfrac[k] / sumx
            return 0, localfrac
        else:
            # fractions summed to zero
            msg = [Color.PURPLE, "fractions summed to zero.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def wt(self) -> npt.NDArray[np.double]:
        """Get species molecular masses."""
        """
        Returns
        -------
            wt: 1-D double array, dimension = number_species
                species molecular masses [gm/mole]

        """
        return self._wt

    # alias
    species_molar_weight = wt

    @property
    def wtm(self) -> float:
        """Get mean molar mass of the gas mixture."""
        """
        Returns
        -------
            wtm: double
                mean molecular mass of the mixture [gm/mol]

        """
        mwt = 0.0e0
        if self._x_set == 1:
            # mole fractions are given
            # remove negative values and normalize fractions
            ierr, x = Mixture.normalize(frac=self._molefrac)
            if ierr == 0:
                # compute mean molar mass
                for k in range(self._kk):
                    mwt += x[k] * self._wt[k]

            return mwt

        elif self._y_set == 1:
            # mass fractions are given
            # remove negative values and normalize fractions
            ierr, y = Mixture.normalize(frac=self._massfrac)
            if ierr == 0:
                # compute mean molar mass
                for k in range(self._kk):
                    mwt += y[k] / self._wt[k]

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
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    # alias
    mean_molar_weight = wtm

    def __xtoy(self) -> tuple[int, npt.NDArray[np.double]]:
        """Convert mole fraction to mass fraction."""
        """
        Returns
        -------
            error_code: integer
                error code
            y: 1-D double array, dimension = number_species
                mass fractions

        """
        # compute mean molar mass
        mwt = self.wtm
        if mwt > 0.0e0:
            # remove negative values and normalize fractions
            ierr, y = Mixture.normalize(frac=self._molefrac)
            if ierr == 0:
                # convert mole fractions to mass fractions
                for k in range(self._kk):
                    y[k] = y[k] * self._wt[k] / mwt
                return 0, y
            else:
                return ierr, self._molefrac
        else:
            # zero mean molar mass
            msg = [Color.PURPLE, "mean molar mass = 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def __ytox(self) -> tuple[int, npt.NDArray[np.double]]:
        """Convert mass fraction to mole fraction."""
        """
        Returns
        -------
            error_code: integer
                error code
            x: 1-D double array, dimensdion = number_species
                mole fractions

        """
        # compute mean molar mass
        mwt = self.wtm
        if mwt > 0.0e0:
            # remove negative values and normalize fractions
            ierr, x = Mixture.normalize(frac=self._massfrac)
            if ierr == 0:
                # convert mass fractions to mole fractions
                for k in range(self._kk):
                    x[k] = x[k] * mwt / self._wt[k]
                return 0, x
            else:
                return ierr, self._massfrac
        else:
            # zero mean molar mass
            msg = [Color.PURPLE, "mean molar mass = 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @staticmethod
    def mean_molar_mass(
        frac: npt.NDArray[np.double], wt: npt.NDArray[np.double], mode: str
    ) -> float:
        """Get mean molar mass of the gas mixture."""
        """
        Parameters
        ----------
            frac: 1-D double array, dimensdion = number_species
                mixture composition in 'mass' or mole fraction
                as indicated by mode
            wt: 1-D double array, dimensdion = number_species
                species molar mass [gm/mol]
            mode: string, {'mole', 'mass'}
                flag indicates the frac array is 'mass' or 'mole' fractions

        Returns
        -------
            mwt: double
                mean molar mass [gm/mol]

        """
        # initialization
        mwt = 0.0e0
        # check sizes
        kgas = len(frac)
        k = len(wt)
        if k != kgas:
            # mismatch input arrays
            msg = [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if mode.lower() == "mole":
            # mole fractions are given
            # remove negative values and normalize fractions
            ierr, x = Mixture.normalize(frac=frac)
            if ierr == 0:
                # compute mean molar mass
                for k in range(kgas):
                    mwt += x[k] * wt[k]

            return mwt

        elif mode.lower() == "mass":
            # mass fractions are given
            # remove negative values and normalize fractions
            ierr, y = Mixture.normalize(frac=frac)
            # compute mean molar mass
            if ierr == 0:
                for k in range(kgas):
                    mwt += y[k] / wt[k]

            if mwt > 0.0:
                return 1.0e0 / mwt
            else:
                # zero mean molar mass
                return mwt
        else:
            # fractions summed to zero
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @staticmethod
    def mole_fraction_to_mass_fraction(
        molefrac: npt.NDArray[np.double], wt: npt.NDArray[np.double]
    ) -> npt.NDArray[np.double]:
        """Convert mole fraction to mass fraction."""
        """
        Parameters
        ----------
            molefrac: 1-D double array, dimension = number_species
                mixture composition in mole fractions
            wt: 1-D double array, dimension = number_species
                species molar mass [gm/mol]

        Returns
        -------
            massfrac: 1-D double array, dimension = number_species
                mass fractions

        """
        # check size
        kgas = len(molefrac)
        k = len(wt)
        if k != kgas:
            # mismatch input arrays
            msg = [
                Color.PURPLE,
                "mole fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # compute mean molar mass
        mwt = Mixture.mean_molar_mass(frac=molefrac, wt=wt, mode="mole")
        if mwt > 0.0e0:
            # remove negative values and normalize fractions
            ierr, massfrac = Mixture.normalize(frac=molefrac)
            if ierr == 0:
                # convert mole fractions to mass fractions
                for k in range(kgas):
                    massfrac[k] = massfrac[k] * wt[k] / mwt

            return massfrac
        else:
            # zero mean molar mass
            msg = [Color.PURPLE, "mean molar mass = 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @staticmethod
    def mass_fraction_to_mole_fraction(
        massfrac: npt.NDArray[np.double], wt: npt.NDArray[np.double]
    ) -> npt.NDArray[np.double]:
        """Convert mass fraction to mole fraction."""
        """
        Parameters
        ----------
            massfrac: 1-D double array, dimension = number_species
                mixture composition in mass fractions
            wt: 1-D double array, dimension = number_species
                species molar mass [gm/mol]

        Returns
        -------
            molefrac: 1-D double array, dimension = number_species
                mole fractions

        """
        # check size
        kgas = len(massfrac)
        k = len(wt)
        if k != kgas:
            msg = [
                Color.PURPLE,
                "mass fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # compute mean molar mass
        mwt = Mixture.mean_molar_mass(frac=massfrac, wt=wt, mode="mass")
        if mwt > 0.0e0:
            # remove negative values and normalize fractions
            ierr, molefrac = Mixture.normalize(frac=massfrac)
            if ierr == 0:
                # convert mass fractions to mole fractions
                for k in range(kgas):
                    molefrac[k] = molefrac[k] * mwt / wt[k]

            return molefrac
        else:
            # zero mean molar mass
            return massfrac

    @staticmethod
    def mass_fraction_to_concentration(
        chemid: int,
        p: float,
        t: float,
        massfrac: npt.NDArray[np.double],
        wt: npt.NDArray[np.double],
    ) -> npt.NDArray[np.double]:
        """Convert mass fractions to molar concentrations."""
        """
        Parameters
        ----------
            chemid: integer
                chemistry set index associated with the mixture
            p: double
                pressure [dynes/cm2]
            t: double
                temperature [K]
            massfrac: 1-D double array, dimension = number_species
                mixture compisition in mass fractions
            wt: 1-D double array, dimension = number_species
                species molecular masses [gm/mole]

        Returns
        -------
            c: 1-D double array, dimension = number_species
                molar concentrations [mole/cm3]

        """
        # check size
        kgas = len(massfrac)
        k = len(wt)
        if k != kgas:
            msg = [
                Color.PURPLE,
                "mass fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # compute density
        den = Mixture.density(chemid, p, t, frac=massfrac, wt=wt, mode="mass")
        if den > 0.0e0:
            # remove negative values and normalize fractions
            ierr, c = Mixture.normalize(frac=massfrac)
            if ierr == 0:
                # convert mass fractions to mole fractions
                for k in range(kgas):
                    c[k] = c[k] * den / wt[k]
            return c
        else:
            # zero mean molar mass
            return massfrac

    @staticmethod
    def mole_fraction_to_concentration(
        chemid: int,
        p: float,
        t: float,
        molefrac: npt.NDArray[np.double],
        wt: npt.NDArray[np.double],
    ) -> npt.NDArray[np.double]:
        """Convert mole fractions to molar concentrations."""
        """
        Parameters
        ----------
            chemid: integer
                chemistry set index associated with the mixture
            p: double
                pressure [dynes/cm2]
            t: double
                temperature [K]
            molefrac: 1-D double array, dimension = number_species
                mixture compisition in mole fractions
            wt: 1-D double array, dimension = number_species
                species molecular masses [gm/mole]

        Returns
        -------
            c: 1-D double array, dimension = number_species
                molar concentrations [mole/cm3]

        """
        # check size
        kgas = len(molefrac)
        k = len(wt)
        if k != kgas:
            msg = [
                Color.PURPLE,
                "molr fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # compute mean molar mass
        mwt = Mixture.mean_molar_mass(frac=molefrac, wt=wt, mode="mole")
        # compute density
        den = Mixture.density(chemid, p, t, frac=molefrac, wt=wt, mode="mole")
        if mwt * den > 0.0e0:
            # remove negative values and normalize fractions
            ierr, c = Mixture.normalize(frac=molefrac)
            if ierr == 0:
                # convert mass fractions to mole fractions
                fac = den / mwt
                for k in range(kgas):
                    c[k] *= fac
            return c
        else:
            # zero mean molar mass
            return molefrac

    def list_composition(self, mode: str, option: str = " ", bound: float = 0.0e0):
        """List the mixture composition."""
        """
        Parameters
        ----------
            mode: string, {'mole', 'mass'}
                flag indicates the fractions returned are
                'mass' or 'mole' fractions
            option: string, {'all, ' '}, default = 'all'
                flag indicates to list 'all' species or just the species
                with non-zero fraction
            bound: double
                minimum fraction value for the species to be printed

        """
        #
        if option.lower() == "all":
            # list all species
            if mode.lower() == "mass":
                print(f"listing mixture composition in {mode} fractions")
                for k in range(self._kk):
                    print(f"{self._specieslist[k]:18} :  {self.y[k]:e}")
            elif mode.lower() == "mole":
                print(f"listing mixture composition in {mode} fractions")
                for k in range(self._kk):
                    print(f"{self._specieslist[k]:18} :  {self.x[k]:e}")
            else:
                msg = [
                    Color.PURPLE,
                    'must specify output "mole" or "mass" fractions.',
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        else:
            # list non-zero components
            if mode.lower() == "mass":
                print(f"listing mixture composition in {mode} fractions")
                for k in range(self._kk):
                    if self.y[k] > np.max([bound, 0.0e0]):
                        print(f"{self._specieslist[k]:18} :  {self.y[k]:e}")
            elif mode.lower() == "mole":
                print(f"listing mixture composition in {mode} fractions")
                for k in range(self._kk):
                    if self.x[k] > np.max([bound, 0.0e0]):
                        print(f"{self._specieslist[k]:18} :  {self.x[k]:e}")
            else:
                msg = [
                    Color.PURPLE,
                    'must specify output "mole" or "mass" fractions.',
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    @staticmethod
    def density(
        chemid: int,
        p: float,
        t: float,
        frac: npt.NDArray[np.double],
        wt: npt.NDArray[np.double],
        mode: str,
    ) -> float:
        """Get mass density from the given mixture condition."""
        """Get mass density from the given mixture condition:
        pressure, temperature, and species composition.

        Parameters
        ----------
            chemid: integer
                chemistry set index associated with the mixture
            p: double
                mixture pressure in [dynes/cm2]
            t: double
                mixture temperature in [K]
            frac: 1-D double array, dimension = number_species
                mixture composition given by either mass or mole fractions
                as specified by mode
            wt: 1-D double array, dimension = number_species
                molar masses of the species in the mixture in [gm/mol]
            mode: string, {'mole', 'mass'}
                flag indicates the frac array is 'mass' or 'mole' fractions

        Returns
        -------
            den: double
                mass density in [gm/cm3]

        """
        # check inputs
        if chemid < 0:
            msg = [Color.PURPLE, "invalid chemistry.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if p <= 0.0 or (p * t) <= 0.0:
            msg = [
                Color.PURPLE,
                "invalid pressure and/or temperature value(s).",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            msg = [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # initialization
        den_c = c_double(0.0)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.mole_fraction_to_mass_fraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            ierr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            msg = [
                Color.PURPLE,
                'must specify "mole" or "mass" fractions given.',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemid)
        pp = c_double(p)  # pressure scalar
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # compute mass density from mass fraction
        ierr = ck_wrapper.chemkin.KINGetMassDensity(chemset_index, tt, pp, yy, den_c)
        if ierr == 0:
            return den_c.value
        else:
            # failed to compute mixture density
            msg = [Color.PURPLE, "failed to compute mixture density.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def rho(self) -> float:
        """Get mixture mass density."""
        """
        Get mixture mass density.

        Returns
        -------
            den: double
                mixture density [gm/cm3]

        """
        # initialization
        den = 0.0e0
        # check pressure
        if self._p_set == 0:
            msg = [
                Color.PURPLE,
                "mixture pressure [dynes/cm2] is not provided.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        #
        if self._x_set == 1:
            # mixture mole fraction given
            den = Mixture.density(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._wt,
                mode="mole",
            )
            return den
        elif self._y_set == 1:
            # mixture mass fraction given
            den = Mixture.density(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._wt,
                mode="mass",
            )
            return den
        else:
            # mixture composition is not provided
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @staticmethod
    def mixture_specific_heat(
        chemid: int,
        p: float,
        t: float,
        frac: npt.NDArray[np.double],
        wt: npt.NDArray[np.double],
        mode: str,
    ) -> float:
        """Get mixture specific heat capacity at constant pressure."""
        """Get mixture specific heat capacity from the given mixture condition:
        pressure, temperature, and species composition.

        Parameters
        ----------
            chemid: integer
                chemistry set index associated with the mixture
            p: double
                mixture pressure in [dynes/cm2]
            t: double
                mixture temperature in [K]
            frac: 1-D double array, dimension = number_species
                mixture composition given by either mass or mole fractions
                as specified by mode
            wt: 1-D double array, dimension = number_species
                molar masses of the species in the mixture in [gm/mol]
            mode: string, {'mole', 'mass'}
                flag indicates the frac array is 'mass' or 'mole' fractions

        Returns
        -------
            CpB: double
                mixture specific heat capacity [erg/mol-K]

        """
        # check inputs
        if chemid < 0:
            msg = [Color.PURPLE, "invalid chemistry.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if t <= 10.0:
            msg = [Color.PURPLE, "invalid temperature value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            msg = [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # initialization
        cpb_c = c_double(0.0)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.mole_fraction_to_mass_fraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            ierr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            msg = [
                Color.PURPLE,
                'must specify "mole" or "mass" fractions given.',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemid)
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # real-gas
        if check_realgas_status(chemid):
            # real-gas cubic EOS is active, set current pressure that is
            # required by the chemkin real-gas module
            set_current_pressure(chemid, p)
        # compute mass density from mass fraction
        ierr = ck_wrapper.chemkin.KINGetGasMixtureSpecificHeat(
            chemset_index, tt, yy, cpb_c
        )
        # compute mean molar mass
        mwt = Mixture.mean_molar_mass(frac=y, wt=wt, mode="mass")
        if ierr == 0:
            return cpb_c.value * mwt
        else:
            # failed to compute mixture specific heat
            msg = [
                Color.PURPLE,
                "failed to compute mixture specific heat capacity.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @staticmethod
    def mixture_enthalpy(
        chemid: int,
        p: float,
        t: float,
        frac: npt.NDArray[np.double],
        wt: npt.NDArray[np.double],
        mode: str,
    ) -> float:
        """Get mixture enthalpy."""
        """Get mixture enthalpy from the given mixture condition:
        pressure, temperature, and species composition.

        Parameters
        ----------
            chemid: integer
                chemistry set index associated with the mixture
            p: double
                mixture pressure in [dynes/cm2]
            t: double
                mixture temperature in [K]
            frac: 1-D double array, dimension = number_species
                mixture composition given by either mass or mole fractions
                as specified by mode
            wt: 1-D double array, dimension = number_species
                molar masses of the species in the mixture in [gm/mol]
            mode: string, {'mole', 'mass'}
                flag indicates the frac array is 'mass' or 'mole' fractions

        Returns
        -------
            H: double
                mixture enthalpy [erg/mol]

        """
        # check inputs
        if chemid < 0:
            msg = [Color.PURPLE, "invalid chemistry.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if t <= 10.0:
            msg = [Color.PURPLE, "invalid temperature value.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            msg = [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # initialization
        h_c = c_double(0.0)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.mole_fraction_to_mass_fraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            ierr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            msg = [
                Color.PURPLE,
                'must specify "mole" or "mass" fractions given.',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemid)
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # real-gas
        if check_realgas_status(chemid):
            # real-gas cubic EOS is active, set current pressure that is required
            # by the chemkin real-gas module
            set_current_pressure(chemid, p)
        # compute enthalpy from mass fraction
        ierr = ck_wrapper.chemkin.KINGetGasMixtureEnthalpy(chemset_index, tt, yy, h_c)
        # compute mean molar mass
        mwt = Mixture.mean_molar_mass(frac=y, wt=wt, mode="mass")
        if ierr == 0:
            return h_c.value * mwt
        else:
            # failed to compute mixture enthalpy
            msg = [Color.PURPLE, "failed to compute mixture enthalpy.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @staticmethod
    def rate_of_production(
        chemid: int,
        p: float,
        t: float,
        frac: npt.NDArray[np.double],
        wt: npt.NDArray[np.double],
        mode: str,
    ) -> npt.NDArray[np.double]:
        """Get species molar rate of production."""
        """Get species molar rate of production from the given mixture condition:
        pressure, temperature, and species composition.

        Parameters
        ----------
            chemid: integer
                chemistry set index associated with the mixture
            p: double
                mixture pressure in [dynes/cm2]
            t: double
                mixture temperature in [K]
            frac: 1-D double array, dimension = number_species
                mixture composition given by either mass or mole fractions
                as specified by mode
            wt: 1-D double array, dimension = number_species
                molar masses of the species in the mixture in [gm/mol]
            mode: string, {'mole', 'mass'}
                flag indicates the frac array is 'mass' or 'mole' fractions

        Returns
        -------
            rop: 1-D double array, dimension = number_species
                species molar rate of production in [mol/cm3-sec]

        """
        # check inputs
        if chemid < 0:
            msg = [Color.PURPLE, "invalid chemistry.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if p <= 0.0 or (p * t) <= 0.0:
            msg = [
                Color.PURPLE,
                "invalid pressure and/or temperature value(s).",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            msg = [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # initialization
        rop = np.zeros(kgas, dtype=np.double)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.mole_fraction_to_mass_fraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            ierr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            msg = [
                Color.PURPLE,
                'must specify "mole" or "mass" fractions given.',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemid)
        pp = c_double(p)  # pressure scalar
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # compute mass density from mass fraction
        ierr = ck_wrapper.chemkin.KINGetGasROP(chemset_index, tt, pp, yy, rop)
        if ierr == 0:
            return rop
        else:
            # failed to compute species molar rates of production
            msg = [
                Color.PURPLE,
                "failed to compute species molar rates of production.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @staticmethod
    def reaction_rates(
        chemid: int,
        numbreaction: int,
        p: float,
        t: float,
        frac: npt.NDArray[np.double],
        wt: npt.NDArray[np.double],
        mode: str,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get the molar rates of the gas reactions."""
        """Get molar rates of the gas reactions from the given mixture condition:
        pressure, temperature, and species composition.

        Parameters
        ----------
            chemid: integer
                chemistry set index associated with the mixture
            numbreaction: integer
                number of gas reactions associated with the chemistry set
            p: double
                mixture pressure in [dynes/cm2]
            t: double
                mixture temperature in [K]
            frac: 1-D double array, dimension = number_species
                mixture composition given by either mass or mole fractions
                as specified by mode
            wt: 1-D double array, dimension = number_species
                molar masses of the species in the mixture in [gm/mol]
            mode: string, {'mole', 'mass'}
                flag indicates the frac array is 'mass' or 'mole' fractions

        Returns
        -------
            k_forward: 1-D double array, dimension = numbreaction
                forward molar rates of the reactions in [mol/cm3-sec]
            k_reverse: 1-D double array, dimension = numbreaction
                reverse molar rates of the reactions in [mol/cm3-sec]

        """
        # check inputs
        if chemid < 0:
            msg = [Color.PURPLE, "invalid chemistry.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if p <= 0.0 or (p * t) <= 0.0:
            msg = [
                Color.PURPLE,
                "invalid pressure and/or temperature value(s).",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # number species
        kgas = len(frac)
        if kgas != len(wt):
            msg = [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # initialization
        k_forward = np.zeros(numbreaction, dtype=np.double)
        k_reverse = np.zeros_like(k_forward, dtype=np.double)
        if mode.lower() == "mole":
            # convert mole fraction to mass fraction and normalize
            y = Mixture.mole_fraction_to_mass_fraction(molefrac=frac, wt=wt)
        elif mode.lower() == "mass":
            # normalize mass fractions
            ierr, y = Mixture.normalize(frac=frac)
        else:
            # fraction type not given or incorrect
            msg = [
                Color.PURPLE,
                'must specify "mole" or "mass" fractions given.',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chemid)
        pp = c_double(p)  # pressure scalar
        tt = c_double(t)  # temperature scalar
        yy = np.ctypeslib.as_array(y)  # mass fraction array
        # compute mass density from mass fraction
        ierr = ck_wrapper.chemkin.KINGetGasReactionRates(
            chemset_index, tt, pp, yy, k_forward, k_reverse
        )
        if ierr == 0:
            return k_forward, k_reverse
        else:
            # failed to compute reaction rates
            msg = [
                Color.PURPLE,
                "failed to compute reaction rates,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def find_equilibrium(self):
        """Create the equilibrium state mixture corresponding to mixture itself."""
        """
        Create the equilibrium state mixture corresponding to mixture itself with
        both pressure and temperature fioxed.

        Returns
        -------
            eq_state: Mixture object
                gas mixture at the equilibrium state

        """
        # initialization a Mixture object by duplication
        eq_state = copy.deepcopy(self)
        # reset mass/mole fractions
        eq_state._x_set = 0
        eq_state._molefrac[:] = 0.0e0
        eq_state._y_set = 0
        eq_state._massfrac[:] = 0.0e0
        # compute the equilibrium state (mass fraction for now)
        eqvars, eq_state._massfrac = calculate_equilibrium(
            self._chemset_index.value,
            p=eq_state.pressure,
            t=eq_state.temperature,
            frac=self.y,
            wt=self._wt,
            mode_in="mass",
            mode_out="mass",
        )
        if np.sum(eq_state._massfrac, dtype=np.double) > 0.0e0:
            eq_state._y_set = 1
        return eq_state

    def hml(self) -> float:
        """Get enthalpy of the mixture."""
        """
        Get enthalpy of the mixture.

        Returns
        -------
            hml: double
                mixture enthalpy [erg/mol]

        """
        # initialization
        hml = 0.0e0
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if self._x_set == 1:
            # mixture mole fraction given
            hml = Mixture.mixture_enthalpy(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._wt,
                mode="mole",
            )
            return hml
        elif self._y_set == 1:
            # mixture mass fraction given
            hml = Mixture.mixture_enthalpy(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._wt,
                mode="mass",
            )
            return hml
        else:
            # mixture composition is not provided
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def cpbl(self) -> float:
        """Get specific heat capacity of the mixture."""
        """
        Get specific heat capacity of the mixture.

        Returns
        -------
            cpbl: double
                mixture specific heat capacity [erg/mol-K]

        """
        # initialization
        cpbl = 0.0e0
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        #
        if self._x_set == 1:
            # mixture mole fraction given
            cpbl = Mixture.mixture_specific_heat(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._wt,
                mode="mole",
            )
            return cpbl
        elif self._y_set == 1:
            # mixture mass fraction given
            cpbl = Mixture.mixture_specific_heat(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._wt,
                mode="mass",
            )
            return cpbl
        else:
            # mixture composition is not provided
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def rop(self) -> npt.NDArray[np.double]:
        """Get species rate of productions."""
        """Get species molar rate of production from the given mixture condition:
        pressure, temperature, and species compositions.

        Returns
        -------
            rop: 1-D double array, dimension = number_species
                species molar rate of production in [mol/cm3-sec]

        """
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check pressure
        if self._p_set == 0:
            msg = [
                Color.PURPLE,
                "mixture pressure [dynes/cm2] is not provided.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        #
        if self._x_set == 1:
            # mixture mole fraction given
            rop = Mixture.rate_of_production(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._wt,
                mode="mole",
            )
            return rop
        elif self._y_set == 1:
            # mixture mass fraction given
            rop = Mixture.rate_of_production(
                self._chemset_index.value,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._wt,
                mode="mass",
            )
            return rop
        else:
            # mixture composition is not provided
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def rxn_rates(self) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get molar reaction rates."""
        """Get molar rates of the gas reactions from the given mixture condition:
        pressure, temperature, and species composition.

        Returns
        -------
            k_forward: 1-D double array, dimension = numbreaction
                forward molar rates of the reactions in [mol/cm3-sec]
            k_reverse: 1-D double array, dimension = numbreaction
                reverse molar rates of the reactions in [mol/cm3-sec]

        """
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check pressure
        if self._p_set == 0:
            msg = [
                Color.PURPLE,
                "mixture pressure [dynes/cm2] is not provided.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # initialization
        k_forward = np.zeros(self._ii_gas, dtype=np.double)
        k_reverse = np.zeros_like(k_forward, dtype=np.double)
        #
        if self._x_set == 1:
            # mixture mole fraction given
            k_forward, k_reverse = Mixture.reaction_rates(
                self._chemset_index.value,
                numbreaction=self._ii_gas,
                p=self._press,
                t=self._temp,
                frac=self._molefrac,
                wt=self._wt,
                mode="mole",
            )
            return k_forward, k_reverse
        elif self._y_set == 1:
            # mixture mass fraction given
            k_forward, k_reverse = Mixture.reaction_rates(
                self._chemset_index.value,
                numbreaction=self._ii_gas,
                p=self._press,
                t=self._temp,
                frac=self._massfrac,
                wt=self._wt,
                mode="mass",
            )
            return k_forward, k_reverse
        else:
            # mixture composition is not provided
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def species_cp(self) -> npt.NDArray[np.double]:
        """Get species specific heat capacity at constant pressure."""
        """
        Get species specific heat capacity at constant pressure.

        Returns
        -------
            Cp: 1-D double array, dimension = number_species
                species specific heat capacities at constant pressure [ergs/mol-K]

        """
        tt = c_double(self.temperature)
        cp = np.zeros(self._kk, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetGasSpecificHeat(self._chemset_index, tt, cp)
        if ierr == 0:
            # convert [ergs/g-K] to [ergs/mol-K]
            cp *= self._wt
        else:
            # failed to compute specific heats
            msg = [
                Color.PURPLE,
                "failed to compute species specific heat capacities.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return cp

    def species_h(self) -> npt.NDArray[np.double]:
        """Get species enthalpy."""
        """
        Get species enthalpy.

        Returns
        -------
            h: 1-D double array, dimension = number_species
                species enthalpy [ergs/mol]

        """
        tt = c_double(self.temperature)
        h = np.zeros(self._kk, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetGasSpeciesEnthalpy(self._chemset_index, tt, h)
        if ierr == 0:
            # convert [ergs/gm] to [ergs/mol]
            h *= self._wt
        else:
            # failed to compute enthalpies
            msg = [Color.PURPLE, "failed to compute species enthalpies.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return h

    def species_visc(self) -> npt.NDArray[np.double]:
        """Get species viscosity."""
        """
        Get species viscosity.

        Returns
        -------
            visc: : 1-D double array, dimension = number_species
                species viscosity [gm/cm-sec]

        """
        if self.transport_data != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        tt = c_double(self.temperature)
        visc = np.zeros(self._kk, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetViscosity(self._chemset_index, tt, visc)
        if ierr != 0:
            # failed to compute viscosities
            msg = [Color.PURPLE, "failed to compute species viscosities.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return visc

    def species_cond(self) -> npt.NDArray[np.double]:
        """Get species conductivity."""
        """
        Get species conductivity.

        Returns
        -------
            cond: 1-D double array, dimension = number_species
                species conductivity [ergs/cm-K-sec]

        """
        if self.transport_data != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        tt = c_double(self.temperature)
        cond = np.zeros(self._kk, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetConductivity(self._chemset_index, tt, cond)
        if ierr != 0:
            # failed to compute conductivities
            msg = [Color.PURPLE, "failed to compute species conductivities.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return cond

    def species_diffusion_coeffs(self) -> npt.NDArray[np.double]:
        """Get species diffusion coefficients."""
        """
        Get species diffusion coefficients.

        Returns
        -------
            diffusioncoeffs: 2-D double array,
            dimension = [number_species, number_species]
                species diffusion coefficients [cm2/sec]

        """
        if self.transport_data != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        pp = c_double(self.pressure)
        tt = c_double(self.temperature)
        dim = (self._kk, self._kk)
        diffusioncoeffs = np.zeros(dim, dtype=np.double, order="F")
        ierr = ck_wrapper.chemkin.KINGetDiffusionCoeffs(
            self._chemset_index, pp, tt, diffusioncoeffs
        )
        if ierr != 0:
            # failed to compute diffusion coefficients
            msg = [
                Color.PURPLE,
                "failed to compute species diffusion coefficients.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return diffusioncoeffs

    def mixture_viscosity(self) -> float:
        """Get viscosity of the gas mixture."""
        """
        Get viscosity of the gas mixture.

        Returns
        -------
            visc: double
                mixture viscosity [gm/cm-sec]

        """
        if self.transport_data != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # initialization
        visc = c_double(0.0e0)
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get mixture viscosity
        tt = c_double(self.temperature)
        ierr = ck_wrapper.chemkin.KINGetMixtureViscosity(
            self._chemset_index, tt, self.y, visc
        )
        if ierr != 0:
            # error message
            msg = [Color.PURPLE, "failed to compute mixture viscosity.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # mixture viscosity in gm/cm-sec
        return visc.value

    def mixture_conductivity(self) -> float:
        """Get conductivity of the gas mixture."""
        """
        Get conductivity of the gas mixture.

        Returns
        -------
            cond: double
                mixture conductivity [erg/cm-K-sec]

        """
        if self.transport_data != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # initialization
        cond = c_double(0.0e0)
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get mixture viscosity
        tt = c_double(self.temperature)
        ierr = ck_wrapper.chemkin.KINGetMixtureConductivity(
            self._chemset_index, tt, self.y, cond
        )
        if ierr != 0:
            # error message
            msg = [Color.PURPLE, "failed to compute mixture conductivity.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # mixture conductivity in ergs/cm-K-sec
        return cond.value

    def mixture_diffusion_coeffs(self) -> npt.NDArray[np.double]:
        """Get mixture-averaged species diffusion coefficients of the gas mixture."""
        """
        Get mixture-averaged species diffusion coefficients of the gas mixture.

        Returns
        -------
            diffusioncoeffs: 1-D double array, dimension = number_species
                mixture-averaged diffusion coefficients [cm2/sec]

        """
        if self.transport_data != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # initialization
        diffusioncoeffs = np.zeros(self._kk, dtype=np.double)
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check pressure
        if self._p_set == 0:
            msg = [
                Color.PURPLE,
                "mixture pressure [dynes/cm2] is not provided.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get mixture viscosity
        pp = c_double(self.pressure)
        tt = c_double(self.temperature)
        ierr = ck_wrapper.chemkin.KINGetMixtureDiffusionCoeffs(
            self._chemset_index, pp, tt, self.y, diffusioncoeffs
        )
        if ierr != 0:
            # error message
            msg = [
                Color.PURPLE,
                "failed to compute",
                "mixture-averaged diffusion coefficients.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        # mixture-averaged diffusion coefficients in cm2/sec
        return diffusioncoeffs

    def mixture_binary_diffusion_coeffs(self) -> npt.NDArray[np.double]:
        """Get multi-component species binary diffusion coefficients."""
        """
        Get multi-component species binary diffusion coefficients of
        the gas mixture.

        Returns
        -------
            binarydiffusioncoeffs: 2-D double array,
            dimension = [number_species, number_species]
                binary diffusion coefficients [cm2/sec]

        """
        if self.transport_data != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # initialization
        dim = (self._kk, self._kk)
        binarydiffusioncoeffs = np.zeros(dim, dtype=np.double, order="F")
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check pressure
        if self._p_set == 0:
            msg = [
                Color.PURPLE,
                "mixture pressure [dynes/cm2] is not provided.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get mixture viscosity
        pp = c_double(self.pressure)
        tt = c_double(self.temperature)
        ierr = ck_wrapper.chemkin.KINGetOrdinaryDiffusionCoeffs(
            self._chemset_index, pp, tt, self.y, binarydiffusioncoeffs
        )
        if ierr != 0:
            # error message
            msg = [
                Color.PURPLE,
                "failed to compute",
                "multi-component binary diffusion coefficients.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # mixture multi-component binary diffusion coefficients in cm2/sec
        return binarydiffusioncoeffs

    def mixture_thermal_diffusion_coeffs(self) -> npt.NDArray[np.double]:
        """Get thermal diffusivity of the gas mixture."""
        """
        Get thermal diffusivity of the gas mixture.

        Returns
        -------
            thermaldiffusioncoeffs: 1-D double array,
            dimension = number_species
                thermal diffusivity [gm/cm-sec]

        """
        if self.transport_data != 1:
            msg = [Color.PURPLE, "no transport data processed.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # initialization
        thermaldiffusioncoeffs = np.zeros(self._kk, dtype=np.double)
        cond = c_double(0.0e0)  # mixture thermal conductivity
        # check temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature [K] is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check pressure
        if self._p_set == 0:
            msg = [
                Color.PURPLE,
                "mixture pressure [dynes/cm2] is not provided.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get mixture viscosity
        pp = c_double(self.pressure)
        tt = c_double(self.temperature)
        ierr = ck_wrapper.chemkin.KINGetThermalDiffusionCoeffs(
            self._chemset_index, pp, tt, self.y, thermaldiffusioncoeffs, cond
        )
        if ierr != 0:
            # error message
            msg = [
                Color.PURPLE,
                "failed to compute",
                "mixture thermal diffusion coefficients.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # mixture thermal diffusion coefficients in gm/cm-sec
        return thermaldiffusioncoeffs

    def vol_hrr(self) -> float:
        """Get volumetric heat release rate."""
        """
        Get volumetric heat release rate.

        Returns
        -------
            vol_HRR: double
                volumetric heat release rate [ergs/cm3-sec]

        """
        vol_hrr = 0.0e0
        # get species enthalpy
        tt = c_double(self.temperature)
        h = np.zeros(self._kk, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetGasSpeciesEnthalpy(self._chemset_index, tt, h)
        if ierr == 0:
            # convert H from ergs/gm to ergs/mol
            h *= self._wt
        else:
            msg = [
                Color.PURPLE,
                "failed to compute volumetric heatrelease rate.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get species molar rate of production mol/cm3-sec
        rop = self.rop()
        # volumetric heat release rate = SUM(H_k * ROP_k)  ergs/cm3-sec
        vol_hrr = np.dot(h, rop)
        return vol_hrr

    def rop_mass(self) -> npt.NDArray[np.double]:
        """Get species mass rates of production."""
        """
        Get species mass rates of production.

        Returns
        -------
            rop_mass: 1-D double array, dimension = number_species
                mass rates of production [gm/cm3-sec]

        """
        # get species molar rate of production mol/cm3-sec
        rop = self.rop()
        # species mass rate of production = ROP_k * WT_k
        rop_mass = rop * self._wt
        return rop_mass

    def list_rop(
        self, threshold: float = 0.0
    ) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.double]]:
        """List species molar production rate."""
        """List information about species molar production rate
        in descending order.

        Parameters
        ----------
            threshold: double, optional, default = 0.0
                minimum absolute ROP value to be printed

        Returns
        -------
            order: 1-D integer array, dimension = number_species
                sorted species index
            sorted_rop: 1-D double array, dimension = number_species
                sorted ROP values [gm/cm3-sec]

        """
        # get species molar rate of production mol/cm3-sec
        rop = self.rop()
        # create a copy of non-zero rates
        temp_rop = np.zeros_like(rop, dtype=np.double)
        temp_order = np.zeros_like(rop, dtype=np.int32)
        # include reactions with non-zero rate only
        count = 0
        for i in range(len(rop)):
            if abs(rop[i]) > threshold:
                # non-zero entries
                temp_rop[count] = rop[i]
                temp_order[count] = i
                count += 1

        # sort on the temporary array in descending order
        sorted_rop = np.flip(np.sort(temp_rop[:count]))
        # find species index
        new_order = np.zeros_like(sorted_rop, dtype=np.int32)
        for i in range(len(sorted_rop)):
            count, species_id = where_element_in_array_1d(temp_rop, sorted_rop[i])
            # get the species index corresponding to the first occurrence
            # in case there are multiple species having the same rate
            new_order[i] = temp_order[species_id[0]]
        # print out the list of species with its ROP value in descending order
        if verbose():
            print("non-zero species molar rate of production ")
            print("=" * 50)
            print(" order    species symbol     rate of production")
            print("                             [mol/cm3-sec]")
            for i in range(len(new_order)):
                print(
                    f" {i + 1:-2d} {self._specieslist[new_order[i]]:>16}"
                    f"              {sorted_rop[i]: e}"
                )
        return new_order, sorted_rop

    def list_rop_mass(
        self, threshold: float = 0.0
    ) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.double]]:
        """List species mass rate of production in descending order."""
        """List information about species mass rate of production
        in descending order.

        Parameters
        ----------
            threshold: double, optional, default = 0.0
                minimum absolute mass ROP value to be printed

        Returns
        -------
            order: 1-D integer array, dimension = number_species
                sorted species index
            sorted_rop_mass: 1-D double array, dimension = number_species
                sorted mass ROP values [gm/cm3-sec]

        """
        # get species mass rate of production gm/cm3-sec
        rop_mass = self.rop_mass()
        # create a copy of non-zero rates
        temp_rop = np.zeros_like(rop_mass, dtype=np.double)
        temp_order = np.zeros_like(rop_mass, dtype=np.int32)
        # include reactions with non-zero rate only
        count = 0
        for i in range(len(rop_mass)):
            if abs(rop_mass[i]) > threshold:
                # non-zero entries
                temp_rop[count] = rop_mass[i]
                temp_order[count] = i
                count += 1

        # sort on the temporary array in descending order
        sorted_rop = np.flip(np.sort(temp_rop[:count]))
        # find species index
        new_order = np.zeros_like(sorted_rop, dtype=np.int32)
        for i in range(len(sorted_rop)):
            count, species_id = where_element_in_array_1d(temp_rop, sorted_rop[i])
            # get the species index corresponding to the first occurrence
            # in case there are multiple species having the same rate
            new_order[i] = temp_order[species_id[0]]
        # print out the list of species with its ROP value in descending order
        if verbose():
            print("non-zero species mass rate of production ")
            print("=" * 50)
            print(" order    species symbol     rate of production")
            print("                             [gm/cm3-sec]")
            for i in range(len(new_order)):
                print(
                    f" {i + 1:-2d} {self._specieslist[new_order[i]]:>16}"
                    f"              {sorted_rop[i]: e}"
                )
        return new_order, sorted_rop

    def list_reaction_rates(
        self, threshold: float = 0.0
    ) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.double]]:
        """List information about reaction rate in descending order."""
        """
        Parameters
        ----------
            threshold: double, optional, default = 0.0
                minimum absolute reaction rate value to be printed

        Returns
        -------
            order: 1-D integer array, dimension = numb_reactions
                sorted reaction index
            sorted_rxn_rates: 1-D double array, dimension = numb_reactions
                sorted reaction rate values [mol/cm3-sec]

        """
        # molar rates of reactions
        rf, rr = self.rxn_rates()
        # create a copy of non-zero rates
        temp_net_rr = np.zeros_like(rf, dtype=np.double)
        temp_order = np.zeros_like(rr, dtype=np.int32)
        # include reactions with non-zero rate only
        count = 0
        for i in range(len(rf)):
            net_rr = rf[i] - rr[i]
            if abs(net_rr) > threshold:
                # non-zero entries
                temp_net_rr[count] = net_rr
                temp_order[count] = i
                count += 1

        # sort on the temporary array in descending order
        sorted_rr = np.flip(np.sort(temp_net_rr[:count]))
        # find species index
        new_order = np.zeros_like(sorted_rr, dtype=np.int32)
        new_rr = copy.deepcopy(sorted_rr)
        for i in range(len(sorted_rr)):
            # find the instances of this reaction rate
            count, rxn_id = where_element_in_array_1d(temp_net_rr, new_rr[i])
            # get the reaction number corresponding to the first occurrence
            # in case there are multiple reactions having the same rate
            new_order[i] = temp_order[rxn_id[0]]
            # remove this instance from the reaction rate array
            new_rr[i] = 0.0e0
        # print out the list of reaction with its net reaction rate value
        # in descending order
        if verbose():
            print("non-zero molar rates of reaction ")
            print("=" * 50)
            print(" order    reaction number    molar rate of reaction")
            print("                             [mol/cm3-sec]")
            for i in range(len(new_order)):
                print(
                    f" {i + 1:-2d}          {new_order[i] + 1:-4d}"
                    f"              {sorted_rr[i]: e}"
                )
        return new_order, sorted_rr

    def x_by_equivalence_ratio(
        self,
        chemistryset: Chemistry,
        fuel_molefrac: npt.NDArray[np.double],
        oxid_molefrac: npt.NDArray[np.double],
        add_molefrac: npt.NDArray[np.double],
        products: list[str],
        equivalenceratio: float,
        threshold: float = 1.0e-10,
    ) -> int:
        """Set mole fractions using equivalence ratio."""
        """Specify the mixture molar composition by providing the equivalence ratio,
        the mole fractions of the fuel mixture, the oxidizer mixture, and
        the additives mixture, and the list of the complete combustion product species.

        Parameters
        ----------
            chemistryset: Chemistry object
                the chemistry set used to create the mixtures
            fuel_molefrac: 1-D double array, dimension = number_species
                mole fractions of the fuel mixture
            oxid_molefrac: 1-D double array, dimension = number_species
                mole fractions of the oxidizer mixture
            add_molefrac: 1-D double array, dimension = number_species
                mole fractions of the additives mixture
            products: list of string
                list of the complete combustion species symbols
            equivalenceratio: double
                equivalence ratio of the final mixture (double scalar)
            threshold: double, optional, default = 1.0e-10
                minimum species fraction value to be included in
                the stoichiometric coefficient calculation

        Returns
        -------
            Error status: integer

        """
        # check chemistry set
        if not isinstance(chemistryset, Chemistry):
            msg = [
                Color.PURPLE,
                "the first argument must be a Chemistry object.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 1
        # number of gas species in the mechanism
        kspecies = chemistryset.kk
        # find fuel mole array size
        kfuel = len(fuel_molefrac)
        # find oxidizer mole array size
        koxid = len(oxid_molefrac)
        # find additives mole array size
        kadd = len(add_molefrac)
        # check species number consistency
        ierr = 0
        if kspecies != kfuel:
            msg = [
                Color.PURPLE,
                "the fuel mole fraction array must have size",
                str(kspecies),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr += 1
        if kspecies != koxid:
            msg = [
                Color.PURPLE,
                "the oxidizer mole fraction array must have size",
                str(kspecies),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr += 1
        if kspecies != kadd:
            msg = [
                Color.PURPLE,
                "the additive mole fraction array must have size",
                str(kspecies),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr += 1
        if ierr > 0:
            return 2
        # check equivalence ratio value
        if equivalenceratio <= 0.0e0:
            msg = [Color.PURPLE, "the equivalence ratio must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 3
        # check product species
        kprod = len(products)
        if kprod == 0:
            msg = [
                Color.PURPLE,
                "complete combustion products must be provided.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
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
            prod_index[j] = chemistryset.get_specindex(s)
            j += 1
        # find the stoichiometric coefficients assuming complete combustion
        alpha, nu = calculate_stoichiometrics(
            chemistryset, fuel_molefrac, oxid_molefrac, prod_index
        )
        if alpha <= 0.0e0 or nu[0] == 0:
            msg = [
                Color.PURPLE,
                "failed to find the stoichiometric coefficients.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
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
            self._x_set = 1
            self._massfrac[:] = 0.0e0
            self._y_set = 0
            return 0
        else:
            msg = [
                Color.PURPLE,
                "failed to find the stoichiometric coefficients.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 6

    def y_by_equivalence_ratio(
        self,
        chemistryset: Chemistry,
        fuel_massfrac: npt.NDArray[np.double],
        oxid_massfrac: npt.NDArray[np.double],
        add_massfrac: npt.NDArray[np.double],
        products: list[str],
        equivalenceratio: float,
        threshold: float = 1.0e-10,
    ) -> int:
        """Set mass fractions using equivalence ratio."""
        """Specify the mixture molar composition by providing the equivalence ratio,
        the mole fractions of the fuel mixture, the oxidizer mixture, and
        the additives mixture, and the list of the complete combustion product species.

        Parameters
        ----------
            chemistryset: Chemistry object
                the chemistry set used to create the mixtures
            fuel_massfrac: 1-D double array, dimension = number_species
                mass fractions of the fuel mixture
            oxid_massfrac: 1-D double array, dimension = number_species
                mass fractions of the oxidizer mixture
            add_massfrac: 1-D double array, dimension = number_species
                mass fractions of the additives mixture
            products: list of string
                list of the complete combustion species symbols
            equivalenceratio: double
                equivalence ratio of the final mixture (double scalar)
            threshold: double, optional, default = 1.0e-10
                minimum species fraction value to be included in
                he stoichiometric coefficient calculation

        Returns
        -------
            Error status: integer

        """
        # check chemistry set
        if not isinstance(chemistryset, Chemistry):
            msg = [
                Color.PURPLE,
                "the first argument must be a Chemistry object.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert mass fractions to mole fractions
        fuel_molefrac = Mixture.mass_fraction_to_mole_fraction(
            massfrac=fuel_massfrac, wt=chemistryset.WT
        )
        oxid_molefrac = Mixture.mass_fraction_to_mole_fraction(
            massfrac=oxid_massfrac, wt=chemistryset.WT
        )
        add_molefrac = Mixture.mass_fraction_to_mole_fraction(
            massfrac=add_massfrac, wt=chemistryset.WT
        )
        # find the final mixture mole fractions and set the flags
        ierr = self.x_by_equivalence_ratio(
            chemistryset,
            fuel_molefrac,
            oxid_molefrac,
            add_molefrac,
            products,
            equivalenceratio,
            threshold,
        )
        return ierr

    def get_egr_mole_fraction(
        self, egr_ratio: float, threshold: float = 1.0e-8
    ) -> npt.NDArray[np.double]:
        """Compute the EGR composition in mole fraction."""
        """
        Compute the EGR composition in mole fraction corresponding to
        this unburned mixture.

        Parameters
        ----------
            egr_ratio: double
                exhaust gas recirculation (EGR) molar ratio
            threshold: double, optional, default = 1.0e-8
                minimum species fraaction value to be included in the EGR stream

        Returns
        -------
            egr_molefrac: 1-D double array, dimension = number_species
                EGR stream compostion in mole fractions

        """
        # create burned mixture
        burned = self.find_equilibrium()
        # compute the EGR stream mole fractions
        klength = len(burned.x)
        egr_molefrac = np.zeros(klength, dtype=np.double)
        for i in range(klength):
            if burned.x[i] > threshold:
                egr_molefrac[i] = egr_ratio * burned.x[i]
        del burned
        return egr_molefrac

    def validate(self) -> int:
        """Check whether the mixture is fully defined."""
        """
        Check whether the mixture is fully defined before being
        used by other methods.

        Returns
        -------
            Error status: integer

        """
        ierr = 0
        # check mixture temperature
        if self._t_set == 0:
            msg = [Color.PURPLE, "mixture temperature is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr = 1
        if self._p_set == 0:
            msg = [Color.PURPLE, "mixture pressure is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr = 2
        if self._x_set == 0 and self._y_set == 0:
            msg = [Color.PURPLE, "mixture composition is not provided.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr = 3
        return ierr

    def use_realgas_cubic_eos(self):
        """Turn ON the real-gas cubic EOS."""
        """Turn ON the real-gas cubic EOS to compute mixture properties
        if the mechanism contains necessary data.
        """
        if self._eos.value < 1:
            # no real gas EOS data in the mechanism
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return
        # check real-gas EOS status
        iflag = c_int(0)
        ierr = ck_wrapper.chemkin.KINRealGas_UseCubicEOS(self._chemset_index, iflag)
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to activate the real-gas EOS model",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if iflag.value == 0:
            msg = [
                Color.YELLOW,
                "real-gas cubic EOS model",
                Chemistry.realgas_cueos[self._eos.value],
                "is turned ON.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            # set default mixing rule to Van der Waals
            mixingrule = 0
            # set default mixing rule to Van der Waals
            self.set_realgas_mixing_rule(mixingrule)
            self.userealgas = True
        else:
            self.userealgas = False

    def use_idealgas_law(self):
        """Turn on the ideal gas law to compute mixture properties."""
        if self._eos.value < 1:
            # no real gas EOS data in the mechanism
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = False
            return
        # check real-gas EOS status
        iflag = c_int(0)
        ierr = ck_wrapper.chemkin.KINRealGas_UseIdealGasLaw(self._chemset_index, iflag)
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to activate the ideal gas law,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if iflag.value == 0:
            msg = [Color.YELLOW, "the ideal gas law is turned ON.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = False

    def set_realgas_mixing_rule(self, rule: int = 0):
        """Set the mixing rule for calculating the real-gas mixture properties."""
        """Set the mixing rule to be used for calculating
        the real-gas mixture properties.

        Parameters
        ----------
            rule: integer, optional, default = 0
                mixing rule:
                    0 for the Van der Waals mixing rule;
                    1 for the critical properties mixing rule (integer scalar)

        """
        if self._eos.value < 1:
            # no real gas EOS data in the mechanism
            msg = [Color.YELLOW, "mechanism is for ideal gas law only.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            return
        # set default mixing rule to Van der Waals
        mixingrule = c_int(rule + 1)
        iflag = c_int(0)
        ierr = ck_wrapper.chemkin.KINRealGas_SetMixingRule(
            self._chemset_index, mixingrule, iflag
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to activate the real-gas mixing rule,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if iflag.value == 2:
            # real-gas cubic EOS is turned OFF
            msg = [Color.YELLOW, "the ideal gas law is in use.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = False
        elif iflag.value != 0:
            msg = [
                Color.PURPLE,
                "fail to set up the real-gas mixing rule,",
                "error code =",
                str(iflag.value),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            msg = [
                Color.YELLOW,
                "the real-gas cubic EOS is activated,",
                "set the mixing rule to",
                '"' + Chemistry.realgas_mixing_rules[rule] + '"',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            self.userealgas = True


# mixture mixing
def isothermal_mixing(
    recipe: list[tuple[Mixture, float]], mode: str, finaltemperature: float
) -> Mixture:
    """Mixing multiple gas mixtures at gioven temperature."""
    """Find the resulting gas mixture properties from mixing a number of
    gas mixtures at the given mixture temperature.

    Parameters
    ----------
        recipe: list of tuples, [(Mixture object, fraction), ... ]
            non-zero mixture composition corresponding to the given mole/mass
            fraction array
        mode: string, {'mass', 'mole'}, default = 'mole'
            indicting the fractions given in the recipe are in 'mole' or 'mass'
            ratios
        finaltemperature: double
            temperature of the resulting gas mixture after mixing

    Returns
    -------
        finalmixture: Mixture object
            the resulting gas mixture after mixing

    """
    # check number of mixtures
    numb_mixture = len(recipe)
    numb_species = 0
    chem_index_check = -1
    # create the final mixture object
    finalmixture = copy.deepcopy(recipe[0][0])
    if numb_mixture == 0:
        # nothing there
        msg = [
            Color.PURPLE,
            "the mixing recipe is empty.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check types
    if isinstance(recipe[0][0], Mixture):
        numb_species = recipe[0][0]._kk
        # reset the compositions
        finalmixture._y_set = 0
        finalmixture._massfrac[:] = 0.0e0
        finalmixture._x_set = 0
        finalmixture._molefrac[:] = 0.0e0
        # set the chemistry set index
        chem_index_check = finalmixture.chemid
        if chem_index_check < 0:
            msg = [
                Color.PURPLE,
                "Mixture object #0 is not associated to any chemistry set.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check if the chemstry set is active
        if not check_active_chemistryset(chem_index_check):
            msg = [
                Color.PURPLE,
                "the Chemistry Set associated with",
                "the Mixture is not currently active.\n",
                Color.SPACEx6,
                "activate Chemistry Set using the 'active()' method.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check chemistry sets
        count = 0
        ierrs = 0
        for r in recipe:
            if r[0].chemid != chem_index_check:
                msg = [
                    Color.PURPLE,
                    "Mixture #",
                    str(count),
                    "belongs to a different Chemistry Set.\n",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierrs += 1
            count += 1
        if ierrs != 0:
            exit()
    else:
        # incorrect object type
        # delete the finalmixture object
        del finalmixture
        msg = [
            Color.PURPLE,
            "the first component of the recipe tuple",
            "must be a Chemkin Mixture object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check given final mixture temperature
    if finaltemperature <= 10.0:
        # final mixture temperature is not provided
        # delete the finalmixture object
        del finalmixture
        msg = [
            Color.PURPLE,
            "temperature of the final mixture must be provided.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()

    # initialization
    mixfrac = np.zeros(numb_mixture, dtype=np.double)
    mixfrac_sum = 0.0e0
    count = 0
    for s, v in recipe:
        # check object type
        if not isinstance(s, Mixture):
            msg = [
                Color.PURPLE,
                "the first component of the recipe tuple",
                "must be a Chemkin Mixture object.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check ratio value
        if v <= 0.0e0:
            msg = [
                Color.PURPLE,
                "the second component of the recipe tuple",
                "must be a positive number for the mole/mass ratio.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check chemistry set consistency
        if s.chemid != chem_index_check:
            msg = [
                Color.PURPLE,
                "Mixture #",
                str(count),
                "has inconsistent Chemistry setup.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get the mixture's mean molar mass g/mol
        mwt = s.wtm
        # find the composition of the final mixture
        speciesfrac_sum = 0.0e0
        if mode.lower() == "mole":
            # molar ratios are given
            # mole ratio
            mixfrac[count] = v
            # compute the composition of the final mixture
            for k in range(numb_species):
                finalmixture._molefrac[k] += s.x[k] * v
                speciesfrac_sum += finalmixture._molefrac[k]
            if speciesfrac_sum > 0.0e0:
                finalmixture._x_set = 1
        else:
            # assume mass ratios are given by default
            # mass ratio
            mixfrac[count] = v / mwt
            # compute the composition of the final mixture
            for k in range(numb_species):
                finalmixture._molefrac[k] += s.x[k] * mixfrac[count]
                speciesfrac_sum += finalmixture._molefrac[k]
            if speciesfrac_sum > 0.0e0:
                finalmixture._x_set = 1
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


def adiabatic_mixing(recipe: list[tuple[Mixture, float]], mode: str) -> Mixture:
    """Mixing multiple gas mixtures adiabatically."""
    """Find the resulting gas mixture properties from mixing a number of gas mixtures
    with constant total enthalpy.

    Parameters
    ----------
        recipe: list of tuples, [(Mixture object, fraction), ... ]
            non-zero mixture composition corresponding to the given mole/mass
            fraction array
        mode: string, {'mass', 'mole'}, default = 'mole'
            indicting the fractions given in the recipe are in 'mole' or 'mass'
            ratios

    Returns
    -------
        finalmixture: Mixture object
            the resulting gas mixture after mixing

    """
    # check number of mixtures
    numb_mixture = len(recipe)
    numb_species = 0
    chem_index_check = -1
    # create the final mixture object
    finalmixture = copy.deepcopy(recipe[0][0])
    if numb_mixture == 0:
        # nothing there
        msg = [Color.PURPLE, "the mixing recipe is empty.", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check types
    if isinstance(recipe[0][0], Mixture):
        numb_species = recipe[0][0]._kk
        # reset the compositions
        finalmixture._y_set = 0
        finalmixture._massfrac[:] = 0.0e0
        finalmixture._x_set = 0
        finalmixture._molefrac[:] = 0.0e0
        # set the chemistry set index
        chem_index_check = finalmixture.chemid
        if chem_index_check < 0:
            msg = [
                Color.PURPLE,
                "Mixture object #0 is not associated to any chemistry set.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check if the chemstry set is active
        if not check_active_chemistryset(chem_index_check):
            msg = [
                Color.PURPLE,
                "the Chemistry Set associated with",
                "the Mixture is not currently active.\n",
                Color.SPACEx6,
                "activate Chemistry Set using the 'active()' method.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check chemistry sets
        count = 0
        ierrs = 0
        for r in recipe:
            if r[0].chemid != chem_index_check:
                msg = [
                    Color.PURPLE,
                    "Mixture #",
                    str(count),
                    "belongs to a different Chemistry Set.\n",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierrs += 1
            count += 1
        if ierrs != 0:
            exit()
    else:
        # incorrect object type
        # delete the finalmixture object
        del finalmixture
        msg = [
            Color.PURPLE,
            "the first component of the recipe tuple",
            "must be a Chemkin Mixture object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()

    # initialization
    mixfrac = np.zeros(numb_mixture, dtype=np.double)
    mixfrac_sum = 0.0e0
    mix_h = 0.0e0
    count = 0
    for s, v in recipe:
        # check object type
        if not isinstance(s, Mixture):
            msg = [
                Color.PURPLE,
                "the first component of the recipe tuple",
                "must be a Chemkin Mixture object.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check ratio value
        if v <= 0.0e0:
            msg = [
                Color.PURPLE,
                "the second component of the recipe tuple",
                "must be a positive number for the mole/mass ratio.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check chemistry set consistency
        if s.chemid != chem_index_check:
            msg = [
                Color.PURPLE,
                "Mixture #",
                str(count),
                "has inconsistent Chemistry setup.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # get the mixture's mean molar mass g/mol
        mwt = s.wtm
        # find the composition of the final mixture
        speciesfrac_sum = 0.0e0
        if mode.lower() == "mole":
            # molar ratios are given
            # mole ratio
            mixfrac[count] = v
            # compute the composition of the final mixture
            for k in range(numb_species):
                finalmixture._molefrac[k] += s.x[k] * v
                speciesfrac_sum += finalmixture._molefrac[k]
            if speciesfrac_sum > 0.0e0:
                finalmixture._x_set = 1
        else:
            # assume mass ratios are given by default
            # mass ratio
            mixfrac[count] = v / mwt
            # compute the composition of the final mixture
            for k in range(numb_species):
                finalmixture._molefrac[k] += s.x[k] * mixfrac[count]
                speciesfrac_sum += finalmixture._molefrac[k]
            if speciesfrac_sum > 0.0e0:
                finalmixture._x_set = 1

        # compute the final mixture's enthalpy ergs/mol
        mix_h += s.hml() * mixfrac[count]
        count += 1

    # normalize the mixing mole ratios
    mixfrac_sum = np.sum(mixfrac)
    mixfrac /= mixfrac_sum
    # normalize the mole fractions of the final mixture
    finalmixture._molefrac /= mixfrac_sum
    # normalize the total mixture enthalpy ergs/mol
    # (= the enthalpy of the final mixture)
    mix_h /= mixfrac_sum
    # compute temperature of the final mixture from the mixture enthalpy
    # set the guessed temperature
    t_guessed = 0.0e0
    ierr = cal_mixture_temperature_from_enthalpy(
        mixture=finalmixture, h_mixture=mix_h, guesstemperature=t_guessed
    )
    if ierr != 0:
        msg = [
            Color.PURPLE,
            "failed to compute the final mixture temperature,",
            "error code =",
            str(ierr),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    if verbose():
        print(f"final mixture temperature = {finalmixture.temperature}[K]")
    return finalmixture


def cal_mixture_temperature_from_enthalpy(
    mixture: Mixture,
    h_mixture: float,
    guesstemperature: float = 0.0,
) -> int:
    """Compute the mixture temperature from the given mixture enthalpy."""
    """The solved mixture temperature is stored as the temperature attribute of
    the given gas mixture (i.e., as mixture.temperature)

    Parameters
    ----------
        mixture: Mixture object
            gas mixture of interest
        h_mixture: double
            mixture enthalpy of the given gas mixture [erg/mol]
        guesstemperature: double, optional
            a guessed value for the mixture temperature at the start of
            the iteration process

    Returns
    -------
        error code: integer

    """
    # check argument
    if not isinstance(mixture, Mixture):
        msg = [Color.PURPLE, "the first argument must be a Mixture object.", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # make a copy of the mixture object
    localmixture = copy.deepcopy(mixture)
    # set converge tolerance
    tolerance = 0.1  # accurate to 0.1 K
    # iteration count limit
    maxcount = 200
    count = 0
    ierr = 0
    dt = 1.0e3
    # set guessed temperature value if given
    if guesstemperature > 0.0e0:
        localmixture.temperature = guesstemperature
    # solve for the temperature by using the Newton's method
    while True:
        # function: H(T) = h_mixture
        # compute value at T = localmixture.temperature
        f = localmixture.hml() - h_mixture
        # compute slope at T = localmixture.temperature
        df = localmixture.cpbl()
        try:
            # compute correction
            dt = f / df
        except ZeroDivisionError:
            # diverge
            msg = [Color.PURPLE, "search diverged.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierr = 1
            break

        if abs(dt) <= tolerance:
            # search converges
            break
        # update temperature T
        localmixture.temperature -= dt
        count += 1

    if count >= maxcount:
        # not converging within count limit
        msg = [
            Color.PURPLE,
            "failed to reach the desired tolerance within",
            str(maxcount),
            "iterations\n",
            Color.SPACEx6,
            "the final temperature tolerance =",
            str(abs(dt)),
            "[K].",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        ierr = 2
    # update temperature
    if ierr != 1:
        mixture.temperature = localmixture.temperature
    # print(f'** iteration count = {count:d}')
    del localmixture
    return ierr


def interpolate_mixtures(
    mixtureleft: Mixture, mixtureright: Mixture, ratio: float
) -> Mixture:
    """Get Mixture by interpolation."""
    """Create a new mixture object by interpolating the two mixture objects
    with a specific weight ratio.

    ::
        mixture_new = (1 - ratio) * mixtureleft + ratio * mixtureright

    Parameters
    ----------
        mixtureleft: Mixture object
            mixture A to be mixed
        mixtureright: Mixture object
            mixture B to be mixed
        ratio: double
            the weight parameters for interpolation, 0 <= ratio <= 1

    Returns
    -------
        mixturenew: Mixture object
            the resulting gas mixture

    """
    # check mixtures
    if not isinstance(mixtureleft, Mixture):
        msg = [
            Color.PURPLE,
            "the first argument must be a Mixture object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    ierr = mixtureleft.validate()
    if ierr != 0:
        msg = [
            Color.PURPLE,
            "the first mixture is not fully defined.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    #
    if not isinstance(mixtureright, Mixture):
        msg = [
            Color.PURPLE,
            "the second argument must be a Mixture object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    ierr = mixtureright.validate()
    if ierr != 0:
        msg = [
            Color.PURPLE,
            "the second mixture is not fully defined.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check ratio
    if ratio < 0.0e0 or ratio > 1.0e0:
        msg = [
            Color.PURPLE,
            "the weight ratio must be 0 <= and <= 1.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check chemistry sets
    if mixtureright.chemid != mixtureleft.chemid:
        msg = [
            Color.PURPLE,
            "the Mixtures belong to different Chemistry Sets.\n",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check if the chemstry set is active
    if not check_active_chemistryset(mixtureright.chemid):
        msg = [
            Color.PURPLE,
            "the Chemistry Set associated with",
            "the Mixture is not currently active.\n",
            Color.SPACEx6,
            "activate Chemistry Set using the 'active()' method.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
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
    fracleft = mixtureleft.y
    fracright = mixtureright.y
    frac = np.zeros(len(fracleft), dtype=np.double)
    frac = ratiom * fracleft + ratio * fracright
    mixturenew.y = frac
    # clean up
    del frac
    #
    return mixturenew


def compare_mixtures(
    mixture_a: Mixture,
    mixture_b: Mixture,
    atol: float = 1.0e-10,
    rtol: float = 1.0e-3,
    mode: str = "mass",
) -> tuple[bool, float, float]:
    """Compare properties of mixture B against those of mixture A."""
    """The mixture properties include pressure [atm], temperature [K],
    and species mass/mole fractions. When the differences in the property values
    satisfy both the absolute and the relative tolerances, this method will
    return "True", that is, mixture B is essentially identical to mixture A;
    otherwise, "False" will be returned.

    Parameters
    ----------
        mixture_a: Mixture object
            mixture A, the target mixture
        mixture_b: Mixture object
            mixture B, the sample mixture
        atol: double, default = 1.0e-10
            the absolute tolerance for the max property differences
        rtol: double, default = 1.0e-3
            the relative tolerance for the max property differences
        mode: string {"mass", "mole"}, default = "mass"
            compare species "mass" or "mole" fractions

    Returns
    -------
        issame: boolean
            the equivalence of the two mixtures
        atol_max: double
            the max absolute difference value
        rtol_max: double
            the max relative difference value

    """
    # check mixtures
    if not isinstance(mixture_a, Mixture):
        msg = [
            Color.PURPLE,
            "the first argument must be a Mixture object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    ierr = mixture_a.validate()
    if ierr != 0:
        msg = [
            Color.PURPLE,
            "the first mixture is not fully defined.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    #
    if not isinstance(mixture_b, Mixture):
        msg = [
            Color.PURPLE,
            "the second argument must be a Mixture object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    ierr = mixture_b.validate()
    if ierr != 0:
        msg = [
            Color.PURPLE,
            "the second mixture is not fully defined.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check chemistry sets
    if mixture_a.chemid != mixture_b.chemid:
        msg = [
            Color.PURPLE,
            "the Mixtures belong to different Chemistry Sets.\n",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check if the chemstry set is active
    if not check_active_chemistryset(mixture_a.chemid):
        msg = [
            Color.PURPLE,
            "the Chemistry Set associated with the Mixture is not currently active.\n",
            Color.SPACEx6,
            "activate Chemistry Set using the 'active()' method.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # compare mixture pressure
    pres_diff = abs(mixture_a.pressure - mixture_b.pressure)
    # find relative difference
    pres_var = pres_diff / mixture_a.pressure
    # convert the difference to [atm]
    pres_diff /= P_ATM
    # check tolerances
    issame = pres_diff <= atol
    issame = issame or pres_var <= rtol
    diff_max = pres_diff
    var_max = pres_var
    if not issame:
        print(f"pressure difference: {pres_diff}   {pres_var}")
    # compare mixture temperature
    temp_diff = abs(mixture_a.temperature - mixture_b.temperature)
    # find relative difference
    temp_var = temp_diff / mixture_a.temperature
    # check tolerances
    issame = issame or temp_diff <= atol
    issame = issame or temp_var <= rtol
    diff_max = max(diff_max, temp_diff)
    var_max = max(var_max, temp_var)
    if not issame:
        print(f"temperature difference: {temp_diff}   {temp_var}")
    # compare composition
    spec_index_count = 0
    spec_index_max = []
    spec_diff_max: list[float] = []
    spec_var_max: list[float] = []
    if mode == "mole":
        # compare mole fractions
        k = 0
        while k < mixture_a._kk:
            frac = mixture_a.x[k]
            spec_diff = abs(frac - mixture_b.x[k])
            if np.isclose(frac, 0.0, atol=atol):
                found = issame or spec_diff <= atol
                spec_var = 0.0
            else:
                spec_var = spec_diff / frac
                found = spec_diff <= atol
                found = found or spec_var <= rtol
            #
            if found:
                spec_index_max.append(k)
                spec_diff_max.append(spec_diff)
                spec_var_max.append(spec_var)
                spec_index_count += 1
            k += 1
    else:
        # compare mass fractions
        k = 0
        while k < mixture_a._kk:
            frac = mixture_a.y[k]
            spec_diff = abs(frac - mixture_b.y[k])
            if np.isclose(frac, 0.0, atol=atol):
                found = issame or spec_diff <= atol
                spec_var = 0.0
            else:
                spec_var = spec_diff / frac
                found = spec_diff <= atol
                found = found or spec_var <= rtol
            #
            if not found:
                spec_index_max.append(k)
                spec_diff_max.append(spec_diff)
                spec_var_max.append(spec_var)
                spec_index_count += 1
            k += 1
    # check tolerances
    if spec_index_count > 0:
        issame = False
        print("composition differences:")
        count = 0
        for k in spec_index_max:
            print(f"species {mixture_a._specieslist[k]}")
            print(f"   difference: {spec_diff_max[count]}   {spec_var_max[count]}")
            count += 1
        diff_spec = max(spec_diff_max)
        diff_max = max(diff_max, diff_spec)
        var_spec = max(spec_var_max)
        var_max = max(var_max, var_spec)
        print(f"spec value {diff_spec}   {var_spec}")
        print(f"max value {diff_max}   {var_max}")
    #
    return issame, diff_max, var_max


# equilibrium
#
def calculate_equilibrium(
    chemid: int,
    p: float,
    t: float,
    frac: npt.NDArray[np.double],
    wt: npt.NDArray[np.double],
    mode_in: str,
    mode_out: str,
    eq_option: int = 1,
    use_realgas: int = 0,
) -> tuple[list[float], npt.NDArray[np.double]]:
    """Get the equilibrium mixture composition."""
    """Get the equilibrium mixture composition corresponding to
    the given initial mixture composition at the given pressure
    and temperature.

    Parameters
    ----------
        chemid: integer
            chemistry set index associated with the mixture
        p: double
            initial mixture pressure in [dynes/cm2]
        t: double
            initial mixture temperature in [K]
        frac: 1-D double array
            initial mixture composition given by either mass or mole fractions
            as specified by mode_in
        wt: 1-D double arrays
            molar masses of the species in the mixture in [gm/mol]
        mode_in: string, {'mass', 'mole'}, default = 'mole'
            flag indicates the frac array is 'mass' or 'mole' fractions
        mode_out: string, {'mass', 'mole'}, default = 'mole'
            flag to indicate the returning composition is in 'mole' or 'mass' fraction
        eq_option: integer, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
            equilibrium type (see below)

            ::
                1.  SPECIFIED T AND P
                2.  SPECIFIED T AND V
                3.  SPECIFIED T AND S
                4.  SPECIFIED P AND V
                5.  SPECIFIED P AND H
                6.  SPECIFIED P AND S
                7.  SPECIFIED V AND U
                8.  SPECIFIED V AND H
                9.  SPECIFIED V AND S
                10. CHAPMAN-JOUGUET DETONATION
        use_realgas: integer, {0, 1}
            option to turned ON/OFF (1/0) the real-gas cubic EOS if available

    Returns
    -------
        state_variables_equilibrium: list of doubles,
            equilibrium pressure [dynes/cm2],
            equilibrium temperature [K],
            speed of sound at equilibrium [cm/sec],
            detonation wave speed [cm/sec].
            Note: if Chapmen-Jouguet option is not used,
            both speed of sound and detonation wave speed are set to 0.0
        equilibrium composition: 1-D double array
            given in fractions indicated by the parameter mode_out

    """
    # find the equilibrium composition at the mixture pressure and temperature
    # check inputs
    if chemid < 0:
        msg = [Color.PURPLE, "invalid chemistry.", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check if the chemstry set is active
    if not check_active_chemistryset(chemid):
        msg = [
            Color.PURPLE,
            "the Chemistry Set associated with the Mixture is not currently active.\n",
            Color.SPACEx6,
            "activate Chemistry Set using the 'active()' method.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()

    if p <= 0.0 or (p * t) <= 0.0:
        msg = [
            Color.PURPLE,
            "invalid pressure and/or temperature value(s).",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()

    # number species
    kgas = len(frac)
    if kgas != len(wt):
        msg = [
            Color.PURPLE,
            mode_in,
            "fraction and molar mass arrays",
            "must have the same size =",
            str(kgas),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()

    # initialization
    xx_eq = np.zeros(kgas, dtype=np.double)
    #
    if mode_in.lower() == "mole":
        # normalize mass fractions
        ierr, x = Mixture.normalize(frac=frac)
    elif mode_in.lower() == "mass":
        # convert mass fraction to mole fraction and normalize
        x = Mixture.mass_fraction_to_mole_fraction(massfrac=frac, wt=wt)
    else:
        # fraction type not given or incorrect
        msg = [
            Color.PURPLE,
            'must specify "mole" or "mass" fractions given.',
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check equilibrium calculation option
    if eq_option in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10):
        eq_option_c = c_int(eq_option)
    else:
        # set to constant T-P option by default
        eq_option_c = c_int(1)

    # check real gas option
    if use_realgas == 1:
        i_realgas = c_int(1)
        set_current_pressure(chemid, pressure=p)
    else:
        # use ideal gas law by default
        i_realgas = c_int(0)

    # convert parameters to c pointers
    _chemset_index = ctypes.c_int(chemid)
    pp = c_double(p)  # pressure scalar
    tt = c_double(t)  # temperature scalar
    xx = np.ctypeslib.as_array(x)  # mole fraction array
    #
    pp_eq = c_double(p)
    tt_eq = c_double(t)
    detonationwavespeed = c_double(0.0e0)
    soundspeed_eq = c_double(0.0e0)
    # perform gas-phase equilibrium calculationk
    if not check_chemistryset(_chemset_index.value):
        # need to initialize Chemkin-CFD-API
        msg = [Color.YELLOW, "initializing Chemkin", "...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)

        ierr = ck_wrapper.chemkin.KINInitialize(_chemset_index, c_int(0))
        if ierr != 0:
            msg = [
                Color.RED,
                "Chemkin-CFD-API initialization failed;",
                "code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        else:
            chemistryset_initialized(_chemset_index.value)
    else:
        ierr = 0

    ierr = ck_wrapper.chemkin.KINCalculateEqGasWithOption(
        _chemset_index,
        eq_option_c,
        i_realgas,
        pp,
        tt,
        xx,
        pp_eq,
        tt_eq,
        soundspeed_eq,
        detonationwavespeed,
        xx_eq,
    )

    if ierr == 0:
        # process solution
        if eq_option_c.value == 10 and verbose():
            # CHAPMAN-JOUGUET DETONATION
            print(
                f"** detonation wave speed = {detonationwavespeed.value / 1.0e2} "
                "[m/sec]"
            )
            print(
                f"** speed of sound at final state = {soundspeed_eq.value / 1.0e2} "
                "[m/sec]"
            )

        if mode_out.lower() == "mass":
            # convert mass fraction to mole fraction and normalize
            y_eq = Mixture.mole_fraction_to_mass_fraction(molefrac=xx_eq, wt=wt)
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
            ierr, x_eq = Mixture.normalize(frac=xx_eq)
            statevars = [
                pp_eq.value,
                tt_eq.value,
                soundspeed_eq.value,
                detonationwavespeed.value,
            ]
            return statevars, x_eq

    else:
        msg = [Color.PURPLE, "failed to find the equilibrium state.", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()


def equilibrium(mixture: Mixture, opt: int = 1) -> Mixture:
    """Find the equilibrium state mixture corresponding to the given mixture."""
    """
    Find the equilibrium state mixture corresponding to the given mixture
    with given constraints.

    Parameters
    ----------
        mixture: Mixture object
            initial gas mixture
        opt: integer, {1, 2, 4, 5, 7, 8}
            equilibrium type

            ::
                1.  SPECIFIED T AND P
                2.  SPECIFIED T AND V
                3.  SPECIFIED T AND S (*)
                4.  SPECIFIED P AND V
                5.  SPECIFIED P AND H
                6.  SPECIFIED P AND S (*)
                7.  SPECIFIED V AND U
                8.  SPECIFIED V AND H
                9.  SPECIFIED V AND S (*)
                10. CHAPMAN-JOUGUET DETONATION (*)

                (*) indicates the options are not available

    Returns
    -------
        finalmixture: Mixture object
            gas mixture at the equilibrium state

    """
    # check argument
    if not isinstance(mixture, Mixture):
        msg = [
            Color.PURPLE,
            "the first argument must be a Mixture object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check if the chemstry set is active
    if not check_active_chemistryset(mixture.chemid):
        msg = [
            Color.PURPLE,
            "the Chemistry Set associated with the Mixture is not currently active.\n",
            Color.SPACEx6,
            "activate Chemistry Set using the 'active()' method.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # initialization a Mixture object by duplication
    eq_state = copy.deepcopy(mixture)
    # reset mass/mole fractions
    eq_state._x_set = 0
    eq_state._molefrac[:] = 0.0e0
    eq_state._y_set = 0
    eq_state._massfrac[:] = 0.0e0
    # check option
    if opt in [3, 6, 9, 10]:
        msg = [
            Color.PURPLE,
            "equilibrium option",
            str(opt),
            "is not available.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    else:
        option = opt

    if mixture.eos == 0:
        userealgas = 0
    else:
        userealgas = mixture.userealgas
    # compute the equilibrium state (mass fraction for now)
    eqvars, eq_state._massfrac = calculate_equilibrium(
        mixture.chemid,
        p=eq_state.pressure,
        t=eq_state.temperature,
        frac=mixture.y,
        wt=mixture._wt,
        mode_in="mass",
        mode_out="mass",
        eq_option=option,
        use_realgas=userealgas,
    )
    if np.sum(eq_state._massfrac, dtype=np.double) > 0.0e0:
        eq_state._y_set = 1
    eq_state.pressure = eqvars[0]
    eq_state.temperature = eqvars[1]
    return eq_state


def detonation(mixture: Mixture) -> tuple[list[float], Mixture]:
    """Find the Chapman-Jouguet state mixture and detonation wave speed."""
    """Find the Chapman-Jouguet state mixture and detonation wave speed
    corresponding to the given mixture.

    Parameters
    ----------
        mixture: Mixture object
            initial gas mixture

    Returns
    -------
        speed_values: list of doubles
            speed of sound [cm/sec],
            detonation wave speed [cm/sec].
        finalmixture: Mixture object
            gas mixture at the equilibrium state

    """
    # check argument
    if not isinstance(mixture, Mixture):
        msg = [Color.PURPLE, "the argument must be a Mixture object.", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check if the chemstry set is active
    if not check_active_chemistryset(mixture.chemid):
        msg = [
            Color.PURPLE,
            "the Chemistry Set associated with the Mixture is not currently active.\n",
            Color.SPACEx6,
            "activate Chemistry Set using the 'active()' method.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # initialization a Mixture object by duplication
    eq_state = copy.deepcopy(mixture)
    # reset mass/mole fractions
    eq_state._x_set = 0
    eq_state._molefrac[:] = 0.0e0
    eq_state._y_set = 0
    eq_state._massfrac[:] = 0.0e0
    # use the C-J option
    option = 10
    if mixture.eos == 0:
        userealgas = 0
    else:
        userealgas = mixture.userealgas
    # compute the equilibrium state (mass fraction for now)
    eqvars, eq_state._massfrac = calculate_equilibrium(
        mixture.chemid,
        p=eq_state.pressure,
        t=eq_state.temperature,
        frac=mixture.y,
        wt=mixture._wt,
        mode_in="mass",
        mode_out="mass",
        eq_option=option,
        use_realgas=userealgas,
    )
    if np.sum(eq_state._massfrac, dtype=np.double) > 0.0e0:
        eq_state._y_set = 1
    eq_state.pressure = eqvars[0]
    eq_state.temperature = eqvars[1]
    return [eqvars[2], eqvars[3]], eq_state

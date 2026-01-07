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

"""Constants used by Chemkin utilities and models."""

# == Chemkin module global parameters
# -- DO NOT MODIFY without asking Chemkin development team members
BOLTZMANN = 1.3806504e-16  # BOLTZMANN constant [ergs/K] (double scalar)
AVOGADRO = 6.02214179e23  # AVOGADRO number [1/mole] (double scalar)
P_ATM = 1.01325e06  # atmospheric pressure [dynes/cm2] (double scalar)
P_TORRS = P_ATM / 760.0  # 1 torr in [dynes/cm2] (double scalar)
ERGS_PER_JOULE = 1.0e7  # ergs per joule [ergs/J] (double scalar)
JOULES_PER_CALORIE = 4.184e0  # joules per calorie [J/cal] (double scalar)
ERGS_PER_CALORIE = (
    JOULES_PER_CALORIE * ERGS_PER_JOULE
)  # ergs per calorie [erg/cal] (double scalar)
ERGS_PER_EV = 1.602176487e-12  # ergs per eV [erg/volt] (double scalar)
EV_PER_K = ERGS_PER_EV / BOLTZMANN  # eV per K [volt/K] (double scalar)
R_GAS = BOLTZMANN * AVOGADRO  # universal gas constant R [ergs/mol-K] (double scalar)
R_GAS_CAL = (
    R_GAS * 1.0e-7 / JOULES_PER_CALORIE
)  # universal gas constant R [cal/mol-K] (double scalar)
# == end of global constants


class Air:
    """define the "air" composition in PyChemkin with a fixed mixture "recipe"."""

    """
    A "recipe" is a list of tuples of ("species symbol", fraction) to define a
    gas mixture in PyChemkin.
    This class uses the upper case symbols for oxygen and nitrogen.
    """

    @staticmethod
    def x(cap: str = "U") -> list[tuple[str, float]]:
        """Return the 'air' composition in mole fractions."""
        """
        Return the 'air' composition in mole fractions.

        Parameters
        ----------
            cap: string, {"U", "L"}
                indicating the upper or the lower case species symbols to be used
        """
        if cap.upper() == "L":
            # Return the 'air' composition using
            # the lower case symbols for oxygen and nitrogen
            return [("o2", 0.21), ("n2", 0.79)]
        else:
            # Return the 'air' composition using
            # the upper case symbols for oxygen and nitrogen
            return [("O2", 0.21), ("N2", 0.79)]

    @staticmethod
    def y(cap: str = "U") -> list[tuple[str, float]]:
        """Return the 'air' composition in mass fractions."""
        """
        Return the 'air' composition in mass fractions.

        Parameters
        ----------
            cap: string, {"U", "L"}
                indicating the upper or the lower case species symbols to be used
        """
        if cap.upper() == "L":
            # Return the 'air' composition using
            # the lower case symbols for oxygen and nitrogen
            return [("o2", 0.23), ("n2", 0.77)]
        else:
            # Return the 'air' composition using
            # the upper case symbols for oxygen and nitrogen
            return [("O2", 0.23), ("N2", 0.77)]


def water_heat_vaporization(temperature: float) -> float:
    """Get the heat if vporization of water at the given temperature [K]."""
    """
    Parameters
    ----------
        temperature: double
            water temperature [K]

    Returns
    -------
        enthalpy: double
            enthalpy of vaporization of water at the given temperature [erg/g-water]

    """
    # critical temperature of water [K]
    tc = 647.096
    # reduce temperature
    tr = temperature / tc
    # check value
    if temperature < 1.5e2:
        print("temperature value is too low.")
        exit()
    elif tr >= 1.0:
        print(
            f"temperature value is above the critical temperature of wwater {tc} [K]."
        )
        exit()
    #
    a = (5.1546e7, 0.28402, -0.15843, 0.2375)
    tr_m1 = 1.0 - tr
    #
    index = 1
    fac = 1.0
    exponent = 0.0
    while index < 4:
        exponent += a[index] * fac
        fac *= tr
        index += 1
    h = tr_m1**exponent
    # heat of vaporization [J/kmol]
    h *= a[0]
    # convert to [erg/g-water]
    h_erg = h * 1.0e-3 / 18.0 * 1.0e7
    return h_erg

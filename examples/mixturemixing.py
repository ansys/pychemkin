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

"""
.. _ref_mixing_mixtures:

===================================
Combining gas mixtures in PyChemkin
===================================
PyChemkin provides a set of basic mixture utilities enabling the creation of new mixtures from the existing ones.
The mixing method, as its name suggested, let you combine two mixtures at constant pressure with specific constraints.
The ``adiabatic_mixing`` method combines the two mixtures while keeping the overall enthalpy constant. The temperature of
the combined mixture, in this case, is determined by the conservation of the overall enthalpy of the two *"parent"* mixtures.
The ``isothermal_mixing``, on the other hand, simply combines the two mixtures and determines the composition of the final mixture
according to the mole or mass ratios specified. You must assign any temperature value of the combined mixture since the
enthalpy conservation is not considered by this method.

This tutorial will show you the usage of these two ``Mixture`` *mixing* methods, ``isothermal_mixing`` and ``adiabatic_mixing``,
and the difference in the temperature of the combined/final mixture between the two mixing methods. You will first create a
fuel (CH\ :sub:`4`\ ) mixture and an air (O\ :sub:`2`\ +N\ :sub:`2`\ ) mixture and then make a fuel-air mixture by mixing
them *isothermally* (that is, without considering energy conservation) with a given air-to-fuel mass ratio. Afterwards, you will
dilute the fuel-air mixture with argon (AR) *adiabatically* with the molar/volumetric ratio specified.
"""

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin.logger import logger

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)


#####################################
# Create a ``Chemistry Set`` instance
# ===================================
# Load the GRI 3.0 mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"/reaction/data"*.

# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
# transport data not needed


###################################
# Pre-process the ``Chemistry Set``
# =================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()


####################################################################
# Set up gas mixtures based on the species in this ``Chemistry Set``
# ==================================================================
# Use the *equivalence ratio method* to set up the combustible mixture
# so that you can easily change the mixture composition by assigning a
# different *equivalence ratio* value.


#########################
# Create the fuel mixture
# =======================
# The fuel mixture consists of 100% methane.
#
# .. note::
#   Mixture pressures are not specified here because it is not required by the calculations.
#   The mixing process takes place at fixed pressure by assumption, that is, the mixtures are at the same pressure.

fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.X = [("CH4", 1.0)]
fuel.temperature = 300.0


#########################
# Create the air mixture
# =======================
# "Air" is a mixture of oxygen and nitrogen.

air = ck.Mixture(MyGasMech)
# set mole fraction
air.X = [("O2", 0.21), ("N2", 0.79)]
air.temperature = 300.0


#####################################
# Create a fuel-air mixture by mixing
# ===================================
# Mix mixtures ``fuel`` and ``air`` created in the previous steps by the ``isothermal_mixing`` method.
# The *mixing formula* of the two mixtures are defined by the ``mixture_recipe`` with mass ratio:
# ``fuel:air=1.00:17.19``. You must assign the *temperature* of the new mixture ``premixed``
# through the use of parameter ``finaltemperature=300``. Set ``mode="mass"`` because
# the ratios given in the ``mixture_recipe`` are mass ratios.

# mix the fuel and the air with an air-fuel ratio of 17.19 (almost stoichiometric)
mixture_recipe = [(fuel, 1.0), (air, 17.19)]
# create the new mixture (the air-fuel ratio is by mass)
premixed = ck.isothermal_mixing(
    recipe=mixture_recipe, mode="mass", finaltemperature=300.0
)


###############################################
# Display the molar composition of ``premixed``
# =============================================
# Use the ``list_composition`` with ``mode="mole"`` to list the mole fractions.
# The molar composition should resemble the stoichiometric methane-air mixture.

# list the molar composition
premixed.list_composition(mode="mole")
print()


#######################################
# Create a diluent mixture with pure AR
# =====================================
# Now create an argon mixture and use it to dilute ``premixed`` later. Set the mixture
# to a higher temperature of 600[K].

ar = ck.Mixture(MyGasMech)
# species composition
ar.X = [("AR", 1.0)]
# mixture temperature
ar.temperature = 600.0


##############################################################
# Create a new AR-diluted fuel-air mixture by adiabatic mixing
# ============================================================
# Dilute ``premixed`` by mixing 30% ``ar`` by volume adiabatically. You do not
# need to provide the final mixture temperature to ``adiabatic_mixing`` because it is
# determined by the enthalpy conservation. Parameter ``mode="mole"`` is set to indicate
# the ratios in ``dilute_recipe`` are molar ratios.

# create the mixing recipe
dilute_recipe = [(premixed, 0.7), (ar, 0.3)]
# create the diluted mixture
diluted = ck.adiabatic_mixing(recipe=dilute_recipe, mode="mole")


####################################
# Display information of ``diluted``
# ==================================
# Use the ``list_composition`` method with ``mode="mole"`` to display the molar composition.
# also display the temperatures of the three mixtures involved in the *adiabatic* mixing process
# for verification. The ``diluted`` temperature should sit in between those of ``premixed`` and ``ar``.

# list molar composition
diluted.list_composition(mode="mole")
# show the mixture temperatures
print(f"the diluted mixture temperature is  {diluted.temperature:f} [K]")
print(f"the ar mixture temperature is       {ar.temperature:f} [K]")
print(f"the premixed mixture temperature is {premixed.temperature:f} [K]")

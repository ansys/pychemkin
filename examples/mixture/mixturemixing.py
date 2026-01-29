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

r""".. _ref_mixing_mixtures:

====================
Combine gas mixtures
====================

PyChemkin provides a set of basic mixture utilities enabling the creation of
new mixtures from existing ones. The mixing methods let you combine two mixtures
at constant pressure with specific constraints.

- The ``adiabatic_mixing()`` method combines two mixtures while keeping
  the overall enthalpy constant. The temperature of the combined mixture is
  determined by the conservation of the overall enthalpy of the two parent mixtures.

- The ``isothermal_mixing()`` method simply combines two mixtures and determines
  the composition of the final mixture according to the mole or mass ratios specified.
  Because this method does not consider the enthalpy conservation, you must assign
  the temperature value of the combined mixture.

This example shows how to use these two mixing methods and understand the differences
in the temperatures of their combined mixtures. It first creates a fuel
(CH\ :sub:`4`\ ) mixture and an air (O\ :sub:`2`\ +N\ :sub:`2`\ ) mixture and then
makes a fuel-air mixture by mixing them *isothermally* (that is, without considering
the enthalpy conservation) with a given air-to-fuel mass ratio. The example
then dilutes the fuel-air mixture with argon (AR) *adiabatically* with
the molar/volumetric ratio specified.
"""

# sphinx_gallery_thumbnail_path = '_static/merge_operation.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================
# Import the PyChemkin packages, check the working directory, and start
# the logger in verbose mode.

from pathlib import Path

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)


########################
# Create a chemistry set
# ======================
# Load the GRI 3.0 mechanism and its associated data files,
# which come with the standard Ansys Chemkin installation
# in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
# transport data is not needed


##############################
# Preprocess the chemistry set
# ============================
# Preprocess the mechanism files to prepare the chemistry set.

# preprocess the mechanism files
ierror = MyGasMech.preprocess()


#####################
# Set up gas mixtures
# ===================
# Set up gas mixtures based on the species in this chemistry set.
# Use the equivalence ratio method to set up the combustible mixture
# so that you can easily change the mixture composition by assigning a
# different equivalence ratio value.


#######################
# Create a fuel mixture
# =====================
# Create a fuel mixture of 100% methane.
#
# .. note::
#   Mixture pressures are not specified here because they are not required
#   by the calculations. The mixing process assumes a fixed pressure,
#   meaning that the mixtures are at the same pressure.

fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.x = [("CH4", 1.0)]
fuel.temperature = 300.0


#######################
# Create an air mixture
# =====================
# Create an air mixture of oxygen and nitrogen.

air = ck.Mixture(MyGasMech)
# set mole fraction
air.x = [("O2", 0.21), ("N2", 0.79)]
air.temperature = 300.0


#####################################
# Create a fuel-air mixture by mixing
# ===================================
# Use the ``isothermal_mixing()`` method to mix the ``fuel`` and ``air`` mixtures
# created earlier. Define the mixing formula using ``mixture_recipe`` with
# this mass ratio: ``fuel:air=1.00:17.19``. Use the ``finaltemperature`` parameter
# to set the temperature of the new ``premixed`` mixture to 300 [K].
# Set ``mode="mass"`` because the ratios given in ``mixture_recipe`` are
# mass ratios.

# mix the fuel and the air with an air-fuel ratio of 17.19 (almost stoichiometric)
mixture_recipe = [(fuel, 1.0), (air, 17.19)]
# create the new mixture (the air-fuel ratio is by mass)
premixed = ck.isothermal_mixing(
    recipe=mixture_recipe, mode="mass", finaltemperature=300.0
)


#######################################################
# Display the molar composition of the premixed mixture
# =====================================================
# Use the ``list_composition()`` method with ``mode="mole"`` to
# list the mole fractions. The molar composition should resemble
# the stoichiometric methane-air mixture.

# list the molar composition
premixed.list_composition(mode="mole")
print()


##########################################
# Create a diluent mixture with pure argon
# ========================================
# Create an argon mixture to dilute the premixed mixture later. Set the mixture
# temperature to a temperature of 600 [K].

ar = ck.Mixture(MyGasMech)
# species composition
ar.x = [("AR", 1.0)]
# mixture temperature
ar.temperature = 600.0


############################################################
# Dilute the fuel-air mixture with argon by adiabatic mixing
# ==========================================================
# Dilute the premixed mixture by mixing 30% argon by volume adiabatically.
# The ``adiabatic_mixing()`` method determines the final mixture temperature
# based on the enthalpy conservation. Set ``mode="mole"`` to indicate that the ratios
# in ``dilute_recipe`` are molar ratios.

# create the mixing recipe
dilute_recipe = [(premixed, 0.7), (ar, 0.3)]
# create the diluted mixture
diluted = ck.adiabatic_mixing(recipe=dilute_recipe, mode="mole")


#############################################
# Display information for the diluted mixture
# ===========================================
# Use the ``list_composition()`` method with ``mode="mole"`` to display
# the molar composition. Also display the temperatures of the three mixtures
# involved in the adiabatic mixing process for verification. The temperature of
# the diluted mixture should sit in between those of the ``premixed`` and ``ar``
# mixtures.

# list molar composition
diluted.list_composition(mode="mole")
# show the mixture temperatures
print(f"The diluted mixture temperature is  {diluted.temperature:f} [K].")
print(f"The ar mixture temperature is       {ar.temperature:f} [K].")
print(f"The premixed mixture temperature is {premixed.temperature:f} [K].")

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
.. _ref_equilibrium_composition:

=================================================================================
Estimate the steady-state NO emission level of a complete burned fuel-air mixture
=================================================================================

This example shows how to quickly estimate the steady-state NO level
that is formed by the combustion of a fuel-air mixture at a given temperature.
Without needing any reaction, the NO concentration gets to its steady state
(or the maximum level) when the product mixture from the fuel-air combustion
reaches the equilibrium state at the given temperature. To find the equilibrium state of
the fresh fuel-air mixture, the ``equilibrium()`` method is used with the ``fixed pressure``
and ``fixed temperature`` options.

This example explores the influence of temperature on the predicted adiabatic flame temperature.
Knowing that nitrogen oxides (NOx) are stable at high temperatures (> 2000 [K]),
you can expect that the steady-state NO level should increase sharply when the temperature
gets high enough.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_equilibrium_composition.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True

########################
# Create a chemistry set
# ======================
# The first mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
# transport data is not needed

##############################
# Preprocess the chemistry set
# ============================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

#####################
# Set up gas mixtures
# ===================
# Set up gas mixtures based on the species in this chemistry set.
# PyChemkin provides a few methods for creating a gas mixture. Here,
# the ``isothermal_mixing()`` method is used to create a ``fuel-air`` mixture
# by mixing the ``fuel`` and ``air`` mixtures with a predetermined
# air/fuel rate.

#######################
# Create a fuel mixture
# =====================
# Create a fuel mixture of 80% methane and 20% hydrogen.

fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.X = [("CH4", 0.8), ("H2", 0.2)]
fuel.temperature = 300.0
fuel.pressure = ck.Patm  # 1 atm

#######################
# Create an air mixture
# =====================
# Create an air mixture of oxygen and nitrogen.

air = ck.Mixture(MyGasMech)
# set mass fraction
air.Y = [("O2", 0.23), ("N2", 0.77)]
air.temperature = 300.0
air.pressure = ck.Patm  # 1 atm


#####################################
# Create a fuel-air mixture by mixing
# ===================================
# Use the ``isothermal_mixing()`` method to mix the ``fuel`` and ``air`` mixtures created earlier.
# Define the mixing formula using ``mixture_recipe`` with this mass ratio: ``fuel:air=1.00:17.19``.
# Use the ``finaltemperature`` parameter to set the temperature of the
# new ``premixed`` mixture to 300 [K]. Set ``mode="mass"`` because
# the ratios given in ``mixture_recipe`` are mass ratios.

mixture_recipe = [(fuel, 1.0), (air, 17.19)]
# create the new mixture (the air-fuel ratio is by mass)
premixed = ck.isothermal_mixing(
    recipe=mixture_recipe, mode="mass", finaltemperature=300.0
)

############################################################
# Find the equilibrium composition at different temperatures
# ==========================================================
# Perform the equilibrium calculation with fixed pressure and
# fixed temperature. The NO mole fraction of the equilibrium state
# at each temperature is stored in an array. The gas temperature
# is increased from 500 to 2480 [K].

# NO species index
NO_index = MyGasMech.get_specindex("NO")

# set up plotting temperatures
Temp = 500.0
dTemp = 20.0
points = 100
# curve
T = np.zeros(points, dtype=np.double)
NO = np.zeros_like(T, dtype=np.double)
# start the temperature loop
for k in range(points):
    # reset mixture temperature
    premixed.temperature = Temp
    # find the equilibrium state mixture at the given mixture temperature and pressure
    eqstate = ck.equilibrium(premixed, opt=1)
    #
    NO[k] = eqstate.X[NO_index] * 1.0e6  # convert to ppm
    T[k] = Temp
    Temp += dTemp


##################
# Plot the results
# ================
# Plot the equilibrium NO mole fractions versus the temperatures
# of the fuel-air mixtures.

plt.plot(T, NO, "bs--", markersize=3, markevery=4)
plt.xlabel("Temperature [K]")
plt.ylabel("NO [ppm]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_equilibrium_composition.png", bbox_inches="tight")

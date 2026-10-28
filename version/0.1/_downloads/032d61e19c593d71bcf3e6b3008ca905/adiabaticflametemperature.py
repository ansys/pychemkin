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

r""".. _ref_adiabatic_flame_temperature:

=========================================================
Estimate the adiabatic flame temperature of a gas mixture
=========================================================

This example shows how to find the equilibrium state of a mixture.
It uses the ``equilibrium()`` method with the ``constant pressure and enthalpy`` option
to estimate the adiabatic flame temperature of a methane-oxygen mixture. This example
also explores the influence of the equivalence ratio on the predicted
adiabatic flame temperature.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_adiabatic_flame_temperature.png'

###############################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
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
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir

# create a chemistry set based on the GRI 3.0 methane combustion mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
# skip the transport data file where is not needed by this example

##############################
# Preprocess the chemistry set
# ============================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()

#####################
# Set up gas mixtures
# ===================
# Set up gas mixtures based on the species in this chemistry set.
# PyChemkin has a few methods for creating a gas mixture.
# Here, the equivalence ratio method is used to set up the combustible mixture
# so that you can easily change the mixture composition by assigning a
# different equivalence ratio value.

# create an "oxid" mixture associated with the 'MyGasMech' chemistry set
oxid = ck.Mixture(MyGasMech)
# use a "recipe" to set the mole fractions of the mixture
# the "oxid" mixture consists of 100% O2
oxid.x = [("O2", 1.0)]
oxid.temperature = 295.15  # [K]
oxid.pressure = ck.P_ATM  # 1 atm

# create the "fuel" mixture
fuel = ck.Mixture(MyGasMech)
# set the "fuel" molar composition to 100% CH4
fuel.x = [("CH4", 1.0)]
fuel.temperature = oxid.temperature
fuel.pressure = oxid.pressure

# create the final fuel-oxidizer mixture
mixture = ck.Mixture(MyGasMech)
mixture.pressure = oxid.pressure
mixture.temperature = oxid.temperature

# the use of equivalence ratio requires the definition of the complete combustion
# products of the fuel-oxidizer pair:
# CH4 + 2O2 => CO2 + 2H2O
products = ["CO2", "H2O"]

# create an array to specify the composition of the additives to
# the fuel-oxidizer mixture For example, a diluent such as argon or
# helium might be added to the fuel-oxidizer mixture Use an all-zero array
# if there is no additive
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)

############################
# Set up the parameter study
# ==========================
# Set up a parameter study to find out the impact of the equivalence ratio
# on the adiabatic flame temperature of the fuel-oxidizer mixture.
# The equivalence ratio varies from 0.5 to 1.6 with an increment of 0.1.

points = 12
deq = 0.1
equiv_ini = 0.5

# create the solution arrays as double arrays
t = np.zeros(points, dtype=np.double)
equiv = np.zeros_like(t, dtype=np.double)

#########################
# Run the parameter study
# =======================
# Use the ``equilibrium()`` method to estimate the adiabatic flame temperature of
# the fuel-oxidizer mixture. Choose the option corresponding to constant pressure
# and constant enthalpy for the equilibrium calculation because you are finding
# the adiabatic temperature.
#
# The ``equilibrium()`` method returns a mixture object representing the mixture
# at the equilibrium state. Use the ``temperature()`` method to get the
# equilibrium temperature. To see all available options for this method, use
# the ``ck.help(topic="equilibrium")`` method.
#
# This example uses the ``x_by_equivalence_ratio()`` method to set
# the fuel-oxidizer composition with the given equivalence ratios
# because the composition of both the fuel and oxidizer mixtures is specified
# in mole fractions.

for i in range(points):
    # set the current mixture equivalence ratio
    equiv_current = equiv_ini

    # create the fuel-oxidizer mixture with the given equivalence ratio
    ierror = mixture.x_by_equivalence_ratio(
        MyGasMech, fuel.x, oxid.x, add_frac, products, equivalenceratio=equiv_current
    )
    # check fuel-oxidizer mixture creation status
    if ierror != 0:
        print("Error: Failed to create the fuel-oxidizer mixture.")
        print(f"       Equivalence ratio = {equiv_current}.")
        exit()

    # use "equilibrium()" method to calculate the gas mixture at the equilibrium state
    # Option #5 ("opt=5") corresponds to the constant pressure and constant-enthalpy
    # constraints of the equilibrium state
    # "EQ_mixture" is the gas mixture at the equilibrium state
    EQ_mixture = ck.equilibrium(mixture, opt=5)

    # save the results to the solution arrays
    # use "temperature" to obtain the temperature of the equilibrium state
    t[i] = EQ_mixture.temperature
    equiv[i] = equiv_current
    equiv_ini = equiv_ini + deq

##########################################
# Plot the result from the parameter study
# ========================================
# When you plot the result from the parameter study, the adiabatic flame temperature
# should exhibit a peak at around the stoichiometric, that is, equivalence ratio = 1.

# plot equilibrium/adiabatic temperatures against mixture equivalence ratios
plt.plot(equiv, t, "bs--")
# set up axis labels
plt.xlabel("Equivalence ratio")
plt.ylabel("Temperature [K]")
# display or save the plot
if interactive:
    plt.show()
else:
    plt.savefig("plot_adiabatic_flame_temperature.png", bbox_inches="tight")

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

r""".. _ref_detonation_wave:

=========================================================
Calculate the detonation wave speed of a real-gas mixture
=========================================================
Use the ``detonation()`` method on a combustible mixture to compute the
Chapman-Jouguet (C-J) state and detonation wave speed. This example shows
how to predict the detonation wave speeds of a natural gas-air mixture
at various initial pressures and compare the *ideal-gas* and the *real-gas* results
against the experimental data.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_detonation.png'

################################################
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
# The ``C2 NOx`` mechanism comes with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.
# This mechanism includes information about the *Soave* cubic
# Equation of State (EOS) for the real-gas applications. The PyChemkin preprocessor
# indicates the availability of the real-gas model in the chemistry set processed.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on C2 NOx using an alternative method
MyMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files individually
# Because this mechanism file contains all the necessary thermodynamic and
# transport data, you do not need to specify thermodynamic and
# transport data files.
MyMech.chemfile = str(mechanism_dir / "C2_NOx_SRK.inp")

#####################################
# Preprocess the C2 NOx chemistry set
# ===================================
# During preprocessing, you should see this printed:
# ``real-gas cubic EOS 'Soave' is available``.
# Because no transport data file is provided and the ``preprocess_transportdata()``
# method is not used, transport property methods are not available in this project.

# preprocess the mechanism files
ierror = MyMech.preprocess()

######################################################################
# Set up gas mixtures based on the species in the C2 NOx chemistry set
# ====================================================================
# Create gas mixtures named ``fuel`` (natural gas) and ``air`` based on
# ``MyMech``. Then, use these two mixtures to form the combustible ``premixed``
# mixture for the detonation calculations. Use the ``x_by_equivalence_ratio()`` method
# to set the equivalence ratio of the fuel-air mixture to 1.

# create the fuel mixture
fuel = ck.Mixture(MyMech)
# set mole fraction
fuel.x = [("CH4", 0.8), ("C2H6", 0.2)]
fuel.temperature = 290.0
fuel.pressure = 40.0 * ck.P_ATM
# create the air mixture
air = ck.Mixture(MyMech)
# set mass fraction
air.x = [("O2", 0.21), ("N2", 0.79)]
air.temperature = fuel.temperature
air.pressure = fuel.pressure
# create the initial mixture
# create the premixed mixture to be defined by equivalence ratio
premixed = ck.Mixture(MyMech)
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture.
# Can also create an additives mixture here.
add_frac = np.zeros(MyMech.kk, dtype=np.double)  # no additives: all zeros

ierror = premixed.x_by_equivalence_ratio(
    MyMech, fuel.x, air.x, add_frac, products, equivalenceratio=1.0
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
    print("Error: Failed to create the premixed mixture.")
    exit()

#######################################################
# Display the molar composition of the premixed mixture
# =====================================================
# List the composition of the premixed mixture for verification.
premixed.list_composition(mode="mole")

################################
# Run the detonation calculation
# ==============================
# Use the ``detonation()`` method to find the C-J state and detonation
# wave speed of the fuel-air mixture with the initial mixture pressure
# increasing from 40 to 80 [atm]. The initial mixture temperature is kept at 290 [K].
#
# The ``detonation()`` method returns two objects: a ``speed`` tuple containing the
# speed of sound and the detonation wave speed at the C-J state. The ``CJState``
# mixture contains the mixture properties at the C-J state. For instance,
# you can use ``CJState.pressure`` to get the mixture pressure at the C-J state.
#
# .. note::
#
#    - By default, Chemkin variables are in cgs units.
#    - You can enter the ``ansys.chemkin.core.help("equilibrium")`` command
#      at the Python prompt to see the input and output parameters of
#      the ``detonation()`` method.
#

#########################
# Run the parameter study
# =======================
# Set up the parameter study of the detonation wave speed with respect
# to the initial pressure. The predicted detonation wave speed values are saved
# in the ``Det`` array. The experimental data is stored in the ``Det_data`` array.
# By default, the *ideal-gas law* is assumed. You can use the
# ``use_realgas_cubicEOS()`` method to turn on the *real-gas model* if the mechanism
# contains the real-gas parameters in the EOS block. Use the ``use_idealgas_law()``
# method to reactivate the ideal-gas law assumption.
points = 5
dpres = 10.0 * ck.P_ATM
pres = fuel.pressure
p = np.zeros(points, dtype=np.double)
det = np.zeros_like(p, dtype=np.double)
premixed.pressure = pres
premixed.temperature = fuel.temperature

# start of pressure loop
for i in range(points):
    # compute the C-J state corresponding to the initial mixture
    speed, cj_state = ck.detonation(premixed)
    # update plot data
    # convert pressure to atm
    p[i] = pres / ck.P_ATM
    # convert speed to m/sec
    det[i] = speed[1] / 1.0e2
    # update pressure value
    pres += dpres
    premixed.pressure = pres

# create plot for ideal-gas results
plt.plot(p, det, "bo--", label="ideal gas", markersize=5, fillstyle="none")

##################################
# Switch to the real-gas EOS model
# ================================
# Use the ``use_realgas_cubicEOS()`` method to turn on the real-gas EOS model.
# You can enter ``ansys.chemkin.core.help("real gas")`` to see usage information
# on real-gas models or ``ansys.chemkin.core.help("manuals")`` to access
# the online *Chemkin Theory* manual for descriptions of the real-gas EOS models.
#
# .. note::
#   By default, the *Van der Waals* mixing rule is applied to evaluate
#   thermodynamic properties of a real-gas mixture. You can use the
#   ``set_realgas_mixing_rule()`` method to switch to a different mixing rule.

# turn on real-gas cubic equation of state
premixed.use_realgas_cubic_eos()
# set mixture mixing rule to Van der Waals (default)
# premixed.set_realgas_mixing_rule(rule=0)
# restart the calculation with real-gas EOS
premixed.pressure = fuel.pressure
pres = fuel.pressure
p[:] = 0.0e0
det[:] = 0.0e0
# set verbose mode to false to turn off extra printouts
ck.set_verbose(False)
# start of pressure loop
for i in range(points):
    # compute the C-J state corresponding to the initial mixture
    speed, cj_state = ck.detonation(premixed)
    # update plot data
    p[i] = pres / ck.P_ATM
    det[i] = speed[1] / 1.0e2
    # update pressure value
    pres += dpres
    premixed.pressure = pres

# stop Chemkin
ck.done()
# create plot for real-gas results
plt.plot(p, det, "r^-", label="real gas", markersize=5, fillstyle="none")
# plot data
p_data = [44.1, 50.6, 67.2, 80.8]
det_data = [1950.0, 1970.0, 2000.0, 2020.0]
plt.plot(p_data, det_data, "gD:", label="data", markersize=4)

##########################################
# Plot the result from the parameter study
# ========================================
# You should see that the ideal-gas assumption fails to show any noticeable
# pressure influence on the detonation wave speeds. Because of the relatively high
# pressures in this study, you can observe significant differences in the predicted
# detonation wave speeds between the ideal-gas and real-gas models.
plt.legend(loc="upper left")
plt.xlabel("Pressure [atm]")
plt.ylabel("Detonation wave speed [m/sec]")
plt.suptitle("Natural Gas/Air Detonation", fontsize=16)
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_detonation.png", bbox_inches="tight")

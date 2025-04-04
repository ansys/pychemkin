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
.. _ref_detonation_wave:

===========================================================
Calculating the detonation wave speed of a real gas mixture
===========================================================
PyChemkin ``Mixture`` object offers the ``detonation`` method that serves as a convenient tool
to compute the *Chapman-Jouguet state* and the *detonation wave speed* of a combustible mixture.

This tutorial utilizes the ``detonation`` method to predict the detonation wave speeds of a
natural gas-air mixture at various initial pressures and compare the *ideal-gas* and the
*real-gas* results against the experimental data.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_detonation.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

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
# interactive = False: save plot as a png file
global interactive
interactive = True

#####################################
# Create a ``Chemistry Set`` instance
# ===================================
# The 'C2 NOx' mechanism is from the default *"/reaction/data"* directory.
# This mechanism also includes information about the *Soave* cubic
# Equation of State (EOS) for the real-gas applications. PyChemkin preprocessor
# will indicate the availability of the real-gas model in the ``Chemistry Set`` processed.

# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on C2_NOx using an alternative method
MyMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files individually
# this mechanism file contains all the necessary thermodynamic and transport data
# therefore no need to specify the therm and the tran data files
MyMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")

##########################################
# Pre-process the C2 NOx ``Chemistry Set``
# ========================================
# you should see the print-out *"real-gas cubic EOS 'Soave' is available"* during
# the preprocess. Since no transport data file is provided nor the ``preprocess_transportdata``
# method is used, the transport property methods will *not* be available in this project.

# preprocess the mechanism files
iError = MyMech.preprocess()

##########################################################################
# Set up gas mixtures based on the species in the C2 NOx ``Chemistry Set``
# ========================================================================
# Create a gas mixture instances ``fuel`` (natural gas) and ``air`` based on
# ``MyMech``. Then use these two mixtures to form the combustible ``premixed``
# mixture for the detonation calculations. the ``X_by_Equivalence_Ratio`` method
# is used to set the *equivalence ratio* of the fuel-air mixture to *1*.

# create the fuel mixture
fuel = ck.Mixture(MyMech)
# set mole fraction
fuel.X = [("CH4", 0.8), ("C2H6", 0.2)]
fuel.temperature = 290.0
fuel.pressure = 40.0 * ck.Patm
# create the air mixture
air = ck.Mixture(MyMech)
# set mass fraction
air.X = [("O2", 0.21), ("N2", 0.79)]
air.temperature = fuel.temperature
air.pressure = fuel.pressure
# create the initial mixture
# create the premixed mixture to be defined by equivalence ratio
premixed = ck.Mixture(MyMech)
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
add_frac = np.zeros(MyMech.KK, dtype=np.double)  # no additives: all zeros

iError = premixed.X_by_Equivalence_Ratio(
    MyMech, fuel.X, air.X, add_frac, products, equivalenceratio=1.0
)
# check fuel-oxidizer mixture creation status
if iError != 0:
    print("Error: failed to create the premixed mixture!")
    exit()

###############################################
# Display the molar composition of ``premixed``
# =============================================
# list the composition of the premixed mixture for verification.
premixed.list_composition(mode="mole")

########################################
# Perform the ``detonation`` calculation
# ======================================
# Find the *Chapman-Jouguet state* (C-J state) and the *detonation wave speed* of the
# fuel-air mixture by utilizing the ``detonation`` method to find the C-J state and
# the detonation wave speed of the fuel-air mixture with the initial mixture pressure
# increasing from 40 to 80 [atm]. The initial mixture temperature is kept at 290 [K].
#
# The ``detonation`` method will return two objects: a ``speed`` *tuple* containing the
# *speed of sound* and the *detonation wave speed* at the C-J state; the ``CJState``
# ``Mixture`` object containing the mixture properties at the C-J state. For instance,
# you can get the mixture pressure at the C-J state using the method ``CJState.pressure``.
#
# .. note::
#   By default *Chemkin* variables are in the **cgs units**.
#
#
# .. note::
#   You can check out the input/output parameters of the ``detonation`` method by issuing
#   command ``annsys.chemkin.help("equilibrium")`` at the python prompt.
#

#########################
# Run the parameter study
# =======================
# set up the parameter study of detonation wave speed with respect to the initial pressure.
# The predicted *detonation wave speed* values are saved in the ``Det`` array, and the
# experimental data are stored in the ``Det_data`` array. By default, the *ideal gas law*
# is assumed. You may use the ``use_realgas_cubicEOS`` method to turn *ON* the
# *real gas model* if the mechanism contains the real-gas parameters in the "EOS" block.
# Use ``use_idealgas_law`` method to reactivate the ideal gas law assumption.
points = 5
dpres = 10.0 * ck.Patm
pres = fuel.pressure
P = np.zeros(points, dtype=np.double)
Det = np.zeros_like(P, dtype=np.double)
premixed.pressure = pres
premixed.temperature = fuel.temperature

# start of pressure loop
for i in range(points):
    # compute the C-J state corresponding to the initial mixture
    speed, CJstate = ck.detonation(premixed)
    # update plot data
    # convert pressure to atm
    P[i] = pres / ck.Patm
    # convert speed to m/sec
    Det[i] = speed[1] / 1.0e2
    # update pressure value
    pres += dpres
    premixed.pressure = pres

# create plot for ideal gas results
plt.plot(P, Det, "bo--", label="ideal gas", markersize=5, fillstyle="none")

##################################
# Switch to the real gas EOS model
# ================================
# Use the ``use_realgas_cubicEOS`` method to turn ON the real-gas EOS model. For more
# information either type ``ansys.chemkin.help("real gas")`` for the real-gas model
# usage or type ``ansys.chemkin.help("manuals")`` to access the on-line **Chemkin Theory**
# manual for descriptions of the real-gas EOS models.
#
# .. note::
#   By default the *Van der Waals* mixing rule is applied to evaluate thermodynamic properties
#   of a real gas mixture. You can use ``set_realgas_mixing_rule`` to switch to a different
#   mixing rule.

# turn on real-gas cubic equation of state
premixed.use_realgas_cubicEOS()
# set mixture mixing rule to Van der Waals (default)
# premixed.set_realgas_mixing_rule(rule=0)
# restart the calculation with real-gas EOS
premixed.pressure = fuel.pressure
pres = fuel.pressure
P[:] = 0.0e0
Det[:] = 0.0e0
# set verbose mode to false to turn OFF extra printouts
ck.set_verbose(False)
# start of pressure loop
for i in range(points):
    # compute the C-J state corresponding to the initial mixture
    speed, CJstate = ck.detonation(premixed)
    # update plot data
    P[i] = pres / ck.Patm
    Det[i] = speed[1] / 1.0e2
    # update pressure value
    pres += dpres
    premixed.pressure = pres

# stop Chemkin
ck.done()
# create plot for real gas results
plt.plot(P, Det, "r^-", label="real gas", markersize=5, fillstyle="none")
# plot data
P_data = [44.1, 50.6, 67.2, 80.8]
Det_data = [1950.0, 1970.0, 2000.0, 2020.0]
plt.plot(P_data, Det_data, "gD:", label="data", markersize=4)

###########################################
# Plot the result from this parameter study
# =========================================
# You should see that the ideal gas assumption fails to show any noticeable
# pressure influence on the detonation wave speeds. Because of the relatively high
# pressures in this study, significant differences in the predicted detonation wave speeds
# between the ideal gas and the real-gas models are observed.
plt.legend(loc="upper left")
plt.xlabel("Pressure [atm]")
plt.ylabel("Detonation wave speed [m/sec]")
plt.suptitle("Natural Gas/Air Detonation", fontsize=16)
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_detonation.png", bbox_inches="tight")

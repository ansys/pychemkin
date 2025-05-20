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
.. _ref_cooling_vapor:

===========================
Explore cooling water vapor
===========================

How does the volume of water vapor evolves when it is cooled from a temperature above the boiling point
to a temperature that is just above the freezing point at constant pressure? How fast does the vapor
volume drop? This example uses the ``GivenPressureBatchReactor_FixedTemperature`` model to explore
the different behaviors between an *ideal gas* water vapor and its *real gas* counterpart.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_vapor_condensation.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color

# chemkin batch reactor model (transient)
from ansys.chemkin.batchreactors.batchreactor import (
    GivenPressureBatchReactor_FixedTemperature,
)
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True

########################
# Create a chemistry set
# ======================
# To compare the different behaviors of the water vapor under the ideal gas and real gas
# assumptions, you must use a *real gas EOS-enabled gas mechanism. The 'C2 NOx' mechanism that
# includes information about the *Soave* cubic Equation of State (EOS) is suitable for this endeavor.
# Therefore, the ``MyMech`` chemistry set is created from this gas phase mechanism.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on C2_NOx using an alternative method
MyMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files individually
# this mechanism file contains all the necessary thermodynamic and transport data
# thus, there is no need to specify thermodynamic and transport data files
MyMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")

#####################################
# Preprocess the C2 NOx chemistry set
# ===================================

# preprocess the mechanism file
iError = MyMech.preprocess()
if iError == 0:
    print(Color.GREEN + ">>> Preprocessing succeeded.", end=Color.END)
else:
    print(Color.RED + ">>> Preprocessing failed.", end=Color.END)
    exit()

####################################
# Set up the air/water vapor mixture
# ==================================
# Create the mixture of air and water vapor inside the imaginary strong but elastic container.
# The initial gas temperature is set to 500 [K] and at a constant pressure of 100 [atm]. At this
# temperature and pressure, the mixture should be entirely in gas form.
mist = ck.Mixture(MyMech)
# set mole fraction
mist.X = [("H2O", 2.0), ("O2", 1.0), ("N2", 3.76)]
mist.temperature = 500.0  # [K]
mist.pressure = 100.0 * ck.Patm

######################################################################
# Create the reactor ``tank`` to perform the vapor cooling simulation
# ====================================================================
# Use the ``GivenPressureBatchReactor_FixedTemperature()`` method to create
# a constant-pressure batch reactor (with a given temperature). Use the ``mist`` mixture
# that you just created to set the initial gas condition inside the ``tank``
# reactor.
#
tank = GivenPressureBatchReactor_FixedTemperature(mist, label="tank")

############################################
# Set up additional reactor model parameters
# ==========================================
# You must provide reactor parameters, solver controls, and output instructions
# before running the simulations. For a batch reactor, the initial volume and
# simulation end time are required inputs.

# verify initial gas composition inside the reactor
tank.list_composition(mode="mole")

# set other reactor properties
tank.volume = 10.0  # cm3
tank.time = 0.5  # sec

#########################################################
# Set the gas temperature profile of the ``tank`` reactor
# =======================================================
# Create a time-temperature profile by using two arrays. Use the ``set_temperature_profile()``
# method to add the profile to the reactor model.

# number of profile data points
npoints = 3
# position array of the profile data
x = np.zeros(npoints, dtype=np.double)
# value array of the profile data
TPROprofile = np.zeros_like(x, dtype=np.double)
# set tank temperature data points
x = [0.0, 0.2, 2.0]  # [sec]
TPROprofile = [500.0, 275.0, 275.0]  # [K]
# set the temperature profile
tank.set_temperature_profile(x, TPROprofile)

####################################
# Switch on the real-gas EOS model
# ==================================
# Use the ``use_realgas_cubicEOS()`` method to turn on the real-gas EOS model. For more
# information, type either ``ansys.chemkin.help("real gas")`` for information on real-gas model
# usage or ``ansys.chemkin.help("manuals")`` to access the online **Chemkin Theory**
# manual for descriptions of the real-gas EOS models.
#
# .. note::
#   By default the *Van der Waals* mixing rule is applied to evaluate thermodynamic properties
#   of a real-gas mixture. You can use ``set_realgas_mixing_rule`` to switch to a different
#   mixing rule.
#
tank.userealgasEOS(mode=True)

####################
# Set output options
# ==================

# set timestep between saving solution
tank.timestep_for_saving_solution = 0.01

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver related methods, for example,
# ``tolerances``.

# set tolerances in tuple: (absolute tolerance, relative tolerance)
tank.tolerances = (1.0e-10, 1.0e-8)
# get solver parameters
ATOL, RTOL = tank.tolerances
print(f"default absolute tolerance = {ATOL}")
print(f"default relative tolerance = {RTOL}")
# turn on the force non-negative solutions option in the solver
tank.force_nonnegative = True

##########################################################
# Run the vapor cooling simulation with the *real gas EOS*
# ========================================================
# Run the CONP reactor model with given temperature profile with the *real gas EOS*
# model switched *ON*.
runstatus = tank.run()

# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)

##########################
# Postprocess the solution
# ========================
# The postprocessing step parses the solution and packages the solution values at each
# time point into a mixture. There are two ways to access the solution profiles:
#
# - The raw solution profiles (value as a function of distance) are available for distance,
#   temperature, pressure, volume, and species mass fractions.
#
#  -The mixtures permit the use of all property and rate utilities to extract
#   information such as viscosity, density, and mole fractions.
#
# You can use the ``get_solution_variable_profile()`` method to get the raw solution profiles. You
# can get solution mixtures using either the ``get_solution_mixture_at_index()`` method for the
# solution mixture at the given saved location or the ``get_solution_mixture()`` method for the
# solution mixture at the given distance. (In this case, the mixture is constructed by interpolation.)
#
# .. note::
#   Use the ``getnumbersolutionpoints()`` method to get the size of the solution profiles before
#   creating the arrays.
#

tank.process_solution()
# get the number of solution time points
solutionpoints = tank.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile = tank.get_solution_variable_profile("time")
# get the temperature profile
tempprofile = tank.get_solution_variable_profile("temperature")
# get the volume profile
volprofile = tank.get_solution_variable_profile("volume")
# create array for mixture density
denprofile = np.zeros_like(timeprofile, dtype=np.double)
# create array for mixture enthalpy
Hprofile = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tank.get_solution_mixture_at_index(solution_index=i)
    # get mixture density profile
    denprofile[i] = solutionmixture.RHO
    # get mixture enthalpy profile
    Hprofile[i] = solutionmixture.HML() / ck.ergs_per_joule * 1.0e-3

#################################
# Turn off the real gas EOS model
# ===============================
# Use the ``use_realgas_cubicEOS()`` method to turn off the real gas EOS model.
# Alternatively, use the ``use_idealgas_law()`` method to turn on the ideal
# gas law.
tank.userealgasEOS(mode=False)

#########################################################
# Run the vapor cooling simulation with the ideal gas law
# =======================================================
# Run the CONP reactor model with given temperature profile with the
# *ideal gas law* turned back on.
runstatus = tank.run()

# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)

####################################
# Postprocess the ideal gas solution
# ==================================
# postprocess the solutions
tank.process_solution()
# get the number of solution time points
solutionpoints = tank.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile_IG = tank.get_solution_variable_profile("time")
# get the volume profile
volprofile_IG = tank.get_solution_variable_profile("volume")
# create array for mixture density
denprofile_IG = np.zeros_like(timeprofile, dtype=np.double)
# create array for mixture enthalpy
Hprofile_IG = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tank.get_solution_mixture_at_index(solution_index=i)
    # get mixture density profile
    denprofile_IG[i] = solutionmixture.RHO
    # get mixture enthalpy profile
    Hprofile_IG[i] = solutionmixture.HML() / ck.ergs_per_joule * 1.0e-3

ck.done()

################################
# Plot the RCM solution profiles
# ==============================
# You should observe that the mixture volume obtained by the real gas model is
# noticeably lower than the ideal gas volume. When water vapor is cooled below
# the boiling point, formation of liquid water is expected due to condensation.
# The ideal gas law assumes that the mixture is always in gas phase and is unable to
# address the phase change phenomenon. The real gas EOS, on the other hand, can
# capture the formation of the liquid water as indicated by the sharper rise of the
# mixture density during the cooling process. The real gas mixture also has a lower
# enthalpy level than that of the ideal gas mixture. The enthalpy differences
# should largely represent the heat of vaporization of water at the temperature.

# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
thispres = str(mist.pressure / ck.Patm)
thistitle = "Cooling Vapor + Air at " + thispres + " atm"
plt.suptitle(thistitle, fontsize=16)
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, volprofile, "b-", label="real gas")
plt.plot(timeprofile_IG, volprofile_IG, "b--", label="ideal gas")
plt.legend(loc="upper right")
plt.ylabel("Volume [cm3]")
plt.subplot(223)
plt.plot(timeprofile, Hprofile, "g-", label="real gas")
plt.plot(timeprofile_IG, Hprofile_IG, "g--", label="ideal gas")
plt.legend(loc="upper right")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Enthalpy [kJ/mole]")
plt.subplot(224)
plt.plot(timeprofile, denprofile, "m-", label="real gas")
plt.plot(timeprofile_IG, denprofile_IG, "m--", label="ideal gas")
plt.legend(loc="upper left")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Density [g/cm3]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_vapor_condensation.png", bbox_inches="tight")

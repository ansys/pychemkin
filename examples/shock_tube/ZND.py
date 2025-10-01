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
.. _ref_shock_ZND_reactor:

==========================================================
Perform ZND analysis of a combustible hydrogen mixture
==========================================================

The ZND (Zeldovich-von Neumann-Deoring) model is commonly applied to study the
evolution of the gas mixture behind a detonation wave.

The ZND model is a transient reactor, and the solutions describe the state of the gas mixture
(particle) as it moving away from the detonation front in the shock tube reactor.

Use the ``ZNDCalculator()`` method to create a ZND shock tube reactor.
The shock tube reactor models, such as the IncidentShcok, the ReflectedShock and the ZNDCalculator models,
are initiated by a stream, which is simply a mixture with the addition of the shock wave velocity.
You already know how to create a stream if you know how to create
a mixture. You can specify the shock wave velocity using the combination of the ``velocity`` method of the
initial gas stream and the ``location`` parameter when you instantiate the ``incidentShock`` or
the ``ReflectedShock`` object.

This example shows how to use the Chemkin ZND shock tube model to study the
gas mixture evolution behind an incident detonation wave front. The speed of the incident shock wave
is estimated by the ZND model automatically by finding the stable detonation wave speed corresponding to the given
gas mixture before the incident shock front. Once the solution is obtained, you can extract the induction
zone length behind the wave front and utilize this information to gain insights about the characteristic size of
the cellular structure in multi-dimensional flow as well as the stability of the cellular structure.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_shock_ZND_reactor.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger

# chemkin plug flow reactor model
from ansys.chemkin.shock.shocktubereactors import ZNDCalculator as ZND
from ansys.chemkin.utilities import find_file
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
# The mechanism used here is the MFL 2021 hydrogen mechanism.
# The mechanism comes with the standard Ansys Chemkin
# installation in the ``/reaction/data/ModelFuelLibrary/Skeletal`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(
    ck.ansys_dir, "reaction", "data", "ModelFuelLibrary", "Skeletal"
)
mechanism_dir = data_dir
# create a chemistry set based on the MFL 2021 hydrogen mechanism
MyGasMech = ck.Chemistry(label="hydrogen")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = find_file(mechanism_dir, "Hydrogen_chem_MFL", "inp")

#######################################
# Preprocess the hydrogen chemistry set
# =====================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

###################################
# Instantiate and set up the stream
# =================================
# Create a combustible gas stream by assigning the mole fractions of the
# species. Here, several gas streams are created with different
# recipes consisting of H\ :sub:`2`\ , O\ :sub:`2` and N\ :sub:`2`\ .
# You can pick any one of the streams to instantiate the ``ZNDCalculator``.
# You can compare the "induction lengths" by different streams to see
# how the combustion intensity (reversely proportion to the amount of
# N\ :sub:`2` dilution) impacts the size of the induction zone behind
# the detonation wave front.
#
streamA = Stream(MyGasMech)
# initial gas state before the incident shock front
# 1 [atm]
streamA.pressure = ck.Patm
# 300 [K]
streamA.temperature = 300.0
# 50% hydrogen + 50% oxygen by volume
streamA.X = [("h2", 1.0), ("o2", 1.0)]
#
streamB = Stream(MyGasMech)
# a diluted initial gas state before the incident shock front
# 1 [atm]
streamB.pressure = ck.Patm
# 300 [K]
streamB.temperature = 300.0
# the 50% hydrogen + 50% oxygen mixture diluted by 40% nitrogen by volume
streamB.X = [("h2", 0.3), ("o2", 0.3), ("n2", 0.4)]
#
streamC = Stream(MyGasMech)
# a diluted initial gas state before the incident shock front
# 1 [atm]
streamC.pressure = ck.Patm
# 300 [K]
streamC.temperature = 300.0
# the 50% hydrogen + 50% oxygen mixture diluted by 40% nitrogen by volume
streamC.X = [("h2", 0.2), ("o2", 0.2), ("n2", 0.6)]

######################################
# Create the shock tube reactor object
# ====================================
# Use the ``ZNDCalculator()`` method to create an incident shock reactor.
# The required input parameter is the stream representing the state of the initial
# gas mixture before the incident shock front. In this case, ``streamA``
# is used to initialize the ZND model ``ZNDIncident``.

ZNDIncident = ZND(streamC, label="ZND")

############################################
# Set up additional reactor model parameters
# ==========================================
# For the ZND calculator model, the required reactor parameters is the total simulation time [sec].
# The initial gas mixture conditions are defined by the stream when the ``ZNDIncident`` is instantiated.

# set total simulation time (particle time) [sec]
ZNDIncident.time = 3.0e-6

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver related methods, for example,
# ``tolerances``.

# tolerances are given in tuple: (absolute tolerance, relative tolerance)
ZNDIncident.tolerances = (1.0e-10, 1.0e-6)

######################################################
# Run the ZND analysis
# ====================================================
# Use the ``run()`` method to start the ZND analysis.
#
# .. note ::
#   You can use two ``time`` calls (one before the run and one after the run) to
#   get the simulation run time (wall time).
#

# set the start wall time
start_time = time.time()
# run the ZND model
runstatus = ZNDIncident.run()
# compute the total runtime
runtime = time.time() - start_time
# check run status
if runstatus != 0:
    # Run failed.
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run succeeded.
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
print(f"Total simulation duration: {runtime * 1.0e3} [msec]")

##########################
# Postprocess the solution
# ========================
# The postprocessing step parses the solution and packages the solution values at each
# time point into a mixture. There are two ways to access the solution profiles:
#
# - The raw solution profiles (value as a function of time) are available for distance,
#   temperature, pressure, velocity, density, and species mass fractions.
#
# - The mixtures permit the use of all property and rate utilities to extract
#   information such as viscosity, density, total thermicity, speed of sound, Mach number,
#   and mole fractions.
#
# You can use the ``get_solution_variable_profile()`` method to get the raw solution profiles. You
# can get solution mixtures using either the ``get_solution_stream_at_index()`` method for the
# solution mixture at the given saved location or the ``get_solution_stream()`` method for the
# solution mixture at the given distance. (In this case, the mixture is constructed by interpolation.)
#
# You can get ZND analysis information about the induction zone length by using the
# ``get_inductionlength_size()`` and the ``get_induction_lengths()`` methods. The induction zone length
# can be used to derive the characteristic cell size and the stability of the multi-dimensional cellular
# wave structure.

# postprocess the solution profiles
ZNDIncident.process_solution()

# get the number of solution time points
solutionpoints = ZNDIncident.get_solution_size()
print(f"number of solution points = {solutionpoints}")
# get information about the induction length
numb_induct_length = ZNDIncident.get_inductionlength_size()
print(f"number of induction length = {numb_induct_length}")
if numb_induct_length > 0:
    induct_legth, sigmamax = ZNDIncident.get_induction_lengths(numb_induct_length)
    for i in range(numb_induct_length):
        print(f"Induction length # {i+1}:")
        print(f"        Induction length = {induct_legth[i] * 1.0e1} [mm].")
        print(f"    Local max thermicity = {sigmamax[i]} [1/sec].\n")

# estimated detonation wave cell structure parameters
cell_width, chi = ZNDIncident.calculate_cell_width_ng()
print(f"Instability parameter Chi = {chi}")  # stable = Chi > 1
print(f"Estimated cell size = {cell_width} [cm].")

# get the grid profile [cm]
xprofile = ZNDIncident.get_solution_variable_profile("distance")
# get the temperature profile [K]
tempprofile = ZNDIncident.get_solution_variable_profile("temperature")
# get the pressure profile [dynes/cm2]
pressprofile = ZNDIncident.get_solution_variable_profile("pressure")
# get the velocity profile [cm/sec]
velprofile = ZNDIncident.get_solution_variable_profile("velocity")
# get the total thermiocity profile [1/sec]
sigmaprofile = ZNDIncident.get_solution_variable_profile("thermicity")

# create arrays for the gas Mach number profile
Machprofile = np.zeros_like(xprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = ZNDIncident.get_solution_stream_at_index(solution_index=i)
    # gas speed of sound [cm/sec]
    soundspeed = solutionmixture.sound_speed()
    # gas Mach number
    Machprofile[i] = velprofile[i] / soundspeed
    # convert pressure from [dynes/cm2] to [bar]
    pressprofile[i] /= 1.0e6

############################
# Plot the solution profiles
# ==========================
# Plot the temperature, pressure, total thermicity, and the gas Mach number
# profiles as a function of distance behind the wave front.
#

plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("ZND Analysis", fontsize=16)
plt.subplot(221)
plt.plot(xprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(xprofile, pressprofile, "b-")
plt.ylabel("Pressure [bar]")
plt.subplot(223)
plt.plot(xprofile, Machprofile, "g-")
plt.xlabel("distance behind shock [cm]")
plt.ylabel("Mach Number [-]")
plt.subplot(224)
plt.plot(xprofile, sigmaprofile, "m-")
plt.xlabel("distance behind shock [cm]")
plt.ylabel("Total Thermicity [1/sec]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_shock_ZND_reactor.png", bbox_inches="tight")

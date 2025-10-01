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
.. _ref_reflected_shock_reactor:

==============================================================
Dissociation of nitrous oxide behind a reflected shock wave
==============================================================

Shock tube experiments are commonly used to study reaction paths and to measure reaction rates
at elevated temperatures. You can apply the incident shock reactor model ``IncidentShock()``
to validate the reaction mechanism or kinetic parameters derived from such experiments.

The shock tube reactor models, such as the IncidentShcok, the ReflectedShock and the ZNDCalculator models,
are initiated by a stream, which is simply a mixture with the addition of the shock wave velocity.
You already know how to create a stream if you know how to create
a mixture. You can specify the shock wave velocity using the combination of the ``velocity`` method of the
initial gas stream and the ``location`` parameter when you instantiate the ``incidentShock`` or
the ``ReflectedShock`` object.

The gas condition behind a reflected shock is hot and uniform, and is favorable for studying chemical
kinetics at high temperatures. This example utilizes the N\ :sub:`2`\ O dissociation mechanism proposed
by Baber and Dean to simulate the evolution of N\ :sub:`2`\ O behind a reflected shock.
The temperature profile and the mass fraction profiles of major species N\ :sub:`2`\ O, NO,
and O\ :sub:`2` in the hot zone behind the reflected shock will be plotted as a function of time.

Reference:
S.C. Baber and A.M. Dean, International Journal of Chemical Kinetics, vol. 7, pp. 381-398 (1975)
"""

# sphinx_gallery_thumbnail_path = '_static/plot_reflected_shock.png'

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
from ansys.chemkin.shock.shocktubereactors import ReflectedShock
import matplotlib.pyplot as plt  # plotting

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True

######################################################
# Create the nitrous oxide dissociation mechanism file
# ====================================================
# Create a new mechanism input file 'n2o_dissociation_chem.inp' that
# contains the reactions to describe N2O dissociation.
# This file is saved to the working directory ``current_dir``.
#
mymechfile = os.path.join(current_dir, "n2o_dissociation_chem.inp")
m = open(mymechfile, "w")
# the mechanism contains only the necessary species (oxygen, nitrogen, nitrous oxide, nitric oxide, and major byproducts)
# decalre elements
m.write("ELEMENT O N AR END\n")
# declare species
m.write("SPECIES\n")
m.write("O N O2 N2 NO N2O NO2 AR\n")
m.write("END\n")
# write reactions for N2O dissociation
# Reference:
# S.C. Baber and A.M. Dean, International Journal of Chemical Kinetics, vol. 7, pp. 381-398 (1975)
# S.C. Baber and A.M. Dean, Journal of Chemical Physics, vol. 60, Mo. 1, pp. 307-313 (1974)
m.write("REACTIONS\n")
m.write("N2O+M=N2+O+M            7.83E+14  0.0    56845.32\n")
m.write("N2O+O=NO+NO             1.15E+13  0.0    25078.82\n")
m.write("N2O+O=N2+O2             1.15E+13  0.0    25078.82\n")
m.write("O+NO=N+O2               1.55E+09  1.0    38640.\n")
m.write("N+NO=N2+O               3.10E+13  0.0      334.\n")
m.write("NO+O2=O+NO2             1.99E+12  0.0    47740.\n")
m.write("END\n")
# close the mechnaism file
m.close()

########################
# Create a chemistry set
# ======================
# The mechanism used here is the air dissociation mechanism.
# The mechanism will be created in situ in the working directory.
# The thermodynamic data file is the standard one that comes with
# the Ansys Chemkin installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the N2O dissociation mechanism
MyGasMech = ck.Chemistry(label="N2O_dissociation")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = mymechfile
MyGasMech.thermfile = os.path.join(
    data_dir,
    "therm.dat",
)

#######################################
# Preprocess the hydrogen chemistry set
# =====================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

###################################
# Instantiate and set up the stream
# =================================
# Create a diluted air stream by assigning the mole fractions of the
# species.
n2o_mixture = Stream(MyGasMech)
# initial gas state after the reflected shock wave
# 2290 [K]
n2o_mixture.temperature = 2290.0
# AR + N2O ixture composition based on the experiment setup
n2o_mixture.X = [("AR", 0.9899), ("N2O", 0.0101)]
# calculate the gas pressure after the reflected shock from the temperature and density
# density [g/cm3]
state3_den = 0.00015937
n2o_mixture.set_pressure_by_density(rho=state3_den)
# check the pressure value
print(
    f"Gas pressure after the reflected shock = {n2o_mixture.pressure / ck.Patm} [atm]."
)
#
######################################
# Create the shock tube reactor object
# ====================================
# Use the ``ReflectedShock()`` method to create an incident shock reactor.
# The ``ReflectedShock()`` method has two required input parameters.
# The first parameter is the "stream" representing the state of the initial
# gas mixture. In this case, it is the ``n2o_mixture``. The second required parameter
# is the location of the initial gas stream relative to the incident shock front. A value of '1'
# indicate the initial gas stream is before the incident shock front, and a value of '2' indicates
# the gas stream is behind the reflected shock. When the gas condition before the incident shock
# are provided (location = 1), the gas velocity (same as the incident shock velocity)
# must be given by the ``velocity`` method for ``ReflectedShock()`` model.
#

# instantiate the shock tube reactor
# the location '2' means the gas stream is after the reflected shock front
hottube = ReflectedShock(n2o_mixture, location=2, label="reflected_shock")

############################################
# Set up additional reactor model parameters
# ==========================================
# For the reflected shock model, the required reactor parameters is the total simulation time [sec].
# The initial gas mixture conditions are defined by the stream when the ``ReflectedShock`` is instantiated.

# set total simulation time (particle time) [sec]
hottube.time = 3.0e-4

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver related methods, for example,
# ``tolerances``.

# tolerances are given in tuple: (absolute tolerance, relative tolerance)
hottube.tolerances = (1.0e-10, 1.0e-6)

######################################################
# Run the N2O dissociation simulation
# ====================================================
# Use the ``run()`` method to start the reflected shock simulation.
#
# .. note ::
#   You can use two ``time`` calls (one before the run and one after the run) to
#   get the simulation run time (wall time).
#

# set the start wall time
start_time = time.time()
# run the reflected shock tube reactor model
runstatus = hottube.run()
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

# postprocess the solution profiles
hottube.process_solution()

# get the number of solution time points
solutionpoints = hottube.get_solution_size()
print(f"number of solution points = {solutionpoints}")

# get the time profile [sec]
timeprofile = hottube.get_solution_variable_profile("time")
# convert to [msec]
timeprofile *= 1.0e3
# get the temperature profile [K]
tempprofile = hottube.get_solution_variable_profile("temperature")
# get the NO mass fraction profile [-]
NOprofile = hottube.get_solution_variable_profile("NO")
# get the N2O mass fraction profile [-]
N2Oprofile = hottube.get_solution_variable_profile("N2O")
# get the O2 mass fraction profile [-]
O2profile = hottube.get_solution_variable_profile("O2")

# clean up
if os.path.exists(mymechfile):
    os.remove(mymechfile)

############################
# Plot the solution profiles
# ==========================
# Plot the temperature and the NO, N2O, and O2 species
# mass fractions profiles as a function of time.
#

plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("Behind Reflected Shock", fontsize=16)
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, NOprofile, "b-")
plt.ylabel("NO Mass Fraction [-]")
plt.subplot(223)
plt.plot(timeprofile, N2Oprofile, "g-")
plt.xlabel("time [msec]")
plt.ylabel("N2O Mass Fraction [-]")
plt.subplot(224)
plt.plot(timeprofile, O2profile, "m-")
plt.xlabel("time [msec]")
plt.ylabel("O2 Mass Fraction [-]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_reflected_shock.png", bbox_inches="tight")

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

r""".. _ref_rcm:

====================================
Simulate a rapid compression machine
====================================

Ansys Chemkin offers some idealized reactor models commonly used for studying chemical
processes and for developing reaction mechanisms. The *batch reactor* is
a transient 0-D numerical portrayal of the *closed homogeneous/perfectly mixed*
gas-phase reactor. There are two basic types of batch reactor models:

- **constrained-pressure**
- **constrained-volume**

You can choose either to specify the reactor temperature (as a fixed value or
by a piecewise-linear profile) or to solve the energy conservation equation
for each reactor type. In total, you get four variations out of
the base batch reactor model.

**Rapid Compression Machine (RCM)** is often employed to study fuel auto-ignition
at high temperature and high-pressure conditions that are compatible to
the engine-operating environments. The fuel-air mixture inside the RCM chamber is
at relatively low pressure and temperature initially. The gas mixture is then
suddenly compressed causing both the pressure and the temperature of the mixture
to rise rapidly. The reactor/chamber pressure is monitored to identify the onset
of auto-ignition after the compression stopped. This example models the RCM as a
``GivenVolumeBatchReactorEnergyConservation``, and the compression process
is simulated by a predetermined time-volume profile.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_RCM_solution.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

from pathlib import Path

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.core.batchreactors.batchreactor import (
    GivenVolumeBatchReactorEnergyConservation,
)
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
# The mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on the GRI 3.0 mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
MyGasMech.tranfile = str(mechanism_dir / "grimech30_transport.dat")

##############################
# Preprocess the chemistry set
# ============================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()

################################################################
# Set up gas mixtures based on the species in this chemistry set
# ==============================================================
# Use the *equivalence ratio method* so that you can easily set up
# the premixed fuel-oxidizer mixture composition by assigning an
# equivalence ratio value.

# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.x = [("CH4", 1.0)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = 5.0 * ck.P_ATM
fuelmixture.temperature = 1500.0

# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.x = [("O2", 0.21), ("N2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = 5.0 * ck.P_ATM
air.temperature = 1500.0

# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros

# create the premixed mixture to be defined
premixed = ck.Mixture(MyGasMech)

ierror = premixed.x_by_equivalence_ratio(
    MyGasMech, fuelmixture.x, air.x, add_frac, products, equivalenceratio=0.7
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
    print("Error: Failed to create the fuel-oxidizer mixture.")
    exit()

# list the composition of the premixed mixture for verification
premixed.list_composition(mode="mole")

# set mixture temperature and pressure
# (equivalent to setting the initial temperature and pressure of the reactor)
premixed.temperature = 800.0
premixed.pressure = 3.0 * ck.P_ATM

######################################
# Set up the rapid-compression machine
# ====================================
# Create the rapid-compression machine as an instance of the
# ``GivenVolumeBatchReactorEnergyConservation`` object because the reactor volume is
# assigned as a function of time. The batch reactors must be associated with a
# mixture that implicitly links the chemistry set (gas-phase mechanism and properties)
# to the batch reactor. Additionally, it also defines the initial reactor conditions
# (pressure, temperature, volume, and gas composition).

# create a constant volume batch reactor (with energy equation)
MyCONV = GivenVolumeBatchReactorEnergyConservation(premixed, label="RCM")
# show initial gas composition inside the reactor
MyCONV.list_composition(mode="mole")

############################################
# Set up additional reactor model parameters
# ==========================================
# You must provide reactor parameters, solver controls, and output instructions
# before running the simulations. For a batch reactor, the initial volume and the
# simulation end time are required inputs.
#
# .. note::
#   You can reset the initial reactor temperature by using the
#   ``MyCONV.temperature = 800.0`` method. In the run output, you see
#   a warning message about the change.
#

# set other reactor properties
# set initial reactor volume [cm3]
MyCONV.volume = 10.0
# simulation end time [sec]
MyCONV.time = 0.1

########################
# Set the volume profile
# ======================
# Create a time-volume profile by using two arrays. Use
# the ``set_volume_profile()`` method to add the profile to
# the reactor model. The profile data overrides the initial volume
# value set earlier with the ``volume()`` method.

# number of profile data points
npoints = 3
# position array of the profile data
x = np.zeros(npoints, dtype=np.double)
# value array of the profile data
vol_profile = np.zeros_like(x, dtype=np.double)
# set reactor volume data points
x = [0.0, 0.01, 2.0]  # [sec]
vol_profile = [10.0, 4.0, 4.0]  # [cm3]

####################
# Set output options
# ==================
# You can turn on the adaptive solution saving to resolve the steep variations
# in the solution profile. Here additional solution data points are saved for every
# **100 [K]** change in gas temperature. The ``set_ignition_delay()`` method must be
# included for the reactor model to report the ignition delay times after
# the simulation is done. If ``method="T_inflection"`` is set, the reactor model
# treats the inflection points in the predicted gas temperature profile as
# the indication of an auto-ignition. You can choose a different
# auto-ignition definition.
#
# .. note::
#   Type ``ansys.chemkin.core.show_ignition_definitions()`` to get\
#   the list of all available ignition delay time definitions in Chemkin.
#
# .. note::
#   By default, time intervals for both print and save solution are **1/100**
#   of the simulation end time. In this case :math:`dt=time/100=0.001`\ .
#   You can change them to different values.
#

# set the volume profile
MyCONV.set_volume_profile(x, vol_profile)
# output controls
# set timestep between saving solution
MyCONV.timestep_for_saving_solution = 0.01
# turn ON adaptive solution saving
MyCONV.adaptive_solution_saving(mode=True, value_change=100, target="TEMPERATURE")
# specify the ignition definitions
MyCONV.set_ignition_delay(method="T_inflection")

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances.

# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyCONV.tolerances = (1.0e-10, 1.0e-8)
# get solver parameters
atol, rtol = MyCONV.tolerances
print(f"Dfault absolute tolerance = {atol}.")
print(f"Default relative tolerance = {rtol}.")
# turn on the force non-negative solutions option in the solver
MyCONV.force_nonnegative = True
# show solver option
print(f"Timestep between solution printing: {MyCONV.timestep_for_printing_solution}.")
# show timestep between printing solution
print(f"Forced non-negative solution values: {MyCONV.force_nonnegative}.")

#########################################
# Display the added parameters (keywords)
# =======================================
# You can use the ``showkeywordinputlines()`` method to verify the preceding parameters
# are correctly assigned to the reactor model.

# show the additional keywords given by user
MyCONV.showkeywordinputlines()

####################
# Run the simulation
# ==================
# Use the ``run()`` method to start the RCM simulation.

# run the CONV reactor model
runstatus = MyCONV.run()
# check run status
if runstatus != 0:
    # Run failed.
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()

# Run succeeded.
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)

###############################################
# Get the ignition delay time from the solution
# =============================================
# Use the ``get_ignition_delay()`` method to extract the ignition delay time after
# the run is completed.
#
# .. note::
#   You need to deduct the initial compression time = 0.01 [sec] to get
#   the *actual* ignition delay time.
#

# get ignition delay time (need to deduct the initial compression time = 0.01 [sec])
delaytime = MyCONV.get_ignition_delay() - 0.01 * 1.0e3
print(f"Ignition delay time = {delaytime} [msec].")

##########################
# Postprocess the solution
# ========================
# The postprocessing step parses the solution and package the solution values at each
# time point into a mixture. There are two ways to access the solution profiles:
#
# - The raw solution profiles (value as a function of distance) are available
#   for distance, temperature, pressure, volume, and species mass fractions.
#
#  -The mixtures permit the use of all property and rate utilities to extract
#   information such as viscosity, density, and mole fractions.
#
# You can use the ``get_solution_variable_profile()`` method to get
# the raw solution profiles. You can get solution mixtures using either the
# ``get_solution_mixture_at_index()`` method for the solution mixture at
# the given saved location or the ``get_solution_mixture()`` method for
# the solution mixture at the given distance. (In this case, the mixture is
# constructed by interpolation.)
#
# .. note::
#   Use the ``getnumbersolutionpoints()`` method to get the size of
#   the solution profiles before creating the arrays.
#

# postprocess the solutions
MyCONV.process_solution()
# get the number of solution time points
solutionpoints = MyCONV.getnumbersolutionpoints()
print(f"Number of solution points = {solutionpoints}.")

# easily access raw solution profiles
# get the time profile
timeprofile = MyCONV.get_solution_variable_profile("time")
# get the temperature profile
tempprofile = MyCONV.get_solution_variable_profile("temperature")
# get the volume profile
volprofile = MyCONV.get_solution_variable_profile("volume")

# more involved postprocessing using mixtures
# reactor mass
massprofile = np.zeros_like(timeprofile, dtype=np.double)
# create arrays for CH4 mole fraction, CH4 ROP, and mixture viscosity
ch4_profile = np.zeros_like(timeprofile, dtype=np.double)
ch4_rop_profile = np.zeros_like(timeprofile, dtype=np.double)
viscprofile = np.zeros_like(timeprofile, dtype=np.double)
current_rop = np.zeros(MyGasMech.kk, dtype=np.double)
# find CH4 species index
ch4_index = MyGasMech.get_specindex("CH4")

# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = MyCONV.get_solution_mixture_at_index(solution_index=i)
    # get gas density [g/cm3]
    den = solutionmixture.rho
    # reactor mass [g]
    massprofile[i] = den * volprofile[i]
    # get CH4 mole fraction profile
    ch4_profile[i] = solutionmixture.x[ch4_index]
    # get CH4 ROP profile
    current_rop = solutionmixture.rop()
    ch4_rop_profile[i] = current_rop[ch4_index]
    # get mixture vicosity profile
    viscprofile[i] = solutionmixture.mixture_viscosity()

################################
# Validate the simulation result
# ==============================
# Since the RCM is a closed reactor, the total gas mass inside RCM must be
# kept constant. You can verify it by computing the maximum mass deviation
# in the solution profile. If you find the mass variations are too large to
# be acceptable, you can set smaller tolerance values and/or adjust
# some solver parameters such as the maximum solver time step size
# and re-run the simulation.
del_mass = np.zeros_like(timeprofile, dtype=np.double)
mass0 = massprofile[0]
for i in range(solutionpoints):
    del_mass[i] = abs(massprofile[i] - mass0)
#
print(f">>> Maximum magnitude of reactor mass deviation = {np.max(del_mass)} [g].")

################################
# Plot the RCM solution profiles
# ==============================

# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, ch4_profile, "b-")
plt.ylabel("CH4 Mole Fraction")
plt.subplot(223)
plt.plot(timeprofile, ch4_rop_profile, "g-")
plt.xlabel("time [sec]")
plt.ylabel("CH4 Production Rate [mol/cm3-sec]")
plt.subplot(224)
plt.plot(timeprofile, viscprofile, "m-")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Viscosity [g/cm-sec]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_RCM_solution.png", bbox_inches="tight")

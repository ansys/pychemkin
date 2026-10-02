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

r""".. _ref_closed_homogeneous:

===========================================================
Simulate hydrogen combustion in a constant-pressure reactor
===========================================================

Ansys Chemkin offers some idealized reactor models commonly used for studying chemical
processes and for developing reaction mechanisms. The batch reactor is a transient 0-D
numerical portrayal of the closed homogeneous/perfectly mixed gas-phase reactor.
There are two basic types of batch reactor models:

- **constrained-pressure**
- **constrained-volume**

You can choose either to specify the reactor temperature (as a fixed value or by a
piecewise-linear profile) or to solve the energy conservation equation
for each reactor type. In total, you get four variations out of
the base batch reactor model.

This example models the ignition of a stoichiometric hydrogen-air mixture
(\ :math:`\phi = 1`\ ) in a balloon (constant-pressure) reactor. The reactor is created
as an instance of the ``GivenPressureBatchReactorEnergyConservation`` object.
The initial gas mixture in the batch reactor is set by the properties of
the hydrogen-air mixture. You get the ignition delay time (if auto-ignition of
the hydrogen-air mixture occurs during the simulation time) and plot the predicted
profiles of gas temperature, gas density, H\ :sub:`2`\ O mole fraction, and
the net production rate of H\ :sub:`2`\ O.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_close_homogeneous.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.core.batchreactors.batchreactor import (
    GivenPressureBatchReactorEnergyConservation,
)
from ansys.chemkin.core.logger import logger
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

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
# create a chemistry set based on the diesel 14-components mechanism
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
# Compose the species molar ratios of the fuel-air mixture in a
# recipe list and use the ``fuelmixture.x()`` method to set the mixture
# composition (mole fractions) directly.
#
# Since you are going to use the ``fuelmixture`` mixture to instantiate
# the reactor object later, setting the mixture pressure and temperature
# is equivalent to setting the initial temperature and pressure of the
# batch reactor.

# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition (mole ratios)
fuelmixture.x = [("H2", 2.0), ("N2", 3.76), ("O2", 1.0)]
# setting the mixture pressure and temperature is equivalent to setting
# the initial temperature and pressure of the reactor in this case
fuelmixture.pressure = ck.P_ATM
fuelmixture.temperature = 1000

# list the composition of the premixed mixture
# this serves as the baseline for verification later
fuelmixture.list_composition(mode="mole")

#################################################################
# Set up a constant-pressure batch reactor (with energy equation)
# ===============================================================
# Create the constant-pressure batch reactor as an instance of the
# ``GivenPressureBatchReactorEnergyConservation`` object
# because the reactor pressure is kept constant
# (or assigned as a function of time). The batch reactors must be
# associated with a mixture, which implicitly links the chemistry set
# (gas-phase mechanism and properties) to the batch reactor. Additionally,
# it defines the initial conditions (pressure, temperature, volume,
# and gas composition) of the batch reactor.
MyCONP = GivenPressureBatchReactorEnergyConservation(fuelmixture, label="tran")

##############################
# List the mixture composition
# ============================
# List the initial gas composition inside the reactor for verification. You can
# use the composition printout from the ``fuelmixture`` mixture to verify that
# the initial gas composition inside ``MyCONP`` is set correctly, that is,
# ``MyCONP`` has been initialized by the ``fuelmixture`` mixture.
MyCONP.list_composition(mode="mole")

############################################
# Set up additional reactor model parameters
# ==========================================
# Before you can run the simulation, you must provide reactor parameters,
# solver controls, and output instructions. For a batch reactor,
# the initial volume and the simulation end time are required inputs.

# set other reactor properties
# reactor volume [cm3]
MyCONP.volume = 1
MyCONP.temperature = 1000
# simulation end time [sec]
MyCONP.time = 0.0005

####################
# Set output options
# ==================
# You can turn on adaptive solution saving to resolve the steep variations
# in the solution profile. Here, additional solution data points are saved
# for every 20 internal solver steps. You must include the ``set_ignition_delay()``
# method for the reactor model to report the ignition delay times after
# the simulation is done. If ``method="T_rise"`` is set, the reactor model considers
# the gas is auto-ignited when the predicted gas temperature goes above
# the initial temperature by the amount indicated by the parameter ``val=400``.
# You can choose a different auto-ignition definition.
#
# .. note::
#
#   - Type ``ansys.chemkin.core.show_ignition_definitions()`` to get the list of
#     all available ignition delay time definitions in Chemkin.
#
#   - By default, time intervals for both print and save solution are 1/100 of the
#     simulation end time, which in this example is :math:`dt=time/100=0.001`\ .
#     You can change them to different values.
#

# turn on adaptive solution saving
MyCONP.adaptive_solution_saving(mode=True, steps=20)
# specify the ignition definitions
ck.show_ignition_definitions()
MyCONP.set_ignition_delay(method="T_rise", val=400)

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances.

# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyCONP.tolerances = (1.0e-20, 1.0e-8)

# get solver parameters
ATOL, RTOL = MyCONP.tolerances
print(f"Default absolute tolerance = {ATOL}.")
print(f"Default relative tolerance = {RTOL}.")
# turn on the force non-negative solutions option in the solver
MyCONP.force_nonnegative = True
# show solver option
print(f"Timestep between solution printing: {MyCONP.timestep_for_printing_solution}")
# show timestep between printing solution
print(f"Forced non-negative solution values: {MyCONP.force_nonnegative}")

#########################################
# Display the added parameters (keywords)
# =======================================
# Use the ``showkeywordinputlines()`` method to verify that the preceding
# parameters are correctly assigned to the reactor model.
MyCONP.showkeywordinputlines()

####################
# Run the simulation
# ==================
# Use the ``run()`` method to start the batch reactor simulation.
runstatus = MyCONP.run()

# check run status
if runstatus != 0:
    # Run failed.
    print(Color.RED + ">>> Run failed. <<<", end="\n" + Color.END)
    exit()
# Run succeeded.
print(Color.GREEN + ">>> Run completed. <<<", end="\n" + Color.END)

###############################################
# Get the ignition delay time from the solution
# =============================================
# Use the ``get_ignition_delay()`` method to extract the ignition delay time after
# the run is completed.
delaytime = MyCONP.get_ignition_delay()
print(f"Ignition delay time = {delaytime} [msec].")

##########################
# Postprocess the solution
# ========================
# The postprocessing step parses the solution and packages the solution values at each
# time point into a mixture object. There are two ways to access the solution profiles:
#
# - The raw solution profiles (value as a function of time) are available for time,
#   temperature, pressure, volume, and species mass fractions.
#
# - The mixture objects that permit the use of all property and rate utilities
#   to extract information such as viscosity, density, species production rates,
#   and mole fractions.
#
# To obtain the raw solution profiles, you can use
# the ``get_solution_variable_profile()`` method. To obtain
# the solution mixture objects, you can use either
# the ``get_solution_mixture_at_index()`` method for the solution mixture at
# a given time point or the ``get_solution_mixture()`` method for the solution mixture
# at a given time. (In this case, the mixture is constructed by # interpolation.)
#
# .. note::
#   Use the ``getnumbersolutionpoints()`` method to get the size of
#   the solution profiles before creating the arrays.
#
MyCONP.process_solution()
# get the number of solution time points
solutionpoints = MyCONP.getnumbersolutionpoints()
print(f"Number of solution points = {solutionpoints}.")

# get the time profile
timeprofile = MyCONP.get_solution_variable_profile("time")
# get the temperature profile
tempprofile = MyCONP.get_solution_variable_profile("temperature")
# more involving postprocessing by using mixtures
# create arrays for H2O mole fraction, H2O ROP, and mixture density
h2o_profile = np.zeros_like(timeprofile, dtype=np.double)
h2o_rop_profile = np.zeros_like(timeprofile, dtype=np.double)
denprofile = np.zeros_like(timeprofile, dtype=np.double)
current_rop = np.zeros(MyGasMech.kk, dtype=np.double)
# find H2O species index
h2o_index = MyGasMech.get_specindex("H2O")
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = MyCONP.get_solution_mixture_at_index(solution_index=i)
    # get gas density [g/cm3]
    denprofile[i] = solutionmixture.rho
    # reactor mass [g]
    # get H2O mole fraction profile
    h2o_profile[i] = solutionmixture.x[h2o_index]
    # get H2O ROP profile
    current_rop = solutionmixture.rop()
    h2o_rop_profile[i] = current_rop[h2o_index]

############################
# Plot the solution profiles
# ==========================

# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, h2o_profile, "b-")
plt.ylabel("H2O Mole Fraction")
plt.subplot(223)
plt.plot(timeprofile, denprofile, "m-")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Density [g/cm3]")
plt.subplot(224)
plt.plot(timeprofile, h2o_rop_profile, "g-")
plt.xlabel("time [sec]")
plt.ylabel("H2O Production Rate [mol/cm3-sec]")
# display the plots
if interactive:
    plt.show()
else:
    plt.savefig("plot_close_homogeneous.png", bbox_inches="tight")

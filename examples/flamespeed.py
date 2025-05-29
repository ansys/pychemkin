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
.. _ref_hydrogen_flame_speed:

========================================================================
Compute the laminar flame speed of diluted hydrogen at low pressure
========================================================================

The *freely propagating premixed flame* model can be utilized to obtain the
laminar *flame speed* of a gas mixture at given initial temperature and pressure.
The premixed flame model will calculate the temperature and species composition profiles
across the flame, and the laminar flame speed will be derived from the mass flow rate of
the gas mixture, a constant across the entire calculation domain because of the mass
conservation.

This tutorial demonstrates the application of the "freely propagating" premixed flame model
to calculate the laminar flame speed of an N\ :sub:`2` diluted hydrogen-air mixture at
low pressure. Since the transport processes are critical for flame calculations, the
*transport data* must be included in the mechanism data and pre-processed.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_hydrogen_premixed_flame.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.premixedflames.premixedflame import FreelyPropagating as FlameSpeed
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

#####################################
# Create a chemistry set
# ===================================
# The 'C2 NOx' mechanism is from the default *"/reaction/data"* directory.
# This mechanism also includes information about the *Soave* cubic
# Equation of State (EOS) for the real-gas applications. PyChemkin preprocessor
# will indicate the availability of the real-gas model in the ``Chemistry Set`` processed.
#
# .. note::
#   The transport data *must* be included and pre-processed because the transport processes,
#   *convection and diffusion*, are important to sustain the flame structure.
#
# .. note::
#   Use the ``preprocess_transportdata`` method if the transport data is embedded in the
#   gas-phase mechanism file.
#

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the C2 NOx mechanism
MyGasMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")

# direct the preprocessor to include the transport properties
# only when the mechanism file contains all the transport data
MyGasMech.preprocess_transportdata()

############################################
# Pre-process the ``Chemistry Set``
# ==========================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()
if iError != 0:
    print("Error: failed to preprocess the mechanism!")
    print(f"       error code = {iError}")
    exit()

#######################################################################################
# Set up the (H\ :sub:`2` + N\ :sub:`2`\ )-air mixture for the flame speed calculation
# =====================================================================================
# Instantiate an ``Stream`` object ``premixed`` for the inlet gas mixture.
# The ``Stream`` object is a ``Mixture`` object with the addition of the
# *inlet flow rate*. You can specify the inlet gas properties the same way you
# set up a ``Mixture``. Here the ``X_by_Equivalence_Ratio`` method is used.
# You create the ``fuel`` and the ``air`` mixtures first. Then define the
# *complete combustion product species* and provide the *additives* composition
# if applicable. And finally you can simply set ``equivalenceratio=1`` to create
# the stoichiometric hydrogen-air mixture. The estimated inlet mass flow rate can be
# assigned by the ``mass_flowrate`` method.

# create the fuel mixture
fuel = ck.Mixture(MyGasMech)
# set fuel composition: hydrogen diluted by nitrogen
fuel.X = [("H2", 0.7), ("N2", 0.3)]
# setting pressure and temperature is not required in this case
fuel.pressure = 0.0125 * ck.Patm
fuel.temperature = 300.0  # inlet temperature

# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.X = ck.Air.X()
# setting pressure and temperature is not required in this case
air.pressure = fuel.pressure
air.temperature = fuel.temperature

# create the fuel-air Stream for the premixed flame speed calculation
premixed = Stream(MyGasMech, label="premixed")
# products from the complete combustion of the fuel mixture and air
products = ["H2O", "N2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros
# mean equivalence ratio
equiv = 1.0
iError = premixed.X_by_Equivalence_Ratio(
    MyGasMech, fuel.X, air.X, add_frac, products, equivalenceratio=equiv
)
# check fuel-oxidizer mixture creation status
if iError != 0:
    print("Error: failed to create the hydrogen-air mixture!")
    exit()

# setting inlet pressure [dynes/cm2]
premixed.pressure = fuel.pressure
# set inlet/unburnt gas temperature [K]
premixed.temperature = fuel.temperature
# set estimated value of the flame speed [cm/sec]
premixed.velocity = 65.0

##################################################
# Instantiate the laminar speed calculator
# ================================================
# Set up the *freely propagating premixed flame* model by using the ``Stream``
# object representing the premixed fuel-oxidizer mixture (with the estimated
# mass flow rate value). When the flame speed is expected to be very small
# (< 10 [cm/sec]) or very large (> 300 [cm/sec]), it might be beneficial to
# set the ``mass_flowrate`` of this inlet to the mass flux [g/cm\ :sup:`2`\ -sec]
# based on the estimated flame speed value. There are many options and parameters
# related to the treatment of the species boundary condition, the transport properties.
# All the available options and parameters are described in the *Chemkin* *Input* manual.
#
# .. note::
#   The ``Stream`` parameter used to instantiate a ``FlameSpeed`` object is
#   the properties of the unburned fuel-oxidizer mixture of which the *laminar
#   flame speed* will be determined.
#

flamespeedcalculator = FlameSpeed(premixed, label="premixed_hydrogen")

###############################################
# Set up initial mesh and grid adaption options
# =============================================
# The premixed flame models provides several methods to set up the initial
# mesh.  Here a uniform mesh of 11 grid points is used at the start of the simulation.
# The flame models would add more grid points to where they are needed as determined by
# the solution quality parameters specified by the ``set_solution_quality`` method.
#
# The ``end_poistion`` is a required input as it defines the length of the calculation domain.
# Typically, the length of the calculation domain is between 1 to 10 [cm]. For low pressure
# conditions, the flame thickness becomes wider and a larger calculation domain is required.
#
# .. note::
#   There are three methods to set up the initial mesh for the premixed flame calculations:
#
#   1. ``use_TPRO_grids`` method (default) to use the grid points in the estimate temperature profile.
#
#   2. ``set_numb_grid_points`` method to create a uniform mesh of the given number of grid points.
#
#   3. ``set_grid_profile`` method to specify the initial grid point profile.
#

# set the initial mesh to 35 uniformly distributed grid points
flamespeedcalculator.set_numb_grid_points(35)
# set the maximum total number of grid points allowed in the calculation (optional)
flamespeedcalculator.set_max_grid_points(150)
# define the calculation domain [cm]
flamespeedcalculator.end_position = 40.0
# maximum number of grid points can be added during each grid adaption event (optional)
flamespeedcalculator.set_max_adaptive_points(20)
# set the maximum values of the grdient and the curvature of the solution profiles (optional)
flamespeedcalculator.set_solution_quality(gradient=0.1, curvature=0.2)

#################################
# Set transport property options
# ===============================
# *Chemkin* offers three methods to compute the *mixture* properties: the *mixture averaged*
# method, the *multi-component* method, and the constant *Lewis number* method. When the system
# pressure is not too low, the *mixture averaged* method should be adequate. The *multi-component*
# method, although it is slightly more accurate, would make the simulation time longer and harder
# to converge. Using the constant *Lewis number* method implies that all the species would have
# the same transport properties.
#
# Include the thermal diffusion effect, when there are large amount of light species
# (molecular weight < 5.0).
#

# use the mixture-averaged formulism to evaluate the mixture transport properties
flamespeedcalculator.use_mixture_averaged_transport()
# include the thermal diffusion effect (because the unburned mixture has hydrogen (molecular weight < 5.0))
flamespeedcalculator.use_thermal_diffusion(mode=True)

#########################################
# Set species composition boundary option
# =======================================
# There two types of boundary condition treatments for the species composition available
# from the premixed flame models: *"comp"* and *"flux"*. You can find the descriptions of
# these two treatments in the *Chemkin* *Input* manual.

# specific the species composition boundary treatment ("comp" or "flux")
# use "COMP" to keep the inlet species mass fraction values the same as the "given inlet stream".
flamespeedcalculator.set_species_boundary_types(mode="comp")

############################
# Set solver parameters
# ==========================
# The steady state solver parameters for the premixed flame model are optional because all the
# solver parameters have their own default values. Change the solver parameters when the premixed
# flame simulation does not converge with the default settings.
#
# .. note::
#   The ``FreelyPropagating`` flame speed calculator has an option to automatically generate
#   an estimate temperature profile that might improve the convergence performance. Use the
#   ``automatic_temperature_profile_estimate`` method to turn this option on.
#

# reset the tolerances in the steady-state solver (optional)
flamespeedcalculator.steady_state_tolerances = (1.0e-9, 1.0e-6)
flamespeedcalculator.timestepping_tolerances = (1.0e-6, 1.0e-4)
# reset the gas species floor value in the steady-state solver (optional)
flamespeedcalculator.set_species_floor(-1.0e-4)
# skip the fixed-temperature step (optional)
flamespeedcalculator.skip_fix_T_solution(mode=True)
# reduce the Jacobian age during the pseudo time stepping phase
flamespeedcalculator.set_pseudo_Jacobian_age(10)

####################################
# Run the premixed flame calculation
# ==================================
# Use the ``run`` method to run the freely propagating premixed flame (flame speed) model.
# will solve the reactors one by one in the order they are added to the network. After the
# premixed flame calculation concludes successfully, use the ``process_solution`` method to
# post-process the solutions. The predicted laminar flame speed can be obtained by using the
# ``get_flame_speed`` method. You can create other property profiles by looping through the
# solution ``Streams`` with proper ``Mixture`` methods.
#
# .. note::
#   When the inlet stream condition is close to the flammability limit, the flame speed
#   calculation might fail. Remember that the reaction mechanism (reaction rates, thermodynamic
#   properties, and transport properties) and the reactor model are *models* that contain assumptions
#   and uncertainties.
#

# set the start wall time
start_time = time.time()

status = flamespeedcalculator.run()
if status != 0:
    print(Color.RED + "failed to calculate the laminar flame speed!" + Color.END)
    exit()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"total simulation duration: {runtime} [sec]")
print()

############################################
# Post-process the premixed flame results
# ==========================================
# The post-processing step will parse the solution and package the solution values at each
# time point into a ``Stream`` object. There are two ways to access the solution profiles:
#
#   1. the "raw" solution profiles (value as a function of time) are available for "distance",
#   "temperature", and species "mass fractions";
#
#   2. the ``Stream`` objects that permit the use of all property and rate utilities to extract
#   information such as viscosity, density, and mole fractions.
#
# The "raw" solution profiles can be obtained by using the ``get_solution_variable_profile`` method. The
# solution ``Stream`` objects are accessed via either the ``get_solution_stream_at_grid`` for the
# solution stream at the given *grid point* or the ``get_solution_stream`` for the solution stream
# at the given *location* (in this case, the "stream" is constructed by interpolation).
#
# .. note::
#   Use the ``get_solution_size`` to get the number of grid pints in the solution profiles before
#   creating the arrays.
#
# .. note::
#   The ``mass_flowrate`` from the solution ``Streams`` is actually the *mass flux* [g/cm\ :sup:`2`\ -sec].
#   It can used to derive the velocity at the corresponding location by dividing it by the local gas mixture
#   density [g/cm\ :sup:`3`\ ]. Here the ``get_flame_speed`` is used to retrieve the predict laminar flame
#   speed value.
#

# post-process the solutions
flamespeedcalculator.process_solution()

# print the predicted laminar flame speed
print(
    f"the predicted laminar flame speed = {flamespeedcalculator.get_flame_speed()} [cm/sec]"
)
# get the number of solution grid points
solutionpoints = flamespeedcalculator.get_solution_size()
print(f"number of solution points = {solutionpoints}")
# get the grid profile
mesh = flamespeedcalculator.get_solution_variable_profile("distance")
# get the temperature profile
tempprofile = flamespeedcalculator.get_solution_variable_profile("temperature")
# get HO2 mass fraction profile
HO2profile = flamespeedcalculator.get_solution_variable_profile("HO2")

# create arrays for mixture density, mixture viscosity, and mixture specific heat capacity
denprofile = np.zeros_like(mesh, dtype=np.double)
viscprofile = np.zeros_like(mesh, dtype=np.double)
# loop over all solution grid points
for i in range(solutionpoints):
    # get the stream at the grid point
    solutionstream = flamespeedcalculator.get_solution_stream_at_grid(grid_index=i)
    # get gas density [g/cm3]
    denprofile[i] = solutionstream.RHO
    # get mixture viscosity profile [g/cm-sec] or [Poise]
    viscprofile[i] = solutionstream.mixture_viscosity() * 1.0e2

###########################################
# Plot the premixed flame solution profiles
# =========================================
# Plot the solution profiles of the premixed flame.
#
# .. note ::
#   You can get profiles of the thermodynamic and the transport properties
#   by applying ``Mixture`` utility methods to the solution ``Stream``.
#
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(mesh, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(mesh, denprofile, "b-")
plt.ylabel("Mixture Density [g/cm3]")
plt.subplot(223)
plt.plot(mesh, HO2profile, "g-")
plt.xlabel("Distance [cm]")
plt.ylabel("HO2 Mass Fraction")
plt.subplot(224)
plt.plot(mesh, viscprofile, "m-")
plt.xlabel("Distance [cm]")
plt.ylabel("Mixture Viscosity [cP]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_hydrogen_premixed_flame.png", bbox_inches="tight")

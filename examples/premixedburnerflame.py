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
.. _ref_premixed_burner_stabilized_flame:

========================================================================
Simulate a burner stabilized premixed flame
========================================================================

The *burner stabilized premixed flame* model is frequently used as a numerical tool
to develop and to validate combustion mechanisms in the context of a burner stabilized flat
flame. The premixed burner stabilized flat flame experiments are important because, in addition
to the reaction pathways and rate parameters, they provide opportunities to learn about
the impact of the transport processes on the combustion chemistry and flame structure. The
premixed burner stabilized flame model calculates the temperature and the species concentrations
along the burner centerline, and the results will be compared against the measurements. The
burner stabilized flame model assumes that the flow is laminar and the temperature at the burner
rim is constant (in order to anchor the flame).

This tutorial demonstrates the application of the "burner stabilized" premixed flame model
to calculate the species profiles of a C\ :sub:`2`\ H\ :sub:`4`\ -O\ :sub:`2`\ -AR mixture along
the burner centerline with the measured temperature profile data. Since the transport
processes are critical for flame calculations, the *transport data* must be included in the mechanism
data and pre-processed.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_premixed_burner_stabilised_flame.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger

# chemkin perfectly-stirred reactor (PSR) model (steady-state)
from ansys.chemkin.premixedflames.premixedflame import (
    BurnedStabilized_GivenTemperature as Burner,
)
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

####################################################
# Create an instance of the first  ``Chemistry Set``
# ==================================================
# The mechanism loaded is the GRI 3.0 mechanism for methane combustion.
# The mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"/reaction/data"*.
#
# .. note::
#   The transport data *must* be included and pre-processed because the transport processes,
#   *convection and diffusion*, are important to sustain the flame structure.
#

# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# inclusion of the full file path is recommended
chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

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
# Set up the (C\ :sub:`2`\ H\ :sub:`4` + O\ :sub:`2`\ + AR ) mixture to the flat burner
# =====================================================================================
# Instantiate an ``Stream`` object ``premixed`` for the inlet gas mixture.
# The ``Stream`` object is a ``Mixture`` object with the addition of the
# *inlet flow rate*. You can specify the inlet gas properties the same way you
# set up a ``Mixture``. You can simply use a composition *"recipe"* to create
# the C\ :sub:`2`\ H\ :sub:`4` + O\ :sub:`2`\ + AR mixture. The burner inlet mass flow
# rate can be assigned by the ``mass_flowrate`` method.
#
# .. note::
#   By default, the burner flow area is set to unity (1.0 [cm\ :sup:`2`\ ]) because the flame
#   models use the "mass flux" [g/cm\ :sup:`2`\ -sec] in the calculations instead of
#   the "flow rate" [g/sec]. If the mass flow rate is used, the actual burner cross-sectional
#   flow area must be provided by using the ``Stream`` method ``flowarea``.
#

# create the fuel-air Stream for the premixed flame speed calculation
premixed = Stream(MyGasMech, label="premixed")
# set inlet premixed stream molar composition
premixed.X = [("C2H4", 0.163), ("O2", 0.237), ("AR", 0.6)]
# setting inlet pressure [dynes/cm2]
premixed.pressure = 1.0 * ck.Patm
# set inlet/unburnt gas temperature [K]
premixed.temperature = 300.0  # inlet temperature

# set the burner outlet cross sectional flow area to unity (1 [cm2])
# because the flame models use the "mass flux" in the calculation instead of the mass flow rate
premixed.flowarea = 1.0
# set value for the inlet mass flow rate [g/sec]
premixed.velocity = 0.0147

##################################################
# Instantiate the laminar flat flame burner
# ================================================
# Set up the *premixed burner stabilized flame* model by using the ``Stream``
# object representing the premixed fuel-oxidizer mixture.There are many options
# and parameters related to the treatment of the species boundary condition and
# the transport properties. All the available options and parameters are described
# in the *Chemkin* *Input* manual.
#
# .. note::
#   The ``Stream`` parameter used to instantiate a ``Burner`` object is
#   the properties of the unburned fuel-oxidizer mixture leaving the burner outlet.
#

flatflame = Burner(premixed, label="ethylene flame")

############################################################
# Set up the temperature profile along the burner centerline
# ==========================================================
# Use the ``set_temperature_profile`` method to set up the temperature
# profile data. The temperature profile along the burner centerline is
# typically obtained from the corresponding experiments.
#

TPRO_data_points = 21
gridpoints = np.zeros(TPRO_data_points, dtype=np.double)
temp_data = np.zeros_like(gridpoints, dtype=np.double)
gridpoints = [
    0.0,
    0.01,
    0.02,
    0.0345643,
    0.0602659,
    0.100148,
    0.118759,
    0.135598,
    0.15421,
    0.179911,
    0.211817,
    0.233087,
    0.260561,
    0.293353,
    0.348301,
    0.436928,
    0.517578,
    0.7161,
    0.935894,
    1.16632,
    1.2,
]
temp_data = [
    300.0,
    450.0,
    600.0,
    828.498,
    1104.4,
    1496.7,
    1634.65,
    1714.85,
    1768.33,
    1801.12,
    1807.2,
    1803.78,
    1799.51,
    1788.35,
    1778.09,
    1767.01,
    1760.23,
    1744.13,
    1737.55,
    1730.99,
    1729.31,
]
flatflame.set_temperature_profile(x=gridpoints, temp=temp_data)

###############################################
# Set up initial mesh and grid adaption options
# =============================================
# The premixed flame models provides several methods to set up the initial
# mesh.  In this case, the temperature profile is given by the experimental data.
# The ``use_TPRO_grid`` method will let the flame model recycle the grid points of
# the temperature profile as the initial mesh at the start of the calculation.
# The flame models would add more grid points to where they are needed as determined by
# the solution quality parameters specified by the ``set_solution_quality`` method.
#
# The ``end_poistion`` is a required input as it defines the length of the calculation domain.
# Typically, the length of the calculation domain is between 1 to 10 [cm].
#

# set the maximum total number of grid points allowed in the calculation (optional)
flatflame.set_max_grid_points(250)
# define the calculation domain [cm]
flatflame.end_position = 1.2
# maximum number of grid points can be added during each grid adaption event (optional)
flatflame.set_max_adaptive_points(10)
# set the maximum values of the grdient and the curvature of the solution profiles (optional)
flatflame.set_solution_quality(gradient=0.1, curvature=0.1)

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

# use the mixture average formulism to evaluate the mixture transport properties
flatflame.use_mixture_averaged_transport()

#########################################
# Set species composition boundary option
# =======================================
# There two types of boundary condition treatments for the species composition available
# from the premixed flame models: *"comp"* and *"flux"*. You can find the descriptions of
# these two treatments in the *Chemkin* *Input* manual.

# specific the species composition boundary treatment ("comp" or "flux")
# use "FLUX" to ensure that the "net" species mass fluxes are zero at the burner outlet.
flatflame.set_species_boundary_types(mode="flux")

############################
# Set solver parameters
# ==========================
# The steady state solver parameters for the premixed flame model are optional because all the
# solver parameters have their own default values. Change the solver parameters when the premixed
# flame simulation does not converge with the default settings.
#

# reset the tolerances in the steady-state solver (optional)
flatflame.steady_state_tolerances = (1.0e-9, 1.0e-4)
# reset the gas species floor value in the steady-state solver (optional)
flatflame.set_species_floor(-1.0e-3)

####################################
# Run the premixed flame calculation
# ==================================
# Use the ``run`` method to run the freely propagating premixed flame (flame speed) model.
# will solve the reactors one by one in the order they are added to the network. After the
# premixed flame calculation concludes successfully, use the ``process_solution`` method to
# post-process the solutions. You can create other property profiles by looping through the
# solution ``Streams`` with proper ``Mixture`` methods.
#

# set the start wall time
start_time = time.time()

status = flatflame.run()
if status != 0:
    print(Color.RED + "failed to solve the reactor network!" + Color.END)
    exit()

# get the number of solution grid points
solutionpoints = flatflame.get_solution_size()
print(f"number of solution points = {solutionpoints}")

#################################
# Refine the solution profiles
# ===============================
# When the simulation is hard to converge in one go, you can use a number
# of *"continuation"* to gradually get to the solution of the intended condition.
# For example, you can get to the solution of a very fuel-lean mixture by starting with
# a slightly fuel-rich mixture flame and use the ``continuation`` method to run the more
# difficult fuel-lean case from the converged solution.
#

# tightening the converge criteria
flatflame.set_solution_quality(gradient=0.05, curvature=0.1)

# re-start the flame simulation by continuation
flatflame.continuation()

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
#   density [g/cm\ :sup:`3`\ ].
#

# post-process the solutions
flatflame.process_solution()

# get the number of solution grid points
solutionpoints = flatflame.get_solution_size()
print(f"number of final solution points = {solutionpoints}")
# get the grid profile
mesh = flatflame.get_solution_variable_profile("distance")
# get the temperature profile
tempprofile = flatflame.get_solution_variable_profile("temperature")
# get OH mass fraction profile
OHprofile = flatflame.get_solution_variable_profile("OH")

# create arrays for mixture conductivity, and mixture specific heat capacity
Cpprofile = np.zeros_like(mesh, dtype=np.double)
condprofile = np.zeros_like(mesh, dtype=np.double)
# loop over all solution grid points
for i in range(solutionpoints):
    # get the stream at the grid point
    solutionstream = flatflame.get_solution_stream_at_grid(grid_index=i)
    # get mixture specific heat capacity profile [erg/mole-K]
    Cpprofile[i] = solutionstream.CPBL() / ck.ergs_per_joule * 1.0e-3
    # get thermal conductivity profile [ergs/cm-K-sec]
    condprofile[i] = solutionstream.mixture_conductivity() * 1.0e-5

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
plt.plot(mesh, Cpprofile, "b-")
plt.ylabel("Mixture Cp [kJ/mole]")
plt.subplot(223)
plt.plot(mesh, OHprofile, "g-")
plt.xlabel("Distance [cm]")
plt.ylabel("OH Mass Fraction")
plt.subplot(224)
plt.plot(mesh, condprofile, "m-")
plt.xlabel("Distance [cm]")
plt.ylabel("Mixture conductivity [W/m-K]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_premixed_burner_stabilised_flame.png", bbox_inches="tight")

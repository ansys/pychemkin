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
.. _ref_opposed_flow_flame:

========================================================================
Study flame interactions in an opposed-flow flame configuration
========================================================================

The *opposed-flow flame* model is frequently used as a numerical tool
to study the flame structures under fluid dynamic stress without the interferences from the wall.
The opposed-flow flame experiments provide insights on the impact of the transport processes
(mainly diffusion) on the combustion chemistry and flame structure. They also reveal
how the flame structure would response to the strain rate imposed by the mean flow field (or to
some extent by the turbulence).

 .. figure:: opposed_flow_flame.png
   :scale: 80 %
   :alt: opposed-flow flame configuration

Typically, the "fuel" and the "oxidizer" are introduced separately
from the two opposing inlets, and a non-premixed (diffusion) flame would settle around the stagnating plan
in the middle of the separation gap. Since the axial velocity near the stagnating plan is small, the flame is
sustained mainly by the balance of combustion chemistry and the species/heat diffusion. The flame is stretched
(strained) due to the radial velocity so changing the flow rates at the inlets will vary the stretching/strain rate
applied to the flame. It is possible to establish more than one flame in the separation gap. For example, this project
forms two flames by introducing a fuel-rich mixture from the "fuel" inlet. Alternatively, replacing the "oxidizer" with
a fuel-lean mixture or replacing both inlet streams with premixed mixtures will also yields two flames. Having a fuel-rich
mixture for the "fuel" and a fuel-lean mixture for the "oxidizer" might create three flames: one rich premixed flame, one
lean premixed flame and one diffusion flame between them.

The flame model calculates the temperature, velocities, and the species concentrations
along the centerline between the two opposing nozzles, and the results can be presented graphically.
The mixture fraction profile is also available. The opposed-flow flame model evaluates the mixture fraction
by using Bilger's elemental fraction formulation, therefore, the mixture fraction could be negative or have
values greater than 1. Since there is no commonly recognized strain rate definition, you can derive your own strain rate
from the velocity solution.

This tutorial demonstrates the application of the axisymmetric opposed-flow flame model
to study the interactions (species and heat) between two strained flames, a fuel-rich premixed flame
and a diffusion flame, situated in the space between the two inlets. The computed temperature profile does not appear to
show two distinct flame zones. However, by plotting profile of radicals such as NO\ :sub:`2`\ , the locations of the two
flame fronts are clearly marked. The species and temperature gradients between the two flames indicate the two flames
support each other by exchanging reactants and heat. If you plot the hydrogen cyanide (HCN) profile, you would see that
the rich premixed flame can produce NOx through the Fenimore route (prompt NOx); while the oxygen atom (O) profile would
further reveal that, while both flames form thermal NOx (Zeldovich route) , most of the thermal NOx is likely emitted from
the diffusion flame. The nitrogen dioxide (NO\ :sub:`2`\ ) formed in the flames will convert to nitric oxide (NO) in
the hot zone sandwiched by the two flames.

Since the transport processes are critical for flame calculations, the *transport data* must be included in the mechanism
data and pre-processed.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_opposed_flow_flame.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color

# Chemkin 1-D opposed-flow flame model (steady-state)
from ansys.chemkin.diffusionflames.opposedflowflame import OpposedFlame as Flame
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting

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

##########################################
# Create an instance of the Chemistry Set
# ========================================
# The mechanism loaded is the GRI 3.0 mechanism for methane combustion.
# The mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"/reaction/data"*.
#
# .. note::
#   The transport data *must* be included and preprocessed because the transport processes,
#   *convection and diffusion*, are important to sustain the flame structure.
#

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# including the full file path is recommended
chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

##############################
# Preprocess the Chemistry Set
# ============================

# preprocess the mechanism files
iError = MyGasMech.preprocess()
if iError != 0:
    print("Error: Failed to preprocess the mechanism!")
    print(f"       Error code = {iError}")
    exit()

#########################################################################
# Set up the opposed-flow inlet streams for the dual flame simulation
# =======================================================================
# The opposed-flow flame has two opposing inlet streams separated by a small gap.
# Conventionally the "fuel" inlet is located at x = 0, and the "oxidizer" inlet is at the
# opposite end of the separation (in reality the nozzles are arranged vertically).
# A strained non-premixed flame (or diffusion flame) can be established in the gap between
# the inlets by tuning the velocities and the mixture properties of the two inlet streams.
#
# This stream is a mixture with the addition of the
# *inlet flow rate*. You can specify the inlet gas properties the same way you
# set up a mixture. Here the recipe is used to set up the species
# mole fractions of the fuel mixture. For convenience, the pre-defined ``Air()`` in PyChemkin
# is used to set the composition of the ``air`` mixture.
# The inlet velocity is assigned by the ``velocity`` method.

# create the "fuel" stream
fuel = Stream(MyGasMech, label="FUEL")
# set fuel composition: fuel-rich methane-air mixture (equivalence ratio ~ 1.55)
fuel.X = [("CH4", 0.14001807), ("O2", 0.18066847), ("N2", 0.67931346)]
# system pressure [dynes/cm2]
fuel.pressure = ck.Patm
# fuel temperature [K]
fuel.temperature = 300.0
# fuel inlet velocity [cm/sec]
fuel.velocity = 16.0

# create the oxidizer mixture: air
air = Stream(MyGasMech, label="OXID")
air.X = ck.Air.X()
# oxidizer pressure (same as the fuel stream)
air.pressure = fuel.pressure
# oxidizer temperature (same as the fuel temperature)
air.temperature = fuel.temperature
# oxidizer inlet velocity [cm/sec]
air.velocity = 16.0

##################################################
# Instantiate the opposed-flow flame model
# ================================================
# Set up the *opposed-flow flame* model by using the stream
# representing the "fuel" mixture at the origin. The "oxidizer" inlet
# is added to the ``OpposedFlame`` object later by using the ``set_oxidizer_inlet`` method.
# There are many options and parameters related to the treatment of the species boundary
# condition, the transport properties. All the available options and parameters are described
# in the *Chemkin Input* manual.
#
# .. note::
#   The parameter used to instantiate the ``OpposedFlame`` is
#   the stream representing the "fuel" inlet at x = 0.
#

dual_flame = Flame(fuel, label="opposed_flame")

#################################
# Configure the opposed flame
# ===============================
# Use the ``set_oxidizer_inlet`` method to specify the "oxidizer"
# stream at the opposite end of the separation gap. The distance
# between the two opposing inlets is defined by the ``end_position`` method.
# The ``end_poistion`` is a required input as it defines the length of the calculation domain.
# Typically, the length of the calculation domain is less than 10 [cm].
#
# The opposed-flow flame model will set up a *flame zone* to improve the convergence
# performance. The center location of this "estimated" *flame zone* can be given by the
# ``set_reaction_zone_center``, and the width of the *flame zone* is defined by the
# ``set_reaction_zone_width``. The temperature and the gas species profiles are assumed to
# vary linearly from the inlets and plateaued inside the *flame zone*. The gas temperature
# value in the *flame zone* is specified by the ``set_max_flame_temperature()`` method
# (or 2200.0 [K] by default). The gas composition in the *flame zone* is set by the
# equilibrium composition of the mixture of the "fuel" and the "oxidizer" streams at the
# given "maximum flame temperature".
#
# Alternatively, a different temperature profile can be set up by using the ``setprofile()`` method
# to replace the default "plateau" profile.
#

# add the "oxidizer" inlet
dual_flame.set_oxidizer_inlet(air)

# define the gap between the two opposing inlets (calculation domain) [cm]
dual_flame.end_position = 1.5
# set up of the "flame zone" to establish the guessed species profiles
# flame zone center location [cm]
dual_flame.set_reaction_zone_center(0.75)
# flame zone width [cm]
dual_flame.set_reaction_zone_width(0.5)
# set the estimated maximum gas temperature [K]
dual_flame.set_max_flame_temperature(2200.0)

###############################################
# Set up initial mesh and grid adaption options
# =============================================
# The opposed-flow flame models provides several methods to set up the initial
# mesh.  Here a uniform mesh of 26 grid points is used at the start of the simulation.
# The flame models would add more grid points to where they are needed as determined by
# the solution quality parameters specified by the ``set_solution_quality()`` method.
#
# .. note::
#   There are two ways to set up the initial mesh for the opposed-flow flame calculations:
#
#   1. ``set_numb_grid_points`` method to create a uniform mesh of the given number of grid points.
#
#   2. ``set_grid_profile`` method to specify the initial grid point profile.
#

# set the initial mesh to 26 uniformly distributed grid points
dual_flame.set_numb_grid_points(26)
# set the maximum total number of grid points allowed in the calculation (optional)
dual_flame.set_max_grid_points(250)
# maximum number of grid points can be added during each grid adaption event (optional)
dual_flame.set_max_adaptive_points(5)
# set the maximum values of the grdient and the curvature of the solution profiles (optional)
dual_flame.set_solution_quality(gradient=0.1, curvature=0.3)

#################################
# Set transport property options
# ===============================
# Ansys Chemkin offers three methods for computing mixture properties:
#
# - **Mixture averaged**
# - **Multi-component**
# - **Constant Lewis number**
#
# When the system pressure is not too low, the mixture averaged method should be adequate.
# The multi-component method, although it is slightly more accurate, makes the simulation time longer
# and is harder to converge. Using the constant Lewis number method implies that all the species
# would have the same transport properties. Include the thermal diffusion effect, when there are large
# amount of light species (molecular weight < 5.0).
#

# use the mixture averaged formulism to evaluate the mixture transport properties
dual_flame.use_mixture_averaged_transport()
# do NOT include the thermal diffusion effect
dual_flame.use_thermal_diffusion(mode=False)

#########################################
# Set species composition boundary option
# =======================================
# There two types of boundary condition treatments for the species composition available
# from the premixed flame models: ``comp`` and ``flux``. You can find the descriptions of
# these two treatments in the *Chemkin Input* manual.

# specific the species composition boundary treatment ('comp' or 'flux')
# use 'flux' to keep the net species mass fluxes the same as given by the "inlet streams".
dual_flame.set_species_boundary_types(mode="flux")

############################
# Set solver parameters
# ==========================
# The steady state solver parameters for the opposed-flow flame model are optional because
# all the solver parameters have their own default values. Change the solver parameters when the
# flame simulation does not converge with the default settings.
#

# reset the tolerances in the steady-state solver (optional)
dual_flame.steady_state_tolerances = (1.0e-9, 1.0e-5)
dual_flame.time_stepping_tolerances = (1.0e-6, 1.0e-4)

##########################################
# Run the opposed-flow flame simulation
# ========================================
# Use the ``run()`` method to run the opposed-flow flame model.
# After the calculation concludes successfully, use the ``process_solution()`` method
# to postprocess the solutions. You can create other property profiles by looping through the
# solution streams by using proper ``Mixture`` methods.
#

# set the start wall time
start_time = time.time()

status = dual_flame.run()
if status != 0:
    print(
        Color.RED
        + "Failed to get a converged solution of the opposed flow flame!"
        + Color.END
    )
    exit()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"Total simulation duration: {runtime} [sec].")
print()

#############################################
# Postprocess the opposed-flow flame results
# ===========================================
# The post-processing step will parse the solution and package the solution values at each
# time point into a streams. There are two ways to access the solution profiles:
#
#   1. the raw solution profiles (value as a function of time) are available for "distance",
#   "temperature", and species "mass fractions";
#
#   2. the streams that permit the use of all property and rate utilities to extract
#   information such as viscosity, density, and mole fractions.
#
# You can use the ``get_solution_variable_profile()`` method to get the raw solution profiles. You
# solution streams are accessed using either the ``get_solution_stream_at_grid()`` method for the
# solution stream at the given grid point or the ``get_solution_stream()`` method for the
# solution stream at the given location. (In this case, the stream is constructed by interpolation.)
#
# .. note::
#   - Use the ``get_solution_size()`` to get the number of grid pints in the solution profiles before
#     creating the arrays.
#   - The ``mass_flowrate`` from the solution streams is actually the *mass flux* [g/cm\ :sup:`2`\ -sec].
#     It can used to derive the velocity at the corresponding location by dividing it by the local gas mixture
#     density [g/cm\ :sup:`3`\ ].
#   - In addition to the usual raw solution profiles such as "distance", "temperature" and species mass fractions,
#     the opposed-flow flame model provides solution profiles for the velocities ("axial_velocity" and
#     "radial_velocity_gradient") and the mixture fraction ("mixture_fraction").
#   - The existence of two separated flames can be by shown by plotting the heat release rate
#     or the concentrations of radical species found typically in hydrocarbon flames.
#

# postprocess the solutions
dual_flame.process_solution()

# get the number of solution grid points
solutionpoints = dual_flame.get_solution_size()
print(f"Number of solution points = {solutionpoints}.")
# get the grid profile
mesh = dual_flame.get_solution_variable_profile("distance")
# get the temperature profile
tempprofile = dual_flame.get_solution_variable_profile("temperature")
# get the axial velocity profile
velprofile = dual_flame.get_solution_variable_profile("axial_velocity")
# get the mixture fraction profile
mfprofile = dual_flame.get_solution_variable_profile("mixture_fraction")
# get NO2 mass fraction profile
NO2profile = dual_flame.get_solution_variable_profile("NO2")

###############################################
# Plot the opposed-flow flame solution profiles
# =============================================
# Plot the solution profiles of the opposed-flow flame.
#
# .. note ::
#   You can get profiles of the thermodynamic and the transport properties
#   by applying ``Mixture`` utility methods to the solution streams.
#
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(mesh, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(mesh, velprofile, "b-")
plt.ylabel("Axial Velocity [cm/sec]")
plt.subplot(223)
plt.plot(mesh, NO2profile, "g-")
plt.xlabel("Distance [cm]")
plt.ylabel("NO2 Mass Fraction")
plt.subplot(224)
plt.plot(mesh, mfprofile, "m-")
plt.xlabel("Distance [cm]")
plt.ylabel("Mixture Fraction [-]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_opposed_flow_flame.png", bbox_inches="tight")

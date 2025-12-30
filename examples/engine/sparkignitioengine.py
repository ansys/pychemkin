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

r""".. _ref_sparkignition_engine:

================================
Simulate a spark ignition engine
================================

Ansys chemkin offers some idealized internal combustion (IC) engine models commonly
used for fuel combustion and engine performance research.
The Chemkin IC engine model is a specialized transient 0-D closed gas-phase reactor
that mainly performs combustion simulation between the intake valve closing (IVC)
and the exhaust valve opening (EVO), that is, when the engine cylinder
resembles a closed chamber. The cylinder volume is derived from the piston motion
as a function of the engine crank angle (CA) and engine parameters such as
engine speed (RPM) and stroke. The energy equation is always solved.
There are several wall heat transfer models specifically designed
for engine simulations.

.. note ::
    For additional information on Chemkin IC engine models, use the
    ``ansys.chemkin.core.manuals()`` method to view the online **Theory** manual.

The Chemkin spark ignition (SI) engine model offers a simple way to simulate
the chemical kinetics taking place in the spark ignition engine.
The Chemkin SI engine model does not predict the fuel mass burning rate profile.
On the contrary, it requires the burning rate profile as input in the form of
the Wiebe function parameters, the burn profile anchor points,
or normalized profile data. The main uses of the Chemkin SI engine model are
conducting parameter studies of the burning rate profile on engine performance,
emissions, and the onset of engine knock due to end gas autoignition.

This example shows how to set up and run the simplest Chemkin IC engine model,
the SI engine model. In addition to the basic engine parameters,
many engine model-specific features such as the exhaust gas recirculation and
the wall heat transfer can be included in the engine simulation.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_spark_ignition_engine.png'

###############################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# chemkin spark ignition (SI) engine model (transient)
from ansys.chemkin.core.engines.SI import SIengine
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory:" + current_dir)
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
# For PRF, the encrypted 14-component gasoline mechanism,
# ``gasoline_14comp_WBencrypted.inp``, is used. The chemistry set
# is named ``gasoline``.
#
# .. note::
#   Because this gasoline mechanism does not come with any transport data,
#   you do not need to provide a transport data file.
#

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on the gasoline 14 components mechanism
MyGasMech = ck.Chemistry(label="Gasoline")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "gasoline_14comp_WBencrypt.inp")

#######################################
# Preprocess the gasoline chemistry set
# =====================================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()
print("Mechanism information:")
print(f"Number of gas species = {MyGasMech.kk:d}.")
print(f"Number of gas reactions = {MyGasMech.ii_gas:d}.")

################################################
# Set up the stoichiometric gasoline-air mixture
# ==============================================
# You must set up the stoichiometric gasoline-air mixture for the subsequent
# SI engine calculations. Here the ``x_by_equivalence_ratio()``method is used.
# You create the ``fuel`` and the ``air`` mixtures first. You then define the
# *complete combustion product species* and provide the *additives* composition
# if applicable. Finally, you simply set ``equivalenceratio=1`` to create
# the stoichiometric gasoline-air mixture.
#
# For PRF 90 gasoline, the recipe is ``[("ic8h18", 0.9), ("nc7h16", 0.1)]``.

# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.x = [("ic8h18", 0.9), ("nc7h16", 0.1)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = ck.P_ATM
fuelmixture.temperature = 353.0

# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.x = [("o2", 0.21), ("n2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = ck.P_ATM
air.temperature = 353.0

# products from the complete combustion of the fuel mixture and air
products = ["co2", "h2o", "n2"]
# species mole fractions of added/inert mixture.
# You can also create an additives mixture here.
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros

# create the unburned fuel-air mixture
fresh = ck.Mixture(MyGasMech)

# mean equivalence ratio
equiv = 1.0
ierror = fresh.x_by_equivalence_ratio(
    MyGasMech, fuelmixture.x, air.x, add_frac, products, equivalenceratio=equiv
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
    print("Error: Failed to create the fuel-oxidizer mixture.")
    exit()

##############################
# List the mixture composition
# ============================
# List the composition of the premixed mixture for verification.
fresh.list_composition(mode="mole")

##########################################################
# Specify pressure and temperature of the fuel-air mixture
# ========================================================
# Since you are going to use ``fresh`` fuel-air mixture to instantiate
# the engine object later, setting the mixture pressure and temperature
# is equivalent to setting the initial temperature and pressure of the
# engine cylinder.
fresh.temperature = fuelmixture.temperature
fresh.pressure = fuelmixture.pressure

###########################################
# Add EGR to the fresh fuel-air mixture
# =========================================
# Many engines have the configuration for exhaust gas recirculation (EGR). Chemkin
# engine models let you add the EGR mixture to the fresh fuel-air mixture entering
# the cylinder. If the engine you are modeling has EGR, you should have the EGR ratio,
# which is generally the volume ratio between the EGR mixture and
# the fresh fuel-air ratio. However, because you know nothing about the composition
# of the exhaust gas, you cannot simply combine these two mixtures. In this case,
# you use the ``get_egr_mole_fraction()`` method to estimate the major components of
# the exhaust gas from the combustion of the fresh fuel-air mixture.
# The ``threshold=1.0e-8`` parameter tells the method to ignore any species with
# a mole fraction below the threshold value. Once you have
# the EGR mixture composition, use the ``x_by_equivalence_ratio()`` method
# a second time to re-create the fuel-air mixture ``fresh`` with the original
# ``fuelmixture`` and ``air`` mixtures, along with the EGR composition that
# you just got as the *additives*.
egr_ratio = 0.3
# compute the EGR stream composition in mole fractions
add_frac = fresh.get_egr_mole_fraction(egr_ratio, threshold=1.0e-8)
# recreate the initial mixture with EGR
ierror = fresh.x_by_equivalence_ratio(
    MyGasMech,
    fuelmixture.x,
    air.x,
    add_frac,
    products,
    equivalenceratio=equiv,
    threshold=1.0e-8,
)
# list the composition of the fuel+air+EGR mixture for verification
fresh.list_composition(mode="mole", bound=1.0e-8)

##############################
# Set up the SI engine reactor
# ============================
# Use the ``SIengine()`` method to create an SI engine named ``MyEngine``
# and make the new ``fresh`` mixture as the initial in-cylinder gas mixture
# at the IVC. The Chemkin SI engine model consists of two*zones: the unburned
# zone and the burned zone. At the IVC, the unburned zone is filled with the
# ``fresh`` mixture and occupies the entire engine cylinder. On the other hand,
# the burned zone is empty and has zero volume/mass.
MyEngine = SIengine(reactor_condition=fresh)
# show initial gas composition inside the reactor
MyEngine.list_composition(mode="mole", bound=1.0e-8)

################################
# Set up basic engine parameters
# ==============================
# Set the required engine parameters as shown in the following code. These
# engine parameters are used to describe the cylinder volume during the
# simulation. The ``starting_ca`` argument should be the crank angle corresponding
# to the cylinder IVC. The ``ending_ca`` argument is typically the EVC crank angle.

# cylinder bore diameter [cm]
MyEngine.bore = 8.5
# engine stroke [cm]
MyEngine.stroke = 10.82
# connecting rod length [cm]
MyEngine.connecting_rod_length = 17.853
# compression ratio [-]
MyEngine.compression_ratio = 12
# engine speed [RPM]
MyEngine.rpm = 600

# set other parameters
# simulation start CA [degree]
MyEngine.starting_ca = -120.2
# simulation end CA [degree]
MyEngine.ending_ca = 139.8
# list the engine parameters
MyEngine.list_engine_parameters()
print(f"Engine displacement volume = {MyEngine.get_displacement_volume()} [cm3].")
print(f"Engine clearance volume = = {MyEngine.get_clearance_volume()} [cm3].")

##########################################
# Set up the SI engine specific parameters
# ========================================

#########################################
# Set up the mass burned fraction profile
# =======================================
# The burning rate profile is required because the SI engine model uses it to
# determine the mass transfer rate from the unburned zone to the burned zone
# during the simulation. You can use one of three methods to specify
# the mass burning fraction profile for the SI engine simulation:
#
# - Wiebe function
#
#   You must provide two sets of parameters:
#
#   - ``set_burn_timing``: Start of combustion crank angle (CA) and burn duration.
#   - ``wiebe_parameters``: Wiebe function parameters.
#
#  - Burn profile anchor points
#
#    You must provide one set of parameters: ``set_burn_anchor_points``.
#
#  - Burned mass fraction profile data
#
#     You must provide two sets of parameters:
#
#     - ``set_burn_timing``: Start of combustion CA and burn duration.
#     - ``set_mass_burned_profile``: Normalized burn rate profile.
#

# start of combustion CA
MyEngine.set_burn_timing(soc=-14.5, duration=45.6)
MyEngine.wiebe_parameters(n=4.0, b=7.0)
# set minimum zonal mass [g]
MyEngine.set_minimum_zone_mass(minmass=1.0e-5)

########################################
# Set up engine wall heat transfer model
# ======================================
# By default, the engine cylinder is adiabatic. You must set up a
# wall heat transfer model to include the heat loss effects in your
# engine simulation. Chemkin supports three widely used engine wall
# heat transfer models. The models and their parameters follow:
#
# - ``dimensionless``: [<a> <b> <c> <Twall>]
# - ``dimensional``: [<a> <b> <c> <Twall>]
# - ``hohenburg``: [<a> <b> <c> <d> <e> <Twall>]
#
# There is also the in-cylinder gas velocity correlation
# (the Woschni correlation) that is associated with the engine
# wall heat transfer models. Here are the parameters of the Woschni correlation:
#
# ``[<C11> <C12> <C2> <swirl ratio>]``
#
# You can also specify the surface areas of the piston head and the cylinder head
# for more precision heat transfer wall area. By default, both the piston head and
# the cylinder head surfaces are flat.
heattransferparameters = [0.1, 0.8, 0.0]
# set cylinder wall temperature [K]
Twall = 434.0
MyEngine.set_wall_heat_transfer("dimensionless", heattransferparameters, Twall)
# in-cylinder gas velocity correlation parameter (Woschni)
# [<C11> <C12> <C2> <swirl ratio>]
GVparameters = [2.28, 0.318, 0.324, 0.0]
MyEngine.set_gas_velocity_correlation(GVparameters)
# set piston head top surface area [cm2]
MyEngine.set_piston_head_area(area=56.75)
# set cylinder clearance surface area [cm2]
MyEngine.set_cylinder_head_area(area=56.75)

####################
# Set output options
# ==================
# You can turn on the adaptive solution saving to resolve the steep variations
# in the solution profile. Here additional solution data point are saved for
# every 20 solver internal steps. Since ignition in SI engine is controlled
# by the spark timing, that is, the start of combustion CA of
# the chemkin SI engine model, the ignition delay time is not reported.
#
# .. note::
#   By default, time/crank angle intervals for both print and save solution
#   are 1/100 of the simulation duration, which in this case is
#   :math:`dCA=(EVO-IVC)/100=2.58`\ . You can make the model report more frequently
#   by using the ``ca_step_for_saving_solution()`` or
#   ``ca_step_for_printing_solution()`` method to set different interval values
#   in the CA.
#

# set the number of crank angles between saving solution
MyEngine.ca_step_for_saving_solution = 0.5
# set the number of crank angles between printing solution
MyEngine.ca_step_for_printing_solution = 10.0
# turn ON adaptive solution saving
MyEngine.adaptive_solution_saving(mode=True, steps=20)

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances.

# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyEngine.tolerances = (1.0e-15, 1.0e-6)
# get solver parameters
atol, rtol = MyEngine.tolerances
print(f"Default absolute tolerance = {atol}.")
print(f"Default relative tolerance = {rtol}.")
# turn on the force non-negative solutions option in the solver
MyEngine.force_nonnegative = True
# show solver option
# show the number of crank angles between printing solution
print(
    f"Crank angles between solution printing: {MyEngine.ca_step_for_printing_solution}"
)
# show other transient solver setup
print(f"Forced non-negative solution values: {MyEngine.force_nonnegative}")

#########################################
# Display the added parameters (keywords)
# =======================================
# Use the ``showkeywordinputlines()`` method to verify the preceding parameters
# are correctly assigned to the engine model.
#
MyEngine.showkeywordinputlines()
# set the start wall time
start_time = time.time()

####################
# Run the simulation
# ==================
# Use the ``run()`` method to start the SI engine simulation. The simulation
# might take more than 30 seconds because the gasoline mechanism contains a lot
# of species and reactions.
runstatus = MyEngine.run()
# compute the total runtime
runtime = time.time() - start_time
# check run status
if runstatus != 0:
    # Run failed.
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
print(f"Total simulation duration: {runtime} [sec]")

##########################
# Postprocess the solution
# ========================
# The postprocessing step parses the solution and packages the solution values at each
# time point into a mixture. There are two ways to access the solution profiles:
#
# - The raw solution profiles (value as a function of time) are available for time,
#   temperature, pressure, volume, and species mass fractions.
#
# - The mixtures permit the use of all property and rate utilities to extract
#   information such as viscosity, density, and mole fractions.
#
# You can use the ``get_solution_variable_profile()`` method to get
# the raw solution profiles. You can get solution mixtures using either the
# ``get_solution_mixture_at_index()`` method for the solution mixture at
# a given time point or the ``get_solution_mixture()`` method for the solution mixture
# at a given time. (In this case, the mixture is constructed by interpolation.)
#
# .. note ::
#
#   - For engine models, use the ``process_engine_solution()`` method to postprocess
#     the solutions.
#   - Use the ``getnumbersolutionpoints()`` method to get the size of
#     the solution profiles before creating the arrays.
#   - Use the ``get_ca()`` method to convert the time values reported in the solution
#     to crank angles.
#
#

####################################################
# Postprocess the solution profiles in selected zone
# ==================================================
# The solution of the SI engine model contains the results of
# the *unburned zones* and *burned zone*, along with the cylinder averaged results.
# That is, there are three solution records. The unburned zone is designated
# as ``zone 1``, and the burned zone is designated as ``zone 2``. To process
# the result of the burned zone, you use ``zone_id=2`` when you call
# the engine postprocessor with the ``process_engine_solution()`` method.
# If you omit the ``zone_id`` parameter, the cylinder averaged results are
# postprocessed by default.
#
# .. note ::
#   Because the ``process_engine_solution`` method can process only one set
#   of result at a time (one zonal result or the cylinder averaged result),
#   you must postprocess the zones one by one to obtain all solution data
#   of the multi-zone simulation.
#
unburnedzone = 1
burnedzone = 2
zonestrings = ["Unburned Zone", "Burned Zone"]
thiszone = burnedzone
MyEngine.process_engine_solution(zone_id=thiszone)
plottitle = zonestrings[thiszone - 1] + " Solution"
# get the number of solution time points
solutionpoints = MyEngine.getnumbersolutionpoints()
print(f"Number of solution points = {solutionpoints}.")
# get the time profile
timeprofile = MyEngine.get_solution_variable_profile("time")
# convert time to crank angle
ca_profile = np.zeros_like(timeprofile, dtype=np.double)
count = 0
for t in timeprofile:
    ca_profile[count] = MyEngine.get_ca(timeprofile[count])
    count += 1
# get the cylinder pressure profile
presprofile = MyEngine.get_solution_variable_profile("pressure")
presprofile *= 1.0e-6
# get zonal temperature profile [K]
tempprofile = MyEngine.get_solution_variable_profile("temperature")
# get the zonal volume profile
volprofile = MyEngine.get_solution_variable_profile("volume")
# create arrays for zonal CO mole fraction
# find CO species index
co_index = MyGasMech.get_specindex("co")
co_profile = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the zonal mixture at the time point
    solutionmixture = MyEngine.get_solution_mixture_at_index(solution_index=i)
    # get zonal CO mole fraction
    co_profile[i] = solutionmixture.x[co_index]

# postprocess cylinder-averged solution
MyEngine.process_average_engine_solution()
# get the cylinder volume profile
cylindervolprofile = MyEngine.get_solution_variable_profile("volume")
# create arrays for cylinder-averaged mixture temperature
cylindertempprofile = MyEngine.get_solution_variable_profile("temperature")

###################################
# Plot the engine solution profiles
# =================================
# Plot the zonal and the cylinder averaged profiles from
# the SI engine simulation.
#
# .. note ::
#   You can get profiles of the thermodynamic and the transport properties
#   by applying ``Mixture`` utility methods to the solution mixtures.
#
plt.subplots(2, 2, sharex="col", figsize=(13, 6.5))
plt.suptitle(plottitle, fontsize=16)
plt.subplot(221)
plt.plot(ca_profile, presprofile, "r-")
plt.ylabel("Pressure [bar]")
plt.subplot(222)
plt.plot(ca_profile, volprofile, "b-")
plt.plot(ca_profile, cylindervolprofile, "b--")
plt.ylabel("Volume [cm3]")
plt.legend(["Zone", "Cylinder"], loc="lower right")
plt.subplot(223)
plt.plot(ca_profile, tempprofile, "g-")
plt.plot(ca_profile, cylindertempprofile, "g--")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("Mixture Temperature [K]")
plt.legend(["Zone", "Averaged"], loc="upper left")
plt.subplot(224)
plt.plot(ca_profile, co_profile, "m-")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("CO Mole Fraction")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_spark_ignition_engine.png", bbox_inches="tight")

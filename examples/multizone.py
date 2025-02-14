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
.. _ref_multizone_engine:

==================================
Simulate a multi-zone HCCI engine
==================================

**Ansys chemkin** offers some idealized *internal combustion engine*(IC engine) models commonly used
for fuel combustion and engine performance researches. The *Chemkin* ``IC engine model`` is a specialized
transient 0-D *closed* gas-phase reactor that mainly perform combustion simulation between the
intake valve closing (IVC) and the exhaust valve opening (EVO), that is, when the engine cylinder resembles
a closed chamber. The cylinder volume is derived from the piston motion as a function of the engine
crank angle (CA) and engine parameters such as engine speed (RPM) and stroke. The energy equation
is always solved and there are several wall heat transfer models specifically designed for engine
simulations.

.. note ::
    For additional information about the *Chemkin* IC engine models, please use the
    ``ansys.chemkin.manuals()`` method to check out the online **Theory** manual.

The multi-zone homogeneous charged compression ignition (HCCI) model is mainly intended to address
the temperature variation inside the cylinder caused by the wall heat transfer and the imperfect
mixing of the incylinder gas mixture. In addition, the multi-zone HCCI engine model allows the
introduction of non-uniform temperature and/or the equivalence ratio distribution of the gas mixture
at the IVC.

This tutorial models will describe the basic procedures to set up and run the chemkin multi-zone HCCI engine
model. Many engine model specific features such as the basic engine parameters, the "exhaust gas recirculation",
and the "wall heat transfer" will be applied to the engine simulation. You will also learn the steps to create
initial distributions of zone size, temperature, and composition.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_multizone_HCCI_engine.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color

# chemkin homogeneous charge compression ignition (HCCI) engine model (transient)
from ansys.chemkin.engines.HCCI import HCCIengine
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
# The mechanism is the GRI 3.0 mechanism for methane combustion.
# The mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"/reaction/data"*.

# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
MyGasMech.tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")

###################################
# Pre-process the ``Chemistry Set``
# =================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

#############################
# Set up the fuel-air mixture
# ============================
# You need to set up the fuel-air mixture inside the engine cylinder
# right after the intake valve is closed. Here the ``X_by_Equivalence_Ratio``
# method is used. You create the ``fuelmixture`` and the ``air`` mixtures first.
# Then define the *complete combustion product species* and provide the
# *additives* composition if there is any. And finally you can simply set
# the value of ``equivalenceratio`` to create the fuel-air mixture. In this case,
# the fuel mixture consists of methane, ethane, and propane as the simulated "natural gas".
# Typically the HCCI engines run on lean fuel-air mixtures so the equivalence ratio is
# set to 0.8.

# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.X = [("CH4", 0.9), ("C3H8", 0.05), ("C2H6", 0.05)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = 1.5 * ck.Patm
fuelmixture.temperature = 400.0
# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.X = [("O2", 0.21), ("N2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = 1.5 * ck.Patm
air.temperature = 400.0
# create the unburned fuel-air mixture
fresh = ck.Mixture(MyGasMech)
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros
# mean equivalence ratio
equiv = 0.8
iError = fresh.X_by_Equivalence_Ratio(
    MyGasMech, fuelmixture.X, air.X, add_frac, products, equivalenceratio=equiv
)
# check fuel-oxidizer mixture creation status
if iError != 0:
    print("Error: failed to create the Fuel-Oxidizer mixture!")
    exit()

# list the composition of the unburned fuel-air mixture
fresh.list_composition(mode="mole")

##########################################################
# Specify pressure and temperature of the fuel-air mixture
# ========================================================
# Since you are going to use ``fresh`` to instantiate the engine object later,
# setting the mixture pressure and temperature is equivalent to setting
# the initial temperature and pressure of the engine cylinder.
fresh.temperature = 447.0
fresh.pressure = 1.065 * ck.Patm

###########################################
# Add EGR to the ``fresh`` fuel-air mixture
# =========================================
# Many engines have the configuration for exhaust gas recirculation (EGR). Chemkin
# engine models allow you to add the EGR mixture to the fresh fuel-air mixture entered
# the cylinder. If the engine you are modeling has EGR, you should have the EGR ratio which
# in most cases, is the volume ratio between the EGR mixture and the fresh fuel-air ratio.
# However, you know nothing about the composition of the exhaust gas so you cannot simply "mix"
# these two mixtures. In this case, you can use the ``get_EGR_mole_fraction`` method to estimate
# the major components of the exhaust gas from the combustion of the fresh fuel-air mixture. The
# parameter ``threshold=1.0e-8`` tells the method to ignore any species with mole fraction below
# the threshold value. Once you have the EGR mixture composition, use the ``X_by_Equivalence_Ratio``
# method the second time to **recreate** the fuel-air mixture ``fresh`` with the original
# ``fuelmixture`` and ``air`` mixtures **AND** with the EGR composition you just got as the
# *"additives"*.
EGRratio = 0.3
# compute the EGR stream composition in mole fractions
add_frac = fresh.get_EGR_mole_fraction(EGRratio, threshold=1.0e-8)
# recreate the initial mixture with EGR
iError = fresh.X_by_Equivalence_Ratio(
    MyGasMech,
    fuelmixture.X,
    air.X,
    add_frac,
    products,
    equivalenceratio=equiv,
    threshold=1.0e-8,
)

# list the composition of the fuel+air+EGR mixture for verification
fresh.list_composition(mode="mole", bound=1.0e-8)

################################
# Set up the HCCI engine reactor
# ==============================
# Create a smulti-zone HCCI engine object ``MyMZEngine`` using the
# ``HCCIengine`` method and make the *new* ``fresh`` mixture as the initial
# incylinder gas mixture at IVC. Set the parameter ``nzones`` to the number
# of zones in your *multi-zone* HCCI engine model.

# create a 5-zones HCCI engine object
numbzones = 5
MyMZEngine = HCCIengine(reactor_condition=fresh, nzones=numbzones)
# show initial gas composition inside the reactor
MyMZEngine.list_composition(mode="mole", bound=1.0e-8)

#####################################
# Set up basic engine parameters
# ===================================
# Set the required engine parameters as listed in the code below. These
# engine parameters are used to describe the cylinder volume during the
# simulation. The ``starting_CA`` should be the crank angle corresponding
# to the cylinder IVC and the ``ending_CA`` is typically the EVC crank angle.

# cylinder bore diameter [cm]
MyMZEngine.bore = 12.065
# engine stroke [cm]
MyMZEngine.stroke = 14.005
# connecting rod length [cm]
MyMZEngine.connecting_rod_length = 26.0093
# compression ratio [-]
MyMZEngine.compression_ratio = 16.5
# engine speed [RPM]
MyMZEngine.RPM = 1000

# set other parameters
# simulation start CA [degree]
MyMZEngine.starting_CA = -142.0
# simulation end CA [degree]
MyMZEngine.ending_CA = 116.0

# list the engine parameters
MyMZEngine.list_engine_parameters()
print(f"engine displacement volume {MyMZEngine.get_displacement_volume()} [cm3]")
print(f"engine clearance volume {MyMZEngine.get_clearance_volume()} [cm3]")
print(f"number of zone(s) = {MyMZEngine.get_number_of_zones()}")

########################################
# Set up engine wall heat transfer model
# ======================================
# By default, the engine cylinder is adiabatic. You must set up a
# wall heat transfer model to include the heat loss effects in your
# engine simulation. Chemkin support three widely used engine wall
# heat transfer models, and the models and their parameters are
#
#   |     "dimensionless": [<a> <b> <c> <Twall>]
#   |     "dimensional": [<a> <b> <c> <Twall>]
#   |     "hohenburg": [<a> <b> <c> <d> <e> <Twall>]
#
# There is also the incylinder gas velocity correlation
# (the Woschni correlation) that is associated with the engine
# wall heat transfer models. The parameters of the Woschni correlation
# are
#
#   |     [<C11> <C12> <C2> <swirl ratio>]
#
# You can also specify the surface areas of the piston head and the cylinder head
# for more precision heat transfer wall area. By default, both the piston head and
# the cylinder head surfaces are "flat".

heattransferparameters = [0.035, 0.71, 0.0]
# set cylinder wall temperature [K]
Twall = 400.0
MyMZEngine.set_wall_heat_transfer("dimensionless", heattransferparameters, Twall)
# incylinder gas velocity correlation parameter (Woschni)
# [<C11> <C12> <C2> <swirl ratio>]
GVparameters = [2.28, 0.308, 3.24, 0.0]
MyMZEngine.set_gas_velocity_correlation(GVparameters)
# set piston head top surface area [cm2]
MyMZEngine.set_piston_head_area(area=124.75)
# set cylinder clearance surface area [cm2]
MyMZEngine.set_cylinder_head_area(area=123.5)

######################
# Set zonal properties
# ====================
# By default, all zones in the multi-zone HCCI engine model will have the same
# properties. You can artificially stratify the temperature and/or the equivalence
# ratio distribution in the cylinder at the IVC by utilizing the ``set_zonal``
# methods of the ``HCCI`` object.

# zonal temperatures [K]
ztemperature = [447.5, 447.5, 447, 447, 447]
MyMZEngine.set_zonal_temperature(zonetemp=ztemperature)
# zonal volume fractions
zvolumefrac = [0.3, 0.25, 0.2, 0.2, 0.05]
MyMZEngine.set_zonal_volume_fraction(zonevol=zvolumefrac)
# wall heat transfer area fractions
zHTarea = [0.0, 0.15, 0.2, 0.25, 0.4]
MyMZEngine.set_zonal_heat_transfer_area_fraction(zonearea=zHTarea)
# zonal equivalence ratios
zphi = [equiv, equiv, equiv, equiv, equiv]
MyMZEngine.set_zonal_equivalence_ratio(zonephi=zphi)
# zonal EGR ratios
zEGRR = [0.3, 0.3, 0.3, 0.35, 0.35]
MyMZEngine.set_zonal_EGR_ratio(zoneegr=zEGRR)
# set fuel "molar" composition
MyMZEngine.define_fuel_composition([("CH4", 0.9), ("C3H8", 0.05), ("C2H6", 0.05)])
# set oxidizer "molar' composition
MyMZEngine.define_oxid_composition([("O2", 0.21), ("N2", 0.79)])
# set products
MyMZEngine.define_product_composition(["CO2", "H2O", "N2"])
# set EGR composition in mole fractions
zadd = [add_frac, add_frac, add_frac, add_frac, add_frac]
MyMZEngine.define_additive_fractions(addfrac=zadd)

####################
# Set output options
# ==================
# You can turn on the *adaptive solution saving* to resolve the steep variations in the solution
# profile. Here additional solution data point will be saved for every **20** solver internal steps.
# The ``set_ignition_delay`` method must be included for the engine model to
# report the *ignition delay crank angle* after the simulation is done. If ``method="T_inflection"`` is
# set, the reactor model will treat the *inflection points* in the predicted gas temperature profile
# as the indication of an auto-ignition. You can choose a different auto-ignition definition.
#
# .. note::
#   Type ``ansys.chemkin.show_ignition_definitions()`` to get the list of all available ignition
#   delay time definitions in Chemkin.
#
# .. note ::
#   The multi-zone HCCI engine model will report the ignition delay time for each zone as well as
#   the "cylinder averaged" ignition delay time derived from the "cylinder averaged" temperature
#   profile.
#
# .. note::
#   By default, time/crank angle intervals for both print and save solution are **1/100** of the
#   *simulation duration*. In this case :math:`dCA=(EVO-IVC)/100=2.58`\ . You can make the model
#   to report more frequently by using the ``CAstep_for_saving_solution`` or the
#   ``CAstep_for_printing_solution`` method to set different interval values in CA.
#

# set the number of crank angles between saving solution
MyMZEngine.CAstep_for_saving_solution = 0.5
# set the number of crank angles between printing solution
MyMZEngine.CAstep_for_printing_solution = 10.0
# turn ON adaptive solution saving
MyMZEngine.adaptive_solution_saving(mode=True, steps=20)
# specify the ignition definitions
MyMZEngine.set_ignition_delay(method="T_inflection")

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver related methods, for example,
# ``tolerances``.

# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyMZEngine.tolerances = (1.0e-12, 1.0e-10)
# get solver parameters
ATOL, RTOL = MyMZEngine.tolerances
print(f"default absolute tolerance = {ATOL}")
print(f"default relative tolerance = {RTOL}")
# turn on the force non-negative solutions option in the solver
MyMZEngine.force_nonnegative = True
# show solver and output options
# show the number of crank angles between printng solution
print(
    f"crank angles between solution printing: {MyMZEngine.CAstep_for_printing_solution}"
)
# show other transient solver setup
print(f"forced non-negative solution values: {MyMZEngine.force_nonnegative}")

#########################################
# Display the added parameters (keywords)
# =======================================
# You can verify the parameters specified above are correctly assigned to the engine model by
# using the ``showkeywordinputlines`` method.
MyMZEngine.showkeywordinputlines()

####################
# Run the simulation
# ==================
# Use the ``run`` method to start the multi-zone HCCI engine simulation.
runstatus = MyMZEngine.run()
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> RUN COMPLETED <<<", end=Color.END)

######################################################
# Get the ignition delay crank angle from the solution
# ====================================================
# Use the ``get_ignition_delay`` method to extract the "cylinder averaged"
# ignition delay crank angle (CA) after the run is completed.

# get ignition delay "time"
delayCA = MyMZEngine.get_ignition_delay()
print(f"ignition delay CA = {delayCA} [degree]")

###################################
# Get the heat release crank angles
# =================================
# The engine models also report the crank angles when the accumulated heat
# release reaches 10%, 50%, and 90% of the total heat release. Use the
# ``get_engine_heat_release_CAs`` method to extract these heat release
# crank angles (CA).

HR10, HR50, HR90 = MyMZEngine.get_engine_heat_release_CAs()
print("Engine Heat Release Information")
print(f"10% heat release CA = {HR10} [degree]")
print(f"50% heat release CA = {HR50} [degree]")
print(f"90% heat release CA = {HR90} [degree]\n")

###########################
# Post-process the solution
# =========================
# The post-processing step will parse the solution and package the solution values at each
# time point into a ``Mixture`` object. There are two ways to access the solution profiles:
#
#   1. the "raw" solution profiles (value as a function of time) are available for "time",
#   "temperature", "pressure" , "volume", and species "mass fractions";
#
#   2. the ``Mixture`` objects that permit the use of all property and rate utilities to extract
#   information such as viscosity, density, and mole fractions.
#
# The "raw" solution profiles can be obtained by using the ``get_solution_variable_profile`` method. The
# solution ``Mixture`` objects are accessed via either the ``get_solution_mixture_at_index`` for the
# solution mixture at the given *time point* or the ``get_solution_mixture`` for the solution mixture
# at the given *time* (in this case, the "mixture" is constructed by interpolation).
#
# .. note ::
#   For engine models, use the ``process_engine_solution`` to post-process the solutions.
#
# .. note::
#   Use the ``getnumbersolutionpoints`` to get the size of the solution profiles before creating the
#   arrays.
#
# .. note ::
#   Use the ``get_CA`` method to convert the time values reported in the solution to crank angles.
#

######################################################
# Post-process the solution profiles in selected zone
# ====================================================
# The solution of the multi-zone HCCI engine model contains the results of
# the *individual zones* **plus** the *"cylinder averaged"* results. That is,
# if there are n zones in the multi-zone engine model, there will be (n+1) solution
# records: n zonal results and the cylinder averaged results. To process the result
# of the zone number :math:`j`\ , :math:`(1 \leq j \leq n)`\ , set the parameter
# value of ``zoneID`` to :math:`j` when you call the engine post-processor
# ``process_engine_solution``. Otherwise, the cylinder averaged results will be
# post-processed by default, that is, when the ``zoneID`` parameter is omitted.
#
# .. note ::
#   The ``process_engine_solution`` method can process only one set of result at
#   a time (one zonal result or the cylinder averaged result) so you need to
#   post-process the zones one by one to obtain all solution data of the multi-zone
#   simulation.
#
thiszone = 1
MyMZEngine.process_engine_solution(zoneID=thiszone)
plottitle = "Zone " + str(thiszone) + " Solution"
# get the number of solution time points
solutionpoints = MyMZEngine.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile = MyMZEngine.get_solution_variable_profile("time")
# convert time to crank angle
CAprofile = np.zeros_like(timeprofile, dtype=np.double)
count = 0
for t in timeprofile:
    CAprofile[count] = MyMZEngine.get_CA(timeprofile[count])
    count += 1
# get the cylinder pressure profile
presprofile = MyMZEngine.get_solution_variable_profile("pressure")
presprofile *= 1.0e-6
# get the zonal volume profile
volprofile = MyMZEngine.get_solution_variable_profile("volume")
# create arrays for zonal mixture density and mixture specific heat capacity
denprofile = np.zeros_like(timeprofile, dtype=np.double)
viscprofile = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the zonal mixture at the time point
    solutionmixture = MyMZEngine.get_solution_mixture_at_index(solution_index=i)
    # get zonal gas density [g/cm3]
    denprofile[i] = solutionmixture.RHO
    # get zonal mixture viscosity profile [g/cm-sec] or [Poise]
    viscprofile[i] = solutionmixture.mixture_viscosity() * 1.0e2

# post-process cylinder-averged solution
# do NOT set the zoneID parameter
MyMZEngine.process_average_engine_solution()
# get the cylinder volume profile
cylindervolprofile = MyMZEngine.get_solution_variable_profile("volume")
# create arrays for cylinder-averaged mixture density
cylinderdenprofile = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the zonal mixture at the time point
    solutionmixture = MyMZEngine.get_solution_mixture_at_index(solution_index=i)
    # get zonal gas density [g/cm3]
    cylinderdenprofile[i] = solutionmixture.RHO

###################################
# Plot the engine solution profiles
# =================================
# Plot the zonal and the cylinder averaged profiles from
# the multi-zone HCCI engine simulation.
#
# .. note ::
#   You can get profiles of the thermodynamic and the transport properties
#   by applying ``Mixture`` utility methods to the solution ``Mixtures``.
#
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle(plottitle, fontsize=16)
plt.subplot(221)
plt.plot(CAprofile, presprofile, "r-")
plt.ylabel("Pressure [bar]")
plt.subplot(222)
plt.plot(CAprofile, volprofile, "b-")
plt.plot(CAprofile, cylindervolprofile, "b--")
plt.ylabel("Volume [cm3]")
plt.legend(["Zone", "Cylinder"], loc="upper right")
plt.subplot(223)
plt.plot(CAprofile, denprofile, "g-")
plt.plot(CAprofile, cylinderdenprofile, "g--")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("Mixture Density [g/cm3]")
plt.legend(["Zone", "Averaged"], loc="upper left")
plt.subplot(224)
plt.plot(CAprofile, viscprofile, "m-")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("Mixture Viscosity [cP]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_multizone_HCCI_engine.png", bbox_inches="tight")

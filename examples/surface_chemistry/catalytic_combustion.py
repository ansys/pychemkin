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

r""".. _ref_PFR_catalytic_combustion:

===========================================
Simulating CH4 catalytic Combustion process
===========================================
Catalytic combustion was once an promising technology
of ultra low NOx emission for gas power generation. The catalytic combustion process
converts the gas-phase reactants to the gas-phase products on the catalyst surface
while releasing heat, it consists of three major steps: (1) small fuel species
and oxygen get adsorbed onto the catalyst surface, (2) the absorbed species react on
the catalyst surface, and (3) the products such as carbon dioxide and water get
released from the surface into the gas stream. Catalytic combustion can exhibit
something similar to the gas-phase "ignition" phenomenon called "light-off" when
the net heat release from the surface reactions suddenly spikes. This is where the
catalytic combustion process transitions from being kinetic-limited to being
transport-limited. Hence, this is why it is appropriate to use the PFR model in this
project. By going through the "surface route", the combustion process can
avoid triggering the NOx formation reaction pathways, and the "flame" is anchored
at the "light-off" location on the catalyst surface.

Catalytic combustion applications might utilize a two PFRs in series configuration for
large fuel species. The first PFR is for pure gas-phase reactions called the
pre-oxidizer, its main purpose is to partially break down the large hydrocarbon fuel
species for the downstream catalytic combustor. The second PFR allows both gas-phase
and surface reactions where the broken down fuel species are oxidized on the catalyst
surface to generate heat.

This project attempts to obtain the solution profiles in order to catch the surface
light-off in the downstream catalytic combustor. Two ``Chemistry Set`` are created
in this project. The first ``mech_no_surface`` chemistry set is for the first PFR,
the ``pre_burner``, and it is a modified GRI 3.0 combustion mechanism with the addition
of the platinum element (PT). For the downstream ``cat_combustor``, the chemistry set
``mech_catalytic`` includes the surface mechanism "surf_CH4_Cat_Combust_PT.inp" as well
as the same gas mechanism for the ``mech_no_surface`` chemistry set. The gas-phase
solution at the exit of the ``pre_burner`` can be applied to set up the inlet stream
of the ``cat_combustor`` as long as they use the same gas-phase mechanism.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_PFR_catalytic_combustion.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path
import time

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

from ansys.chemkin.core.inlet import Mixture, Stream

from ansys.chemkin.core.logger import logger

# Chemkin PFR model (steady-state)
from ansys.chemkin.core.flowreactors.PFR import PFREnergyConservation as Pfr
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

###########################
# Create the chemistry sets
# =========================
# This example uses the methane catalytic combustion mechanism on platinum that
# comes with the standard Chemkin example,
# ``reactor_network__two_stage_catalytic_combustor``. The reaction mechanism
# consists of two parts: the gas-phase mechanism, ``chem_GRImech3_PT.inp``,
# describing the gas-phase combustion of methane in air; and the surface mechanism,
# ``surf_CH4_Cat_Combust_PT.inp``, describing the reactions lead to the light-off
# (surface ignition) of the adsorbed methane and the adsorbed oxygen (or surface C
# and O atoms) on the platinum catalyst.
#
# The mechanisms are provided in the *examples\data* folder of the local
# *PyChemkin* installation. You must provide the mechanism file names,
# ``chem_GRImech3_PT.inp`` and  ``surf_CH4_Cat_Combust_PT.inp`` with the correct path
# ``example_mech_data`` to the Chemistry Set before pre-processing. The thermodynamic
# data file from the default Ansys Chemkin installation will be used.
#
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
# find the location of this example py file
try:
    script_dir_obj = Path(__file__).parent.resolve()
except NameError:
    script_dir_obj = Path(current_dir)
script_dir = str(script_dir_obj)
# use the relative path to locate the example mechanism folder
example_mech_data = script_dir_obj / ".." / "data"
print(f"example data = {str(example_mech_data)}")


###########################################
# Set up and run the gas phase only reactor
# =========================================
# Create a chemistry set based on the GRI 3.0 methane combustion mechanism
# with no surface mechanism input.

# gas-phase chemistry only
mech_no_surface = ck.Chemistry(label="GRI3")
# set mechanism input files
# including the full file path is recommended
mech_no_surface.chemfile = str(data_dir / "grimech30_chem.inp")
mech_no_surface.thermfile = str(data_dir / "grimech30_thermo.dat")

#############################################
# Preprocess the gas-phase only chemistry set
# ===========================================
# preprocess the mechanism files.
_ = mech_no_surface.preprocess()

###########################################
# Set up the lean premixed fuel-air stream
# =========================================
# Instantiate a stream feed for the inlet gas mixture.
# The stream is a mixture with the addition of the
# inlet flow rate. You specify inlet gas properties in the same way as you
# set up a mixture.
#
# Use the ``x`` method of the ``Stream`` object to specify the natural gas
# composition with a ``recipe``. In this project, the ``sccm`` method,
# [cm3/min] @ the standard condition, is used set the inlet volumetric flow rate.
# The ``air`` composition is copied from the pre-defined ``Air`` object in Pychemkin.
# To create the premixed fuel-air inlet stream with the given equivalence ratio, the
# ``x_by_equivalence_ratio`` method is applied.
#
# Once the inlet stream ``lean_premix`` is created correctly, it is used to instantiate
# the gas-only PFR ``pre-burner``. Provide the necessary input parameters and run the
# PFR. Use the ``get_last_solution_mixture()`` method to get the solution ``Stream``
# ``preburner_exhaust`` at the exit of the ``pre_burner``. ``preburner_exhaust`` will
# be used to define the inlet stream properties of the down stream catalytic combustor
# ``cat_combustor``.

# set the fuel composition and conditions
natural_gas = Mixture(mech_no_surface)
natural_gas.pressure = 2.5 * ck.P_ATM
natural_gas.temperature = 890.0
natural_gas.x = [
    ("CH4", 0.971),
    ("C2H6", 0.0092),
    ("C3H8", 0.0051),
    ("CO2", 0.0053),
    ("N2", 0.0094),
]

# set the air composition and conditions
air = Mixture(mech_no_surface)
air.pressure = natural_gas.pressure
air.temperature = natural_gas.temperature
air.x = ck.Air.x()

# create the lean premixed natural gas-air mixture for the pre-burner
lean_premix = Stream(mech_no_surface)
lean_premix.pressure = natural_gas.pressure
lean_premix.temperature = natural_gas.temperature
# set stream velocity [cm/sec]
lean_premix.velocity = 60.0
# set equivalence ratio
equiv = 0.3
# products from the complete combustion of the natural gas and the air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture.
# You can also create an additives mixture here.
add_frac = np.zeros(mech_no_surface.kk, dtype=np.double)  # no additives: all zeros
# use the equivalence ratio to set the stream composition
ierror = lean_premix.x_by_equivalence_ratio(
    mech_no_surface, natural_gas.x, air.x, add_frac, products, equivalenceratio=equiv
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
    print("Error: Failed to create the fuel-oxidizer mixture.")
    exit()

# instantiate the pre-burner object with the lean premixed natural gas-air stream
pre_burner = Pfr(lean_premix, label="Pre-Burner")
# pre-burner parameters
# reactor length [cm]
pre_burner.length = 100.0
# reactor diameter [cm]
pre_burner.diameter = 12.0
# solver parameters
pre_burner.force_nonnegative = True

# set the start wall time
start_time = time.time()
# run the pre-burner PFR model
runstatus = pre_burner.run()
# check run status
if runstatus != 0:
    # Run failed.
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# Run succeeded.
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)

# post process the homogeneous stage solution
pre_burner.process_solution()
# get the number of solution time points
n_points = pre_burner.getnumbersolutionpoints()

# the exit/last solution stream from the pre-burner becomes the feed stream of the
# downstream catalytic combustor
preburner_exhaust = pre_burner.get_last_solution_mixture()
# list the composition of the unburned fuel-air mixture
print(f"preburner exhaust temperature = {preburner_exhaust.temperature} [K]")
print(f"preburner exhaust mass flow rate = {pre_burner.mass_flowrate} [g/sec]")
# preburner_exhaust.list_composition(mode="mole")

# save the mass flow rate [g/cm]
preburner_exhaust_mflrt = pre_burner.mass_flowrate

######################################
# Set up and run the catalytic reactor
# ====================================
# Create a new chemistry set based on the GRI 3.0 methane combustion mechanism
# with the addition of the catalytic combustion surface mechanism.

mech_catalytic = ck.Chemistry(label="catalytic")
# set mechanism input files
# including the full file path is recommended
local_data_dir = "CH4_Cat_Combust"
mech_catalytic.chemfile = str(
    example_mech_data / local_data_dir / "chem_GRImech3_PT.inp"
)
mech_catalytic.surffile = str(
    example_mech_data / local_data_dir / "surf_CH4_Cat_Combust_PT.inp"
)
mech_catalytic.thermfile = str(data_dir / "grimech30_thermo.dat")
###################################################
# Preprocess the catalytic combustion chemistry set
# =================================================
# preprocess the mechanism files.
_ = mech_catalytic.preprocess()

#####################################################
# Set up the inlet stream for the catalytic combustor
# ===================================================
# Use the exhaust mixture from the ``pre-burner`` to define the inlet
# stream ``cat_feed`` for the catalytic combustor ``cat_combustor``.
# Even though the ``Chemistry Sets`` for the two PFR have the same
# gas phase mechanism, you should avoid applying or copying the
# ``preburner_exhaust`` object to instantiate ``cat_combustor`` because
# ``preburner_exhaust`` is obtained from the ``pre_burner`` which does
# not have any surface chemistry data. You should assign individual properties
# of the ``cat_feed`` from the ``preburner_exhaust``.

cat_feed = Stream(mech_catalytic)

# use the exhaust mixture (the solution) from the gas-phase only pre-burner
# to set up the feed stream properties of the catalytic combustor
cat_feed.x = preburner_exhaust.x
cat_feed.pressure = preburner_exhaust.pressure
cat_feed.temperature = preburner_exhaust.temperature
# the mass flow rate of the feed stream to the catalytic burner
# equals the mass flow rate out of the gas phase only pre-burner
cat_feed.mass_flowrate = preburner_exhaust_mflrt
# list the composition of the unburned fuel-air mixture
print(f"cat burner feed temperature = {cat_feed.temperature} [K]")
print(f"cat burner feed mass flow rate = {cat_feed.mass_flowrate} [g/sec]")
# cat_feed.list_composition(mode="mole")

####################################################################
# Create the PFR to predict the gas composition of the outlet stream
# ===================================================================
# Instantiate the PFR ``cat_combustor`` as a ``PFREnergyConservation``
# object.

cat_combustor = Pfr(cat_feed, label="cat_combustor")

############################################
# Set up additional reactor model parameters
# ==========================================
# Before you can run the simulation, you must provide reactor parameters,
# solver controls, and output instructions. For PFR, you must provide
# either the reactor length, the reactor diameter (or, separately, the flow area
# and the surface area per reactor length). You can also use profiles to assign
# certain reactor conditions, such as the reactor temperature, the reactor pressure,
# and the surface area.

# catalytic combustor parameters
# reactor length [cm]
cat_combustor.length = 10.0
# reactor cross-sectional flow area [cm2]
cat_combustor.flowarea = 5.0
# reactor chemically active surface area per reactor length [cm2/cm]
# cat_combustor.area = 3000.0
area_prof_x = [0.0, 0.5, 1.0, 25.0]
area_prof_a = [2.0e2, 1.0e3, 1.0e4, 1.0e4]
cat_combustor.set_surfacearea_profile(area_prof_x, area_prof_a)
# specify the preactor ressure profile to simulate the pressure drop across
# the reactor length
# x: [cm]
pres_prof_x = [0.0, 25.0]
# pressure: [dynes/cm2]
pres_prof_p = [2.5 * ck.P_ATM, 2.48 * ck.P_ATM]
cat_combustor.set_pressure_profile(pres_prof_x, pres_prof_p)

#####################################
# Set up surface chemistry parameters
# ===================================
# When the ``Chemistry Set`` of the ``Stream`` includes surface chemistry,
# you need to provide some essential surface chemistry related model parameters.
# For example, you need to specify the activate surface area as well as the
# temperature for each material defined in the surface mechanism. You can also provide
# the surface coverage, site fractions and bulk activities, of the materials at the
# PFR inlet to improve convergence performance. Use the ``set_site_fraction()`` and
# the ``set_bulk_activity()`` method to specify estimated site fractions and bulk
# activities at the reactor entrance (x = 0). By default, all site species have the
# same fraction and all bulk species activity is set to 1.
#
# You can use surface chemistry methods such as ``get_material_names()` and
# ``get_site_species_names()`` to obtain the material names and the surface site
# species symbols, respectively.

# get the number of surface material defined in the surface mechanism
n_material = cat_combustor.get_numb_material
# get all surface material names
cat_material_names = cat_combustor.get_material_names()
# set surface area fraction (by default, materials have the same area)
# site species symbols: list[list[str]]
site_names = cat_combustor.get_site_species_names()
# bulk species symbols: list[list[str]]
bulk_names = cat_combustor.get_bulk_species_names()
# list the material names and the site species and the bulk species symbols
# per material
print(f"number of surface materials = {n_material}")
for m in range(n_material):
    print(f"material names: {cat_material_names[m]}")
    print(f"   list of site species:\n     {site_names[m]}")
    print(f"   list of bulk species:\n     {bulk_names[m]}")

# set up surface coverage of all surface materials
site_recipe = [[("O(S)", 1.0), ("PT(S)", 0.0)]]
bulk_recipe = [[("PT(B)", 1.0)]]
# set the surface condition of the surface material
for i, mname in enumerate(cat_material_names):
    # set area fraction [-] (optional: by default, materials have the same area)
    # set material temperature [K] if it is different from the gas temperature
    # (optional: by default, it is the same as the mixture temperature)
    # cat_combustor.set_material_temperature(mname, 1000.0)
    # set reactor surface coverage of the material
    # site fractions (optional: by default, all sites have the same fraction)
    cat_combustor.set_site_fraction(mname, site_recipe[i])
    # bulk activities (optional: by default, all bulk species activities = 1.0)
    cat_combustor.set_bulk_activity(mname, bulk_recipe[i])
# scale surface reaction rates
cat_combustor.surface_ratemultiplier = 1.0

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances.

# set solver the tolerances
cat_combustor.tolerances = (1.0e-8, 1.0e-6)
# solver non-negative mode
cat_combustor.force_nonnegative = False
# transient solver max time step size [cm]
# cat_combustor.set_solver_max_timestep_size(1.0e-2)
# set adaptive solution saving distance
cat_combustor.adaptive_solution_saving(mode=True, steps=20)

###########################################
# Run the catalytic combustor reactor model
# =========================================
# Once all the necessary input parameters are provided, run the ``cat_combustor`` with
# the ``run()`` method. Use the ``process_solution()`` method to retrieve the raw
# solution profiles. If you are interested in the exit solution, use the
# ``get_last_solution_mixture()`` method to get the solution at the exit as a
# ``Stream`` object ``combustor_exhaust``.
#
# .. note ::
#   In Pychemkin, the surface material is referred by its name. Instead of
#   the surface material index, you loop over the surface material names.
#   Use ``get_material_names()`` method to get a list of surface material
#   names defined in the surface mechanism.
#
# .. note ::
#   For a site or a bulk species, use any ``Mixture`` or ``Stream`` created
#   with the ``Chemistry Set`` containing the surface mechanism to get its
#   species index. There are two type of index for surface species, the
#   global index includes all species (the gas plus the site and the bulk of
#   all materials); the local index includes the same type of surface species
#   (site or bulk) of all materials.
#


# get gas-species index
i_co = cat_feed.get_specindex("CO")
i_ch4 = cat_feed.get_specindex("CH4")
i_no = cat_feed.get_specindex("NO")
# get surface species index: <global index> <local index>
# global index includes all species: the gas, the site of all materials, and the bulk
# of all materials
# local index includes the same type of species (site or bulk) of all materials
i_h_s_global, i_h_s_local = cat_feed.get_surf_specindex("H(S)")
i_pt_s_global, i_pt_s_local = cat_feed.get_surf_specindex("PT(S)")
# run the catalytic combustor PFR model
runstatus = cat_combustor.run()
# check run status
if runstatus != 0:
    # Run failed.
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# Run succeeded.
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
# compute the total runtime
runtime = time.time() - start_time
print(f"Total simulation duration: {runtime} [sec].")
# postprocess the solution profiles
cat_combustor.process_solution()
# get the exit/last solution stream from the catalytic combustor
combustor_exhaust = cat_combustor.get_last_solution_mixture()
# display the gas composition at the reactor exit
combustor_exhaust.list_composition(mode="mole")
# get the number of data points in the solution profiles
n_points = cat_combustor.getnumbersolutionpoints()
print(f"number of solution time points = {n_points}")
# store solution profiles
# distance from the reactor entrance [cm]
distance = cat_combustor.get_solution_variable_profile("distance")  #
# temperature solution profile along the reactor length [K]
temp_profile = cat_combustor.get_solution_variable_profile("temperature")
# total heat release rate from the gas-phase reactions [erg/sec]
gashrr_profile = cat_combustor.get_solution_variable_profile("gashrr")
gashrr_profile /= 1.0e3 * ck.ERGS_PER_JOULE
# total heat release rate from all surface reactions [erg/sec]
surfhrr_profile = cat_combustor.get_solution_variable_profile("surfhrr")
surfhrr_profile /= 1.0e3 * ck.ERGS_PER_JOULE
# gas phase CO solution profile [mass fraction]
co_profile = cat_combustor.get_solution_variable_profile("CO")
# gas phase CH4 solution profile [mass fraction]
ch4_profile = cat_combustor.get_solution_variable_profile("CH4")
# gas phase NO solution profile [mass fraction]
no_profile = cat_combustor.get_solution_variable_profile("NO")
# surface PT(S) site fraction solution profile [site fraction]
pt_s_profile = cat_combustor.get_solution_variable_profile("PT(S)")
# surface OH(S) site fraction solution profile [site fraction]
oh_s_profile = cat_combustor.get_solution_variable_profile("OH(S)")

# clean up
ck.done()

################################################
# Plot the catalytic combustor solution profiles
# ==============================================
# Plot the solution profiles along the ``cat_combustor`` length.
# You could observe the surface "light-off" point around x = 1.0 by the spikes
# in the solution profiles. The much bigger peak value of the surface heat release
# rate in comparison to the peak gas-phase heat release rate indicates the
# combustion process takes place at the catalyst surface.
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("PFR Catalytic Combustor Solution", fontsize=16)
plt.subplot(221)
plt.plot(distance, temp_profile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(distance, gashrr_profile, "b-")
plt.plot(distance, surfhrr_profile, "r-")
plt.legend(["Gas Phase", "Surface"], loc="upper right")
plt.ylabel("Total Heat Release Rate [kJ/sec]")
plt.subplot(223)
plt.plot(distance, co_profile, "g-")
plt.plot(distance, ch4_profile, "r--")
plt.plot(distance, no_profile, "b:")
plt.legend(["CO", "CH4", "NO"], loc="upper right")
plt.yscale("log")
plt.xlabel("Distance [cm]")
plt.ylabel("Gas Mass Fraction")
plt.subplot(224)
plt.plot(distance, pt_s_profile, "b-")
plt.plot(distance, oh_s_profile, "g:")
plt.legend(["PT(S)", "OH(S)"], loc="lower right")
plt.yscale("log")
plt.xlabel("Distance [cm]")
plt.ylabel("Surface Site Fraction")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_PFR_catalytic_combustion.png", bbox_inches="tight")

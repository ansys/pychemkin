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

"""CH4 catalytic Combustion."""

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# Chemkin PFR model (steady-state)
from ansys.chemkin.core.flowreactors.PFR import PFREnergyConservation as Pfr
from ansys.chemkin.core.inlet import Mixture, Stream
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
interactive = False

# Create the chemistry sets
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
# find the location of this example py file
try:
    script_dir_obj = Path(__file__).parent.resolve()
except NameError:
    script_dir_obj = Path(current_dir) / ".." / "data"
script_dir = str(script_dir_obj)
# use the relative path to locate the example mechanism folder
example_mech_data = script_dir_obj / ".." / "data"

# gas-phase chemistry only
mech_no_surface = ck.Chemistry(label="GRI3")
# set mechanism input files
# including the full file path is recommended
mech_no_surface.chemfile = str(data_dir / "grimech30_chem.inp")
mech_no_surface.thermfile = str(data_dir / "grimech30_thermo.dat")

# preprocess the mechanism files.
_ = mech_no_surface.preprocess()


# Set up the lean premixed fuel-air stream
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

# Set up and run the catalytic reactor
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

# Preprocess the catalytic combustion chemistry set
_ = mech_catalytic.preprocess()


# Set up the inlet stream for the catalytic combustor
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

# Create the PFR to predict the gas composition of the outlet stream
cat_combustor = Pfr(cat_feed, label="cat_combustor")

# Set up additional reactor model parameters
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

# Set up surface chemistry parameters
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

# Set solver controls
# set solver the tolerances
cat_combustor.tolerances = (1.0e-8, 1.0e-6)
# solver non-negative mode
cat_combustor.force_nonnegative = False
# transient solver max time step size [cm]
# cat_combustor.set_solver_max_timestep_size(1.0e-2)
# set adaptive solution saving distance
cat_combustor.adaptive_solution_saving(mode=True, steps=20)

# Run the catalytic combustor reactor model
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

# Plot the catalytic combustor solution profiles
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

resultfile = Path(current_dir) / "catalytic_combustion.result"
results = {}
results["state-distance"] = distance.tolist()
results["state-temperature"] = temp_profile.tolist()
results["species-CO_mass_fraction"] = co_profile.tolist()
results["species-PT(S)_site_fraction"] = pt_s_profile.tolist()
results["rate-total_surface_heat_release_rate"] = surfhrr_profile.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

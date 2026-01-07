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

"""Test for the spark ignition engine model."""

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
interactive = False

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on the gasoline 14 components mechanism
MyGasMech = ck.Chemistry(label="Gasoline")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "gasoline_14comp_WBencrypt.inp")
# preprocess the mechanism files
ierror = MyGasMech.preprocess()
print("mechanism information:")
print(f"number of gas species = {MyGasMech.kk:d}")
print(f"number of gas reactions = {MyGasMech.ii_gas:d}")
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
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
# create the unburned fuel-air mixture
fresh = ck.Mixture(MyGasMech)
# products from the complete combustion of the fuel mixture and air
products = ["co2", "h2o", "n2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros
# mean equivalence ratio
equiv = 1.0
ierror = fresh.x_by_equivalence_ratio(
    MyGasMech, fuelmixture.x, air.x, add_frac, products, equivalenceratio=equiv
)
if ierror != 0:
    raise RuntimeError
# list the composition of the unburned fuel-air mixture
fresh.list_composition(mode="mole")
# set mixture temperature and pressure
# (equivalent to setting the initial temperature and pressure of the reactor)
fresh.temperature = fuelmixture.temperature
fresh.pressure = fuelmixture.pressure
# set exhaust gas recirculation (EGR) ratio with volume fraction
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
# list the composition of the fuel+air+EGR mixture
fresh.list_composition(mode="mole", bound=1.0e-8)
# SI engine
# create an SI engine object
MyEngine = SIengine(reactor_condition=fresh)
# show initial gas composition inside the reactor
MyEngine.list_composition(mode="mole", bound=1.0e-8)
#
# set engine parameters
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
print(f"engine displacement volume {MyEngine.get_displacement_volume()} [cm3]")
print(f"engine clearance volume {MyEngine.get_clearance_volume()} [cm3]")
# set mass burned fraction profile for the SI engine combustion
# >>> option 1
# use Wiebe function
# start of combustion crank angle
MyEngine.set_burn_timing(soc=-14.5, duration=45.6)
MyEngine.wiebe_parameters(n=4.0, b=7.0)
# >>> option 2
# use anchor points
# MyEngine.set_burn_anchor_points(CA10=5.2, CA50=14.22, CA90=22.01)
# >>> option 3
# use burned mass fraction profile data
# start of combustion crank angle
# MyEngine.set_burn_timing(SOC=-14.5, duration=45.6)
# nMBFdata = 9
# MBangles = np.zeros(nMBFdata, dtype=np.double)
# MBFrac = np.zeros_like(MBangles, dtype=np.double)
# normalized crank angles (CA - SOC) / duration
# MBangles = [0.0, 0.26985, 0.3371, 0.3742, 0.4322, 0.63, 0.704, 0.8, 1.0]
# mass burned fraction
# MBFrac =   [0.0, 0.01, 0.03, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0]
# MyEngine.set_mass_burned_profile(crankangles=MBangles, fractions=MBFrac)
# <<<
# set minimum zonal mass [g]
MyEngine.set_minimum_zone_mass(minmass=1.0e-5)
# wall heat transfer model
# set model parameters
# "dimensionless": [<a> <b> <c> <Twall>]
# "dimensional": [<a> <b> <c> <Twall>]
# "hohenburg": [<a> <b> <c> <d> <e> <Twall>]
heattransferparameters = [0.1, 0.8, 0.0]
# set cylinder wall temperature [K]
Twall = 434.0
MyEngine.set_wall_heat_transfer("dimensionless", heattransferparameters, Twall)
# incylinder gas velocity correlation parameter (Woschni)
# [<C11> <C12> <C2> <swirl ratio>]
gv_parameters = [2.28, 0.318, 0.324, 0.0]
MyEngine.set_gas_velocity_correlation(gv_parameters)
# set piston head top surface area [cm2]
MyEngine.set_piston_head_area(area=56.75)
# set cylinder clearance surface area [cm2]
MyEngine.set_cylinder_head_area(area=56.75)
# output controls
# set the number of crank angles between saving solution
MyEngine.ca_step_for_saving_solution = 0.5
# set the number of crank angles between printing solution
MyEngine.ca_step_for_printing_solution = 10.0
# turn OFF adaptive solution saving
MyEngine.adaptive_solution_saving(mode=False, steps=20)
# turn OFF adaptive solution saving
# MyEngine.adaptive_solution_saving(mode=False)
# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyEngine.tolerances = (1.0e-15, 1.0e-6)
# get solver parameters
atol, rtol = MyEngine.tolerances
print(f"default absolute tolerance = {atol}")
print(f"default relative tolerance = {rtol}")
# turn on the force non-negative solutions option in the solver
MyEngine.force_nonnegative = True
# show solver option
# show the number of crank angles between printng solution
print(
    f"crank angles between solution printing: {MyEngine.ca_step_for_printing_solution}"
)
# show other transient solver setup
print(f"forced non-negative solution values: {MyEngine.force_nonnegative}")
# show the additional keywords given by user
MyEngine.showkeywordinputlines()
# set the start wall time
start_time = time.time()
# run the single-zone HCCI engine model
runstatus = MyEngine.run()
# compute the total runtime
runtime = time.time() - start_time
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
print(f"total simulation duration: {runtime} [sec]")
#
# post-process the solution profiles in selected zone
unburnedzone = 1
burnedzone = 2
zonestrings = ["Unburned Zone", "Burned Zone"]
thiszone = burnedzone
MyEngine.process_engine_solution(zone_id=thiszone)
plottitle = zonestrings[thiszone - 1] + " Solution"
# post-process cylinder-averged solution profiles
# MyMZEngine.process_average_engine_solution()
# plottitle = "Cylinder Averaged Solution"
# get the number of solution time points
solutionpoints = MyEngine.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
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
#
# post-process cylinder-averged solution
MyEngine.process_average_engine_solution()
# get the cylinder volume profile
cylindervolprofile = MyEngine.get_solution_variable_profile("volume")
# create arrays for cylinder-averaged mixture temperature
cylindertempprofile = MyEngine.get_solution_variable_profile("temperature")
#
# plot the profiles
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
    plt.savefig("spark_ignition_engine.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "sparkignitionengine.result"
results = {}
results["state-crank_angle"] = ca_profile.tolist()
results["state-temperature"] = tempprofile.tolist()
results["state-pressure"] = presprofile.tolist()
results["state-volume"] = volprofile.tolist()
results["species-CO_mole_fraction"] = co_profile.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

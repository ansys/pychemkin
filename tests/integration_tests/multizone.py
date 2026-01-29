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

"""Test of the multizone HCCI engine model."""

from pathlib import Path

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# chemkin homonegeous charge compression ignition (HCCI) engine model (transient)
from ansys.chemkin.core.engines.HCCI import HCCIengine
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

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
MyGasMech.tranfile = str(mechanism_dir / "grimech30_transport.dat")
# preprocess the mechanism files
ierror = MyGasMech.preprocess()
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.x = [("CH4", 0.9), ("C3H8", 0.05), ("C2H6", 0.05)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = 1.5 * ck.P_ATM
fuelmixture.temperature = 400.0
# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.x = [("O2", 0.21), ("N2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = 1.5 * ck.P_ATM
air.temperature = 400.0
# create the unburned fuel-air mixture
fresh = ck.Mixture(MyGasMech)
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros
# mean equivalence ratio
equiv = 0.8
ierror = fresh.x_by_equivalence_ratio(
    MyGasMech, fuelmixture.x, air.x, add_frac, products, equivalenceratio=equiv
)
if ierror != 0:
    raise RuntimeError
# list the composition of the unburned fuel-air mixture
fresh.list_composition(mode="mole")
# set mixture temperature and pressure
# (equivalent to setting the initial temperature and pressure of the reactor)
fresh.temperature = 447.0
fresh.pressure = 1.065 * ck.P_ATM
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
# HCCI engine
# create a 5-zones HCCI engine object
numbzones = 5
MyMZEngine = HCCIengine(reactor_condition=fresh, nzones=numbzones)
# show initial gas composition inside the reactor
MyMZEngine.list_composition(mode="mole", bound=1.0e-8)
#
# set engine parameters
# cylinder bore diameter [cm]
MyMZEngine.bore = 12.065
# engine stroke [cm]
MyMZEngine.stroke = 14.005
# connecting rod length [cm]
MyMZEngine.connecting_rod_length = 26.0093
# compression ratio [-]
MyMZEngine.compression_ratio = 16.5
# engine speed [RPM]
MyMZEngine.rpm = 1000
# set other parameters
# simulation start CA [degree]
MyMZEngine.starting_ca = -142.0
# simulation end CA [degree]
MyMZEngine.ending_ca = 116.0
# list the engine parameters
MyMZEngine.list_engine_parameters()
print(f"engine displacement volume {MyMZEngine.get_displacement_volume()} [cm3]")
print(f"engine clearance volume {MyMZEngine.get_clearance_volume()} [cm3]")
print(f"number of zone(s) = {MyMZEngine.get_number_of_zones()}")
# wall heat transfer model
# set model parameters
# "dimensionless": [<a> <b> <c> <Twall>]
# "dimensional": [<a> <b> <c> <Twall>]
# "hohenburg": [<a> <b> <c> <d> <e> <Twall>]
heattransferparameters = [0.035, 0.71, 0.0]
# set cylinder wall temperature [K]
t_wall = 400.0
MyMZEngine.set_wall_heat_transfer("dimensionless", heattransferparameters, t_wall)
# incylinder gas velocity correlation parameter (Woschni)
# [<C11> <C12> <C2> <swirl ratio>]
gv_parameters = [2.28, 0.308, 3.24, 0.0]
MyMZEngine.set_gas_velocity_correlation(gv_parameters)
# set piston head top surface area [cm2]
MyMZEngine.set_piston_head_area(area=124.75)
# set cylinder clearance surface area [cm2]
MyMZEngine.set_cylinder_head_area(area=123.5)
# set zonal properties
# zonal temperatures [K]
ztemperature = [447.5, 447.5, 447, 447, 447]
MyMZEngine.set_zonal_temperature(zonetemp=ztemperature)
# zonal volume fractions
zvolumefrac = [0.3, 0.25, 0.2, 0.2, 0.05]
MyMZEngine.set_zonal_volume_fraction(zonevol=zvolumefrac)
# wall heat transfer area fractions
zone_ht_area = [0.0, 0.15, 0.2, 0.25, 0.4]
MyMZEngine.set_zonal_heat_transfer_area_fraction(zonearea=zone_ht_area)
# zonal equivalence ratios
zphi = [equiv, equiv, equiv, equiv, equiv]
MyMZEngine.set_zonal_equivalence_ratio(zonephi=zphi)
# zonal EGR ratios
zone_egrr = [0.3, 0.3, 0.3, 0.35, 0.35]
MyMZEngine.set_zonal_egr_ratio(zoneegr=zone_egrr)
# set fuel "molar" composition
MyMZEngine.define_fuel_composition([("CH4", 0.9), ("C3H8", 0.05), ("C2H6", 0.05)])
# set oxidizer "molar' composition
MyMZEngine.define_oxid_composition([("O2", 0.21), ("N2", 0.79)])
# set products
MyMZEngine.define_product_composition(["CO2", "H2O", "N2"])
# set EGR composition in mole fractions
zadd = [add_frac, add_frac, add_frac, add_frac, add_frac]
MyMZEngine.define_additive_fractions(addfrac=zadd)
# output controls
# set the number of crank angles between saving solution
MyMZEngine.ca_step_for_saving_solution = 0.5
# set the number of crank angles between printing solution
MyMZEngine.ca_step_for_printing_solution = 10.0
# turn OFF adaptive solution saving
MyMZEngine.adaptive_solution_saving(mode=False, steps=20)
# turn OFF adaptive solution saving
# MyMZEngine.adaptive_solution_saving(mode=False)
# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyMZEngine.tolerances = (1.0e-12, 1.0e-10)
# get solver parameters
atol, rtol = MyMZEngine.tolerances
print(f"default absolute tolerance = {atol}")
print(f"default relative tolerance = {rtol}")
# turn on the force non-negative solutions option in the solver
MyMZEngine.force_nonnegative = True
# specify the ignition definitions
# ck.show_ignition_definitions()
MyMZEngine.set_ignition_delay(method="T_inflection")
# stop the simulation when ignition is detected
# MyMZEngine.stop_after_ignition()
# show solver option
# show the number of crank angles between printng solution
print(
    f"crank angles between solution printing: "
    f"{MyMZEngine.ca_step_for_printing_solution}"
)
# show other transient solver setup
print(f"forced non-negative solution values: {MyMZEngine.force_nonnegative}")
# show the additional keywords given by user
MyMZEngine.showkeywordinputlines()
# run the single-zone HCCI engine model
runstatus = MyMZEngine.run()
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
#
# get ignition delay "time"
delay_ca = MyMZEngine.get_ignition_delay()
print(f"ignition delay CA = {delay_ca} [degree]")
#
# get heat release information
hr10, hr50, hr90 = MyMZEngine.get_engine_heat_release_cas()
print("Engine Heat Release Information")
print(f"10% heat release CA = {hr10} [degree]")
print(f"50% heat release CA = {hr50} [degree]")
print(f"90% heat release CA = {hr90} [degree]\n")
#
# post-process the solution profiles in selected zone
thiszone = 1
MyMZEngine.process_engine_solution(zone_id=thiszone)
plottitle = "Zone " + str(thiszone) + " Solution"
# post-process cylinder-averged solution profiles
# MyMZEngine.process_average_engine_solution()
# plottitle = "Cylinder Averaged Solution"
# get the number of solution time points
solutionpoints = MyMZEngine.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile = MyMZEngine.get_solution_variable_profile("time")
# convert time to crank angle
ca_profile = np.zeros_like(timeprofile, dtype=np.double)
count = 0
for t in timeprofile:
    ca_profile[count] = MyMZEngine.get_ca(timeprofile[count])
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
    denprofile[i] = solutionmixture.rho
    # get zonal mixture viscosity profile [g/cm-sec] or [Poise]
    viscprofile[i] = solutionmixture.mixture_viscosity() * 1.0e2
#
# post-process cylinder-averged solution
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
    cylinderdenprofile[i] = solutionmixture.rho
#
# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle(plottitle, fontsize=16)
plt.subplot(221)
plt.plot(ca_profile, presprofile, "r-")
plt.ylabel("Pressure [bar]")
plt.subplot(222)
plt.plot(ca_profile, volprofile, "b-")
plt.plot(ca_profile, cylindervolprofile, "b--")
plt.ylabel("Volume [cm3]")
plt.legend(["Zone", "Cylinder"], loc="upper right")
plt.subplot(223)
plt.plot(ca_profile, denprofile, "g-")
plt.plot(ca_profile, cylinderdenprofile, "g--")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("Mixture Density [g/cm3]")
plt.legend(["Zone", "Averaged"], loc="upper left")
plt.subplot(224)
plt.plot(ca_profile, viscprofile, "m-")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("Mixture Viscosity [cP]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("multizone_HCCI_engine.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "multizone.result"
results = {}
results["state-crank_angle"] = ca_profile.tolist()
results["state-density"] = denprofile.tolist()
results["state-pressure"] = presprofile.tolist()
results["state-volume"] = volprofile.tolist()
results["state-viscosity"] = viscprofile.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

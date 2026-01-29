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

"""Test of the constrained volume option of the closed reactor."""

from pathlib import Path

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.core.batchreactors.batchreactor import (
    GivenVolumeBatchReactorEnergyConservation,
)
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
# create a chemistry set based on the diesel 14 components mechanism
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
fuelmixture.x = [("CH4", 1.0)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = 5.0 * ck.P_ATM
fuelmixture.temperature = 1500.0
# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.x = [("O2", 0.21), ("N2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = 5.0 * ck.P_ATM
air.temperature = 1500.0
# create the premixed mixture to be defined
premixed = ck.Mixture(MyGasMech)
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros
ierror = premixed.x_by_equivalence_ratio(
    MyGasMech, fuelmixture.x, air.x, add_frac, products, equivalenceratio=0.7
)
if ierror != 0:
    raise RuntimeError
# list the composition of the premixed mixture
premixed.list_composition(mode="mole")
# set mixture temperature and pressure
# (equivalent to setting the initial temperature and pressure of the reactor)
premixed.temperature = 800.0
premixed.pressure = 3.0 * ck.P_ATM
# Rapid Compression Machine
# create a constant volume batch reactor (with energy equation)
#
MyCONV = GivenVolumeBatchReactorEnergyConservation(premixed, label="RCM")
# set the initial reactor temperature (see the warning message in the run output)
# MyCONV.temperature = 800.0  # K
# show initial gas composition inside the reactor
MyCONV.list_composition(mode="mole")
# set other reactor properties
# reactor volume [cm3]
MyCONV.volume = 10.0
# simulation end time [sec]
MyCONV.time = 0.1
# set RCM volume profile (overriding the volume value set earlier)
# number of profile data points
npoints = 3
# position array of the profile data
x = np.zeros(npoints, dtype=np.double)
# value array of the profile data
volprofile = np.zeros_like(x, dtype=np.double)
# set reactor volume data points
x = [0.0, 0.01, 2.0]  # [sec]
volprofile = [10.0, 4.0, 4.0]  # [cm3]
# set the volume profile
MyCONV.set_volume_profile(x, volprofile)
# output controls
# set timestep between saving solution
MyCONV.timestep_for_saving_solution = 0.01
# turn OFF adaptive solution saving
MyCONV.adaptive_solution_saving(mode=False, value_change=100, target="TEMPERATURE")
# turn OFF adaptive solution saving
# MyCONV.adaptive_solution_saving(mode=False)
# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyCONV.tolerances = (1.0e-10, 1.0e-8)
# get solver parameters
atol, rtol = MyCONV.tolerances
print(f"default absolute tolerance = {atol}")
print(f"default relative tolerance = {rtol}")
# turn on the force non-negative solutions option in the solver
MyCONV.force_nonnegative = True
# specify the ignition definitions
# ck.show_ignition_definitions()
MyCONV.set_ignition_delay(method="T_inflection")
# stop the simulation when ignition is detected
# MyCONV.stop_after_ignition()
# show solver option
print(f"timestep between solution printing: {MyCONV.timestep_for_printing_solution}")
# show timestep between printing solution
print(f"forced non-negative solution values: {MyCONV.force_nonnegative}")
# show the additional keywords given by user
MyCONV.showkeywordinputlines()
# run the CONV reactor model
runstatus = MyCONV.run()
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
# get ignition delay time (need to deduct the initial compression time = 0.01 [sec])
delaytime = MyCONV.get_ignition_delay() - 0.01 * 1.0e3
print(f"ignition delay time = {delaytime} [msec]")
# post-process the solutions
MyCONV.process_solution()
# get the number of solution time points
solutionpoints = MyCONV.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile = MyCONV.get_solution_variable_profile("time")
# get the temperature profile
tempprofile = MyCONV.get_solution_variable_profile("temperature")
# get the volume profile
volprofile = MyCONV.get_solution_variable_profile("volume")
# get CH4 mass fraction profile
# CH4massfraction = MyCONV.get_solution_variable_profile("CH4")
#
# more involving post-processing by using Mixtures
#
# mass
massprofile = np.zeros_like(timeprofile, dtype=np.double)
# create arrays for CH4 mole fraction, CH4 ROP, and mixture viscosity
ch4_profile = np.zeros_like(timeprofile, dtype=np.double)
ch4_rop_profile = np.zeros_like(timeprofile, dtype=np.double)
viscprofile = np.zeros_like(timeprofile, dtype=np.double)
current_rop = np.zeros(MyGasMech.kk, dtype=np.double)
# find CH4 species index
ch4_index = MyGasMech.get_specindex("CH4")
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = MyCONV.get_solution_mixture_at_index(solution_index=i)
    # get gas density [g/cm3]
    den = solutionmixture.rho
    # reactor mass [g]
    massprofile[i] = den * volprofile[i]
    # get CH4 mole fraction profile
    ch4_profile[i] = solutionmixture.x[ch4_index]
    # get CH4 ROP profile
    current_rop = solutionmixture.rop()
    ch4_rop_profile[i] = current_rop[ch4_index]
    # get mixture vicosity profile
    viscprofile[i] = solutionmixture.mixture_viscosity()
# validation
del_mass = np.zeros_like(timeprofile, dtype=np.double)
mass0 = massprofile[0]
for i in range(solutionpoints):
    del_mass[i] = abs(massprofile[i] - mass0)
#
print(f">>> maximum magnitude of reactor mass deviation = {np.max(del_mass)} [g]")
# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, ch4_profile, "b-")
plt.ylabel("CH4 Mole Fraction")
plt.subplot(223)
plt.plot(timeprofile, ch4_rop_profile, "g-")
plt.xlabel("time [sec]")
plt.ylabel("CH4 Production Rate [mol/cm3-sec]")
plt.subplot(224)
plt.plot(timeprofile, viscprofile, "m-")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Viscosity [g/cm-sec]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("CONV_solution.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "CONV.result"
results = {}
results["state-time"] = timeprofile.tolist()
results["state-temperature"] = tempprofile.tolist()
results["species-CH4_mole_fraction"] = ch4_profile.tolist()
results["rate-CH4_production_rate"] = ch4_rop_profile.tolist()
results["state-viscocity"] = viscprofile.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

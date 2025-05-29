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
import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.batchreactors.batchreactor import (
    GivenPressureBatchReactor_EnergyConservation,
)
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
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the diesel 14 components mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
MyGasMech.tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.X = [("H2", 2.0), ("N2", 3.76), ("O2", 1.0)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = ck.Patm
fuelmixture.temperature = 1000
# products from the complete combustion of the fuel mixture and air
products = ["H2O", "N2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
# add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros
# iError = premixed.X_by_Equivalence_Ratio(
#     MyGasMech, fuelmixture.X, air.X, add_frac, products, equivalenceratio=0.7
# )
# if iError != 0:
#     raise RuntimeError
# list the composition of the premixed mixture
fuelmixture.list_composition(mode="mole")
# set mixture temperature and pressure (equivalent to setting the initial temperature and pressure of the reactor)

# Rapid Compression Machine
# create a constant volume batch reactor (with energy equation)
#
MyCONV = GivenPressureBatchReactor_EnergyConservation(fuelmixture, label="tran")
# set the initial reactor temperature (see the warning message in the run output)
# MyCONV.temperature = 800.0  # K
# show initial gas composition inside the reactor
MyCONV.list_composition(mode="mole")
# set other reactor properties
# reactor volume [cm3]
MyCONV.volume = 1
MyCONV.temperature = 1000
# simulation end time [sec]
MyCONV.time = 0.0005
# set RCM volume profile (overriding the volume value set earlier)
# # number of profile data points
# npoints = 3
# # position array of the profile data
# x = np.zeros(npoints, dtype=np.double)
# # value array of the profile data
# volprofile = np.zeros_like(x, dtype=np.double)
# # set reactor volume data points
# x = [0.0, 0.01, 2.0]  # [sec]
# volprofile = [10.0, 4.0, 4.0]  # [cm3]
# set the volume profile
# MyCONV.set_volume_profile(x, volprofile)
# output controls
# set timestep between saving solution
# MyCONV.timestep_for_saving_solution = 0.01
# turn OFF adaptive solution saving
MyCONV.adaptive_solution_saving(mode=False, steps=20)
# turn OFF adaptive solution saving
# MyCONV.adaptive_solution_saving(mode=False)
# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyCONV.tolerances = (1.0e-20, 1.0e-8)
# get solver parameters
ATOL, RTOL = MyCONV.tolerances
print(f"default absolute tolerance = {ATOL}")
print(f"default relative tolerance = {RTOL}")
# turn on the force non-negative solutions option in the solver
MyCONV.force_nonnegative = True
# specify the ignition definitions
ck.show_ignition_definitions()
MyCONV.set_ignition_delay(method="T_rise", val=400)
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
    print(Color.RED + ">>> Run failed. <<<", end="\n" + Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end="\n" + Color.END)
# get ignition delay time (need to deduct the initial compression time = 0.01 [sec])
# delaytime = MyCONV.get_ignition_delay()
# print(f"ignition delay time = {delaytime} [msec]")
# post-process the solutions
MyCONV.process_solution()
# get the number of solution time points
solutionpoints = MyCONV.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile = MyCONV.get_solution_variable_profile("time")
# get the temperature profile
tempprofile = MyCONV.get_solution_variable_profile("temperature")
# more involving post-processing by using Mixtures
# create arrays for CH4 mole fraction, CH4 ROP, and mixture viscosity
H2Oprofile = np.zeros_like(timeprofile, dtype=np.double)
H2OROPprofile = np.zeros_like(timeprofile, dtype=np.double)
denprofile = np.zeros_like(timeprofile, dtype=np.double)
CurrentROP = np.zeros(MyGasMech.KK, dtype=np.double)
# find CH4 species index
H2O_index = MyGasMech.get_specindex("H2O")
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = MyCONV.get_solution_mixture_at_index(solution_index=i)
    # get gas density [g/cm3]
    denprofile[i] = solutionmixture.RHO
    # reactor mass [g]
    # get CH4 mole fraction profile
    H2Oprofile[i] = solutionmixture.X[H2O_index]
    # get CH4 ROP profile
    currentROP = solutionmixture.ROP()
    H2OROPprofile[i] = currentROP[H2O_index]
# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, H2Oprofile, "b-")
plt.ylabel("H2O Mole Fraction")
plt.subplot(223)
plt.plot(timeprofile, H2OROPprofile, "g-")
plt.xlabel("time [sec]")
plt.ylabel("H2O Production Rate [mol/cm3-sec]")
plt.subplot(224)
plt.plot(timeprofile, denprofile, "m-")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Density [g/cm3]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("close_homogeneous.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "closed_homogeneous__transient.result")
results = {}
results["state-time"] = timeprofile.tolist()
results["state-temperature"] = tempprofile.tolist()
results["species-H2O_mole_fraction"] = H2Oprofile.tolist()
results["rate-H2O_production_rate"] = H2OROPprofile.tolist()
results["state-density"] = denprofile.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

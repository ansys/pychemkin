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
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color

# chemkin plug flow reactor model
from ansys.chemkin.flowreactors.PFR import PlugFlowReactor_FixedTemperature
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the inlet (mixture + flow rate)
feedstock = Stream(MyGasMech)
# set inlet temperature [K]
feedstock.temperature = 1444.48
# set inlet/PFR pressure [atm]
feedstock.pressure = 0.83 * ck.Patm
# set inlet composition
feedstock.X = [
    ("AR", 0.8433),
    ("CO", 0.0043),
    ("CO2", 0.0429),
    ("H2O", 0.0956),
    ("N2", 0.0031),
    ("NH3", 0.0021),
    ("NO", 0.0012),
    ("O2", 0.0074),
    ("OH", 4.6476e-5),
]
# set inlet velocity [cm/sec]
feedstock.velocity = 26.815
#
# create a plug flow reactor instance
tubereactor = PlugFlowReactor_FixedTemperature(feedstock)
# set PFR diameter [cm]
tubereactor.diameter = 5.8431
# set PFR length [cm]
tubereactor.length = 5.0
print(f"PFR inlet mass flow rate {tubereactor.mass_flowrate} [g/sec]")
print(f"PFR inlet velocity {tubereactor.velocity} [cm/sec]")
# show inlet gas composition of the PFR
print("PFR inlet gas compsition")
tubereactor.list_composition(mode="mole", bound=1.0e-8)
# set distance between saving solution
tubereactor.timestep_for_saving_solution = 0.0005
# turn OFF adaptive solution saving
tubereactor.adaptive_solution_saving(mode=False, steps=100)
# show the additional keywords given by user
tubereactor.showkeywordinputlines()
# set the start wall time
start_time = time.time()
# run the PFR model
runstatus = tubereactor.run()
# compute the total runtime
runtime = time.time() - start_time
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
print(f"total simulation duration: {runtime * 1.0e3} [msec]")
#
# post-process the solution profiles
tubereactor.process_solution()
# get the number of solution time points
solutionpoints = tubereactor.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the grid profile [cm]
xprofile = tubereactor.get_solution_variable_profile("time")
# get the temperature profile [K]
tempprofile = tubereactor.get_solution_variable_profile("temperature")
# get the NO mass fraction profile
YNOprofile = tubereactor.get_solution_variable_profile("NO")
# outlet grid index
xout_index = solutionpoints - 1
print(f"At the reactor outlet: x = {xprofile[xout_index]} [cm]")
print(f"the NO mass fraction = {YNOprofile[xout_index]}")
#
# more involving post-processing by using Mixtures
#
# create arrays for CO, NH3, and NO2 mole fractions
COprofile = np.zeros_like(xprofile, dtype=np.double)
NH3profile = np.zeros_like(xprofile, dtype=np.double)
NO2profile = np.zeros_like(xprofile, dtype=np.double)
velocityprofile = np.zeros_like(xprofile, dtype=np.double)
# find species index
CO_index = MyGasMech.get_specindex("CO2")
NH3_index = MyGasMech.get_specindex("NH3")
NO2_index = MyGasMech.get_specindex("NO2")
# reactor mass flow rate (constant) [g/sec]
massflowrate = tubereactor.mass_flowrate
# reactor cross-section area [cm2]
areaflow = tubereactor.flowarea
print(f"mass flow rate: {massflowrate} [g/sec]\nflow area: {areaflow} [cm2]")
# ratio
ratio = massflowrate / areaflow
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tubereactor.get_solution_mixture_at_index(solution_index=i)
    # get gas density [g/cm3]
    den = solutionmixture.RHO
    # gas velocity [g]
    velocityprofile[i] = ratio / den
    # get CO mole fraction profile
    COprofile[i] = solutionmixture.X[CO_index]
    # get NH3 mole fraction profile
    NH3profile[i] = solutionmixture.X[NH3_index]
    # get NO2 mole fraction profile
    NO2profile[i] = solutionmixture.X[NO2_index]
# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("Constant Temperature Plug-Flow Reactor", fontsize=16)
plt.subplot(221)
plt.plot(xprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(xprofile, COprofile, "b-")
plt.ylabel("CO Mole Fraction")
plt.subplot(223)
plt.plot(xprofile, NO2profile, "g-")
plt.xlabel("distance [cm]")
plt.ylabel("NO2 Mole Fraction")
plt.subplot(224)
plt.plot(xprofile, velocityprofile, "m-")
plt.xlabel("distance [cm]")
plt.ylabel("Gas Velocity [cm/sec]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plug_flow_reactor.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "plugflow.result")
results = {}
results["state-distance"] = xprofile.tolist()
results["state-temperature"] = tempprofile.tolist()
results["state-velocity"] = velocityprofile.tolist()
results["species-CO_mole_fraction"] = COprofile.tolist()
results["species-NO2_mole_fraction"] = NO2profile.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

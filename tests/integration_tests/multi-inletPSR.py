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
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.stirreactors.PSR import PSR_SetVolume_EnergyConservation as PSR
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
data_dir = os.path.join(
    ck.ansys_dir, "reaction", "data", "ModelFuelLibrary", "Skeletal"
)
mechanism_dir = data_dir
# create a chemistry set based on the gasoline 14 components mechanism
MyGasMech = ck.Chemistry(label="hydrogen")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(
    mechanism_dir, "Hydrogen-Ammonia-NOx_chem_MFL2021.inp"
)
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# create the fuel inlet
fuel = Stream(MyGasMech, label="Fuel")
# set fuel composition
fuel.X = [("h2", 0.21), ("n2", 0.79)]
# setting pressure and temperature is not required in this case
fuel.pressure = ck.Patm
fuel.temperature = 450.0  # inlet temperature
# set inlet volumetric flow rate [cm3/sec]
fuel.vol_flowrate = 25.0
# create the oxidizer inlet: air
air = Stream(MyGasMech, label="Oxid")
air.X = [("o2", 0.21), ("n2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = fuel.pressure
air.temperature = fuel.temperature
# set inlet volumetric flow rate [cm3/sec]
air.vol_flowrate = 50.0
# create a PSR with fixed reactor volume and
# with the fuel inlet composition as the estimated reactor condition
combustor = PSR(fuel, label="tincan")
# set the estimated reactor temperature [K]
combustor.temperature = 2000.0
# set the reactor volume (cm3): required for PSR_SetVolume_EnergyConservation model
combustor.volume = 200.0
# add external inlets to the PSR
combustor.set_inlet(fuel)
combustor.set_inlet(air)
# reset the tolerances in the steady-state solver
combustor.steady_state_tolerances = (1.0e-9, 1.0e-6)
combustor.timestepping_tolerances = (1.0e-9, 1.0e-6)
# reset the gas species floor value in the steady-state solver
combustor.set_species_floor(-1.0e-10)
# reactor volume increment
deltaVol = -5
numbruns = 9
# solution arrays
residencetime = np.zeros(numbruns, dtype=np.double)
tempSSsolution = np.zeros_like(residencetime, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all inlet temperature values
for i in range(numbruns):
    # run the PSR model
    runstatus = combustor.run()
    # check run status
    if runstatus != 0:
        # run failed!
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()
    # run success!
    print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
    # post-process the solution profiles
    solnmixture = combustor.process_solution()
    # print the steady-state solution values
    print(f"steady-state temperature = {solnmixture.temperature} [K]")
    # solnmixture.list_composition(mode="mole")
    # store solution values
    # final reactor gas density [g/cm3]
    # density = solnmixture.RHO
    # final reactor mass [g]
    # mass = density * combustor.volume
    # PSR apparent residence time [sec]
    residencetime[i] = combustor.volume / combustor.net_vol_flowrate
    tempSSsolution[i] = solnmixture.temperature
    # update reactor volume
    combustor.volume += deltaVol

# compute the total runtime
runtime = time.time() - start_time
print(f"total simulation duration: {runtime} [sec] over {numbruns} runs")
#
# plot results
plt.plot(residencetime, tempSSsolution, "bo-")
plt.xlabel("Apparent Residence Time [sec]")
plt.ylabel("Exit Gas Temperature [K]")
plt.title("PSR Solution")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("multi_inlet_PSR.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "multi-inletPSR.result")
results = {}
results["state-residence_time"] = residencetime.tolist()
results["state-temperature"] = tempSSsolution.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

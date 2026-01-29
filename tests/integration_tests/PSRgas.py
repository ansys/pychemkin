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

"""Test of steady state PSR model with the residence time option."""

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream  # external gaseous inlet
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.stirreactors.PSR import PSRSetResTimeEnergyConservation as Psr
from ansys.chemkin.core.utilities import find_file

# Chemkin PSR model (steady-state)

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
data_dir = Path(ck.ansys_dir) / "reaction" / "data" / "ModelFuelLibrary" / "Skeletal"
mechanism_dir = data_dir
# create a chemistry set based on the gasoline 14 components mechanism
MyGasMech = ck.Chemistry(label="hydrogen")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = find_file(
    str(mechanism_dir),
    "Hydrogen-Ammonia-NOx_chem_MFL",
    "inp",
)
# preprocess the mechanism files
ierror = MyGasMech.preprocess()
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the fuel mixture
fuel = ck.Mixture(MyGasMech)
# set fuel composition
fuel.x = [("h2", 0.8), ("n2", 0.2)]
# setting pressure and temperature is not required in this case
fuel.pressure = ck.P_ATM
fuel.temperature = 298.0  # inlet temperature
# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.x = [("o2", 0.21), ("n2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = fuel.pressure
air.temperature = fuel.temperature
# create the fuel-oxidizer inlet to the PSR
feed = Stream(MyGasMech, label="feed_1")
# products from the complete combustion of the fuel mixture and air
products = ["h2o", "n2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros
# mean equivalence ratio
equiv = 1.0
ierror = feed.x_by_equivalence_ratio(
    MyGasMech, fuel.x, air.x, add_frac, products, equivalenceratio=equiv
)
if ierror != 0:
    raise RuntimeError
# setting reactor pressure [dynes/cm2]
feed.pressure = fuel.pressure
# set inlet gas temperature [K]
feed.temperature = fuel.temperature
# set inlet mass flow rate [g/sec]
feed.mass_flowrate = 432.0
# create the Jet-Stirred Reactor
# use the inlet gas property as the estimated reactor condition
sphere = Psr(feed, label="PSR_1")
# set the estimated reactor temperature [K]
sphere.temperature = 1700.0
# connect the inlet to the reactor
sphere.set_inlet(feed)
# set PSR residence time (sec): required for PSRSetResTimeEnergyConservation model
sphere.residence_time = 3.0 * 1.0e-5
# reset the tolerances in the steady-state solver
sphere.steady_state_tolerances = (1.0e-9, 1.0e-6)
sphere.timestepping_tolerances = (1.0e-9, 1.0e-6)
# reset the gas species floor value in the steady-state solver
sphere.set_species_floor(-1.0e-10)
# inlet gas equivalence ratio increment
deltaequiv = 0.05
numbruns = 9
# solution arrays
inletequiv = np.zeros(numbruns, dtype=np.double)
temp_solution = np.zeros_like(inletequiv, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all inlet temperature values
for i in range(numbruns):
    # run the PSR model
    runstatus = sphere.run()
    # check run status
    if runstatus != 0:
        # run failed!
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()
    # run success!
    print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
    # post-process the solution profiles
    solnmixture = sphere.process_solution()
    # print the steady-state solution values
    # print(f"steady-state temperature = {solnmixture.temperature} [K]")
    # solnmixture.list_composition(mode="mole")
    # store solution values
    inletequiv[i] = equiv
    temp_solution[i] = solnmixture.temperature
    # update inlet gas equivalence ratio (composition)
    equiv += deltaequiv
    ierror = feed.x_by_equivalence_ratio(
        MyGasMech, fuel.x, air.x, add_frac, products, equivalenceratio=equiv
    )
    if ierror != 0:
        print(f"error encountered with inlet equivalence ratio = {equiv}")
        raise RuntimeError

# compute the total runtime
runtime = time.time() - start_time
print(f"total simulation duration: {runtime} [sec] over {numbruns} runs")
#
# plot results
plt.plot(inletequiv, temp_solution, "b-")
plt.xlabel("Inlet Gas Equivalence Ratio")
plt.ylabel("Reactor Temperature [K]")
plt.title("PSR Solution")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("PSR_gas.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "PSRgas.result"
results = {}
results["state-equivalence_ratio_inlet"] = inletequiv.tolist()
results["state-temperature"] = temp_solution.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

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

"""Test for steady state PSR model."""

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream  # external gaseous inlet
from ansys.chemkin.core.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.core.stirreactors.PSR import PSRSetResTimeFixedTemperature as Psr
from ansys.chemkin.core.utilities import find_file

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
# create a chemistry set based on the hydrogen-ammonia mechanism
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
# create a premixed fuel-oxidizer mixture
# create the fuel-oxidizer inlet to the JSR
feed = Stream(MyGasMech)
# set H2-O2-N2 composition
feed.x = [("h2", 1.1e-2), ("n2", 9.62e-1), ("o2", 2.75e-2)]
# setting reactor pressure [dynes/cm2]
feed.pressure = ck.P_ATM
# set inlet gas temperature [K]
temp = 800.0
feed.temperature = temp
# set inlet mass flow rate [g/sec]
feed.mass_flowrate = 0.11
# create the Jet-Stirred Reactor
# use the inlet gas property as the estimated reactor condition
JSR = Psr(feed, label="JSR")
# connect the inlet to the reactor
JSR.set_inlet(feed)
# set PSR residence time (sec): required for PSRSetResTimeFixedTemperature model
JSR.residence_time = 120.0 * 1.0e-3
# set the number of initial pseudo time steps in the steady-state solver
JSR.set_initial_timesteps(1000)
# inlet gas temperature increment
deltatemp = 25.0
numbruns = 19
# find H2O species index
h2o_index = MyGasMech.get_specindex("h2o")
# solution arrays
inlet_temp = np.zeros(numbruns, dtype=np.double)
h2o_solution = np.zeros_like(inlet_temp, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all inlet temperature values
for i in range(numbruns):
    # run the PSR model
    runstatus = JSR.run()
    # check run status
    if runstatus != 0:
        # run failed!
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()
    # run success!
    print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
    # post-process the solution profiles
    solnmixture = JSR.process_solution()
    # print the steady-state solution values
    # print(f"steady-state temperature = {solnmixture.temperature} [K]")
    # solnmixture.list_composition(mode="mole")
    # store solution values
    inlet_temp[i] = solnmixture.temperature
    h2o_solution[i] = solnmixture.x[h2o_index]
    # update reactor temperature
    temp += deltatemp
    JSR.temperature = temp
# compute the total runtime
runtime = time.time() - start_time
print(f"total simulation duration: {runtime} [sec] over {numbruns} runs")
#
# experimental data
# JSR temperature [K]
TEMP_data = [
    803.0,
    823.0,
    850.0,
    875.0,
    902.0,
    925.0,
    951.0,
    973.0,
    1002.0,
    1023.0,
    1048.0,
]
# measured H2O mole fractions in the JSR exit flow
H2O_data = [
    0.000312,
    0.000313,
    0.000318,
    0.000303,
    0.003292,
    0.006226,
    0.008414,
    0.009745,
    0.010103,
    0.011036,
    0.010503,
]
# plot results
plt.plot(inlet_temp, h2o_solution, "b-", label="prediction")
plt.plot(TEMP_data, H2O_data, "ro", label="data", markersize=4, fillstyle="none")
plt.xlabel("Reactor Temperature [K]")
plt.ylabel("H2O Mole Fraction")
plt.legend(loc="lower right")
plt.title("JSR Solution")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("jet_stirred_reactor.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "jetstirredreactor.result"
results = {}
results["state-temperature_inlet"] = inlet_temp.tolist()
results["species-H2O_mole_fraction"] = h2o_solution.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

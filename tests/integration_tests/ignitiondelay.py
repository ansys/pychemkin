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

"""Ignition delay time test for the closed homogeneous reactor."""

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # chemkin
from ansys.chemkin.core import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.core.batchreactors.batchreactor import (
    GivenPressureBatchReactorEnergyConservation,
)
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(False)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
gasoline = ck.Chemistry(label="gasoline 14comp")
# set mechanism input files
# including the full file path is recommended
gasoline.chemfile = str(mechanism_dir / "gasoline_14comp_WBencrypt.inp")
# preprocess the mechanism files
ierror = gasoline.preprocess()
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the fuel mixture
fuelmixture = ck.Mixture(gasoline)
# set fuel = composition PRF 60
fuelmixture.x = [("ic8h18", 0.6), ("nc7h16", 0.4)]
# setting pressure and temperature
fuelmixture.pressure = 5.0 * ck.P_ATM
fuelmixture.temperature = 1500.0
# create the oxidizer mixture: air
air = ck.Mixture(gasoline)
air.x = [("o2", 0.21), ("n2", 0.79)]
# setting pressure and temperature
air.pressure = 5.0 * ck.P_ATM
air.temperature = 1500.0
# create the premixed mixture to be defined
premixed = ck.Mixture(gasoline)
# products from the complete combustion of the fuel mixture and air
products = ["co2", "h2o", "n2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(gasoline.kk, dtype=np.double)  # no additives: all zeros
ierror = premixed.x_by_equivalence_ratio(
    gasoline, fuelmixture.x, air.x, add_frac, products, equivalenceratio=1.0
)
if ierror != 0:
    raise RuntimeError
# list the composition of the premixed mixture
premixed.list_composition(mode="mole")
# set mixture temperature and pressure
# (equivalent to setting the initial temperature and pressure of the reactor)
premixed.pressure = 40.0 * ck.P_ATM
premixed.temperature = 700.0
#
# create a constant pressure batch reactor (with energy equation)
#
MyCONP = GivenPressureBatchReactorEnergyConservation(premixed, label="CONP")
# show initial gas composition inside the reactor
MyCONP.list_composition(mode="mole")
# set other reactor parameters
# reactor volume [cm3]
MyCONP.volume = 10.0
# simulation end time [sec]
MyCONP.time = 1.0
# output controls
# set timestep between saving solution
MyCONP.timestep_for_saving_solution = 0.001
# change timestep between saving solution
MyCONP.timestep_for_saving_solution = 0.01
# turn ON adaptive solution saving
MyCONP.adaptive_solution_saving(mode=True, value_change=100, target="TEMPERATURE")
# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyCONP.tolerances = (1.0e-10, 1.0e-8)
# set ignition delay
# ck.show_ignition_definitions()
MyCONP.set_ignition_delay(method="T_inflection")
# stop after ignition is detected
# (not recommended for ignition delay time calculations)
# MyCONP.stop_after_ignition()
# show solver option
print(f"timestep between solution printing: {MyCONP.timestep_for_printing_solution}")
# show timestep between printing solution
print(f"forced non-negative solution values: {MyCONP.force_nonnegative}")
#
# loop over initial reactor temperature to create an ignition delay time plot
#
npoints = 20
delta_temp = 20.0
init_temp = premixed.temperature
delaytime = np.zeros(npoints, dtype=np.double)
temp_inv = np.zeros_like(delaytime, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all cases with different initial gas/reactor temperatures
for i in range(npoints):
    # update the initial reactor temperature
    MyCONP.temperature = init_temp  # K
    # show the additional keywords given by user
    # (verify inputs before running the simulation)
    # MyCONP.showkeywordinputlines()
    # run the reactor model
    runstatus = MyCONP.run()
    #
    if runstatus == 0:
        # plot 1/T instead of T
        temp_inv[i] = 1.0e0 / init_temp
        # get ignition delay time
        delaytime[i] = MyCONP.get_ignition_delay()
        print(f"ignition delay time = {delaytime[i]} [msec]")
    else:
        # if get this, most likely the END time is too short
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    init_temp += delta_temp
# compute the total runtime
runtime = time.time() - start_time
print(f"total simulation duration: {runtime} [sec] for {npoints} cases")
# create an ignition delay versus 1/T plot for the PRF fuel
# (should exhibit the NTC region)
plt.rcParams.update({"figure.autolayout": True})
fig, ax1 = plt.subplots()
ax1.semilogy(temp_inv, delaytime, "bs--")
ax1.set_xlabel("1/T [1/K]")
ax1.set_ylabel("Ignition delay time [msec]")


# Create a secondary x-axis for T (=1/(1/T))
def one_over(x):
    """Vectorized 1/x, treating x==0 manually."""
    x = np.array(x, float)
    near_zero = np.isclose(x, 0)
    x[near_zero] = np.inf
    x[~near_zero] = 1 / x[~near_zero]
    return x


# the function "1/x" is its own inverse
inverse = one_over
ax2 = ax1.secondary_xaxis("top", functions=(one_over, inverse))
ax2.set_xlabel("T [K]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("ignition_delay.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "ignitiondelay.result"
results = {}
results["state-temperature_inverse"] = temp_inv.tolist()
results["state-ignition_delay"] = delaytime.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

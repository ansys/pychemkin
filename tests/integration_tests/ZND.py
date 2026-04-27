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

"""ZND analysis test for the shick tube model."""

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream
from ansys.chemkin.core.logger import logger

# chemkin plug flow reactor model
from ansys.chemkin.core.shock.shocktubereactors import ZNDCalculator as Znd
from ansys.chemkin.core.utilities import find_file

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# set numpy printing option
np.set_printoptions(legacy="1.25")

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data" / "ModelFuelLibrary" / "Skeletal"
mechanism_dir = data_dir
# create a chemistry set based on the MFL 2021 hydrogen mechanism
MyGasMech = ck.Chemistry(label="hydrogen")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = find_file(str(mechanism_dir), "Hydrogen_chem_MFL", "inp")

# preprocess the mechanism files
ierror = MyGasMech.preprocess()
if ierror != 0:
    print("Error: Failed to preprocess the mechanism!")
    print(f"       Error code = {ierror}")
    exit()

#
stream_a = Stream(MyGasMech)
# initial gas state before the incident shock front
stream_a = Stream(MyGasMech)
# a diluted initial gas state before the incident shock front
# 1 [atm]
stream_a.pressure = ck.P_ATM
# 300 [K]
stream_a.temperature = 300.0
# the 50% hydrogen + 50% oxygen mixture diluted by 40% nitrogen by volume
stream_a.x = [("h2", 0.2), ("o2", 0.2), ("n2", 0.6)]

# Create the shock tube reactor object
ZNDIncident = Znd(stream_a, label="ZND")
# set total simulation time (particle time) [sec]
ZNDIncident.time = 3.0e-6
# tolerances are given in tuple: (absolute tolerance, relative tolerance)
ZNDIncident.tolerances = (1.0e-10, 1.0e-6)

# set the start wall time
start_time = time.time()
# run the ZND model
runstatus = ZNDIncident.run()
# compute the total runtime
runtime = time.time() - start_time
# check run status
if runstatus != 0:
    # Run failed.
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run succeeded.
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
print(f"Total simulation duration: {runtime * 1.0e3} [msec]")

# postprocess the solution profiles
ZNDIncident.process_solution()

# get the number of solution time points
solutionpoints = ZNDIncident.get_solution_size()
print(f"number of solution points = {solutionpoints}")
# get information about the induction length
numb_induct_length = ZNDIncident.get_inductionlength_size()
print(f"number of induction length = {numb_induct_length}")
if numb_induct_length > 0:
    induct_legth, sigmamax = ZNDIncident.get_induction_lengths(numb_induct_length)
    for i in range(numb_induct_length):
        print(f"Induction length # {i + 1}:")
        print(f"        Induction length = {induct_legth[i] * 1.0e1} [mm].")
        print(f"    Local max thermicity = {sigmamax[i]} [1/sec].\n")

# get the grid profile [cm]
xprofile = ZNDIncident.get_solution_variable_profile("distance")
# get the temperature profile [K]
tempprofile = ZNDIncident.get_solution_variable_profile("temperature")
# get the pressure profile [dynes/cm2]
pressprofile = ZNDIncident.get_solution_variable_profile("pressure")
# get the velocity profile [cm/sec]
velprofile = ZNDIncident.get_solution_variable_profile("velocity")
# get the total thermiocity profile [1/sec]
sigmaprofile = ZNDIncident.get_solution_variable_profile("thermicity")

# create arrays for the gas Mach number profile
machprofile = np.zeros_like(xprofile, dtype=np.double)
# total thermicity profile
sigmaprofile = np.zeros_like(xprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = ZNDIncident.get_solution_stream_at_index(solution_index=i)
    # get total thermicity [1/sec]
    sigmaprofile[i] = solutionmixture.thermicity()
    # gas speed of sound [cm/sec]
    soundspeed = solutionmixture.sound_speed()
    # gas Mach number
    machprofile[i] = velprofile[i] / soundspeed
    # convert pressure from [dynes/cm2] to [bar]
    pressprofile[i] /= 1.0e6

# Plot the solution profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("ZND Analysis", fontsize=16)
plt.subplot(221)
plt.plot(xprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(xprofile, pressprofile, "b-")
plt.ylabel("Pressure [bar]")
plt.subplot(223)
plt.plot(xprofile, machprofile, "g-")
plt.xlabel("distance behind shock [cm]")
plt.ylabel("Mach Number [-]")
plt.subplot(224)
plt.plot(xprofile, sigmaprofile, "m-")
plt.xlabel("distance behind shock [cm]")
plt.ylabel("Total Thermicity [1/sec]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_shock_ZND_reactor.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "ZND.result"
results = {}
results["state-temperature"] = tempprofile.tolist()
results["state-Mach_number"] = machprofile.tolist()
results["species-mixture_thermicity"] = sigmaprofile.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

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
from ansys.chemkin.shock.shocktubereactors import IncidentShock
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# Create a new mechanism input file 'no_hot_air_chem.inp' that
# contains the reactions to describe NO formation in heated air.
# This file is saved to the working directory ``current_dir``.
mymechfile = os.path.join(current_dir, "no_hot_air_chem.inp")
m = open(mymechfile, "w")
# the mechanism contains only the necessary species (oxygen, nitrogen, nitric oxide, and major byproducts)
# decalre elements
m.write("ELEMENT O N AR END\n")
# declare species
m.write("SPECIES\n")
m.write("O2 N2 NO N O AR\n")
m.write("END\n")
# write reactions for N2O dissociation
# Reference:
# M. Camac and R.M. Feinberg, Proceedings of Combustion Institute, vol. 11, p. 137-145 (1967)
m.write("REACTIONS\n")
m.write("N2+O2=NO+NO             9.1E24   -2.5   128500.\n")
m.write("N2+O=NO+N               7.0E13    0.     75000.\n")
m.write("O2+N=NO+O               1.34E10   1.0     7080.\n")
m.write("O2+M=O+O+M              3.62E18  -1.0   118000.\n")
m.write("N2/2/  O2/9/   O/25/\n")
m.write("N2+M=N+N+M              1.92E17  -0.5   224900.\n")
m.write("N2/2.5/  N/0/\n")
m.write("N2+N=N+N+N              4.1E22   -1.5   224900.\n")
m.write("NO+M=N+O+M              4.0E20   -1.5   150000.\n")
m.write("NO/20/  O/20/  N/20/\n")
m.write("END\n")
# close the mechnaism file
m.close()

# Create a chemistry set
# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the N2O dissociation mechanism
MyGasMech = ck.Chemistry(label="NO_from_hot_air")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = mymechfile
MyGasMech.thermfile = os.path.join(
    data_dir,
    "therm.dat",
)

# preprocess the mechanism files
iError = MyGasMech.preprocess()

# Create a diluted air stream by assigning the mole fractions of the
# species.
diluted_air = Stream(MyGasMech)
# initial gas state before the incident shock front
# 5 [torrs]
diluted_air.pressure = 5.0 * ck.Ptorrs
# 296 [K]
diluted_air.temperature = 296.0
# AR diluted air based on the experiment setup
diluted_air.X = [("AR", 0.0093), ("O2", 0.2095), ("N2", 0.7812)]

# Create the shock tube reactor object
# set the incident shock velocity [cm/sec]
diluted_air.velocity = 2.8e5

# instantiate the shock tube reactor
# the location '1' means the gas stream is before the incident shock front
Incident = IncidentShock(diluted_air, location=1, label="incident_shock")

# to use the boundary layer correction, both diameter and viscocity must be given
# shock tube diameter [cm]
Incident.diameter = 3.81
# mixture viscosity [g/cm-sec] at 300 [K]
Incident.set_inlet_viscosity(2.0e-4)

# set total simulation time (particle time) [sec]
Incident.time = 2.0e-3

# tolerances are given in tuple: (absolute tolerance, relative tolerance)
Incident.tolerances = (1.0e-8, 1.0e-4)

# set the start wall time
start_time = time.time()
# run the ZND model
runstatus = Incident.run()
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
Incident.process_solution()

# get the number of solution time points
solutionpoints = Incident.get_solution_size()
print(f"number of solution points = {solutionpoints}")

# get the time profile [sec]
timeprofile = Incident.get_solution_variable_profile("time")
# convert to [msec]
timeprofile *= 1.0e3
# get the temperature profile [K]
tempprofile = Incident.get_solution_variable_profile("temperature")
# get the velocity profile [cm/sec]
velprofile = Incident.get_solution_variable_profile("velocity")
# convert to [m/sec]
velprofile *= 1.0e-2
# get the NO mass fraction profile [-]
NOprofile = Incident.get_solution_variable_profile("NO")
# get the O mass fraction profile [-]
Oprofile = Incident.get_solution_variable_profile("O")

# clean up
if os.path.exists(mymechfile):
    os.remove(mymechfile)

# Plot the solution profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("Incident Shock", fontsize=16)
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, velprofile, "b-")
plt.ylabel("Velocity [m/sec]")
plt.subplot(223)
plt.plot(timeprofile, NOprofile, "g-")
plt.xlabel("time [msec]")
plt.ylabel("NO Mass Fraction [-]")
plt.subplot(224)
plt.plot(timeprofile, Oprofile, "m-")
plt.xlabel("time [msec]")
plt.ylabel("O Mass Fraction [-]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_incident_shock.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "incidentshock.result")
results = {}
results["state-temperature"] = tempprofile.tolist()
results["state-velocity"] = velprofile.tolist()
results["species-NO_mass_fraction"] = NOprofile.tolist()
results["species-O_mass_fraction"] = Oprofile.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

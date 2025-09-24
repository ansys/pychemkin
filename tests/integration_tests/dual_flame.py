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

# Chemkin 1-D opposed-flow flame model (steady-state)
from ansys.chemkin.diffusionflames.opposedflowflame import (
    OpposedFlame as Flame,
)
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
# including the full file path is recommended
chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

# preprocess the mechanism files
iError = MyGasMech.preprocess()
if iError != 0:
    print("Error: Failed to preprocess the mechanism!")
    print(f"       Error code = {iError}")
    exit()

# create the "fuel" stream
fuel = Stream(MyGasMech, label="FUEL")
# set fuel composition: fuel-rich methane-air mixture (equivalence ratio ~ 1.55)
fuel.X = [("CH4", 0.14001807), ("O2", 0.18066847), ("N2", 0.67931346)]
# system pressure [dynes/cm2]
fuel.pressure = ck.Patm
# fuel temperature [K]
fuel.temperature = 300.0
# fuel inlet velocity [cm/sec]
fuel.velocity = 16.0

# create the oxidizer mixture: air
air = Stream(MyGasMech, label="OXID")
air.X = ck.Air.X()
# oxidizer pressure (same as the fuel stream)
air.pressure = fuel.pressure
# oxidizer temperature (same as the fuel temperature)
air.temperature = fuel.temperature
# oxidizer inlet velocity [cm/sec]
air.velocity = 16.0

# instantiate the opposed-flow flame model
dual_flame = Flame(fuel, label="opposed_flame")

# add the "oxidizer" inlet
dual_flame.set_oxidizer_inlet(air)

# define the gap between the two opposing inlets (calculation domain) [cm]
dual_flame.end_position = 1.5
# set up of the "flame zone" to establish the guessed species profiles
# flame zone center location [cm]
# dual_flame.set_reaction_zone_center(0.75)
# flame zone width [cm]
# dual_flame.set_reaction_zone_width(0.5)
# set the estimated maximum gas temperature [K]
dual_flame.set_max_flame_temperature(2200.0)

# set the initial mesh to 26 uniformly distributed grid points
dual_flame.set_numb_grid_points(26)
# set the maximum total number of grid points allowed in the calculation (optional)
dual_flame.set_max_grid_points(250)
# maximum number of grid points can be added during each grid adaption event (optional)
dual_flame.set_max_adaptive_points(5)
# set the maximum values of the grdient and the curvature of the solution profiles (optional)
dual_flame.set_solution_quality(gradient=0.1, curvature=0.3)

# use the mixture averaged formulism to evaluate the mixture transport properties
dual_flame.use_mixture_averaged_transport()
# do NOT include the thermal diffusion effect
dual_flame.use_thermal_diffusion(mode=False)

# specific the species composition boundary treatment ('comp' or 'flux')
# use 'flux' to keep the net species mass fluxes the same as given by the "inlet streams".
dual_flame.set_species_boundary_types(mode="flux")

# reset the tolerances in the steady-state solver (optional)
dual_flame.steady_state_tolerances = (1.0e-9, 1.0e-5)
dual_flame.time_stepping_tolerances = (1.0e-6, 1.0e-4)

# set the start wall time
start_time = time.time()

status = dual_flame.run()
if status != 0:
    print(Color.RED + "Failed to get a converged solution of the opposed flow flame!" + Color.END)
    exit()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"Total simulation duration: {runtime} [sec].")
print()

# postprocess the solutions
dual_flame.process_solution()

# get the number of solution grid points
solutionpoints = dual_flame.get_solution_size()
print(f"Number of solution points = {solutionpoints}.")
# get the grid profile
mesh = dual_flame.get_solution_variable_profile("distance")
# get the temperature profile
tempprofile = dual_flame.get_solution_variable_profile("temperature")
# get the axial velocity profile
velprofile = dual_flame.get_solution_variable_profile("axial_velocity")
# get the mixture fraction profile
mfprofile = dual_flame.get_solution_variable_profile("mixture_fraction")
# get NO2 mass fraction profile
NO2profile = dual_flame.get_solution_variable_profile("NO2")

# plot the opposed-flow flame solution profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(mesh, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(mesh, velprofile, "b-")
plt.ylabel("Axial Velocity [cm/sec]")
plt.subplot(223)
plt.plot(mesh, NO2profile, "g-")
plt.xlabel("Distance [cm]")
plt.ylabel("NO2 Mass Fraction")
plt.subplot(224)
plt.plot(mesh, mfprofile, "m-")
plt.xlabel("Distance [cm]")
plt.ylabel("Mixture Fraction [-]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_opposed_flow_flame.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "dualflame.result")
results = {}
results["state-temperature"] = tempprofile.tolist()
results["state-velocity"] = velprofile.tolist()
results["species-mixture_fraction"] = mfprofile.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

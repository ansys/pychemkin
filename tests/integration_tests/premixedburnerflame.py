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

# Chemkin 1-D premixed freely propagating flame model (steady-state)
from ansys.chemkin.premixedflames.premixedflame import (
    BurnedStabilized_GivenTemperature as Burner,
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
    print("Error: failed to preprocess the mechanism!")
    print(f"       error code = {iError}")
    exit()

# create the fuel-air Stream for the premixed flame speed calculation
premixed = Stream(MyGasMech, label="premixed")
# set inlet premixed stream molar composition
premixed.X = [("C2H4", 0.163), ("O2", 0.237), ("AR", 0.6)]
# setting inlet pressure [dynes/cm2]
premixed.pressure = 1.0 * ck.Patm
# set inlet/unburnt gas temperature [K]
premixed.temperature = 300.0  # inlet temperature

# set the burner outlet cross sectional flow area to unity (1 [cm2])
# because the flame models use the "mass flux" in the calculation instead of the mass flow rate
premixed.flowarea = 1.0
# set value for the inlet mass flow rate [g/sec]
premixed.velocity = 0.0147

# Instantiate the laminar flat flame burner
flatflame = Burner(premixed, label="ethylene flame")

# Set up the temperature profile along the burner centerline
TPRO_data_points = 21
gridpoints = np.zeros(TPRO_data_points, dtype=np.double)
temp_data = np.zeros_like(gridpoints, dtype=np.double)
gridpoints = [
    0.0,
    0.01,
    0.02,
    0.0345643,
    0.0602659,
    0.100148,
    0.118759,
    0.135598,
    0.15421,
    0.179911,
    0.211817,
    0.233087,
    0.260561,
    0.293353,
    0.348301,
    0.436928,
    0.517578,
    0.7161,
    0.935894,
    1.16632,
    1.2,
]
temp_data = [
    300.0,
    450.0,
    600.0,
    828.498,
    1104.4,
    1496.7,
    1634.65,
    1714.85,
    1768.33,
    1801.12,
    1807.2,
    1803.78,
    1799.51,
    1788.35,
    1778.09,
    1767.01,
    1760.23,
    1744.13,
    1737.55,
    1730.99,
    1729.31,
]
flatflame.set_temperature_profile(x=gridpoints, temp=temp_data)

# set the maximum total number of grid points allowed in the calculation (optional)
flatflame.set_max_grid_points(250)
# define the calculation domain [cm]
flatflame.end_position = 1.2
# maximum number of grid points can be added during each grid adaption event (optional)
flatflame.set_max_adaptive_points(10)
# set the maximum values of the grdient and the curvature of the solution profiles (optional)
flatflame.set_solution_quality(gradient=0.1, curvature=0.1)

# use the mixture average formulism to evaluate the mixture transport properties
flatflame.use_mixture_averaged_transport()

# specific the species composition boundary treatment ("comp" or "flux")
# use "FLUX" to ensure that the "net" species mass fluxes are zero at the burner outlet.
flatflame.set_species_boundary_types(mode="flux")

# reset the tolerances in the steady-state solver (optional)
flatflame.steady_state_tolerances = (1.0e-9, 1.0e-4)
# reset the gas species floor value in the steady-state solver (optional)
flatflame.set_species_floor(-1.0e-3)

# set the start wall time
start_time = time.time()

status = flatflame.run()
if status != 0:
    print(Color.RED + "failed to solve the reactor network!" + Color.END)
    exit()

# get the number of solution grid points
solutionpoints = flatflame.get_solution_size()
print(f"number of solution points = {solutionpoints}")

# tightening the converge criteria
flatflame.set_solution_quality(gradient=0.05, curvature=0.1)

# re-start the flame simulation by continuation
flatflame.continuation()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"total simulation duration: {runtime} [sec]")
print()

# postprocess the solutions
flatflame.process_solution()

# get the number of solution grid points
solutionpoints = flatflame.get_solution_size()
print(f"number of final solution points = {solutionpoints}")
# get the grid profile
mesh = flatflame.get_solution_variable_profile("distance")
# get the temperature profile
tempprofile = flatflame.get_solution_variable_profile("temperature")
# get OH mass fraction profile
OHprofile = flatflame.get_solution_variable_profile("OH")

# create arrays for mixture conductivity, and mixture specific heat capacity
Cpprofile = np.zeros_like(mesh, dtype=np.double)
condprofile = np.zeros_like(mesh, dtype=np.double)
# loop over all solution grid points
for i in range(solutionpoints):
    # get the stream at the grid point
    solutionstream = flatflame.get_solution_stream_at_grid(grid_index=i)
    # get mixture specific heat capacity profile [erg/mole-K]
    Cpprofile[i] = solutionstream.CPBL() / ck.ergs_per_joule * 1.0e-3
    # get thermal conductivity profile [ergs/cm-K-sec]
    condprofile[i] = solutionstream.mixture_conductivity() * 1.0e-5

# Plot the premixed flame solution profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(mesh, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(mesh, Cpprofile, "b-")
plt.ylabel("Mixture Cp [kJ/mole]")
plt.subplot(223)
plt.plot(mesh, OHprofile, "g-")
plt.xlabel("Distance [cm]")
plt.ylabel("OH Mass Fraction")
plt.subplot(224)
plt.plot(mesh, condprofile, "m-")
plt.xlabel("Distance [cm]")
plt.ylabel("Mixture conductivity [W/m-K]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_premixed_burner_stabilised_flame.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "burnerstabilized.result")
results = {}
results["state-temperature"] = tempprofile.tolist()
results["state-conductivity"] = condprofile.tolist()
results["species-OH_mass_fraction"] = OHprofile.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

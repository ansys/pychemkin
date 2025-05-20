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
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
# transport data not needed
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# create the fuel mixture
fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.X = [("CH4", 0.8), ("H2", 0.2)]
fuel.temperature = 300.0
fuel.pressure = ck.Patm  # 1 atm
# create the air mixture
air = ck.Mixture(MyGasMech)
# set mass fraction
air.Y = [("O2", 0.23), ("N2", 0.77)]
air.temperature = 300.0
air.pressure = ck.Patm  # 1 atm
# mix the fuel and the air with an air-fuel ratio of 17.19 (almost stoichiometric?)
mixture_recipe = [(fuel, 1.0), (air, 17.19)]
# create the new mixture (the air-fuel ratio is by mass)
premixed = ck.isothermal_mixing(
    recipe=mixture_recipe, mode="mass", finaltemperature=300.0
)
# find the equilibrium composition at different temperature
# and create a NO mole fraction versus temperature plot
# NO species index
NO_index = MyGasMech.get_specindex("NO")
# set up plotting temperatures
Temp = 500.0
dTemp = 20.0
points = 100
# curve
T = np.zeros(points, dtype=np.double)
NO = np.zeros_like(T, dtype=np.double)
# start the temperature loop
for k in range(points):
    # reset mixture temperature
    premixed.temperature = Temp
    # find the equilibrium state mixture at the given mixture temperature and pressure
    eqstate = ck.equilibrium(premixed, opt=1)
    #
    NO[k] = eqstate.X[NO_index] * 1.0e6  # convert to ppm
    T[k] = Temp
    Temp += dTemp
# create plot
plt.plot(T, NO, "bs--", markersize=3, markevery=4)
plt.xlabel("Temperature [K]")
plt.ylabel("NO [ppm]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("equilibrium_composition.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "equilibriumcomposition.result")
results = {}
results["state-temperature"] = T.tolist()
results["species-NO_mole_fraction"] = NO.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

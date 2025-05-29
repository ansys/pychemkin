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

# This is a pychemkin equivalent of equil_test07

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the diesel 14 components mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")

iError = MyGasMech.preprocess()

oxid = ck.Mixture(MyGasMech)
# set mass fraction
oxid.X = [("O2", 1.0)]
oxid.temperature = 295.15
oxid.pressure = ck.Patm  # 1 atm
# mix the fuel and the air with an air-fuel ratio of 17.19 (almost stoichiometric?

fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.X = [("CH4", 1.0)]
fuel.temperature = oxid.temperature
fuel.pressure = oxid.pressure

mixture = ck.Mixture(MyGasMech)
mixture.pressure = oxid.pressure
mixture.temperature = oxid.temperature
products = ["CO2", "H2O"]

points = 12
deq = 0.1
equiv_ini = 0.5

T = np.zeros(points, dtype=np.double)
equiv = np.zeros_like(T, dtype=np.double)

add_frac = np.zeros(MyGasMech.KK, dtype=np.double)
iError = 0
for i in range(points):
    equiv_current = equiv_ini
    iError = mixture.X_by_Equivalence_Ratio(
        MyGasMech, fuel.X, oxid.X, add_frac, products, equivalenceratio=equiv_current
    )
    if iError != 0:
        raise RuntimeError
    result = ck.equilibrium(mixture, opt=5)
    T[i] = result.temperature
    equiv[i] = equiv_current
    equiv_ini = equiv_ini + deq

plt.plot(equiv, T, "bs--")
plt.xlabel("Equivalence ratio")
plt.ylabel("Temperature [K]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("adiabatic_flame_temperature.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "adiabaticflametemperature.result")
results = {}
results["state-equivalence_ratio"] = equiv.tolist()
results["state-temperature"] = T.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

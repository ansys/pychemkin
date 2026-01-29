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

"""Test for calculating gas species properties."""

from pathlib import Path

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.logger import logger

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
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
MyGasMech.tranfile = str(mechanism_dir / "grimech30_transport.dat")
# preprocess the mechanism files
ierror = MyGasMech.preprocess()
# extract element symbols as a list
elelist = MyGasMech.element_symbols
# extract gas species symbols as a list
specieslist = MyGasMech.species_symbols
# list of gas species interested
plotspeclist = ["CH4", "O2", "N2"]
# find elemental compositions of selected species
print(" ")
for s in plotspeclist:
    species_id = MyGasMech.get_specindex(s)
    print("species " + specieslist[species_id])
    print("elemental composition")
    for elem_id in range(MyGasMech.mm):
        num_elem = MyGasMech.species_composition(elem_id, species_id)
        print(f"    {elelist[elem_id]:>4}: {num_elem:2d}")
    print("=" * 10)
print()
#
# plot Cv value at different temperatures for selected gas species
#
plt.figure(figsize=(12, 6))
# temperature increment
dtemp = 20.0
# number of property data points
points = 100
# curve attributes
curvelist = ["g", "b--", "r:"]
# create arrays
# species specific heat capacity at constant volume data
cv = np.zeros(points, dtype=np.double)
# temperature data
t = np.zeros(points, dtype=np.double)
# start of the plotting loop #1
k = 0
# loop over the selected gas species
for s in plotspeclist:
    # starting temperature at 300K
    temp = 300.0
    # loop over temperature data points
    for i in range(points):
        heatcapacity = MyGasMech.species_cv(temp)
        id = MyGasMech.get_specindex(s)
        t[i] = temp
        # convert ergs to joules
        cv[i] = heatcapacity[id] / ck.ERGS_PER_JOULE
        temp += dtemp
    plt.subplot(121)
    plt.plot(t, cv, curvelist[k])
    k += 1
# plot Cv versus temperature
plt.xlabel("Temperature [K]")
plt.ylabel("Cv [J/mol-K]")
plt.legend(plotspeclist, loc="upper left")
# create arrays
# species conductivity
kappa = np.zeros(points, dtype=np.double)
# start of the plotting loop #2
k = 0
# loop over the selected gas species
for s in plotspeclist:
    # starting temperature at 300K
    temp = 300.0
    # loop over temperature data points
    for i in range(points):
        conductivity = MyGasMech.species_cond(temp)
        id = MyGasMech.get_specindex(s)
        t[i] = temp
        # convert ergs to joules
        kappa[i] = conductivity[id] / ck.ERGS_PER_JOULE
        temp += dtemp
    plt.subplot(122)
    plt.plot(t, kappa, curvelist[k])
    k += 1
# plot conductivity versus temperature
plt.xlabel("Temperature [K]")
plt.ylabel("Conductivity [J/cm-K-sec]")
plt.legend(plotspeclist, loc="upper left")
# calculate species binary diffusion coefficients
# at 2 atm and 500K
diffcoef = MyGasMech.species_diffusioncoeffs(2.0 * ck.P_ATM, 500.0)
id1 = MyGasMech.get_specindex(plotspeclist[0])
id2 = MyGasMech.get_specindex(plotspeclist[1])
c = diffcoef[id1][id2]
print(
    f"diffusion coefficient for {plotspeclist[0]}"
    f" against {plotspeclist[1]} is {c:e} [cm2/sec]"
)
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("species_properties.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "speciesproperties.result"
results = {}
results["state-temperature"] = t.tolist()
results["state-Cv"] = cv.tolist()
results["state-conductivity"] = kappa.tolist()
results["state-binary_diffusivity"] = [float(c)]
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

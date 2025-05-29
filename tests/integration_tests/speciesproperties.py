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
MyGasMech.tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# extract element symbols as a list
elelist = MyGasMech.element_symbols
# extract gas species symbols as a list
specieslist = MyGasMech.species_symbols
# list of gas species interested
plotspeclist = ["CH4", "O2", "N2"]
# find elemental compositions of selected species
print(" ")
for s in plotspeclist:
    speciesID = MyGasMech.get_specindex(s)
    print("species " + specieslist[speciesID])
    print("elemental composition")
    for elemID in range(MyGasMech.MM):
        num_elem = MyGasMech.SpeciesComposition(elemID, speciesID)
        print(f"    {elelist[elemID]:>4}: {num_elem:2d}")
    print("=" * 10)
print()
#
# plot Cv value at different temperatures for selected gas species
#
plt.figure(figsize=(12, 6))
# temperature increment
dTemp = 20.0
# number of property data points
points = 100
# curve attributes
curvelist = ["g", "b--", "r:"]
# create arrays
# species specific heat capacity at constant volume data
Cv = np.zeros(points, dtype=np.double)
# temperature data
T = np.zeros(points, dtype=np.double)
# start of the plotting loop #1
k = 0
# loop over the selected gas species
for s in plotspeclist:
    # starting temperature at 300K
    Temp = 300.0
    # loop over temperature data points
    for i in range(points):
        HeatCapacity = MyGasMech.SpeciesCv(Temp)
        ID = MyGasMech.get_specindex(s)
        T[i] = Temp
        # convert ergs to joules
        Cv[i] = HeatCapacity[ID] / ck.ergs_per_joule
        Temp += dTemp
    plt.subplot(121)
    plt.plot(T, Cv, curvelist[k])
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
    Temp = 300.0
    # loop over temperature data points
    for i in range(points):
        conductivity = MyGasMech.SpeciesCond(Temp)
        ID = MyGasMech.get_specindex(s)
        T[i] = Temp
        # convert ergs to joules
        kappa[i] = conductivity[ID] / ck.ergs_per_joule
        Temp += dTemp
    plt.subplot(122)
    plt.plot(T, kappa, curvelist[k])
    k += 1
# plot conductivity versus temperature
plt.xlabel("Temperature [K]")
plt.ylabel("Conductivity [J/cm-K-sec]")
plt.legend(plotspeclist, loc="upper left")
# calculate species binary diffusion coefficients
# at 2 atm and 500K
diffcoef = MyGasMech.SpeciesDiffusionCoeffs(2.0 * ck.Patm, 500.0)
ID1 = MyGasMech.get_specindex(plotspeclist[0])
ID2 = MyGasMech.get_specindex(plotspeclist[1])
c = diffcoef[ID1][ID2]
print(
    f"diffusion coefficient for {plotspeclist[0]} against {plotspeclist[1]} is {c:e} [cm2/sec]"
)
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("species_properties.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "speciesproperties.result")
results = {}
results["state-temperature"] = T.tolist()
results["state-Cv"] = Cv.tolist()
results["state-conductivity"] = kappa.tolist()
results["state-binary_diffusivity"] = [float(c)]
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

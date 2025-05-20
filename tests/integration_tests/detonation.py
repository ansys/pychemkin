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
# create a chemistry set based on C2_NOx using an alternative method
MyMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files individually
# this mechanism file contains all the necessary thermodynamic and transport data
# therefore no need to specify the therm and the tran data files
MyMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")
# preprocess the 2nd mechanism files
iError = MyMech.preprocess()
# create the fuel mixture
fuel = ck.Mixture(MyMech)
# set mole fraction
fuel.X = [("CH4", 0.8), ("C2H6", 0.2)]
fuel.temperature = 290.0
fuel.pressure = 40.0 * ck.Patm
# create the air mixture
air = ck.Mixture(MyMech)
# set mass fraction
air.X = [("O2", 0.21), ("N2", 0.79)]
air.temperature = fuel.temperature
air.pressure = fuel.pressure
# create the initial mixture
# create the premixed mixture to be defined by equivalence ratio
premixed = ck.Mixture(MyMech)
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
add_frac = np.zeros(MyMech.KK, dtype=np.double)  # no additives: all zeros
iError = premixed.X_by_Equivalence_Ratio(
    MyMech, fuel.X, air.X, add_frac, products, equivalenceratio=1.0
)
if iError != 0:
    raise RuntimeError
# list the composition of the premixed mixture
premixed.list_composition(mode="mole")
# set up parameter study of detonation wave speed with respect to pressure
points = 5
dpres = 10.0 * ck.Patm
pres = fuel.pressure
P = np.zeros(points, dtype=np.double)
Det = np.zeros_like(P, dtype=np.double)
premixed.pressure = pres
premixed.temperature = fuel.temperature
# start of pressure loop
for i in range(points):
    # compute the C-J state corresponding to the initial mixture
    speed, CJstate = ck.detonation(premixed)
    # update plot data
    # convert pressure to atm
    P[i] = pres / ck.Patm
    # convert speed to m/sec
    Det[i] = speed[1] / 1.0e2
    # update pressure value
    pres += dpres
    premixed.pressure = pres
# create plot for ideal gas results
plt.plot(P, Det, "bo--", label="ideal gas", markersize=5, fillstyle="none")
#
# turn on real-gas cubic equation of state
premixed.use_realgas_cubicEOS()
# set mixture mixing rule to Van der Waals (default)
# premixed.set_realgas_mixing_rule(rule=0)
# restart the calculation with real-gas EOS
premixed.pressure = fuel.pressure
pres = fuel.pressure
P[:] = 0.0e0
Det[:] = 0.0e0
# set verbose mode to false to turn OFF extra printouts
ck.set_verbose(False)
# start of pressure loop
for i in range(points):
    # compute the C-J state corresponding to the initial mixture
    speed, CJstate = ck.detonation(premixed)
    # update plot data
    P[i] = pres / ck.Patm
    Det[i] = speed[1] / 1.0e2
    # update pressure value
    pres += dpres
    premixed.pressure = pres
# stop Chemkin
ck.done()
# create plot for real gas results
plt.plot(P, Det, "r^-", label="real gas", markersize=5, fillstyle="none")
# plot data
P_data = [44.1, 50.6, 67.2, 80.8]
Det_data = [1950.0, 1970.0, 2000.0, 2020.0]
plt.plot(P_data, Det_data, "gD:", label="data", markersize=4)
#
plt.legend(loc="upper left")
plt.xlabel("Pressure [atm]")
plt.ylabel("Detonation wave speed [m/sec]")
plt.suptitle("Natural Gas/Air Detonation", fontsize=16)
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("detonation.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "detonation.result")
results = {}
results["state-pressure"] = P.tolist()
results["state-detonation_speed_RealGas"] = Det.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

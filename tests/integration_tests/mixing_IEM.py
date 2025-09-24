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
from ansys.chemkin.mixture import Mixture, mixing_by_exchange_with_the_mean
from ansys.chemkin.logger import logger

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

# get species sybols
specieslist = MyGasMech.species_symbols

# create mixture A
mixtureA = Mixture(MyGasMech)
mixtureA.X = ck.Air.X()
mixtureA.temperature = 900.0
mixtureA.pressure = ck.Patm
mixtureA.volume = 1.0

# create mixture B
mixtureB = Mixture(MyGasMech)
mixtureB.X = [("H2", 0.2), ("CO2", 0.4), ("CH4", 0.4)]
mixtureB.temperature = 400.0
mixtureB.pressure = ck.Patm
mixtureB.volume = 1.0

# set IEM mixing parameters
# mixing duration [sec]
dt = 1.0e-2
# IEM model parameter
Cmix = 1.0e0
# characteristic mixing time scale [sec]
tau = 1.0e-2
mixtureA_new, mixtureB_new = mixing_by_exchange_with_the_mean(
    mixtureA, mixtureB, mix_time=dt, mix_param=Cmix, tau=tau
)
# compare mixture properties before and after the mixing
# number of gas species
numb_spec = MyGasMech.KK

print("\n*Mixture A")
print("==Before==")
print(f"Temperature = {mixtureA.temperature} [K].")
for k in range(numb_spec):
    x = mixtureA.X[k]
    if x > 1.0e-8:
        print(f"Species index {k + 1}  {specieslist[k]}: {x}")
print("==After==")
print(f"Temperature = {mixtureA_new.temperature} [K].")
for k in range(numb_spec):
    x = mixtureA_new.X[k]
    if x > 1.0e-8:
        print(f"Species index {k + 1}  {specieslist[k]}: {x}")

print("\n*Mixture B")
print("==Before==")
print(f"Temperature = {mixtureB.temperature} [K].")
for k in range(numb_spec):
    x = mixtureB.X[k]
    if x > 1.0e-8:
        print(f"Species index {k + 1}  {specieslist[k]}: {x}")
print("==After==")
print(f"Temperature = {mixtureB_new.temperature} [K].")
for k in range(numb_spec):
    x = mixtureB_new.X[k]
    if x > 1.0e-8:
        print(f"Species index {k + 1}  {specieslist[k]}: {x}")

# return results for comparisons
resultfile = os.path.join(current_dir, "mixing_IEM.result")
results = {}
results["state-model_parameters"] = [dt, Cmix, tau]
results["state-temperatureA"] = [mixtureA.temperature, mixtureA_new.temperature]
results["state-mole_fractionA"] = mixtureA.X.tolist()
results["state-mole_fractionA_new"] = mixtureA_new.X.tolist()
results["state-temperatureB"] = [mixtureB.temperature, mixtureB_new.temperature]
results["state-mole_fractionB"] = mixtureB.X.tolist()
results["state-mole_fractionB_new"] = mixtureB_new.X.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

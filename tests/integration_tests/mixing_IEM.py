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

"""IEM micro mixing test for the stochastic process utilities."""

from pathlib import Path

import numpy as np

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.mixture import Mixture, mixing_by_exchange_with_the_mean

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

# set numpy printing option
np.set_printoptions(legacy="1.25")

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# including the full file path is recommended
chemfile = str(mechanism_dir / "grimech30_chem.inp")
thermfile = str(mechanism_dir / "grimech30_thermo.dat")
tranfile = str(mechanism_dir / "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

# preprocess the mechanism files
ierror = MyGasMech.preprocess()
if ierror != 0:
    print("Error: failed to preprocess the mechanism!")
    print(f"       error code = {ierror}")
    exit()

# get species symbols
specieslist = MyGasMech.species_symbols

# create mixture A
mixture_a = Mixture(MyGasMech)
mixture_a.x = ck.Air.x()
mixture_a.temperature = 900.0
mixture_a.pressure = ck.P_ATM
mixture_a.volume = 1.0

# create mixture B
mixture_b = Mixture(MyGasMech)
mixture_b.x = [("H2", 0.2), ("CO2", 0.4), ("CH4", 0.4)]
mixture_b.temperature = 400.0
mixture_b.pressure = ck.P_ATM
mixture_b.volume = 1.0

# set IEM mixing parameters
# mixing duration [sec]
dt = 1.0e-2
# IEM model parameter
c_mix = 1.0e0
# characteristic mixing time scale [sec]
tau = 1.0e-2
mixture_a_new, mixture_b_new = mixing_by_exchange_with_the_mean(
    mixture_a, mixture_b, mix_time=dt, mix_param=c_mix, tau=tau
)
# compare mixture properties before and after the mixing
# number of gas species
numb_spec = MyGasMech.kk

print("\n*Mixture A")
print("==Before==")
print(f"Temperature = {mixture_a.temperature} [K].")
for k in range(numb_spec):
    x_sp = mixture_a.x[k]
    if x_sp > 1.0e-8:
        print(f"Species index {k + 1}  {specieslist[k]}: {x_sp}")
print("==After==")
print(f"Temperature = {mixture_a_new.temperature} [K].")
for k in range(numb_spec):
    x_sp = mixture_a_new.x[k]
    if x_sp > 1.0e-8:
        print(f"Species index {k + 1}  {specieslist[k]}: {x_sp}")

print("\n*Mixture B")
print("==Before==")
print(f"Temperature = {mixture_b.temperature} [K].")
for k in range(numb_spec):
    x_sp = mixture_b.x[k]
    if x_sp > 1.0e-8:
        print(f"Species index {k + 1}  {specieslist[k]}: {x_sp}")
print("==After==")
print(f"Temperature = {mixture_b_new.temperature} [K].")
for k in range(numb_spec):
    x_sp = mixture_b_new.x[k]
    if x_sp > 1.0e-8:
        print(f"Species index {k + 1}  {specieslist[k]}: {x_sp}")

# return results for comparisons
resultfile = Path(current_dir) / "mixing_IEM.result"
results = {}
results["state-model_parameters"] = [dt, c_mix, tau]
results["state-temperatureA"] = [mixture_a.temperature, mixture_a_new.temperature]
results["state-mole_fractionA"] = mixture_a.x.tolist()
results["state-mole_fractionA_new"] = mixture_a_new.x.tolist()
results["state-temperatureB"] = [mixture_b.temperature, mixture_b_new.temperature]
results["state-mole_fractionB"] = mixture_b.x.tolist()
results["state-mole_fractionB_new"] = mixture_b_new.x.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

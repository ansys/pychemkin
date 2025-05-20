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
from ansys.chemkin import Color
from ansys.chemkin.logger import logger

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)
# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# set mechanism input files
# including the full file path is recommended
chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
My1stMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

# preprocess the mechanism files
iError = My1stMech.preprocess()
print()
if iError != 0:
    print(f"Preprocessing error encountered. Code = {iError:d}.")
    print(f"see the summary file {My1stMech.summaryfile} for details")
    exit()
else:
    print(Color.GREEN + "Preprocessing succeeded.", end=Color.END)
    print("mechanism information:")
    print(f"number of elements = {My1stMech.MM:d}")
    print(f"number of gas species = {My1stMech.KK:d}")
    print(f"number of gas reactions = {My1stMech.IIGas:d}")

# create a mixture with My1stMech
mymixture1 = ck.Mixture(My1stMech)
# set mixture temperature [K]
mymixture1.temperature = 1000.0
# set mixture pressure [dynes/cm2]
mymixture1.pressure = ck.Patm
# set molar compositions
mymixture1.X = [("CH4", 0.1), ("O2", 0.21), ("N2", 0.79)]
# compute the constrained H-P equilibrium state
ck.help("equilibrium")
equil_mix1_HP = ck.equilibrium(mymixture1, opt=5)
print(f"equilibrium temperature of mymixture1 : {equil_mix1_HP.temperature} [K]")
#
# load the second mechanism
#
# set the 2nd mechanism directory (the default Chemkin mechanism data directory)
mechanism_dir = data_dir
# create a chemistry set based on C2_NOx using an alternative method
My2ndMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files individually
# this mechanism file contains all the necessary thermodynamic and transport data
# therefore no need to specify the therm and the tran data files
My2ndMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")
# instruct the preprocessor to include the transport properties
# only when the mechanism file contains all the transport data
My2ndMech.preprocess_transportdata()
# preprocess the 2nd mechanism files
iError = My2ndMech.preprocess()
print()
if iError != 0:
    print(f"Preprocessing error encountered. Code = {iError:d}.")
    print(f"see the summary file {My2ndMech.summaryfile} for details")
    exit()
else:
    print(Color.GREEN + "Preprocessing succeeded.", end=Color.END)
    print("mechanism information:")
    print(f"number of elements = {My2ndMech.MM:d}")
    print(f"number of gas species = {My2ndMech.KK:d}")
    print(f"number of gas reactions = {My2ndMech.IIGas:d}")

# create the 2nd mixture with the My2ndMech
mymixture2 = ck.Mixture(My2ndMech)
# set mixture temperature [K]
mymixture2.temperature = 500.0
# set mixture pressure [dynes/cm2]
mymixture2.pressure = 2.0 * ck.Patm
# set mixture molar composition
mymixture2.X = [("H2", 0.02), ("O2", 0.2), ("N2", 0.8)]
# compute detonation wave speed with mymixture2
speeds_mix2, CJ_mix2 = ck.detonation(mymixture2)
print(f"detonation mymixture2 temperature: {CJ_mix2.temperature} [K]")
print(f"detonation wave speed = {speeds_mix2[1] / 100.0} [m/sec]")
#
# re-activate My1stMech
My1stMech.activate()
# compute detonation wave speed with mymixture1
speeds_mix1, CJ_mix1 = ck.detonation(mymixture1)
print(f"detonation mymixture1 temperature: {CJ_mix1.temperature} [K]")
print(f"detonation wave speed = {speeds_mix1[1] / 100.0} [m/sec]")

# return results for comparisons
resultfile = os.path.join(current_dir, "multiplemechanisms.result")
results = {}
# results["state-temperature_IdealGas"] = [CJ_mix1.temperature]
# results["state-detonation_speed_IdealGas"] = [speeds_mix1[1] / 100.0]
results["state-temperature_RealGas"] = [CJ_mix2.temperature]
results["state-detonation_speed_RealGas"] = [speeds_mix2[1] / 100.0]
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

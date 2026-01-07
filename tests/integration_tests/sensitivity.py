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

"""Test for the brute force sensitivity analysis."""

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.core.batchreactors.batchreactor import (
    GivenPressureBatchReactorEnergyConservation,
)
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(False)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on the diesel 14 components mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")
# pre-process
ierror = MyGasMech.preprocess()
if ierror == 0:
    print("mechanism information:")
    print(f"number of gas species = {MyGasMech.kk:d}")
    print(f"number of gas reactions = {MyGasMech.ii_gas:d}")
else:
    exit()
#
# create air-fuel mixture with equivalence ratio = 1.1
oxid = ck.Mixture(MyGasMech)
# set mole fraction
oxid.x = [("O2", 1.0), ("N2", 3.76)]
oxid.temperature = 900
oxid.pressure = ck.P_ATM  # 1 atm
fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.x = [("C3H8", 0.1), ("CH4", 0.8), ("H2", 0.1)]
fuel.temperature = oxid.temperature
fuel.pressure = oxid.pressure
mixture = ck.Mixture(MyGasMech)
mixture.pressure = oxid.pressure
mixture.temperature = oxid.temperature
products = ["CO2", "H2O", "N2"]
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)
# create the air-fuel mixture by using the equivalence ratio method
ierror = mixture.x_by_equivalence_ratio(
    MyGasMech, fuel.x, oxid.x, add_frac, products, equivalenceratio=1.1
)
if ierror != 0:
    raise RuntimeError
if ck.verbose():
    mixture.list_composition(mode="mole")
# get the original rate parameters
Afactor, Beta, ActiveEnergy = MyGasMech.get_reaction_parameters()
if ck.verbose():
    for i in range(MyGasMech.ii_gas):
        print(f"reaction: {i + 1}")
        print(f"A  = {Afactor[i]}")
        print(f"B  = {Beta[i]}")
        print(f"Ea = {ActiveEnergy[i]}\n")
        if np.isclose(0.0, Afactor[i], atol=1.0e-15):
            print("reaction pre-exponential factor = 0")
            exit()
#
# compute the ignition delay time
# create a constant pressure batch reactor (with energy equation)
#
MyCONP = GivenPressureBatchReactorEnergyConservation(mixture, label="CONP")
# show initial gas composition inside the reactor
MyCONP.list_composition(mode="mole")
# set other reactor parameters
# reactor volume [cm3]
MyCONP.volume = 10.0
# simulation end time [sec]
MyCONP.time = 2.0
# turn ON adaptive solution saving
MyCONP.adaptive_solution_saving(mode=True, steps=20)
# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyCONP.tolerances = (1.0e-10, 1.0e-8)
# set ignition delay
# ck.show_ignition_definitions()
MyCONP.set_ignition_delay(method="T_inflection")
# set the start wall time
start_time = time.time()
# run the nominal case
runstatus = MyCONP.run()
#
if runstatus == 0:
    # get ignition delay time
    delaytime_org = MyCONP.get_ignition_delay()
    print(f"ignition delay time = {delaytime_org} [msec]")
else:
    # if get this, most likely the END time is too short
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    print("failed to find the ignition delay time of the nominal case")
    exit()
#
# brute force A factor sensitivity coefficients of ignition delay time
# create sensitivity coefficient array
IGsen = np.zeros(MyGasMech.ii_gas, dtype=np.double)
# set perturbation magnitude
perturb = 0.001  # increase by 0.1%
perturbplus1 = 1.0 + perturb
# loop over all reactions
for i in range(MyGasMech.ii_gas):
    Anew = Afactor[i] * perturbplus1
    # actual reaction index
    ireac = i + 1
    # update the A factor
    MyGasMech.set_reaction_afactor(ireac, Anew)
    # run the reactor model
    runstatus = MyCONP.run()
    #
    if runstatus == 0:
        # get ignition delay time
        delaytime = MyCONP.get_ignition_delay()
        print(f"ignition delay time = {delaytime} [msec]")
        # compute d(delaytime)
        IGsen[i] = delaytime - delaytime_org
        # restore the A factor
        MyGasMech.set_reaction_afactor(ireac, Afactor[i])
    else:
        # if get this, most likely the END time is too short
        print(f"trouble finding ignition delay time for raection {ireac}")
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()

# compute the total runtime
runtime = time.time() - start_time
print(f"\ntotal simulation time: {runtime} [sec] over {MyGasMech.ii_gas + 1} runs")
# normalized sensitivity coefficient = d(delaytime) * A[i] / d(A[i])
IGsen /= perturb
# print top n positive and negative ignition delay time sensitivity coefficients
top = 5
# rank the positive coefficients
posindex = np.argpartition(IGsen, -top)[-top:]
poscoeffs = IGsen[posindex]

# rank the negative coefficients
NegIGsen = np.negative(IGsen)
negindex = np.argpartition(NegIGsen, -top)[-top:]
negcoeffs = IGsen[negindex]
# print the top sensitivity coefficients
if ck.verbose():
    print("positive sensitivity coefficients")
    for i in range(top):
        print(f"reaction {posindex[i] + 1}: coefficient = {poscoeffs[i]}")
    print()
    print("negative sensitivity coefficients")
    for i in range(top):
        print(f"reaction {negindex[i] + 1}: coefficient = {negcoeffs[i]}")
#
# create a rate plot
plt.rcParams.update({"figure.autolayout": True, "ytick.color": "blue"})
plt.subplots(2, 1, sharex="col", figsize=(10, 5))
# convert reaction # from integers to strings
rxnstring = []
for i in range(len(posindex)):
    # the array index starting from 0 so the actual reaction # = index + 1
    rxnstring.append(MyGasMech.get_gas_reaction_string(posindex[i] + 1))
# use horizontal bar chart
plt.subplot(211)
plt.barh(rxnstring, poscoeffs, color="orange", height=0.4)
plt.axvline(x=0, c="gray", lw=1)
#
# convert reaction # from integers to strings
rxnstring.clear()
fnegindex = np.flip(negindex)
fnegcoeffs = np.flip(negcoeffs)
for i in range(len(negindex)):
    # the array index starting from 0 so the actual reaction # = index + 1
    rxnstring.append(MyGasMech.get_gas_reaction_string(fnegindex[i] + 1))
plt.subplot(212)
plt.barh(rxnstring, fnegcoeffs, color="orange", height=0.4)
plt.axvline(x=0, c="gray", lw=1)
plt.xlabel("Sensitivity Coefficients")
plt.suptitle("Ignition Delay Time Sensitivity", fontsize=16)
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("sensitivity_analysis.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "sensitivity.result"
results = {}
results["state-index_positive"] = posindex.tolist()
results["rate-sensitivity_positive"] = poscoeffs.tolist()
results["state-index_negative"] = negindex.tolist()
results["rate-sensitivity_negative"] = negcoeffs.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

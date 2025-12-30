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

"""Test for the profile feature of the closed homogeneous reactor model."""

from pathlib import Path

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# chemkin batch reactor model (transient)
from ansys.chemkin.core.batchreactors.batchreactor import (
    GivenPressureBatchReactorFixedTemperature,
)
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on C2_NOx using an alternative method
MyMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files individually
# this mechanism file contains all the necessary thermodynamic and transport data
# therefore no need to specify the therm and the tran data files
MyMech.chemfile = str(mechanism_dir / "C2_NOx_SRK.inp")
# preprocess the 2nd mechanism files
ierror = MyMech.preprocess()
if ierror == 0:
    print(Color.GREEN + ">>> preprocess OK", end=Color.END)
else:
    print(Color.RED + ">>> preprocess failed!", end=Color.END)
    exit()
# create the air+vapor mixture
mist = ck.Mixture(MyMech)
# set mole fraction
mist.x = [("H2O", 2.0), ("O2", 1.0), ("N2", 3.76)]
mist.temperature = 500.0  # [K]
mist.pressure = 100.0 * ck.P_ATM
# set mixture mixing rule to Van der Waals (default)
# mist.set_realgas_mixing_rule(rule=0)
# create a constant pressure batch reactor (with given temperature)
#
tank = GivenPressureBatchReactorFixedTemperature(mist, label="tank")
# show initial gas composition inside the reactor
tank.list_composition(mode="mole")
# set other reactor properties
tank.volume = 10.0  # cm3
tank.time = 0.5  # sec
# turn on real-gas cubic equation of state
tank.userealgas_eos(mode=True)
# output controls
# set timestep between saving solution
tank.timestep_for_saving_solution = 0.01
# set tolerances in tuple: (absolute tolerance, relative tolerance)
tank.tolerances = (1.0e-10, 1.0e-8)
# get solver parameters
atol, rtol = tank.tolerances
print(f"default absolute tolerance = {atol}")
print(f"default relative tolerance = {rtol}")
# turn on the force non-negative solutions option in the solver
tank.force_nonnegative = True
# set tank profile
# number of profile data points
npoints = 3
# position array of the profile data
x = np.zeros(npoints, dtype=np.double)
# value array of the profile data
tpro_profile = np.zeros_like(x, dtype=np.double)
# set tank temperature data points
x = [0.0, 0.2, 2.0]  # [sec]
tpro_profile = [500.0, 275.0, 275.0]  # [K]
# set the temperature profile
tank.set_temperature_profile(x, tpro_profile)
# run the CONP reactor model with given temperature profile
runstatus = tank.run()
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)

# post-process the solutions
tank.process_solution()
# get the number of solution time points
solutionpoints = tank.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile = tank.get_solution_variable_profile("time")
# get the temperature profile
tempprofile = tank.get_solution_variable_profile("temperature")
# get the volume profile
volprofile = tank.get_solution_variable_profile("volume")
# create array for mixture density
denprofile = np.zeros_like(timeprofile, dtype=np.double)
# create array for mixture enthalpy
Hprofile = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tank.get_solution_mixture_at_index(solution_index=i)
    # get mixture density profile
    denprofile[i] = solutionmixture.rho
    # get mixture enthalpy profile
    Hprofile[i] = solutionmixture.hml() / ck.ERGS_PER_JOULE * 1.0e-3

#
# turn off real-gas cubic equation of state
tank.userealgas_eos(mode=False)
# run the CONP reactor model with given temperature profile
runstatus = tank.run()
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
# post-process the solutions
tank.process_solution()
# get the number of solution time points
solutionpoints = tank.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile_idealgas = tank.get_solution_variable_profile("time")
# get the volume profile
volprofile_idealgas = tank.get_solution_variable_profile("volume")
# create array for mixture density
denprofile_idealgas = np.zeros_like(timeprofile, dtype=np.double)
# create array for mixture enthalpy
Hprofile_idealgas = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tank.get_solution_mixture_at_index(solution_index=i)
    # get mixture density profile
    denprofile_idealgas[i] = solutionmixture.rho
    # get mixture enthalpy profile
    Hprofile_idealgas[i] = solutionmixture.hml() / ck.ERGS_PER_JOULE * 1.0e-3

ck.done()
# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
thispres = str(mist.pressure / ck.P_ATM)
thistitle = "Cooling Vapor + Air at " + thispres + " atm"
plt.suptitle(thistitle, fontsize=16)
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, volprofile, "b-", label="real gas")
plt.plot(timeprofile_idealgas, volprofile_idealgas, "b--", label="ideal gas")
plt.legend(loc="upper right")
plt.ylabel("Volume [cm3]")
plt.subplot(223)
plt.plot(timeprofile, Hprofile, "g-", label="real gas")
plt.plot(timeprofile_idealgas, Hprofile_idealgas, "g--", label="ideal gas")
plt.legend(loc="upper right")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Enthalpy [kJ/mole]")
plt.subplot(224)
plt.plot(timeprofile, denprofile, "m-", label="real gas")
plt.plot(timeprofile_idealgas, denprofile_idealgas, "m--", label="ideal gas")
plt.legend(loc="upper left")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Density [g/cm3]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("vapor_condensation.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "vapor.result"
results = {}
results["state-time"] = timeprofile.tolist()
results["state-temperature"] = tempprofile.tolist()
results["state-volume_RealGas"] = volprofile.tolist()
results["state-volume_IdealGas"] = volprofile_idealgas.tolist()
results["state-enthalpy_RealGas"] = Hprofile.tolist()
results["state-enthalpy_IdealGas"] = Hprofile_idealgas.tolist()
results["state-density_RealGas"] = denprofile.tolist()
results["state-density_IdealGas"] = denprofile_idealgas.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

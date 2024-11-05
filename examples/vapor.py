import os

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import chemkin as ck  # Chemkin
from chemkin import Color

# chemkin batch reactor model (transient)
from chemkin.batchreactors.batchreactor import (
    GivenPressureBatchReactor_FixedTemperature,
)

# check working directory
current_dir = os.getcwd()
print("current working directory: " + current_dir)
# set mechanism directory (the default chemkin mechanism data directory)
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
if iError == 0:
    print(Color.GREEN + ">>> preprocess OK", end=Color.END)
else:
    print(Color.RED + ">>> preprocess failed!", end=Color.END)
    exit()
# create the air+vapor mixture
mist = ck.Mixture(MyMech)
# set mole fraction
mist.X = [("H2O", 2.0), ("O2", 1.0), ("N2", 3.76)]
mist.temperature = 500.0  # [K]
mist.pressure = 100.0 * ck.Patm
# set mixture mixing rule to Van der Waals (default)
# mist.setrealgasmixingrule(rule=0)
# create a constant pressure batch reactor (with given temperature)
#
tank = GivenPressureBatchReactor_FixedTemperature(mist, label="tank")
# show initial gas composition inside the reactor
tank.listcomposition(mode="mole")
# set other reactor properties
tank.volume = 10.0  # cm3
tank.time = 0.5  # sec
# turn on real-gas cubic equation of state
tank.userealgasEOS(mode=True)
# output controls
# set timestep between saving solution
tank.timestepforsavingsolution = 0.01
# set tolerance
tank.settolerances(absolute_tolerance=1.0e-10, relative_tolerance=1.0e-8)
# get solver parameters
ATOL, RTOL = tank.tolerances
print(f"default absolute tolerance = {ATOL}")
print(f"default relative tolerance = {RTOL}")
# turn on the force non-negative solutions option in the solver
tank.forcenonnegative = True
# set tank profile
# number of profile data points
npoints = 3
# position array of the profile data
x = np.zeros(npoints, dtype=np.double)
# value array of the profile data
TPROprofile = np.zeros_like(x, dtype=np.double)
# set tank temperature data points
x = [0.0, 0.2, 2.0]  # [sec]
TPROprofile = [500.0, 275.0, 275.0]  # [K]
# set the temperature profile
tank.settemperatureprofile(x, TPROprofile)
# run the CONP reactor model with given temperature profile
runstatus = tank.run()
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> RUN COMPLETED <<<", end=Color.END)

# post-process the solutions
tank.processsolution()
# get the number of solution time points
solutionpoints = tank.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile = tank.getsolutionvariableprofile("time")
# get the temperature profile
tempprofile = tank.getsolutionvariableprofile("temperature")
# get the volume profile
volprofile = tank.getsolutionvariableprofile("volume")
# create array for mixture density
denprofile = np.zeros_like(timeprofile, dtype=np.double)
# create array for mixture enthalpy
Hprofile = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tank.getsolutionmixtureatindex(solution_index=i)
    # get mixture density profile
    denprofile[i] = solutionmixture.RHO
    # get mixture enthalpy profile
    Hprofile[i] = solutionmixture.HML() / ck.ergsperjoule * 1.0e-3

#
# turn off real-gas cubic equation of state
tank.userealgasEOS(mode=False)
# run the CONP reactor model with given temperature profile
runstatus = tank.run()
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> RUN COMPLETED <<<", end=Color.END)
# post-process the solutions
tank.processsolution()
# get the number of solution time points
solutionpoints = tank.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile_IG = tank.getsolutionvariableprofile("time")
# get the volume profile
volprofile_IG = tank.getsolutionvariableprofile("volume")
# create array for mixture density
denprofile_IG = np.zeros_like(timeprofile, dtype=np.double)
# create array for mixture enthalpy
Hprofile_IG = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tank.getsolutionmixtureatindex(solution_index=i)
    # get mixture density profile
    denprofile_IG[i] = solutionmixture.RHO
    # get mixture enthalpy profile
    Hprofile_IG[i] = solutionmixture.HML() / ck.ergsperjoule * 1.0e-3

# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
thispres = str(mist.pressure / ck.Patm)
thistitle = "Cooling Vapor + Air at " + thispres + " atm"
plt.suptitle(thistitle, fontsize=16)
plt.subplot(221)
plt.plot(timeprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(timeprofile, volprofile, "b-", label="real gas")
plt.plot(timeprofile_IG, volprofile_IG, "b--", label="ideal gas")
plt.legend(loc="upper right")
plt.ylabel("Volume [cm3]")
plt.subplot(223)
plt.plot(timeprofile, Hprofile, "g-", label="real gas")
plt.plot(timeprofile_IG, Hprofile_IG, "g--", label="ideal gas")
plt.legend(loc="upper right")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Enthalpy [kJ/mole]")
plt.subplot(224)
plt.plot(timeprofile, denprofile, "m-", label="real gas")
plt.plot(timeprofile_IG, denprofile_IG, "m--", label="ideal gas")
plt.legend(loc="upper left")
plt.xlabel("time [sec]")
plt.ylabel("Mixture Dernsity [g/cm3]")
# display the plots
plt.show()

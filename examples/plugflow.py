import os
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import chemkin as ck  # Chemkin
from chemkin import Color

# chemkin batch reactor models (transient)
from chemkin.flowreactors.PFR import PlugFlowReactor_GivenTemperature
from chemkin.inlet import Inlet

# check working directory
current_dir = os.getcwd()
print("current working directory: " + current_dir)
# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the inlet (mixture + flow rate)
feedstock = Inlet(MyGasMech)
# set inlet temperature [K]
feedstock.temperature = 1444.48
# set inlet/PFR pressure [atm]
feedstock.pressure = 0.83 * ck.Patm
# set inlet composition
feedstock.X = [
    ("AR", 0.8433),
    ("CO", 0.0043),
    ("CO2", 0.0429),
    ("H2O", 0.0956),
    ("N2", 0.0031),
    ("NH3", 0.0021),
    ("NO", 0.0012),
    ("O2", 0.0074),
    ("OH", 4.6476e-5),
]
# set inlet velocity [cm/sec]
feedstock.velocity = 26.815
print(type(feedstock))
#
# create a plug flow reactor instance
tubereactor = PlugFlowReactor_GivenTemperature(feedstock)
# set PFR diameter [cm]
tubereactor.diameter = 5.8431
# set PFR length [cm]
tubereactor.length = 5.0
print(f"PFR inlet mass flow rate {tubereactor.massflowrate} [g/sec]")
print(f"PFR inlet velocity {tubereactor.velocity} [cm/sec]")
# show inlet gas composition of the PFR
print("PFR inlet gas compsition")
tubereactor.listcomposition(mode="mole", bound=1.0e-8)
# set distance between saving solution
tubereactor.timestepforsavingsolution = 0.0005
# turn ON adaptive solution saving
tubereactor.adaptivesolutionsaving(mode=True, steps=100)
# show the additional keywords given by user
tubereactor.showkeywordinputlines()
# set the start wall time
start_time = time.time()
# run the PFR model
runstatus = tubereactor.run()
# compute the total runtime
runtime = time.time() - start_time
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> RUN COMPLETED <<<", end=Color.END)
print(f"total simulation duration: {runtime * 1.0e3} [msec]")
#
# post-process the solution profiles
tubereactor.processsolution()
# get the number of solution time points
solutionpoints = tubereactor.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the grid profile [cm]
xprofile = tubereactor.getsolutionvariableprofile("time")
# get the temperature profile [K]
tempprofile = tubereactor.getsolutionvariableprofile("temperature")
# get the NO mass fraction profile
YNOprofile = tubereactor.getsolutionvariableprofile("NO")
# outlet grid index
xout_index = solutionpoints - 1
print(f"At the reactor outlet: x = {xprofile[xout_index]} [cm]")
print(f"the NO mass fraction = {YNOprofile[xout_index]}")
#
# more involving post-processing by using Mixtures
#
# create arrays for CO, NH3, and NO2 mole fractions
COprofile = np.zeros_like(xprofile, dtype=np.double)
NH3profile = np.zeros_like(xprofile, dtype=np.double)
NO2profile = np.zeros_like(xprofile, dtype=np.double)
velocityprofile = np.zeros_like(xprofile, dtype=np.double)
# find species index
CO_index = MyGasMech.getspecindex("CO2")
NH3_index = MyGasMech.getspecindex("NH3")
NO2_index = MyGasMech.getspecindex("NO2")
# reactor mass flow rate (constant) [g/sec]
massflowrate = tubereactor.massflowrate
# reactor cross-section area [cm2]
areaflow = tubereactor.flowarea
print(f"mass flow {massflowrate} area {areaflow}")
# ratio
ratio = massflowrate / areaflow
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tubereactor.getsolutionmixtureatindex(solution_index=i)
    # get gas density [g/cm3]
    den = solutionmixture.RHO
    # gas velocity [g]
    velocityprofile[i] = ratio / den
    # get CO mole fraction profile
    COprofile[i] = solutionmixture.X[CO_index]
    # get NH3 mole fraction profile
    NH3profile[i] = solutionmixture.X[NH3_index]
    # get NO2 mole fraction profile
    NO2profile[i] = solutionmixture.X[NO2_index]
# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(xprofile, tempprofile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(xprofile, COprofile, "b-")
plt.ylabel("CO Mole Fraction")
plt.subplot(223)
plt.plot(xprofile, NO2profile, "g-")
plt.xlabel("distance [cm]")
plt.ylabel("NO2 Mole Fraction")
plt.subplot(224)
plt.plot(xprofile, velocityprofile, "m-")
plt.xlabel("distance [cm]")
plt.ylabel("Gas Velocity [cm/sec]")
# display the plots
plt.show()

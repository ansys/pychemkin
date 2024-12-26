import os
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import chemkin as ck  # Chemkin
from chemkin import Color

# chemkin batch reactor models (transient)
from chemkin.batchreactors.batchreactor import (
    GivenPressureBatchReactor_EnergyConservation,
)

# check working directory
current_dir = os.getcwd()
print("current working directory: " + current_dir)
# set verbose mode
ck.setverbose(True)
# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
gasoline = ck.Chemistry(label="gasoline 14comp")
# set mechanism input files
# inclusion of the full file path is recommended
gasoline.chemfile = os.path.join(mechanism_dir, "gasoline_14comp_WBencrypt.inp")
# preprocess the mechanism files
iError = gasoline.preprocess()
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the fuel mixture
fuelmixture = ck.Mixture(gasoline)
# set fuel = composition PRF 60
fuelmixture.X = [("ic8h18", 0.6), ("nc7h16", 0.4)]
# setting pressure and temperature
fuelmixture.pressure = 5.0 * ck.Patm
fuelmixture.temperature = 1500.0
# create the oxidizer mixture: air
air = ck.Mixture(gasoline)
air.X = [("o2", 0.21), ("n2", 0.79)]
# setting pressure and temperature
air.pressure = 5.0 * ck.Patm
air.temperature = 1500.0
# create the premixed mixture to be defined
premixed = ck.Mixture(gasoline)
# products from the complete combustion of the fuel mixture and air
products = ["co2", "h2o", "n2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
add_frac = np.zeros(gasoline.KK, dtype=np.double)  # no additives: all zeros
iError = premixed.XbyEquivalenceRatio(
    gasoline, fuelmixture.X, air.X, add_frac, products, equivalenceratio=1.0
)
if iError != 0:
    raise RuntimeError
# list the composition of the premixed mixture
premixed.listcomposition(mode="mole")
# set mixture temperature and pressure (equivalent to setting the initial temperature and pressure of the reactor)
premixed.pressure = 40.0 * ck.Patm
premixed.temperature = 700.0
#
# create a constant pressure batch reactor (with energy equation)
#
MyCONP = GivenPressureBatchReactor_EnergyConservation(premixed, label="CONP")
# show initial gas composition inside the reactor
MyCONP.listcomposition(mode="mole")
# set other reactor parameters
# reactor volume [cm3]
MyCONP.volume = 10.0
# simulation end time [sec]
MyCONP.time = 1.0
# output controls
# set timestep between saving solution
MyCONP.timestepforsavingsolution = 0.001
# change timestep between saving solution
MyCONP.timestepforsavingsolution = 0.01
# turn ON adaptive solution saving
MyCONP.adaptivesolutionsaving(mode=True, value_change=100, target="TEMPERATURE")
# set tolerance
MyCONP.settolerances(absolute_tolerance=1.0e-10, relative_tolerance=1.0e-8)
# set ignition delay
# ck.showignitiondefinition()
MyCONP.setignitiondelay(method="T_inflection")
# stop after ignition is detected (not recommended for ignition delay time calculations)
# MyCONP.stopafterignition()
# show solver option
print(f"timestep between solution printing: {MyCONP.timestepforprintingsolution}")
# show timestep between printing solution
print(f"forced non-negative solution values: {MyCONP.forcenonnegative}")
#
# loop over initial reactor temperature to create an ignition delay time plot
#
npoints = 20
delta_temp = 20.0
init_temp = premixed.temperature
delaytime = np.zeros(npoints, dtype=np.double)
temp_inv = np.zeros_like(delaytime, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all cases with different initial gas/reactor temperatures
for i in range(npoints):
    # update the initial reactor temperature
    MyCONP.temperature = init_temp  # K
    # show the additional keywords given by user (verify inputs before running the simulation)
    # MyCONP.showkeywordinputlines()
    # run the reactor model
    runstatus = MyCONP.run()
    #
    if runstatus == 0:
        # plot 1/T instead of T
        temp_inv[i] = 1.0e0 / init_temp
        # get ignition delay time
        delaytime[i] = MyCONP.getignitiondelay()
        print(f"ignition delay time = {delaytime[i]} [msec]")
    else:
        # if get this, most likely the END time is too short
        print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
    init_temp += delta_temp
# compute the total runtime
runtime = time.time() - start_time
print(f"total simulation duration: {runtime} [sec] for {npoints} cases")
# create an ignition delay versus 1/T plot for the PRF fuel (should exhibit the NTC region)
plt.rcParams.update({"figure.autolayout": True})
fig, ax1 = plt.subplots()
ax1.semilogy(temp_inv, delaytime, "bs--")
ax1.set_xlabel("1/T [1/K]")
ax1.set_ylabel("Ignition delay time [msec]")


# Create a secondary x-axis for T (=1/(1/T))
def one_over(x):
    """Vectorized 1/x, treating x==0 manually"""
    x = np.array(x, float)
    near_zero = np.isclose(x, 0)
    x[near_zero] = np.inf
    x[~near_zero] = 1 / x[~near_zero]
    return x


# the function "1/x" is its own inverse
inverse = one_over
ax2 = ax1.secondary_xaxis("top", functions=(one_over, inverse))
ax2.set_xlabel("T [K]")
plt.show()

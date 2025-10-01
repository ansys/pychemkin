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
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger

# Chemkin 1-D premixed freely propagating flame model (steady-state)
from ansys.chemkin.premixedflames.premixedflame import FreelyPropagating as FlameSpeed
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

# create the fuel mixture
fuel = ck.Mixture(MyGasMech)
# set fuel composition: methane
fuel.X = [("CH4", 1.0)]
# setting pressure and temperature condition for the flame speed calculations
fuel.pressure = 1.0 * ck.Patm
fuel.temperature = 300.0  # inlet temperature

# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.X = ck.Air.X()
# setting pressure and temperature is not required in this case
air.pressure = fuel.pressure
air.temperature = fuel.temperature

# create the fuel-air Stream for the premixed flame speed calculation
premixed = Stream(MyGasMech, label="premixed")
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros

# setting pressure and temperature is not required in this case
premixed.pressure = fuel.pressure
premixed.temperature = fuel.temperature

# set estimated value of the flame mass flux [g/cm2-sec]
premixed.mass_flowrate = 0.4

# equivalence ratio for the first case
phi = 0.6
# create mixture by using the equivalence ratio
iError = premixed.X_by_Equivalence_Ratio(
    MyGasMech, fuel.X, air.X, add_frac, products, equivalenceratio=phi
)
# check fuel-oxidizer mixture creation status
if iError != 0:
    print(
        +"Error: failed to create the methane-air mixture "
        + "for equivalence ratio = "
        + str(phi)
    )
    exit()

# Instantiate the laminar speed calculator
flamespeedcalculator = FlameSpeed(premixed, label="premixed_methane")

# Set up initial mesh and grid adaption options
# set the maximum total number of grid points allowed in the calculation (optional)
# flamespeedcalculator.set_max_grid_points(150)
# define the calculation domain [cm]
flamespeedcalculator.end_position = 1.0

# Run the flame speed parameter study
# total number of parameter cases
points = 21
# equivalence ratio increment
delta_phi = 0.05
# create solution arrays
equival = np.zeros(points, dtype=np.double)
flamespeed = np.zeros_like(equival, dtype=np.double)
# set the start wall time
start_time = time.time()

# start the parameter study runs
for i in range(points):
    # run the flame speed calculation for this equivalence ratio
    status = flamespeedcalculator.run()
    if status != 0:
        print(
            Color.RED
            + "failed to calculate the laminar flame speed"
            + "for equivalence ratio = "
            + str(phi)
            + Color.END
        )
        exit()
    # get flame speed
    # postprocess the solutions
    flamespeedcalculator.process_solution()
    # save data
    equival[i] = phi
    # get flame speed
    flamespeed[i] = flamespeedcalculator.get_flame_speed()
    # print the predicted laminar flame speed
    print(
        f"methane-air equivalence ratio = {phi} :\n"
        + f"the predicted laminar flame speed = {flamespeed[i]} [cm/sec]"
    )
    #
    # update parameter
    phi += delta_phi
    # create mixture by using the equivalence ratio
    iError = premixed.X_by_Equivalence_Ratio(
        MyGasMech, fuel.X, air.X, add_frac, products, equivalenceratio=phi
    )
    # check fuel-oxidizer mixture creation status
    if iError != 0:
        print(
            "Error: failed to create the methane-air mixture ",
            "for equivalence ratio = ",
            str(phi),
        )
        exit()
    # update initial gas composition
    flamespeedcalculator.set_molefractions(premixed.X)

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"total simulation duration: {runtime} [sec]")
print()

# experimental data by Vagelopoulos
# equivalence ratios
data_equiv = [
    0.6126,
    0.6619,
    0.7109,
    0.7533,
    0.8268,
    0.9109,
    0.9826,
    1.0387,
    1.0901,
    1.1321,
    1.1858,
    1.2347,
    1.2695,
    1.325,
    1.3563,
    1.4279,
    1.4977,
]
# methane flame speeds at 1 atm
data_speed = [
    9.4434,
    12.7281,
    17.4088,
    21.0219,
    26.5237,
    33.2573,
    37.0347,
    38.677,
    38.8412,
    37.8558,
    34.9818,
    31.1223,
    25.9489,
    21.3504,
    17.2445,
    12.7281,
    9.7719,
]

# Plot the predicted flame speeds against the experimental data
plt.plot(data_equiv, data_speed, label="data", linestyle="", marker="^", color="blue")
plt.plot(equival, flamespeed, label=MyGasMech.label, linestyle="-", color="blue")
plt.legend()
plt.ylabel("Flame Speed [cm/sec]")
plt.xlabel("Equivalence Ratio")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_flame_speed_table.png", bbox_inches="tight")

# return results for comparisons
resultfile = os.path.join(current_dir, "flamespeedtable.result")
results = {}
results["state-equivalence_ratio"] = equival.tolist()
results["state-flame_speed"] = flamespeed.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

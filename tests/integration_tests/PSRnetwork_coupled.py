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
from ansys.chemkin.stirreactors.PSRcluster import PSRCluster as ERN
from ansys.chemkin.inlet import Mixture
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.stirreactors.PSR import PSR_SetResTime_EnergyConservation as PSR
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
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")

# preprocess the mechanism files
iError = MyGasMech.preprocess()

# Set up gas mixtures based on the species in this chemistry set
# fuel is pure methane
fuel = Mixture(MyGasMech)
fuel.temperature = 650.0  # [K]
fuel.pressure = 10.0 * ck.Patm  # [atm] => [dyne/cm2]
fuel.X = [("CH4", 1.0)]

# air is modeled as a mixture of oxygen and nitrogen
air = Mixture(MyGasMech)
air.temperature = 650.0  # [K]
air.pressure = 10.0 * ck.Patm
air.X = ck.Air.X()  # mole fractions

# primary fuel-air mixture
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture. (You can also create an additives mixture here.)
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros

# create the unburned fuel-air mixture
premixed = Stream(MyGasMech)
# mean equivalence ratio
equiv = 0.6
iError = premixed.X_by_Equivalence_Ratio(
    MyGasMech, fuel.X, air.X, add_frac, products, equivalenceratio=equiv
)
# check fuel-oxidizer mixture creation status
if iError != 0:
    print("Error: Failed to create the fuel-oxidizer mixture.")
    exit()

# list the composition of the unburned fuel-air mixture
premixed.list_composition(mode="mole")

# set premixed fuel-air inlet temperature and mass flow rate
premixed.temperature = fuel.temperature
premixed.pressure = fuel.pressure
premixed.mass_flowrate = 500.0  # [g/sec]

# primary air stream to mix with the primary fuel-air stream
primary_air = Stream(MyGasMech, label="Primary_Air")
primary_air.X = air.X
primary_air.pressure = air.pressure
primary_air.temperature = air.temperature
primary_air.mass_flowrate = 50.0  # [g/sec]

# secondary bypass air stream
secondary_air = Stream(MyGasMech, label="Secondary_Air")
secondary_air.X = air.X
secondary_air.pressure = air.pressure
secondary_air.temperature = 670.0  # [K]
secondary_air.mass_flowrate = 100.0  # [g/sec]

# find the species index
CH4_index = MyGasMech.get_specindex("CH4")
O2_index = MyGasMech.get_specindex("O2")
NO_index = MyGasMech.get_specindex("NO")
CO_index = MyGasMech.get_specindex("CO")

# Define reactors in the reactor network
# PSR #1: mixing zone
mix = PSR(premixed, label="mixing zone")
# use different guess temperature
mix.set_estimate_conditions(option="TP", guess_temp=800.0)
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
mix.residence_time = 0.5 * 1.0e-3
# add external inlets
mix.set_inlet(premixed)
mix.set_inlet(primary_air)

# PSR #2: recirculation zone
recirculation = PSR(premixed, label="recirculation zone")
# use the equilibrium state of the inlet gas mixture as the guessed solution
recirculation.set_estimate_conditions(option="TP", guess_temp=1600.0)
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
recirculation.residence_time = 1.5 * 1.0e-3

# PSR #3: flame zone
flame = PSR(premixed, label="flame zone")
# use the equilibrium state of the inlet gas mixture as the guessed solution
flame.set_estimate_conditions(option="TP", guess_temp=1600.0)
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
flame.residence_time = 1.5 * 1.0e-3
# add external inlet
flame.set_inlet(secondary_air)

# Create the reactor network
# instantiate the PSR network as a coupled reactor network
PSR_list = [mix, recirculation, flame]
PSRcluster = ERN(PSR_list, label="combustor_cluster")

# Define the PSR connectivity
# PSR #1 outlet flow splitting
split_table = [(flame.label, 1.0)]
PSRcluster.set_recycling_stream(mix.label, split_table)
# PSR #2 outlet flow splitting
# PSR #2 does not have an external outlet.
split_table = [(mix.label, 0.15), (flame.label, 0.85)]
PSRcluster.set_recycling_stream(recirculation.label, split_table)
# PSR #3 outlet flow splitting
# part of the outlet flow from PSR #3 exits the reactor network
split_table = [(recirculation.label, 0.2)]
PSRcluster.set_recycling_stream(flame.label, split_table)

# set the start wall time
start_time = time.time()
# solve the reactor network iteratively
status = PSRcluster.run()
if status != 0:
    print(Color.RED + "Failed to solve the reactor network." + Color.END)
    exit()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"Total simulation duration: {runtime} [sec].")
print()

# postprocess the solutions of all PSRs in the cluster
iErr = PSRcluster.process_cluster_solution()

# verify the mass flow rate in and out of the PSR cluster
print(f"net external inlet mass flow rate = {PSRcluster.total_inlet_mass_flow_rate} [g/sec].")
print(f"net outlet mass flow rate = {PSRcluster.get_cluster_outlet_flowrate()} [g/sec].")

# display the reactor solutions
print("=" * 10)
print("reactor/zone")
print("=" * 10)
temp = []
mflr = []
X_CH4 = []
X_CO = []
X_NO = []

for index in range(PSRcluster.numb_PSRs):
    id = index + 1
    name = PSRcluster.get_reactor_label(id)
    # get the solution stream of the reactor
    sstream = PSRcluster.get_reactor_stream(name)
    print(f"Reactor: {name}.")
    print(f"Temperature = {sstream.temperature} [K].")
    print(f"Mass flow rate = {sstream.mass_flowrate} [g/sec].")
    print(f"CH4 = {sstream.X[CH4_index]}.")
    print(f"O2 = {sstream.X[O2_index]}.")
    print(f"CO = {sstream.X[CO_index]}.")
    print(f"NO = {sstream.X[NO_index]}.")
    print("-" * 10)
    # save results to lists for comparisons
    temp.append(sstream.temperature)
    mflr.append(sstream.mass_flowrate)
    X_CH4.append(sstream.X[CH4_index])
    X_CO.append(sstream.X[CO_index])
    X_NO.append(sstream.X[NO_index])

# return results for comparisons
resultfile = os.path.join(current_dir, "PSRnetwork_coupled.result")
results = {}
results["state-temperature"] = temp
results["state-mass_flow_rate"] = mflr
results["species-mole_fraction_CH4"] = X_CH4
results["species-mole_fraction_CO"] = X_CO
results["species-mole_fraction_NO"] = X_NO
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

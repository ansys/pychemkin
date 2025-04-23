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

"""
.. _ref_connected_reactors:

=====================================================================
Use a chain of individual reactors to model a fictional gas combustor
=====================================================================

This tutorial describes the process of setting up and solving a series of linked perfectly-stirred reactors
(PSR) in PyChemkin. This is the simplest reactor network as it does not contain any recycling stream or
outflow splitting.

The PSR chain model of a fictional can combustor is displayed below

 .. figure:: chain_reactor_network.png
   :scale: 80 %
   :alt: the chain reactor network

The *"primary inlet stream"* to the first reactor, the *"combustor"*, is the fuel-lean methane-air mixture
that is formed by mixing the fuel (methane) and the heated air. The exhaust from the *"combustor"* will enter the
second reactor, the "dilution zone"*, where the hot combustion products will be cooled by the introduction of
additional cool air. The cooled and diluted gas mixture in the *"dilution zone"* will then travel to the third
reactor, the *"reburning zone"*. A mixture of fuel (methane) and carbon dioxide is injected to the gas in the
*"reburning zone"* attempting to convert any remaining carbon monoxide or nitric oxide in the exhaust gas to
carbon dioxide or nitrogen, respectively.

In this tutorial, the reactors will be solved one by one from upstream to downstream. Once the solution of the
upstream reactor is obtained, it will be used to set up the external inlet of the immediate downstream reactor.
This process will continue till all reactors in the chain network are solved. Since there is no recycling stream
in this configuration, the entire reactor network can be solved in one sweep.
"""

# sphinx_gallery_thumbnail_path = '_static/chain_reactor_network.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.inlet import adiabatic_mixing_streams
from ansys.chemkin.logger import logger

# chemkin perfectly-stirred reactor (PSR) model (steady-state)
from ansys.chemkin.stirreactors.PSR import PSR_SetResTime_EnergyConservation as PSR

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a png file
global interactive
interactive = True

#####################################
# Create a ``Chemistry Set`` instance
# ===================================
# The mechanism is the GRI 3.0 mechanism for methane combustion.
# The mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"/reaction/data"*.

# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")

############################################
# Pre-process the gasoline ``Chemistry Set``
# ==========================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

####################################################################
# Set up gas mixtures based on the species in this ``Chemistry Set``
# ==================================================================
# Create the "fuel" and the "air" mixtures to initialize the external
# inlet streams. The "fuel" for this case is pure methane.

# fuel is pure methane
fuel = Stream(MyGasMech)
fuel.temperature = 300.0  # [K]
fuel.pressure = 2.1 * ck.Patm  # [atm] => [dyne/cm2]
fuel.X = [("CH4", 1.0)]
fuel.mass_flowrate = 3.275  # [g/sec]

# air is modeled as a mixture of oxygen and nitrogen
air = Stream(MyGasMech)
air.temperature = 550.0  # [K]
air.pressure = 2.1 * ck.Patm
air.X = ck.Air.X()  # use predefined "air" recipe in mole fractions
air.mass_flowrate = 45.0  # [g/sec]

#################################################
# Create external inlet streams from the mixtures
# ===============================================
# Use the ``Stream`` method ``adiabatic_mixing_streams`` to combine
# the ``fuel`` and the ``air`` streams. The final gas temperature should
# land between the temperatures of the two source streams. The mass flow
# rate of the ``premixed`` stream should be the sum of the sources. A simple
# *PyChemkin* composition "recipe" is used to create the ``reburn_fuel`` stream.

# premixed stream for the combustor
premixed = adiabatic_mixing_streams(fuel, air)
print(str(premixed.temperature))
print(str(premixed.mass_flowrate))

# additional fuel injection for the reburning zone
reburn_fuel = Stream(MyGasMech)
reburn_fuel.temperature = 300.0  # [K]
reburn_fuel.pressure = 2.1 * ck.Patm  # [atm] => [dyne/cm2]
reburn_fuel.X = [("CH4", 0.6), ("CO2", 0.4)]
reburn_fuel.mass_flowrate = 0.12  # [g/sec]

# find the species index
CH4_index = MyGasMech.get_specindex("CH4")
O2_index = MyGasMech.get_specindex("O2")
NO_index = MyGasMech.get_specindex("NO")
CO_index = MyGasMech.get_specindex("CO")

#####################################
# Create individual PSR for each zone
# ===================================
# Set up the PSR for each zone one by one with external inlets only. PyChemkin
# requires that the **first reactor/zone must have at least one external inlet**.
# There are three reactors in the network; from upstream to downstream, they
# are ``combustor``, ``dilution zone``, and ``reburning zone``. And all of them
# have one external inlet. For the two downstream reactors, you need to create
# the "through flow" stream from their respective upstream reactor by using the
# solution of the upstream reactor. Get the solution stream from the upstream
# reactor using the ``process_solution`` method, then use the ``set_inlet`` method
# to connect the resulting solution stream to the downstream reactor.
#
# .. note::
#   *PyChemkin* requires that the **first** reactor/zone must have at least
#   **one external inlet**. The rest of the reactors will have at least the
#   "through flow" from the immediate upstream reactor so they do not require
#   an external inlet.
#
# .. note::
#   The ``Stream`` parameter used to instantiate a ``PSR`` object is used to establish
#   the *guessed reactor solution* and will be modified when the network is solved by
#   the ``ERN``.
#
# .. note::
#   The reactors in the "network" must be post-processed individually.
#

# PSR #1: combustor
combustor = PSR(premixed, label="combustor")
# use the equilibrium state of the inlet gas mixture as the guessed solution
combustor.set_estimate_conditions(option="HP")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
combustor.residence_time = 2.0 * 1.0e-3
# add external inlet
combustor.set_inlet(premixed)
# set the start wall time
start_time = time.time()
# run PSR #1
status = combustor.run()
if status != 0:
    print(Color.RED + combustor.label + " run failed!")
    exit()
# post-process the solution profiles
solnstream1 = combustor.process_solution()
solnstream1.label = "PSR1"
print("=" * 40)
print("combustor exit")
print("=" * 40)
print(f"temperature = {solnstream1.temperature} [K]")
print(f"mass flow rate = {solnstream1.mass_flowrate} [g/sec]")
print(f"CH4 = {solnstream1.X[CH4_index]}")
print(f"O2 = {solnstream1.X[O2_index]}")
print(f"CO = {solnstream1.X[CO_index]}")
print(f"NO = {solnstream1.X[NO_index]}")

# PSR #2: cooling
cooling = PSR(solnstream1, label="cooling zone")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
cooling.residence_time = 1.5 * 1.0e-3
# add external inlet
air.mass_flowrate = 62.0  # [g/sec]
cooling.set_inlet(air)
# add the through flow from PSR #1
cooling.set_inlet(solnstream1)
# run PSR #2
status = cooling.run()
if status != 0:
    print(Color.RED + cooling.label + " run failed!")
    exit()
# post-process the solution profiles
solnstream2 = cooling.process_solution()
solnstream2.label = "PSR2"
print()
print("=" * 40)
print("dilution zone exit")
print("=" * 40)
print(f"temperature = {solnstream2.temperature} [K]")
print(f"mass flow rate = {solnstream2.mass_flowrate} [g/sec]")
print(f"CH4 = {solnstream2.X[CH4_index]}")
print(f"O2 = {solnstream2.X[O2_index]}")
print(f"CO = {solnstream2.X[CO_index]}")
print(f"NO = {solnstream2.X[NO_index]}")

# PSR #3: reburn
reburn = PSR(solnstream2, label="reburn zone")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
reburn.residence_time = 3.5 * 1.0e-3
# add external inlet
reburn.set_inlet(reburn_fuel)
# add the through flow from PSR #2
reburn.set_inlet(solnstream2)
# run PSR #3
status = reburn.run()
if status != 0:
    print(Color.RED + reburn.label + " run failed!")
    exit()
# post-process the solution profiles
outflow = reburn.process_solution()
outflow.label = "outflow"
print()
print("=" * 40)
print("outflow")
print("=" * 40)
print(f"temperature = {outflow.temperature} [K]")
print(f"mass flow rate = {outflow.mass_flowrate} [g/sec]")
print(f"CH4 = {outflow.X[CH4_index]}")
print(f"O2 = {outflow.X[O2_index]}")
print(f"CO = {outflow.X[CO_index]}")
print(f"NO = {outflow.X[NO_index]}")

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"total simulation duration: {runtime} [sec]")

# return results for comparisons
resultfile = os.path.join(current_dir, "PSRChain_declustered.result")
results = {}
results["state-temperature"] = [outflow.temperature]
results["state-mass_flow_rate"] = [outflow.mass_flowrate]
results["species-mole_fraction_CH4"] = [outflow.X[CH4_index]]
results["species-mole_fraction_CO"] = [outflow.X[CO_index]]
results["species-mole_fraction_NO"] = [outflow.X[NO_index]]
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

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
.. _ref_chain_reactor_network:

==============================================================
Use a chain reactor network to model a fictional gas combustor
==============================================================

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

In this tutorial, the chain reactor network will be configured and solved by using the ``ReactorNetwork`` module.
The module automatically handles the tasks of running the individual reactors and setting up the inlet to the downstream
reactor.
"""

# sphinx_gallery_thumbnail_path = '_static/chain_reactor_network.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.hybridreactornetwork import ReactorNetwork as ERN
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.inlet import adiabatic_mixing_streams
from ansys.chemkin.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.stirreactors.PSR import PSR_SetResTime_EnergyConservation as PSR

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)

#####################################
# Create a chemistry set
# ===================================
# The mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
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
# Create the ``fuel`` and the ``air`` streams/mixtures before setting up the
# external inlet streams. The "fuel" in this case is pure methane. The main
# inlet stream ``premixed`` to the ``combustor`` is formed by mixing the
# ``fuel`` and the ``air`` streams adiabatically. The fuel to air mass ratio
# is provided implicitly by the *mass flow rates* of the two streams. The
# external inlet to the second reactor, the ``dilution zone``, is simply the
# ``air`` stream with a different mass flow rate (and different temperature
# if desirable). The ``reburn_fuel`` stream to be injected to the downstream
# ``reburning zone`` is a mixture of methane and carbon dioxide.
#
# .. note::
#   PyChemkin has *"air"* redefined as a convenient way to set up the air
#   stream/mixture in the simulations. Use ``ansys.chemkin.Air.X()`` or
#   ``ansys.chemkin.Air.Y()`` when the mechanism uses "O2" and "N2" for
#   oxygen and nitrogen. Use ``ansys.chemkin.air.X()`` or ``ansys.chemkin.air.Y()``
#   when oxygen and nitrogen are represented by "o2" and "n2".
#

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
# use predefined "air" recipe in mole fractions (with upper cased symbols)
air.X = ck.Air.X()
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

# verify the premixed stream properties
print(f"premixed stream temperature = {premixed.temperature} [K]")
print(f"premixed stream mass flow rate = {premixed.mass_flowrate} [g/sec]")

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
# Set up the PSR for each zone one by one with *external inlets only*.
# For PSRs, use the ``set_inlet`` method to add the external inlets to the
# reactor. PFR always requires *one* external inlet when it is instantiated.
#
# There are three reactors in the network; from upstream to downstream, they
# are ``combustor``, ``dilution zone``, and ``reburning zone``. And all of them
# have one external inlet.
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

# PSR #1: combustor
combustor = PSR(premixed, label="combustor")
# use the equilibrium state of the inlet gas mixture as the guessed solution
combustor.set_estimate_conditions(option="HP")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
combustor.residence_time = 2.0 * 1.0e-3
# add external inlet
combustor.set_inlet(premixed)

# PSR #2: dilution zone
dilution = PSR(premixed, label="dilution zone")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
dilution.residence_time = 1.5 * 1.0e-3
# add external inlet
# first, assign the "correct" mass flow rate to the "air" stream
air.mass_flowrate = 62.0  # [g/sec]
dilution.set_inlet(air)

# PSR #3: reburning zone
reburn = PSR(premixed, label="reburning zone")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
reburn.residence_time = 3.5 * 1.0e-3
# add external inlet
reburn.set_inlet(reburn_fuel)

#####################################
# Create the reactor network
# ===================================
# Create a hybrid ``ReactorNetwork`` object ``PSRChain`` and add the reactors one
# by one from upstream to downstream using the ``add_reactor`` method. For a
# simple chain network such as the one used in the current tutorial, you do not
# need to define the connectivity among the reactors. The reactor network model
# will figure out the through flow connections automatically.
#
# .. note::
#   Use the ``show_reactors`` method to get the list of reactors in the network in
#   the order they are added.
#
# .. note::
#   Use the ``remove_reactor`` method to remove an existing reactor from the
#   network by the reactor ``name/label``. Similarly, use the ``clear_connections``
#   to undo the network connectivity.
#
# .. note::
#   The order of the reactor addition is important as it dictates the solution
#   sequence and thus the convergence rate.
#

# instantiate the chain PSR network as a hy
PSRChain = ERN(MyGasMech)

# add the reactors from upstream to downstream
PSRChain.add_reactor(combustor)
PSRChain.add_reactor(dilution)
PSRChain.add_reactor(reburn)

# list the reactors in the network
PSRChain.show_reactors()

###########################
# Solve the reactor network
# =========================
# Use the ``run`` method to solve the entire reactor network. The hybrid ``ReactorNetwork``
# will solve the reactors one by one in the order they are added to the network.
#

# set the start wall time
start_time = time.time()

# solve the reactor network
status = PSRChain.run()
if status != 0:
    print(Color.RED + "failed to solve the reactor network!" + Color.END)
    exit()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"total simulation duration: {runtime} [sec]")

############################################
# Post-process the reactor network solutions
# ==========================================
# There are two ways to process the results from a reactor network. You
# can extract the solution of an individual reactor member as a ``Stream``
# object by using the ``get_reactor_stream`` by the reactor ``name``. Or,
# get the stream properties of a specific network outlet by using the
# ``get_external_stream`` method with the outlet index. Once you get the
# solution as a ``Stream`` object, you can use any ``Stream`` or ``Mixture``
# method to further manipulate the solutions.
#
# .. note::
#   Use the ``number_external_outlets`` method to find out the number of
#   external exists of the reactor network.
#

# get the outlet stream from the reactor network solutions
# find the number of external outlet stream from the reactor network
print(f"number of outlet stream = {PSRChain.number_external_outlets}")

# get the first (and the only) external outlet stream properties
network_outflow = PSRChain.get_external_stream(1)
# set the stream label
network_outflow.label = "outflow"

# print the desired outlet stream properties
print()
print("=" * 10)
print("outflow")
print("=" * 10)
print(f"temperature = {network_outflow.temperature} [K]")
print(f"mass flow rate = {network_outflow.mass_flowrate} [g/sec]")
print(f"CH4 = {network_outflow.X[CH4_index]}")
print(f"O2 = {network_outflow.X[O2_index]}")
print(f"CO = {network_outflow.X[CO_index]}")
print(f"NO = {network_outflow.X[NO_index]}")

# return results for comparisons
resultfile = os.path.join(current_dir, "PSRChain_network.result")
results = {}
results["state-temperature"] = [network_outflow.temperature]
results["state-mass_flow_rate"] = [network_outflow.mass_flowrate]
results["species-mole_fraction_CH4"] = [network_outflow.X[CH4_index]]
results["species-mole_fraction_CO"] = [network_outflow.X[CO_index]]
results["species-mole_fraction_NO"] = [network_outflow.X[NO_index]]
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

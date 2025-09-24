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
.. _ref_reactor_network_heat_exchange:

============================================
A simple reactor network with heat exchange
============================================

This example shows how to set up and solve a series of linked PSRs (perfectly-stirred reactors)
with heat exchange.

Here is a PSR chain model of a fictional gas combustor with heat exchange:

 .. figure:: reactor_network_heat_exchange.png
   :scale: 80 %
   :alt: Chain reactor network with heat exchange

The primary fuel inlet stream to the first reactor, the *pre-heater*, contains pure methane. The outlet
flow from the *pre-heater* enters the second reactor, the *mixer*, where the hot fuel would mix with
the cool air from the external inlet to form a combustible mathen-air mixture. The combustible mixture from
the *mixer* then travel to the third reactor, the *combustor* in which the methane-air mixture would ignite
and burn. The *combustor* is connected to the *pre-heater* by a heat exchanger so that a portion of the heat
generated from the combustion in the *combustor* will be recycled to the *pre-heeater* to heat up
the fuel stream.

This example uses the ``ReactorNetwork`` module to configure and solve this chain reactor network
with heat exchange iteratively. This module automatically handles the tasks of running the individual
reactors and setting up the inlet to the downstream reactor.
"""

# sphinx_gallery_thumbnail_path = '_static/reactor_network_heat_exchange.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

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

########################
# Create a chemistry set
# ======================
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

#######################################
# Preprocess the gasoline chemistry set
# =====================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

################################################################
# Set up gas mixtures based on the species in this chemistry set
# ==============================================================
# Create the ``fuel`` and ``air`` streams before setting up the
# external inlet streams. The fuel in this case is pure methane. The main
# ``premixed`` inlet stream to the combustor is formed by mixing the
# ``fuel`` and ``air`` streams adiabatically. The fuel-to-air mass ratio
# is provided implicitly by the mass flow rates of the two streams. The
# external inlet to the second reactor, the ``dilution zone``, is simply the
# ``air`` stream with a different mass flow rate (and different temperature
# if desirable). The ``reburn_fuel`` stream to inject to the downstream
# reburning zone is a mixture of methane and carbon dioxide.
#
# .. note::
#   PyChemkin has ``air`` predefined as a convenient way to set up the air
#   stream/mixture in simulations. Use the ``ansys.chemkin.Air.X()`` or
#   ``ansys.chemkin.Air.Y()`` method when the mechanism uses "O2" and "N2" for
#   oxygen and nitrogen. Use the ``ansys.chemkin.air.X()`` or ``ansys.chemkin.air.Y()``
#   method when the mechanism uses "o2" and "n2" for oxygen and nitrogen.
#

# fuel is pure methane
fuel = Stream(MyGasMech)
fuel.temperature = 300.0  # [K]
fuel.pressure = ck.Patm  # [atm] => [dyne/cm2]
fuel.X = [("CH4", 1.0)]
fuel.mass_flowrate = 3.275  # [g/sec]

# air is modeled as a mixture of oxygen and nitrogen
air = Stream(MyGasMech)
air.temperature = 300.0  # [K]
air.pressure = ck.Patm
# use predefined "air" recipe in mole fractions (with upper cased symbols)
air.X = ck.Air.X()
air.mass_flowrate = 66.75  # [g/sec]

# prepare a premixed stream as the guessed condition for the combustor
premixed = adiabatic_mixing_streams(fuel, air)

# find the species index
CH4_index = MyGasMech.get_specindex("CH4")
O2_index = MyGasMech.get_specindex("O2")
NO_index = MyGasMech.get_specindex("NO")
CO_index = MyGasMech.get_specindex("CO")

###########################
# Create PSRs for each zone
# =========================
# Set up the PSR for each zone one by one with *external inlets only*.
# For PSR creation, use the ``set_inlet()`` method to add the external inlets to the
# reactor. A PFR always requires one external inlet when it is instantiated.
#
# There are three reactors in the network. From upstream to downstream, they
# are ``combustor``, ``dilution zone``, and ``reburning zone``. All of them
# have one external inlet.
#
# .. note::
#   PyChemkin requires that the first reactor/zone must have at least
#   one external inlet. Because the rest of the reactors have at least the
#   through flow from the immediate upstream reactor, they do not require
#   an external inlet.
#
# .. note::
#   The ``Stream`` parameter used to instantiate a PSR is used to establish
#   the *guessed reactor solution* and is modified when the network is solved by
#   the ``ERN``.
#

# PSR #1: pre-heater
preheater = PSR(air, label="pre-heater")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
preheater.residence_time = 1.5 * 1.0e-3
# add external inlet
preheater.set_inlet(air)

# PSR #2: mixer
mixer = PSR(fuel, label="mixer")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
mixer.residence_time = 0.5 * 1.0e-3
# add external inlet
mixer.set_inlet(fuel)

# PSR #3: combustor
# use the premixed methane-air mixture to set up the combustor
# premixed.temperature = 1800.0
combustor = PSR(premixed, label="combustor")
# use the equilibrium state of the premixed methane-air mixture as the guessed solution
combustor.set_estimate_conditions(option="HP")
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
combustor.residence_time = 4.0 * 1.0e-3

############################
# Create the reactor network
# ==========================
# Create a hybrid reactor network named ``PSRChain`` and use the ``add_reactor()``
# method to add the reactors one by one from upstream to downstream. For a
# simple chain network such as the one used in this example, you do not
# need to define the connectivity among the reactors. The reactor network model
# automatically figures out the through-flow connections.
#
# .. note::
#
#   - Use the ``show_reactors()`` method to get the list of reactors in the network in
#     the order they are added.
#
#   - Use the ``remove_reactor()`` method to remove an existing reactor from the
#     network by the reactor ``name/label``. Similarly, use the ``clear_connections()``
#     method to undo the network connectivity.
#
#   - The order of the reactor addition is important as it dictates the solution
#   sequence and thus the convergence rate.
#

# instantiate the chain PSR network as a hybrid reactor network
PSRChain = ERN(MyGasMech)

# add the reactors from upstream to downstream
PSRChain.add_reactor(preheater)
PSRChain.add_reactor(mixer)
PSRChain.add_reactor(combustor)

# list the reactors in the network
PSRChain.show_reactors()

##################################
# Add the heat transfer connection
# ================================
# Use ``add_heat_exchange()`` method to define a heat exchange connection
# between the *pre-heater* and the *combustor*.

# effecctive heat transfer coefficient of the heat exchanger [cal/cm2-K-sec]
HT_Coeff = 0.0025
# effective surface area of the heat exchanger [cm2]
HT_area = 100.0
PSRChain.add_heat_exchange(preheater.label, combustor.label, HT_Coeff, HT_area)

# reset the network relative tolerance
PSRChain.set_tear_tolerance(1.0e-7)

###########################
# Solve the reactor network
# =========================
# Use the ``run()`` method to solve the entire reactor network. The hybrid reactor network
# solves the reactors one by one in the order that they are added to the network.
#

# set the start wall time
start_time = time.time()

# solve the reactor network
status = PSRChain.run()
if status != 0:
    print(Color.RED + "Failed to solve the reactor network." + Color.END)
    exit()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"Total simulation duration: {runtime} [sec]")

#####################################
# Postprocess reactor network results
# ===================================
# There are two ways to process results from a reactor network. You
# can extract the solution of an individual reactor member as a stream
# by using the ``get_reactor_stream()`` method to get the reactor by its name. Or,
# you can get the stream properties of a specific network outlet by using the
# ``get_external_stream()`` method to get the stream by its outlet index. Once
# you get the solution as a stream, you can use any stream or mixture
# method to further manipulate the solutions.
#
# .. note::
#   Use the ``number_external_outlets()`` method to find out the number of
#   external outlets of the reactor network.
#

# display the final temperatures in reactor #1 (pre-heater) and reactor #2 (mixer)
psr1_mixture = PSRChain.get_reactor_stream(preheater.label)
psr2_mixture = PSRChain.get_reactor_stream(mixer.label)
print(f"temperature in {preheater.label} = {psr1_mixture.temperature} [K].")
print(f"temperature in {mixer.label} = {psr2_mixture.temperature} [K].")

# get the outlet stream from the reactor network solutions
# find the number of external outlet streams from the reactor network
print(f"Number of outlet streams = {PSRChain.number_external_outlets}.")

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

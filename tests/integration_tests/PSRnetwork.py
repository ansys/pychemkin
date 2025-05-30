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
.. _ref_equivalent_reactor_network:

==============================================================================
Simulate a combustor using an equivalent reactor network with stream recycling
==============================================================================

This tutorial describes the process of setting up and solving an equivalent reactor network (ERN) in PyChemkin.

An equivalent reactor network is employed as a reduced-order model to simulate the steady-state
combustion process inside a gas turbine combustor chamber. This reduced-order reactor network model retains the
complexity of the combustion chemistry by sacrificing details of the combustor geometries, the spatial resolution and
the mass and energy transfer processes. An ERN usually comprises of PSRs (perfectly-stirred reactors) and
PFRs (plug-flow reactors). The network configuration and connectivity,
the reactor parameters, and the mass flow rates can be determined from "hot" steady-state CFD simulation results and/or
from observations and measured data of the actual/similar devices. Once a reactor network is calibrated against the
experimental data of a gas combustor, it becomes a handy tool to quickly estimate the emissions from the combustor
when it is subjected to certain variations in the fuel compositions.

The proposed ERN model of a fictional gas turbine combustor is shown below

 .. figure:: combustor_ERN.png
   :scale: 80 %
   :alt: the combustor reactor network

The *"primary fuel"* is mixed with the incoming air to form a *fuel-lean* mixture before entering the chamber through the
primary inlet. Additional air, the *"primary air"*, is introduced to combustion chamber separately through openings
surrounding the primary inlet. The *"secondary air"* will be entrained into the combustion chamber through well placed holes
on the liners at a location slightly downstream from the primary inlet.

The first PSR (reactor #1) represents the *"mixing zone"* around the main injector where the cool *"premixed fuel-air"* stream and
the *"primary air"* stream are preheated by mixing with the hot combustion products from the *"recirculation zone"*. Downstream
from PSR #1, PSR #2, the *"flame zone"*, is where the combustion of the heated fuel-air mixture takes place. The *"secondary air"*
is injected here to cool down the combustion exhaust before it exits the combustion chamber. A portion of the exhaust gas coming
out of PSR #2 does not leave the combustion chamber directly and is diverted to PSR #3, the *"recirculation zone"*. The majority
of the outlet flow from PSR #3 is recirculated back to the *"flame zone"* to sustain the fuel-lean premixed flame there. The rest
of the hot gas from PSR #3 will travel further back to PSR #1, the *"mixing zone"*, to preheat the fuel-lean mixture just entering
the combustion chamber. And finally, the cooled flue gas leaving the chamber in a stream-like manner. Typically, a PFR will be
applied to simulation the out-flow.

The reactors in the network will be solved individually one by one. When there is a *"tear stream"* in the network,
the ERN should be solved iteratively. A *"tear stream"* is usually a "recycle" stream that can serve as a "pivot point"
for solving the recycle network iteratively. The convergence of the ERN is then determined by the absence of the
per-iteration variation of the computed *"tear stream"* properties such as species composition. In the current project,
the recycling streams from PSR #3 to PSR #1 and PSR #2 are *"tear streams"*.
"""

# sphinx_gallery_thumbnail_path = '_static/combustor_ERN.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.hybridreactornetwork import ReactorNetwork as ERN
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
interactive = True

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
# Create the "fuel" and the "air" mixtures to initialize the external
# inlet streams. The "fuel" for this case is pure methane.

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

#################################################
# Create external inlet streams from the mixtures
# ===============================================
# Create the ``fuel`` and the ``air`` streams/mixtures before setting up the
# external inlet streams. The "fuel" in this case is pure methane. The main
# inlet stream ``premixed`` to the ``mixing zone`` is formed by using the
# ``X_by_Equivalence_Ratio`` method. The external inlets to the first reactor,
# the *"mixing zone"*, and the second reactor, the *"flame zone"*, are simply
# ``air`` mixtures with different mass flow rates and temperatures.
#
# .. note::
#   PyChemkin has *"air"* redefined as a convenient way to set up the air
#   stream/mixture in the simulations. Use ``ansys.chemkin.Air.X()`` or
#   ``ansys.chemkin.Air.Y()`` when the mechanism uses "O2" and "N2" for
#   oxygen and nitrogen. Use ``ansys.chemkin.air.X()`` or ``ansys.chemkin.air.Y()``
#   when oxygen and nitrogen are represented by "o2" and "n2".
#

# primary fuel-air mixture
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
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

# primary air stream to be mixed with the primary fuel-air stream
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

##################################################
# Define reactors in the reactor network
# ================================================
# Set up the PSR for each zone one by one with *external inlets only*.
# For PSRs, use the ``set_inlet`` method to add the external inlets to the
# reactor. PFR always requires *one* external inlet when it is instantiated.
#
# There are three reactors in the network; from upstream to downstream, they
# are ``premix zone``, ``flame zone``, and ``recirculation zone``. The
# *"recirculation zone"* does not have external inlet.
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

# PSR #1: mixing zone
mix = PSR(premixed, label="mixing zone")
# use different guess temperature
mix.set_estimate_conditions(option="TP", guess_temp=800.0)
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
mix.residence_time = 0.5 * 1.0e-3
# add external inlets
mix.set_inlet(premixed)
mix.set_inlet(primary_air)


# PSR #2: flame zone
flame = PSR(premixed, label="flame zone")
# use use the equilibrium state of the inlet gas mixture as the guessed solution
flame.set_estimate_conditions(option="TP", guess_temp=1600.0)
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
flame.residence_time = 1.5 * 1.0e-3
# add external inlet
flame.set_inlet(secondary_air)

# PSR #3: recirculation zone
recirculation = PSR(premixed, label="recirculation zone")
# use use the equilibrium state of the inlet gas mixture as the guessed solution
recirculation.set_estimate_conditions(option="TP", guess_temp=1600.0)
# set PSR residence time (sec): required for PSR_SetResTime_EnergyConservation model
recirculation.residence_time = 1.5 * 1.0e-3

#####################################
# Create the reactor network
# ===================================
# Create a hybrid ``ReactorNetwork`` object ``PSRnetwork`` and add the reactors one
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
PSRnetwork = ERN(MyGasMech)

# add the reactors from upstream to downstream
PSRnetwork.add_reactor(mix)
PSRnetwork.add_reactor(flame)
PSRnetwork.add_reactor(recirculation)

# list the reactors in the network
PSRnetwork.show_reactors()

##################################
# Define the PSR connectivity
# ==============================
# Because the current reactor network contains recycling streams, for example, from PSR #3 to
# PSR #1 and PSR #2, the network connectivity must be specified explicitly. This is done by
# using a *"split"* dictionary to define the outflow connections among the PSRs in the network.
# The *key* of the *"split"* dictionary is the label of the *"originate PSR"*. The *value* is a
# list of tuples consisting of the label of the *"target PSR"* and its *"mass flow rate fraction"*.
#
# ::
#
#   { "originate PSR" : [("target PSR1", fraction), ("target PSR2", fraction)],
#     "another originate PSR" : [("target PSR1", fraction), ... ], ... }
#
# For example, the outflow from reactor "PSR1" is split into two streams. 90% of the
# mass flow rate goes to "PSR2" and the rest is diverted to "PSR5". The split dictionary entry for
# the outflow split of "PSR1" becomes
#
# ::
#
#   "PSR1" : [("PSR2", 0.9), ("PSR5", 0.1)]
#
# .. note::
#   The outlet flow from a reactor that is leaving the reactor network must be labeled as "EXIT>>"
#   when you define the outflow splitting.
#

# PSR #1 outlet flow splitting
split_table = [(flame.label, 1.0)]
PSRnetwork.add_outflow_connections(mix.label, split_table)
# PSR #2 outlet flow splitting
# part of the outlet flow from PSR #2 exits the reactor network
split_table = [(recirculation.label, 0.2), ("EXIT>>", 0.8)]
PSRnetwork.add_outflow_connections(flame.label, split_table)
# PSR #3 outlet flow splitting
# PSR #3 is the last reactor of the network, however it does not have external outlet
split_table = [(mix.label, 0.15), (flame.label, 0.85)]
PSRnetwork.add_outflow_connections(recirculation.label, split_table)

##############################################
# Define the tear point for the iteration
# ============================================
# Because the reactor network contains recycling streams, it has to
# solved iteratively by applying the *tear stream* method. The
# *"tear points"* of the network must be explicitly defined using the
# ``add_tearingpoint`` method. In this project, the stream recycling
# occurs at reactor #3 where the outlet flow is sent back to
# the upstream reactors, reactor #1 and reactor #2. Therefore, reactor
# #3, the "recirculation zone", should be set as the only tear point of
# this reactor network. The relative tolerance for the overall
# reactor network convergence can be set by the ``set_tear_tolerance``
# method. The default tolerance is 1.0e-6.
#

# define the tear point reactor as reactor #3
PSRnetwork.add_tearingpoint(recirculation.label)

# reset the network relative tolerance
PSRnetwork.set_tear_tolerance(1.0e-5)

###########################
# Solve the reactor network
# =========================
# Use the ``run`` method to solve the entire reactor network. The hybrid ``ReactorNetwork``
# will solve the reactors one by one in the order they are added to the network.
#
# .. note::
#   Use the ``set_tear_iteration_limit(count)`` method to change the limit on the number of
#   tear loop tear loop iterations can be taken before declaring the failure.
#   The default limit is 200.
#
# .. note::
#   The ``set_relaxation_factor(factor)`` method can be employed to make solution updating at the
#   end of each iteration to be more aggressive (factor > 1.0) or more conservative (factor < 1.0).
#   By default, the relaxation factor is set to 1.
#

# set the start wall time
start_time = time.time()
# solve the reactor network iteratively
status = PSRnetwork.run()
if status != 0:
    print(Color.RED + "failed to solve the reactor network!" + Color.END)
    exit()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"total simulation duration: {runtime} [sec]")
print()

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
print(f"number of outlet stream = {PSRnetwork.number_external_outlets}")

# display the external outlet stream properties
for m in range(PSRnetwork.number_external_outlets):
    n = m + 1
    network_outflow = PSRnetwork.get_external_stream(n)
    # set the stream label
    network_outflow.label = "outflow"

    # print the desired outlet stream properties
    print("=" * 10)
    print(f"outflow # {n}")
    print("=" * 10)
    print(f"temperature = {network_outflow.temperature} [K]")
    print(f"mass flow rate = {network_outflow.mass_flowrate} [g/sec]")
    print(f"CH4 = {network_outflow.X[CH4_index]}")
    print(f"O2 = {network_outflow.X[O2_index]}")
    print(f"CO = {network_outflow.X[CO_index]}")
    print(f"NO = {network_outflow.X[NO_index]}")
    print("-" * 10)

# display the raector solutions
print()
print("=" * 10)
print("reactor/zone")
print("=" * 10)
temp = []
mflr = []
X_CH4 = []
X_CO = []
X_NO = []
for index, stream in PSRnetwork.reactor_solutions.items():
    name = PSRnetwork.get_reactor_label(index)
    print(f"reactor: {name}")
    print(f"temperature = {stream.temperature} [K]")
    print(f"mass flow rate = {stream.mass_flowrate} [g/sec]")
    print(f"CH4 = {stream.X[CH4_index]}")
    print(f"O2 = {stream.X[O2_index]}")
    print(f"CO = {stream.X[CO_index]}")
    print(f"NO = {stream.X[NO_index]}")
    # save results to lists for comparisons
    temp.append(stream.temperature)
    mflr.append(stream.mass_flowrate)
    X_CH4.append(stream.X[CH4_index])
    X_CO.append(stream.X[CO_index])
    X_NO.append(stream.X[NO_index])
    print("-" * 10)

# return results for comparisons
resultfile = os.path.join(current_dir, "PSRnetwork.result")
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

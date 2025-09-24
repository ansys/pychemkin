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
.. _ref_coupled_equivalent_reactor_network:

=====================================================================================
Simulate a combustor using a coupled equivalent reactor network with stream recycling
=====================================================================================

This example shows how to set up and solve a "coupled" ERN (equivalent reactor network) in PyChemkin.

An ERN is employed as a reduced-order model to simulate the steady-state
combustion process inside a gas turbine combustor chamber. This reduced-order reactor
network model retains the complexity of the combustion chemistry by sacrificing details
of the combustor geometries, the spatial resolution, and the mass and energy transfer processes.
An ERN usually comprises PSRs (perfectly-stirred reactors) and PFRs (plug-flow reactors).
The network configuration and connectivity, the reactor parameters, and the mass flow rates
can be determined from "hot" steady-state CFD simulation results and/or from observations and
measured data of the actual/similar devices. Once a reactor network is calibrated against
the experimental data of a gas combustor, it becomes a handy tool for quickly estimating
the emissions from the combustor when it is subjected to certain variations in the fuel compositions.

This figure shows a proposed ERN model of a fictional gas turbine combustor.

 .. figure:: combustor_ERN.png
   :scale: 80 %
   :alt: Combustor reactor network

The *primary fuel* is mixed with the incoming air to form a *fuel-lean* mixture before entering
the chamber through the primary inlet. Additional air, the *primary air*, is introduced to
the combustion chamber separately through openings surrounding the primary inlet.
The *secondary air* is entrained into the combustion chamber through well-placed holes
on the liners at a location slightly downstream from the primary inlet.

The first PSR (reactor #1) represents the *mixing zone* around the main injector
where the cool *premixed fuel-air* stream and the *primary air* stream are preheated
by mixing with the hot combustion products from the *recirculation zone*.

The configurations of a "coupled" ERN (the same as the reactor network when you use the Chemkin GUI)
has one major difference in requirements. The "hybrid" ERN in PyChemkin allows all PSRs to have an external
outlet. However, for the "coupled" reactor network in PyChemkin and in the Chemkin GUI, ONLY the "last"
PSR of the network can have an external outlet. Since the "recirculation zone* does not have any
downstream outflow (through flow), it cannot be the "last" reactor of a "coupled" ERN. Consequently,
even both the "PSRnetwork" and the current "PSRnetwork_coupled" represent exactly the same ERN,
their configurations are different. In the current "coupled" ERN configuration, the *recirculation zone*
becomes PSR #2 and the *flame zone* becomes the last reactor of the ENR, PSR #3.

Downstream from PSR #1, PSR #3, the *flame zone* is where the combustion of
the heated fuel-air mixture takes place. The secondary air is injected here to cool down
the combustion exhaust before it exits the combustion chamber. A portion of the exhaust gas
coming out of PSR #3 does not leave the combustion chamber directly and is diverted to PSR #2,
the *recirculation zone*. The external outlet of the ERN must be the direct downstream from
this reactor.

The majority of the outlet flow from PSR #2 is recirculated back to the flame zone to
sustain the fuel-lean premixed flame there. The rest of the hot gas from PSR #2 travels
further back to PSR #1, the mixing zone, to preheat the fuel-lean mixture just entering
the combustion chamber. Finally, the cooled flue gas leaves the chamber in a stream-like manner.
Typically, a PFR is applied to simulation of the outflow.

The reactors in the "coupled" network are solved simultaneously. Therefore, there is
no need to define any *tear stream*.
"""

# sphinx_gallery_thumbnail_path = '_static/combustor_ERN.png'

################################################
# Import PyChemkin packages and start the logger
# =======================================+======

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
interactive = True

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
# Create the "fuel" and the "air" mixtures to initialize the external
# inlet streams. The fuel for this case is pure methane.

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
# external inlet streams. The fuel in this case is pure methane. The
# ``X_by_Equivalence_Ratio()`` method is used to form the main
# ``premixed`` inlet stream to the *mixing zone* reactor. The external inlets
# to the first reactor, the *mixing zone*, and the second reactor, the *flame zone*, are simply
# ``air`` mixtures with different mass flow rates and temperatures.
#
# .. note::
#   PyChemkin has *air* redefined as a convenient way to set up the air
#   stream/mixture in the simulations. Use the ``ansys.chemkin.Air.X()`` or
#   ``ansys.chemkin.Air.Y()`` method when the mechanism uses ``O2`` and ``N2`` for
#   oxygen and nitrogen. Use the ``ansys.chemkin.air.X()`` or ``ansys.chemkin.air.Y()``
#   method when oxygen and nitrogen are represented by ``o2`` and ``n2``.
#

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

########################################
# Define reactors in the reactor network
# ======================================
# Set up the PSR for each zone one by one with *external inlets only*.
# For PSRs, use the ``set_inlet()`` method to add the external inlets to the
# reactor. A PFR always requires one external inlet when it is instantiated.
# From upstream to downstream, they are ``premix zone``, ``recirculation zone``,
# and ``flame zone``. The ``recirculation zone`` does not have an external inlet.
#
# .. note::
#   PyChemkin requires that the first reactor/zone must have at least
#   one external inlet. Because the rest of the reactors have at least the
#   through-flow from the immediate upstream reactor, they do not require
#   an external inlet.
#
# .. note::
#   The stream parameter used to instantiate a PSR is used to establish
#   the *guessed reactor solution* and is modified when the network is solved by
#   the ERN.
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

############################
# Create the reactor network
# ==========================
# Create a coupled reactor network named ``PSRcluster`` with the list of the PSRs
# in the correct order. For a simple chain network, you do not
# need to define the connectivity among the reactors. The reactor network model
# automatically figures out the through-flow connections.
#
# .. note::
#
#   - The order of the reactor in the list is important as it dictates the solution
#     sequence and thus the convergence rate.
#

# instantiate the PSR network as a coupled reactor network
PSR_list = [mix, recirculation, flame]
PSRcluster = ERN(PSR_list, label="combustor_cluster")

#############################
# Define the PSR connectivity
# ===========================
# Because the current reactor network contains recycling streams, for example, from PSR #3 to
# PSR #1 and PSR #2, you must explicitly specify the network connectivity. To do this, you use
# a *split* dictionary to define the outflow connections among the PSRs in the network.
# The *key* of the split dictionary is the label of the *originate PSR*. The *value* is a
# list of tuples consisting of the label of the *target PSR* and its mass flow rate fraction.
#
# ::
#
#   { "originate PSR" : [("target PSR1", fraction), ("target PSR2", fraction)],
#     "another originate PSR" : [("target PSR1", fraction), ... ], ... }
#
# For example, the outflow from reactor ``PSR1``is split into two streams. 90% of the
# mass flow rate goes to ``PSR2`` and the rest is diverted to ``PSR5``. The split dictionary entry for
# the outflow split of ``PSR1`` is as follows:
#
# ::
#
#   "PSR1" : [("PSR2", 0.9), ("PSR5", 0.1)]
#
# .. note::
#   The outlet flow from a reactor that is leaving the reactor network must be labeled as ``EXIT>>``
#   when you define the outflow splitting.
#

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

###########################
# Solve the reactor network
# =========================
# Use the ``run()`` method to solve the entire reactor network. The ``PSRCluster``
# solves the PSRs simultaneously.
#

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

###########################################
# Postprocess the reactor network solutions
# =========================================
# Since the reactors in the same cluster are solved together,
# you must get use the ``process_cluster_solution()`` method first,
# and extract the solution stream of a specific PSR outlet by using the
# ``get_reactor_stream()`` method with the (1-based) reactor index.
# Once you get the solution as a stream, you can use any stream or mixture
# method to further manipulate the solutions.
#

# postprocess the solutions of all PSRs in the cluster
iErr = PSRcluster.process_cluster_solution()

# verify the mass flow rate in and out of the PSR cluster
print(f"net external inlet mass flow rate = {PSRcluster.total_inlet_mass_flow_rate} [g/sec].")
print(f"net outlet mass flow rate = {PSRcluster.get_cluster_outlet_flowrate()} [g/sec].")

# display the reactor solutions
print("=" * 10)
print("reactor/zone")
print("=" * 10)
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

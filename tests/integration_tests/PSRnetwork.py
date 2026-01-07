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

"""Test of a hybrid reactor network model with stream recycling."""

###############################################
# Import PyChemkin package and start the logger
# =============================================

from pathlib import Path
import time

import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.hybridreactornetwork import ReactorNetwork as Ern
from ansys.chemkin.core.inlet import (
    Mixture,
    Stream,  # external gaseous inlet
)
from ansys.chemkin.core.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.core.stirreactors.PSR import PSRSetResTimeEnergyConservation as Psr

# check working directory
current_dir = str(Path.cwd())
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
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")

############################################
# Pre-process the gasoline ``Chemistry Set``
# ==========================================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()

####################################################################
# Set up gas mixtures based on the species in this ``Chemistry Set``
# ==================================================================
# Create the "fuel" and the "air" mixtures to initialize the external
# inlet streams. The "fuel" for this case is pure methane.

# fuel is pure methane
fuel = Mixture(MyGasMech)
fuel.temperature = 650.0  # [K]
fuel.pressure = 10.0 * ck.P_ATM  # [atm] => [dyne/cm2]
fuel.x = [("CH4", 1.0)]

# air is modeled as a mixture of oxygen and nitrogen
air = Mixture(MyGasMech)
air.temperature = 650.0  # [K]
air.pressure = 10.0 * ck.P_ATM
air.x = ck.Air.x()  # mole fractions

#################################################
# Create external inlet streams from the mixtures
# ===============================================
# Create the ``fuel`` and the ``air`` streams/mixtures before setting up the
# external inlet streams. The "fuel" in this case is pure methane. The main
# inlet stream ``premixed`` to the ``mixing zone`` is formed by using the
# ``x_by_equivalence_ratio`` method. The external inlets to the first reactor,
# the *"mixing zone"*, and the second reactor, the *"flame zone"*, are simply
# ``air`` mixtures with different mass flow rates and temperatures.
#
# .. note::
#   PyChemkin has *"air"* redefined as a convenient way to set up the air
#   stream/mixture in the simulations. Use ``ansys.chemkin.core.Air.x('U')`` or
#   ``ansys.chemkin.core.Air.Y('U')`` when the mechanism uses "O2" and "N2" for
#   oxygen and nitrogen. Use ``ansys.chemkin.core.Air.x('L')`` or
#   ``ansys.chemkin.core.Air.Y('L')`` when oxygen and nitrogen are represented
#   by "o2" and "n2".
#

# primary fuel-air mixture
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros

# create the unburned fuel-air mixture
premixed = Stream(MyGasMech)
# mean equivalence ratio
equiv = 0.6
ierror = premixed.x_by_equivalence_ratio(
    MyGasMech, fuel.x, air.x, add_frac, products, equivalenceratio=equiv
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
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
primary_air.x = air.x
primary_air.pressure = air.pressure
primary_air.temperature = air.temperature
primary_air.mass_flowrate = 50.0  # [g/sec]

# secondary bypass air stream
secondary_air = Stream(MyGasMech, label="Secondary_Air")
secondary_air.x = air.x
secondary_air.pressure = air.pressure
secondary_air.temperature = 670.0  # [K]
secondary_air.mass_flowrate = 100.0  # [g/sec]

# find the species index
ch4_index = MyGasMech.get_specindex("CH4")
o2_index = MyGasMech.get_specindex("O2")
no_index = MyGasMech.get_specindex("NO")
co_index = MyGasMech.get_specindex("CO")

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
mix = Psr(premixed, label="mixing zone")
# use different guess temperature
mix.set_estimate_conditions(option="TP", guess_temp=800.0)
# set PSR residence time (sec): required for PSRSetResTimeEnergyConservation model
mix.residence_time = 0.5 * 1.0e-3
# add external inlets
mix.set_inlet(premixed)
mix.set_inlet(primary_air)


# PSR #2: flame zone
flame = Psr(premixed, label="flame zone")
# use use the equilibrium state of the inlet gas mixture as the guessed solution
flame.set_estimate_conditions(option="TP", guess_temp=1600.0)
# set PSR residence time (sec): required for PSRSetResTimeEnergyConservation model
flame.residence_time = 1.5 * 1.0e-3
# add external inlet
flame.set_inlet(secondary_air)

# PSR #3: recirculation zone
recirculation = Psr(premixed, label="recirculation zone")
# use use the equilibrium state of the inlet gas mixture as the guessed solution
recirculation.set_estimate_conditions(option="TP", guess_temp=1600.0)
# set PSR residence time (sec): required for PSRSetResTimeEnergyConservation model
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
PSRnetwork = Ern(MyGasMech)

# add the reactors from upstream to downstream
PSRnetwork.add_reactor(mix)
PSRnetwork.add_reactor(flame)
PSRnetwork.add_reactor(recirculation)

# list the reactors in the network
PSRnetwork.show_reactors()

##################################
# Define the PSR connectivity
# ==============================
# Because the current reactor network contains recycling streams, for example,
# from PSR #3 to PSR #1 and PSR #2, the network connectivity must be specified
# explicitly. This is done by using a *"split"* dictionary to define
# the outflow connections among the PSRs in the network.
# The *key* of the *"split"* dictionary is the label of the *"originate PSR"*.
# The *value* is a list of tuples consisting of the label of the *"target PSR"*
# and its *"mass flow rate fraction"*.
#
# ::
#
#   { "originate PSR" : [("target PSR1", fraction), ("target PSR2", fraction)],
#     "another originate PSR" : [("target PSR1", fraction), ... ], ... }
#
# For example, the outflow from reactor "PSR1" is split into two streams.
# 90% of the mass flow rate goes to "PSR2" and the rest is diverted to "PSR5".
# The split dictionary entry for the outflow split of "PSR1" becomes
#
# ::
#
#   "PSR1" : [("PSR2", 0.9), ("PSR5", 0.1)]
#
# .. note::
#   The outlet flow from a reactor that is leaving the reactor network must be
#   labeled as "EXIT>>" when you define the outflow splitting.
#

# PSR #1 outlet flow splitting
split_table = [(flame.label, 1.0)]
PSRnetwork.add_outflow_connections(mix.label, split_table)
# PSR #2 outlet flow splitting
# part of the outlet flow from PSR #2 exits the reactor network
split_table = [(recirculation.label, 0.2), ("EXIT>>", 0.8)]
PSRnetwork.add_outflow_connections(flame.label, split_table)
# PSR #3 outlet flow splitting
# PSR #3 is the last reactor of the network,
# however it does not have external outlet
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
# Use the ``run`` method to solve the entire reactor network.
# The hybrid ``ReactorNetwork`` will solve the reactors one by one
# in the order they are added to the network.
#
# .. note::
#   Use the ``set_tear_iteration_limit(count)`` method to change the limit
#   on the number of tear loop tear loop iterations can be taken before
#   declaring the failure. The default limit is 200.
#
# .. note::
#   The ``set_relaxation_factor(factor)`` method can be employed to make
#   solution updating at the end of each iteration to be more aggressive
#   (factor > 1.0) or more conservative (factor < 1.0). By default,
#   the relaxation factor is set to 1.
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
    print(f"CH4 = {network_outflow.x[ch4_index]}")
    print(f"O2 = {network_outflow.x[o2_index]}")
    print(f"CO = {network_outflow.x[co_index]}")
    print(f"NO = {network_outflow.x[no_index]}")
    print("-" * 10)

# display the raector solutions
print()
print("=" * 10)
print("reactor/zone")
print("=" * 10)
temp = []
mflr = []
x_ch4 = []
x_co = []
x_no = []
for index, stream in PSRnetwork.reactor_solutions.items():
    name = PSRnetwork.get_reactor_label(index)
    print(f"reactor: {name}")
    print(f"temperature = {stream.temperature} [K]")
    print(f"mass flow rate = {stream.mass_flowrate} [g/sec]")
    print(f"CH4 = {stream.x[ch4_index]}")
    print(f"O2 = {stream.x[o2_index]}")
    print(f"CO = {stream.x[co_index]}")
    print(f"NO = {stream.x[no_index]}")
    # save results to lists for comparisons
    temp.append(stream.temperature)
    mflr.append(stream.mass_flowrate)
    x_ch4.append(stream.x[ch4_index])
    x_co.append(stream.x[co_index])
    x_no.append(stream.x[no_index])
    print("-" * 10)

# return results for comparisons
resultfile = Path(current_dir) / "PSRnetwork.result"
results = {}
results["state-temperature"] = temp
results["state-mass_flow_rate"] = mflr
results["species-mole_fraction_CH4"] = x_ch4
results["species-mole_fraction_CO"] = x_co
results["species-mole_fraction_NO"] = x_no
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

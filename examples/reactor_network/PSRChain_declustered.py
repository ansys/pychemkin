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

r""".. _ref_connected_reactors:

===========================================================
Use a chain of individual reactors to model a gas combustor
===========================================================

This example shows how to set up and solve a series of linked PSRs
(perfectly-stirred reactors). This is the simplest reactor network as
it does not contain any recycling streams or outflow splittings.

Here is the PSR chain model of a fictional gas combustor.

 .. figure:: chain_reactor_network.png
   :scale: 80 %
   :alt: Chain reactor network

The primary inlet stream to the first reactor, the *combustor*, is the fuel-lean
methane-air mixture that is formed by mixing the fuel (methane) and the heated air.
The exhaust from the combustor enters the second reactor, the *dilution zone*,
where the hot combustion products are cooled by the introduction of
additional cool air. The cooled and diluted gas mixture in the *dilution zone*
then travels to the third reactor, the *reburning zone*. A mixture of fuel (methane)
and carbon dioxide is injected to the gas in the reburning zone, attempting to
convert any remaining carbon monoxide or nitric oxide in the exhaust gas to
carbon dioxide or nitrogen, respectively.

This example solves the reactors one by one, from upstream to downstream.
Once the solution of the upstream reactor is obtained, it is used to set up
the external inlet of the immediate downstream reactor. This process continues
until all reactors in the chain network are solved. Since there is no recycling stream
in this configuration, the entire reactor network can be solved in one sweep.
"""

# sphinx_gallery_thumbnail_path = '_static/chain_reactor_network.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path
import time

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream  # external gaseous inlet
from ansys.chemkin.core.inlet import adiabatic_mixing_streams
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

########################
# Create a chemistry set
# ======================
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
MyGasMech.chemfile = str((mechanism_dir / "grimech30_chem.inp").resolve())
MyGasMech.thermfile = str((mechanism_dir / "grimech30_thermo.dat").resolve())

#######################################
# Preprocess the gasoline chemistry set
# =====================================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()

################################################################
# Set up gas mixtures based on the species in this chemistry set
# ==============================================================
# Create the ``fuel`` and ``air`` mixtures to initialize the external
# inlet streams. The fuel for this case is pure methane.

# fuel is pure methane
fuel = Stream(MyGasMech)
fuel.temperature = 300.0  # [K]
fuel.pressure = 2.1 * ck.P_ATM  # [atm] => [dyne/cm2]
fuel.x = [("CH4", 1.0)]
fuel.mass_flowrate = 3.275  # [g/sec]

# air is modeled as a mixture of oxygen and nitrogen
air = Stream(MyGasMech)
air.temperature = 550.0  # [K]
air.pressure = 2.1 * ck.P_ATM
air.x = ck.Air.x()  # use predefined "air" recipe in mole fractions
air.mass_flowrate = 45.0  # [g/sec]

#################################################
# Create external inlet streams from the mixtures
# ===============================================
# Use the stream ``adiabatic_mixing_streams()`` method to combine
# the ``fuel`` and ``air`` streams. The final gas temperature should
# land between the temperatures of the two source streams. The mass flow
# rate of the ``premixed`` stream should be the sum of the sources. A simple
# PyChemkin composition recipe is used to create the ``reburn_fuel`` stream.

# premixed stream for the combustor
premixed = adiabatic_mixing_streams(fuel, air)
print(str(premixed.temperature))
print(str(premixed.mass_flowrate))

# additional fuel injection for the reburning zone
reburn_fuel = Stream(MyGasMech)
reburn_fuel.temperature = 300.0  # [K]
reburn_fuel.pressure = 2.1 * ck.P_ATM  # [atm] => [dyne/cm2]
reburn_fuel.x = [("CH4", 0.6), ("CO2", 0.4)]
reburn_fuel.mass_flowrate = 0.12  # [g/sec]

# find the species index
ch4_index = MyGasMech.get_specindex("CH4")
o2_index = MyGasMech.get_specindex("O2")
no_index = MyGasMech.get_specindex("NO")
co_index = MyGasMech.get_specindex("CO")

###########################
# Create PSRs for each zone
# =========================
# Set up the PSR for each zone one by one with external inlets only. PyChemkin
# requires that the first reactor/zone must have at least one external inlet.
# There are three reactors in the network. From upstream to downstream, they
# are ``combustor``, ``dilution zone``, and ``reburning zone``. All of them
# have one external inlet. For the two downstream reactors, you must create
# the through-flow stream from their respective upstream reactor by using the
# solution of the upstream reactor. Use the ``process_solution()`` method to
# get the solution stream from the upstream reactor. Then, use the ``set_inlet()``
# method to connect the resulting solution stream to the downstream reactor.
#
# .. note::
#
#   - PyChemkin requires that the first reactor/zone must have at least
#     one external inlet. The rest of the reactors have at least the
#     through-flow from the immediate upstream reactor so they do not require
#     an external inlet.
#
#   - The ``Stream`` parameter used to instantiate a PSR object is used to establish
#     the *guessed reactor solution* and is modified when the network is solved by
#     the ``ERN``.
#
#   - The reactors in the network must be postprocessed individually.
#

# PSR #1: combustor
combustor = Psr(premixed, label="combustor")
# use the equilibrium state of the inlet gas mixture as the guessed solution
combustor.set_estimate_conditions(option="HP")
# set PSR residence time (sec): required for PSRSetResTimeEnergyConservation model
combustor.residence_time = 2.0 * 1.0e-3
# add external inlet
combustor.set_inlet(premixed)
# set the start wall time
start_time = time.time()
# run PSR #1
status = combustor.run()
if status != 0:
    print(Color.RED + combustor.label + " Run failed.")
    exit()
# postprocess the solution profiles
solnstream1 = combustor.process_solution()
solnstream1.label = "PSR1"
print("=" * 40)
print("Combustor exited.")
print("=" * 40)
print(f"Temperature = {solnstream1.temperature} [K].")
print(f"Mass flow rate = {solnstream1.mass_flowrate} [g/sec].")
print(f"CH4 = {solnstream1.x[ch4_index]}.")
print(f"O2 = {solnstream1.x[o2_index]}.")
print(f"CO = {solnstream1.x[co_index]}.")
print(f"NO = {solnstream1.x[no_index]}.")

# PSR #2: cooling
cooling = Psr(solnstream1, label="cooling zone")
# set PSR residence time (sec): required for PSRSetResTimeEnergyConservation model
cooling.residence_time = 1.5 * 1.0e-3
# add external inlet
air.mass_flowrate = 62.0  # [g/sec]
cooling.set_inlet(air)
# add the through-flow from PSR #1
cooling.set_inlet(solnstream1)
# run PSR #2
status = cooling.run()
if status != 0:
    print(Color.RED + cooling.label + " Run failed.")
    exit()
# post-process the solution profiles
solnstream2 = cooling.process_solution()
solnstream2.label = "PSR2"
print()
print("=" * 40)
print("Dilution zone exited.")
print("=" * 40)
print(f"Temperature = {solnstream2.temperature} [K].")
print(f"Mass flow rate = {solnstream2.mass_flowrate} [g/sec].")
print(f"CH4 = {solnstream2.x[ch4_index]}.")
print(f"O2 = {solnstream2.x[o2_index]}.")
print(f"CO = {solnstream2.x[co_index]}.")
print(f"NO = {solnstream2.x[no_index]}.")

# PSR #3: reburn
reburn = Psr(solnstream2, label="reburn zone")
# set PSR residence time (sec): required for PSRSetResTimeEnergyConservation model
reburn.residence_time = 3.5 * 1.0e-3
# add external inlet
reburn.set_inlet(reburn_fuel)
# add the through flow from PSR #2
reburn.set_inlet(solnstream2)
# run PSR #3
status = reburn.run()
if status != 0:
    print(Color.RED + reburn.label + " Run failed.")
    exit()
# post-process the solution profiles
outflow = reburn.process_solution()
outflow.label = "outflow"
print()
print("=" * 40)
print("Outflow exited.")
print("=" * 40)
print(f"Temperature = {outflow.temperature} [K].")
print(f"Mass flow rate = {outflow.mass_flowrate} [g/sec].")
print(f"CH4 = {outflow.x[ch4_index]}.")
print(f"O2 = {outflow.x[o2_index]}.")
print(f"CO = {outflow.x[co_index]}.")
print(f"NO = {outflow.x[no_index]}.")

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"Total simulation duration: {runtime} [sec].")

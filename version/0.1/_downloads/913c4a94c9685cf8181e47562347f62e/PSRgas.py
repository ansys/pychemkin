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

r""".. _ref_gas_PSR:

===================================================================
Set up a PSR parameter study for the inlet stream equivalence ratio
===================================================================

Ansys Chemkin offers some idealized reactor models commonly used for studying chemical
processes and for developing reaction mechanisms. The PSR (perfectly stirred reactor)
model is a steady-state 0-D model of the open perfectly mixed gas-phase reactor.
There is no limit on the number of inlets to the PSR. As soon as the inlet gases
enter the reactor, they are thoroughly mixed with the gas mixture inside. The PSR has
only one outlet, and the outlet gas is assumed to be exactly the same as
the gas mixture in the PSR.

There are two basic types of PSR models:

- **constrained-pressure** (or set residence time)
- **constrained-volume**

By default, the PSR model is running under constant pressure. PyChemkin PSR models
always require the connected inlets to be defined, that is, the total inlet flow rate
to the PSR is always known. Therefore, either the residence time or
the reactor volume is needed to satisfy the basic setup of the PSR model.
In this case, you specify the residence time of the PSR, and the PSR model
automatically calculates the reactor volume from the given residence time and
the total inlet volumetric flow rate.

For each type of the PSR, you can choose either to specify the reactor temperature
(as a fixed value or by a piecewise-linear profile) or to solve the energy conservation
equation. In total, you get four variations of the PSR model.

The PSR model is mostly employed in chemical kinetics studies. By controlling
the reactor temperature, pressure, and/or residence time, you can gain knowledge
about the major intermediates of a complex chemical process and postulate possible
reaction pathways. The parameter study in this example shows how the inlet
equivalence ratio impacts the hydrogen combustion process in a
fixed residence time PSR.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_PSR_gas.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream  # external gaseous inlet
from ansys.chemkin.core.logger import logger

# Chemkin PSR model (steady-state)
from ansys.chemkin.core.stirreactors.PSR import PSRSetResTimeEnergyConservation as Psr
from ansys.chemkin.core.utilities import find_file

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
# This example uses the encrypted hydrogen-ammonia mechanism,
# ``Hydrogen-Ammonia-NOx_chem_MFL2021.inp``. This mechanism is developed
# under Chemkin's *Model Fuel Library* (MFL) project.
# Like the rest of the MFL mechanisms, it is in the *ModelFuelLibrary* in the
# ``/reaction/data`` directory of the standard Ansys Chemkin installation.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data" / "ModelFuelLibrary" / "Skeletal"
mechanism_dir = str(data_dir.resolve())
# create a chemistry set based on the gasoline 14-components mechanism
MyGasMech = ck.Chemistry(label="hydrogen")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = find_file(
    mechanism_dir,
    "Hydrogen-Ammonia-NOx_chem_MFL",
    "inp",
)

#######################################
# Preprocess the gasoline chemistry set
# =====================================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()

#####################################
# Set up the H\ :sub:`2`\ -air stream
# ===================================
# Instantiate a stream feed for the inlet gas mixture.
# The stream is a mixture with the addition of the
# inlet flow rate. You specify inlet gas properties in the same way as you
# set up a mixture. Here, the ``x_by_equivalence_ratio()`` method is used.
#
# First create the fuel and air mixtures. Then, define the
# complete combustion product species and provide the additives composition
# if applicable. Finally, set the equivalence ratio to 1 to create
# the stoichiometric hydrogen-air mixture. Use the ``mass_flowrate()`` method
# to assign the inlet mass flow rate. Use a fixed inlet mass flow rate
# of 432 [g/sec] since you are to change the PSR residence time
# in the parameter study.

# create the fuel mixture
fuel = ck.Mixture(MyGasMech)
# set fuel composition: hydrogen diluted by nitrogen
fuel.x = [("h2", 0.8), ("n2", 0.2)]
# setting the pressure and temperature is not required in this case
fuel.pressure = ck.P_ATM
fuel.temperature = 298.0  # inlet temperature

# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.x = [("o2", 0.21), ("n2", 0.79)]
# setting the pressure and temperature is not required in this case
air.pressure = fuel.pressure
air.temperature = fuel.temperature

# create the fuel-oxidizer inlet to the PSR
feed = Stream(MyGasMech, label="feed_1")
# products from the complete combustion of the fuel mixture and air
products = ["h2o", "n2"]
# species mole fractions of added/inert mixture.
# You can also create an additives mixture here.
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros
# mean equivalence ratio
equiv = 1.0
ierror = feed.x_by_equivalence_ratio(
    MyGasMech, fuel.x, air.x, add_frac, products, equivalenceratio=equiv
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
    print("Error: Failed to create the fuel-oxidizer mixture.")
    exit()

# set reactor pressure [dynes/cm2]
feed.pressure = fuel.pressure
# set inlet gas temperature [K]
feed.temperature = fuel.temperature
# set inlet mass flow rate [g/sec]
feed.mass_flowrate = 432.0

####################################################################
# Create the PSR to predict the gas composition of the outlet stream
# ===================================================================
# Use the ``PSRSetResTimeEnergyConservation()`` method to instantiate
# the PSR ``sphere`` object because the goal is to see how
# the residence time affects the hydrogen combustion process.
# The gas property of the inlet feed is applied as the estimated
# reactor condition of the ``sphere`` object by default.
# You can overwrite any estimated reactor conditions by
# using appropriate methods. For example, ``sphere.temperature = 1700.0``
# changes the estimated reactor temperature from 298 to 1700 [K].
# The residence time of the nominal case is set by
# the ``residence_time()`` method.
sphere = Psr(feed, label="PSR_1")

############################################
# Set up additional reactor model parameters
# ==========================================
# Before you can run the simulation, you must provide reactor parameters,
# solver controls, and output instructions. For a steady-state PSR,
# you must provide either the residence time or reactor volume. You can
# also make changes to any estimated reactor conditions if desired.

# reset the estimated reactor temperature [K]
sphere.temperature = 1700.0
# set PSR residence time (sec): required for PSRSetResTimeEnergyConservation model
sphere.residence_time = 3.0 * 1.0e-5

##################################
# Connect the inlet to the reactor
# ================================
# You must connect at least one inlet to the open reactor. Use the ``set_inlet()``
# method to add a stream to the PSR. Inversely, use the ``remove_inlet()`` to
# disconnect an inlet from the PSR.
#
# .. note ::
#   There is no limit on the number of inlets that can be connected to a PSR.
#

# connect the inlet to the reactor
sphere.set_inlet(feed)

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances. Here the tolerances that the steady-state solver
# is to use for the *steady-state search* and the *pseudo time stepping* stages
# are changed. Sometimes, during the iterations, some species mass fractions might
# become negative, causing the solver to report an error and stop. To overcome
# this issue, you can use the ``set_species_floor()`` method to provide a
# small cushion, allowing species mass fractions to go slightly negative by resetting
# the mass fraction floor value.

# reset the tolerances in the steady-state solver
sphere.steady_state_tolerances = (1.0e-9, 1.0e-6)
sphere.timestepping_tolerances = (1.0e-9, 1.0e-6)
# reset the gas species floor value in the steady-state solver
sphere.set_species_floor(-1.0e-10)

#################################################
# Run the inlet equivalence ratio parameter study
# ===============================================
# In the parameter study, the equivalence ratio :math:`\phi` of the inlet gas mixture
# is increased from 1.0 to 1.4. This is done by applying the
# ``x_by_equivalence_ratio()`` method on the feed inlet. Changing the equivalence ratio
# of the inlet gas mixture inevitably has impact on the inlet stream density and
# hence the inlet volumetric flow rate. By using the
# ``PSRSetResTimeEnergyConservation`` model, the PSR residence time remains
# constant for all runs. The effects from any variations of the inlet are
# reflected on the PSR volume.
#
# Use the ``process_solution()`` method to convert the result from each PSR to
# a mixture. You can either overwrite the solution mixture or use a new one for
# each simulation result.

# inlet gas equivalence ratio increment
deltaequiv = 0.05
numbruns = 9
# solution arrays
inletequiv = np.zeros(numbruns, dtype=np.double)
temp_ss_solution = np.zeros_like(inletequiv, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all inlet temperature values
for i in range(numbruns):
    # run the PSR model
    runstatus = sphere.run()
    # check run status
    if runstatus != 0:
        # Run failed.
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()
    # Run succeeded.
    print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
    # postprocess the solution profiles
    solnmixture = sphere.process_solution()
    # print the steady-state solution values
    # print(f"Steady-state temperature = {solnmixture.temperature} [K].")
    # solnmixture.list_composition(mode="mole")
    # store solution values
    inletequiv[i] = equiv
    temp_ss_solution[i] = solnmixture.temperature
    # update inlet gas equivalence ratio (composition)
    equiv += deltaequiv
    ierror = feed.x_by_equivalence_ratio(
        MyGasMech, fuel.x, air.x, add_frac, products, equivalenceratio=equiv
    )
    # check fuel-oxidizer mixture creation status
    if ierror != 0:
        print(f"Error encountered with inlet equivalence ratio = {equiv}.")
        exit()

# compute the total runtime
runtime = time.time() - start_time
print(f"Total simulation duration: {runtime} [sec] over {numbruns} runs.")

##################################
# Plot the parameter study results
# ================================
# Plot the steady-state PSR temperature against the equivalence ratio of the
# inlet H\ :sub:`2`\ -air mixture. You should see that the maximum combustion
# temperature does not correspond to the :math:`\phi = 1` mixture. Instead,
# the temperature peak occurs when the mixture is slightly fuel rich. You can
# run the same parameter study on a different fuel species such as CH\ :sub:`4`
# to see if you observe the same behavior.
plt.plot(inletequiv, temp_ss_solution, "b-")
plt.xlabel("Inlet Gas Equivalence Ratio")
plt.ylabel("Reactor Temperature [K]")
plt.title("PSR Solution")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_PSR_gas.png", bbox_inches="tight")

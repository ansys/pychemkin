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
.. _ref_jet_stirred_reactor:

=============================================
Run PSR calculations for mechanism validation
=============================================

Ansys Chemkin offers some idealized reactor models commonly used for studying chemical
processes and for developing reaction mechanisms. The PSR (perfectly stirred reactor) model is
a steady-state 0-D model of the open perfectly mixed gas-phase reactor. There is no limit on
the number of inlets to the PSR. As soon as the inlet gases enter the reactor, they are
thoroughly mixed with the gas mixture inside. The PSR has only one outlet, and the outlet gas
is assumed to be exactly the same as the gas mixture in the PSR.

There are two basic types of PSR models:

- **constrained-pressure** (or set residence time)
- **constrained-volume**

By default, the PSR is running under constant pressure. In this case, you specify the
residence time of the PSR. For the constrained-volume type of application, you must provide
the PSR volume. You can calculate the residence time from the reactor volume and the total
inlet volumetric flow rate.

For each type of PSR, you either specify the reactor temperature (as a fixed
value or by a piecewise-linear profile) or solve the energy conservation equation. In total,
you get four variations of the PSR model.

The JSR (jet-stirred reactor is mostly employed in chemical kinetics studies. By controlling the reactor temperature, pressure, and/or residence time, you can gain knowledge about the major intermediates of a complex chemical process and postulate possible reaction pathways.

This example shows how to use the PSR model to validate the reaction mechanism against the measured data from hydrogen oxidation experiments.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_jet_stirred_reactor.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger

# chemkin perfectly stirred reactor (PSR) model (steady-state)
from ansys.chemkin.stirreactors.PSR import PSR_SetResTime_FixedTemperature as PSR
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
interactive = True

########################
# Create a chemistry set
# ======================
# This code uses the encrypted hydrogen-ammonia mechanism, ``Hydrogen-Ammonia-NOx_chem_MFL2021.inp``.
# This mechanism is developed under Chemkin's **Model Fuel Library (MFL)** project.
# Like the rest of the MFL mechanisms, it is in the ``ModelFuelLibrary`` in the
# ``/reaction/data`` directory of the standard Ansys Chemkin installation.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(
    ck.ansys_dir, "reaction", "data", "ModelFuelLibrary", "Skeletal"
)
mechanism_dir = data_dir
# create a chemistry set based on the hydrogen-ammonia mechanism
MyGasMech = ck.Chemistry(label="hydrogen")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(
    mechanism_dir, "Hydrogen-Ammonia-NOx_chem_MFL2021.inp"
)

#######################################
# Preprocess the gasoline chemistry set
# =====================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

###########################################################
# Set up the H\ :sub:`2`\ -O\ :sub:`2`\ -N\ :sub:`2` stream
# =========================================================
# Instantiate a stream feed for the inlet gas mixture.
# A stream is a mixture with the addition of the
# inlet flow rate. You set up the inlet gas properties in the same way you
# set up a mixture. Set up the gas properties of the feed as described by the experiment.
#
# Use the ``mass_flowrate()`` method to assign the inlet mass flow rate. The experiments use
# a fixed inlet mass flow rate of 0.11 [g/sec].

# create the fuel-oxidizer inlet to the JSR
feed = Stream(MyGasMech)
# set H2-O2-N2 composition
feed.X = [("h2", 1.1e-2), ("n2", 9.62e-1), ("o2", 2.75e-2)]
# set reactor pressure [dynes/cm2]
feed.pressure = ck.Patm
# set inlet gas temperature [K]
temp = 800.0
feed.temperature = temp
# set inlet mass flow rate [g/sec]
feed.mass_flowrate = 0.11

####################################################################
# Create the PSR to predict the gas composition of the outlet stream
# ==================================================================
# Use the ``PSR_SetResTime_FixedTemperature()`` method to instantiate the JSR
# because both the reactor temperature and residence time are fixed during the experiments.
# The gas property of the inlet feed is applied as the estimated reactor condition
# of the JSR.
JSR = PSR(feed, label="JSR")

###################################
# Connect the inlets to the reactor
# =================================
# You must connect at least one inlet to the open reactor. Use the ``set_inlet()`` method to
# add a stream to the PSR. Inversely, use the ``remove_inlet()`` method to disconnect
# an inlet from the PSR.

# connect the inlet to the reactor
JSR.set_inlet(feed)

############################################
# Set up additional reactor model parameters
# ==========================================
# You must provide reactor parameters, solver controls, and output instructions
# before running the simulations. For the steady-state PSR, you must provide either the residence
# time or the reactor volume.

# set PSR residence time (sec): required for PSR_SetResTime_FixedTemperature model
JSR.residence_time = 120.0 * 1.0e-3

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods, such as
# those for tolerances. Here, particular to this mechanism, a number of *pseudo timesteps* is required
# before attempting to actually search for the steady-state solution. This is done by using the
# steady-state solver control method, ``set_initial_timesteps()``.

# set the number of initial pseudo timesteps in the steady-state solver
JSR.set_initial_timesteps(1000)

######################################################
# Run the parameter study to replicate the experiments
# ====================================================
# Different inlet and reactor temperatures are used in the experiments. The temperature
# value varies from 800 to 1050 [K].
#
# Use the ``process_solution()`` method to convert the result from each PSR run to a
# mixture. You can either overwrite the solution mixture or use a new one for each simulation result.

# inlet gas temperature increment
deltatemp = 25.0
numbruns = 19
# find "h2o" species index
H2Oindex = MyGasMech.get_specindex("h2o")
# solution arrays
inletTemp = np.zeros(numbruns, dtype=np.double)
h2oSSsolution = np.zeros_like(inletTemp, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all inlet temperature values
for i in range(numbruns):
    # run the PSR model
    runstatus = JSR.run()
    # check run status
    if runstatus != 0:
        # Run failed.
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()
    # Run succeeded.
    print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
    # postprocess the solution profiles
    solnmixture = JSR.process_solution()
    # print the steady-state solution values
    # print(f"Steady-state temperature = {solnmixture.temperature} [K].")
    # solnmixture.list_composition(mode="mole")
    # store solution values
    inletTemp[i] = solnmixture.temperature
    h2oSSsolution[i] = solnmixture.X[H2Oindex]
    # update reactor temperature
    temp += deltatemp
    JSR.temperature = temp

# compute the total runtime
runtime = time.time() - start_time
print(f"Total simulation duration: {runtime} [sec] over {numbruns} runs")
#
# experimental data
# JSR temperature [K]
TEMP_data = [
    803.0,
    823.0,
    850.0,
    875.0,
    902.0,
    925.0,
    951.0,
    973.0,
    1002.0,
    1023.0,
    1048.0,
]
# measured H2O mole fractions in the JSR exit flow
H2O_data = [
    0.000312,
    0.000313,
    0.000318,
    0.000303,
    0.003292,
    0.006226,
    0.008414,
    0.009745,
    0.010103,
    0.011036,
    0.010503,
]

###############################
# Plot the ignition delay curve
# =============================
# Compare the predicted H\ :sub:`2`\ O mole fraction to the measurement.

# plot results
plt.plot(inletTemp, h2oSSsolution, "b-", label="prediction")
plt.plot(TEMP_data, H2O_data, "ro", label="data", markersize=4, fillstyle="none")
plt.xlabel("Reactor Temperature [K]")
plt.ylabel("H2O Mole Fraction")
plt.legend(loc="lower right")
plt.title("JSR Solution")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_jet_stirred_reactor.png", bbox_inches="tight")

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

r""".. _ref_ignition_delay:

========================================================
Predict the ignition delay time of a combustible mixture
========================================================

One of the important properties of hydrocarbon fuels is the *ignition delay time*.
For automobiles, the gasoline, a mixtures of many fuel species, is graded by the
*octane number* which is closely related to the fuel mixture's ignition
characteristics. Higher octane number fuel tends to ignite more easily. However,
this does not mean you should blindly put the highest octane number fuel
into your car. Car engines are designed for specific operating conditions,
and you should always use the gasoline grades recommended by the car/engine
manufacturer to prevent damages to the engine.

This tutorial employs the **constant-pressure batch reactor model**
(with the energy equation **ON**) to predict the *auto-ignition delay time* of
a **Primary Reference Fuel (PRF)**. PRF is created for automobile fuel researches
and consists of two major components: n-heptane C\ :sub:`7`\ H\ :sub:`16` and
iso-octane C\ :sub:`8`\ H\ :sub:`18`\ . By definition, n-heptane is assigned with
an octane number of 0, and an octane number of 100 for iso-octane. Thus, a PRF
of 20% n-heptane and 80% iso-octane has an octane number of 80. In this tutorial,
a PRF 60 fuel is used to show the **negative temperature coefficient (NTC)**
behavior in the ignition delay curve. The ignition delay curve is constructed by
collecting the predicted auto-ignition delay times with various
initial gas conditions. Here, the initial gas temperature is increased from
700 to 1080 [K].
"""

# sphinx_gallery_thumbnail_path = '_static/plot_ignition_delay.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # chemkin
from ansys.chemkin.core import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.core.batchreactors.batchreactor import (
    GivenPressureBatchReactorEnergyConservation,
)
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(False)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True

#####################################
# Create a chemistry set
# ===================================
# For PRF, the encrypted 14-component gasoline mechanism,
# ``gasoline_14comp_WBencrypted.inp``, is used. ``gasoline``
# is the name given to this ``Chemistry Set``.
#
# .. note::
#   This gasoline mechanism does not come with any transport data
#   so you do not need to provide the transport data file.
#

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
gasoline = ck.Chemistry(label="gasoline 14comp")
# set mechanism input files
# including the full file path is recommended
gasoline.chemfile = str(mechanism_dir / "gasoline_14comp_WBencrypt.inp")

############################################
# Pre-process the gasoline ``Chemistry Set``
# ==========================================

# preprocess the mechanism files
ierror = gasoline.preprocess()

################################################
# Set up the stoichiometric gasoline-air mixture
# ==============================================
# You need to set up the stoichiometric gasoline-air mixture for the subsequent
# ignition delay time calculations. Here the ``x_by_equivalence_ratio``method is used.
# You create the ``fuel`` and the ``air`` mixtures first. Then define the
# *complete combustion product species* and provide the *additives* composition
# if applicable. And finally you can simply set ``equivalenceratio=1`` to create
# the stoichiometric gasoline-air mixture.
#
# For PRF 60 gasoline, the ``recipe`` is ``[("ic8h18", 0.6), ("nc7h16", 0.4)]``.

# create the fuel mixture
fuelmixture = ck.Mixture(gasoline)
# set fuel = composition of PRF 60
fuelmixture.x = [("ic8h18", 0.6), ("nc7h16", 0.4)]
# setting pressure and temperature
fuelmixture.pressure = 5.0 * ck.P_ATM
fuelmixture.temperature = 1500.0

# create the oxidizer mixture: air
air = ck.Mixture(gasoline)
air.x = [("o2", 0.21), ("n2", 0.79)]
# setting pressure and temperature
air.pressure = 5.0 * ck.P_ATM
air.temperature = 1500.0

# products from the complete combustion of the fuel mixture and air
products = ["co2", "h2o", "n2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(gasoline.kk, dtype=np.double)  # no additives: all zeros

# create the premixed mixture to be defined
premixed = ck.Mixture(gasoline)
# define the composition by the equivalence ratio
ierror = premixed.x_by_equivalence_ratio(
    gasoline, fuelmixture.x, air.x, add_frac, products, equivalenceratio=1.0
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
    print("Error: Failed to create the fuel-oxidizer mixture.")
    exit()

##############################
# List the mixture composition
# ============================
# list the composition of the premixed mixture for verification
premixed.list_composition(mode="mole")
# set mixture temperature and pressure
# (equivalent to setting the initial temperature and pressure of the reactor)
premixed.pressure = 40.0 * ck.P_ATM
premixed.temperature = 700.0

################################################################
# Create the reactor object for ignition delay time calculations
# ==============================================================
# Use the ``GivenPressureBatchReactorEnergyConservation`` method to instantiate a
# *constant pressure batch reactor that also includes the energy equation*.
# This ``ReactorModel`` method requires a ``Mixture`` object as an input parameter.
# This ``Mixture`` input serves two purposes: passing the ``Chemistry Set`` to
# be used by the reactor model and setting the *initial condition*
# (pressure, temperature, and gas composition) of the simulation. You should use
# the ``premixed`` mixture you just created.

MyCONP = GivenPressureBatchReactorEnergyConservation(premixed, label="CONP")

############################################
# Set up additional reactor model parameters
# ==========================================
# *Reactor parameters*, *solver controls*, and *output instructions*
# need to be provided before running the simulations. For a batch reactor,
# the *initial volume* and the *simulation end time* are required inputs.

# verify initial gas composition inside the reactor
MyCONP.list_composition(mode="mole")

# set other reactor parameters
# initial reactor volume [cm3]
MyCONP.volume = 10.0
# simulation end time [sec]
MyCONP.time = 1.0

####################
# Set output options
# ==================
# You can turn on the *adaptive solution saving* to resolve the steep variations
# in the solution profile. Here additional solution data point will be saved for
# every **100 [K]** change in gas **temperature**. The ``set_ignition_delay`` method
# must be included for the reactor model to report the *ignition delay times* after
# the simulation is done. If ``method="T_inflection"`` is set, the reactor model will
# treat the *inflection points* in the predicted gas temperature profile as
# the indication of an auto-ignition. You can choose a different
# auto-ignition definition.
#
# .. note::
#   Type ``ansys.chemkin.core.show_ignition_definitions()`` to get the list of
#   all available ignition delay time definitions in Chemkin.
#
# .. note::
#   By default, time intervals for both print and save solution are **1/100**
#   of the *simulation end time*. In this case :math:`dt=time/100=0.001`\ .
#   You can change them to different values.
#

# set timestep between saving solution
MyCONP.timestep_for_saving_solution = 0.001
# change timestep between saving solution
MyCONP.timestep_for_saving_solution = 0.01
# turn ON adaptive solution saving
MyCONP.adaptive_solution_saving(mode=True, value_change=100, target="TEMPERATURE")
# set ignition delay
MyCONP.set_ignition_delay(method="T_inflection")

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver related methods,
# for example, ``tolerances``.

# tolerances are given in tuple: (absolute tolerance, relative tolerance)
MyCONP.tolerances = (1.0e-10, 1.0e-8)
# stop after ignition is detected
# (not recommended for ignition delay time calculations)
# MyCONP.stop_after_ignition()
# show solver option
print(f"timestep between solution printing: {MyCONP.timestep_for_printing_solution}")
# show timestep between printing solution
print(f"forced non-negative solution values: {MyCONP.force_nonnegative}")

#########################
# Run the parameter study
# =======================
# Construct the ignition delay curve by varying the initial reactor temperature
# in 20 [K] increments from 700 to 1080 [K]. Use the ``get_ignition_delay`` method
# to extract the ignition delay times after finishing each run.

# loop over initial reactor temperature to create an ignition delay time plot
npoints = 20
delta_temp = 20.0
init_temp = premixed.temperature
delaytime = np.zeros(npoints, dtype=np.double)
temp_inv = np.zeros_like(delaytime, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all cases with different initial gas/reactor temperatures
for i in range(npoints):
    # update the initial reactor temperature
    MyCONP.temperature = init_temp  # K
    # show the additional keywords given by user
    # (verify inputs before running the simulation)
    # MyCONP.showkeywordinputlines()
    # run the reactor model
    runstatus = MyCONP.run()
    #
    if runstatus == 0:
        # plot 1/T instead of T
        temp_inv[i] = 1.0e0 / init_temp
        # get ignition delay time
        delaytime[i] = MyCONP.get_ignition_delay()
        print(f"ignition delay time = {delaytime[i]} [msec]")
    else:
        # if get this, most likely the END time is too short
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    init_temp += delta_temp

# compute the total runtime
runtime = time.time() - start_time
print(f"total simulation duration: {runtime} [sec] for {npoints} cases")

###############################
# Plot the ignition delay curve
# =============================
# The ignition delay curve is customarily plotted against the inverse temperature.
# The corresponding temperature values are shown on the top axis.
#
# Intuitively, you expect that the reactivity increases as the temperature increases.
# However, as you can see, in this case, the ignition delay curve does not
# vary monotonically with the temperature. Between 850 and 950 [K],
# the ignition delay time actually increases slightly
# (because the reactivity decreases). This is regarded as the
# **Negative Temperature Coefficient (NTC)** behavior which is often observed during
# low-temperature oxidation of large hydrocarbon fuels.

# create an ignition delay versus 1/T plot for the PRF fuel
# (should exhibit the NTC region)
plt.rcParams.update({"figure.autolayout": True})
fig, ax1 = plt.subplots()
ax1.semilogy(temp_inv, delaytime, "bs--")
ax1.set_xlabel("1/T [1/K]")
ax1.set_ylabel("Ignition delay time [msec]")


# Create a secondary x-axis for T (=1/(1/T))
def one_over(x):
    """Vectorized 1/x, treating x==0 manually."""
    x = np.array(x, float)
    near_zero = np.isclose(x, 0)
    x[near_zero] = np.inf
    x[~near_zero] = 1 / x[~near_zero]
    return x


# the function "1/x" is its own inverse
inverse = one_over
ax2 = ax1.secondary_xaxis("top", functions=(one_over, inverse))
ax2.set_xlabel("T [K]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_ignition_delay.png", bbox_inches="tight")

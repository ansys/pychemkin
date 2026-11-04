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

r""".. _ref_trans_PSR_ignition:

=====================================================
Transient simulation of methane-air ignition in a PSR
=====================================================
This project demonstrates the process of running a transient
perfectly stirred reactor model (transPSR) with multiple inlets.
Normally the PSR model is a 0-D steady-state model with constant
pressure. However, the dynamic responses of
a chemistry domanted process are in interest, performing transient
simulations of the system becomes necessary. In this example, a
transient PSR model with the fuel streamand the air stream introduced
separately. The "mean" reactor residence time is given to cap
the time available for ignition to take place inside the reactor.
The reactor is initially filled with hot air to promote the chemical
reactions. The fuel-air mixture must ignite within the time
duration of 1 "residence time" or the PSR would go cold. The solution
profiles also indicate that it takes several "residence times" for the
PSR to reach steady state.

"""

# sphinx_gallery_thumbnail_path = '_static/plot_trans_PSR_ignition.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path
import time

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream  # external gaseous inlet
from ansys.chemkin.core.logger import logger

# Chemkin 0-D transient PSR model with given reactor residence time
from ansys.chemkin.core.stirreactors.transient_PSR import (
    TransientPSRSetResTimeEnergyConservation as TransPsr,
)
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

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

##########################################
# Create an instance of the Chemistry Set
# ========================================
# The mechanism loaded is the GRI 3.0 mechanism for methane combustion.
# The mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"/reaction/data"*.
#

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# including the full file path is recommended
chemfile = str(mechanism_dir / "grimech30_chem.inp")
thermfile = str(mechanism_dir / "grimech30_thermo.dat")
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, label="GRI 3.0")

##############################
# Preprocess the Chemistry Set
# ============================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()
if ierror != 0:
    print("Error: failed to preprocess the mechanism!")
    print(f"       error code = {ierror}")
    exit()

###########################################
# Set up the fuel and the air inlet streams
# =========================================
# Instantiate a stream feed for the inlet gas mixture.
# The stream is a mixture with the addition of the
# inlet flow rate. You specify inlet gas properties in the same way as you
# set up a mixture.
#
# Use the ``x`` method of the ``Stream`` object to specify the inlet gas
# composition with a ``recipe``. In this project, the ``mass_flowrate`` method,
# [g/sec] is used set the volumetric flow rates of the two inlet streams.
# The ``air`` composition is copied from the pre-defined ``Air`` object in Pychemkin.

# set the fuel stream composition and conditions
fuel = Stream(MyGasMech, label="methane")
fuel.pressure = 2.0 * ck.P_ATM
fuel.temperature = 300.0
fuel.x = [("CH4", 1.0)]
# methane mass flow rate [g/sec]
fuel.mass_flowrate = 0.5

# set the air composition and conditions
air = Stream(MyGasMech, label="air")
air.pressure = fuel.pressure
air.temperature = 350.0
air.x = ck.Air.x()
# air stream mass flow rate to create a stoichiiometric combustion
# condition with methane [g/sec]
air.mass_flowrate = 8.1

#######################################
# Set up the transient PSR bomb reactor
# =====================================
# Instantiate the transient PSR ``fire_bomb`` as a
# ``TransientPSRSetResTimeEnergyConservation`` object.
# Since this is a transient simulation, a initial reactor condition is required.
# In this example, the ``fire_bomb`` is initially filled with ``air``.

fire_bomb = TransPsr(air, label="Fire_Bomb")

############################################
# Set up additional reactor model parameters
# ==========================================
# Before you can run the simulation, you must provide reactor parameters,
# solver controls, and output instructions. For a PSR, you must provide
# either the residence time or the reactor volume. You can also make changes to
# any reactor initial conditions such as the reactor temperature and the gas
# composition. For example, the ``reset_initial_temperature()``` method and
# the ``reset_initial_composition()`` method.
# You can use the solver parameters to improve the convergence performance.

# set PSR residence time (sec): required for the
# ``TransientPSRSetResTimeEnergyConservation`` model
fire_bomb.residence_time = 0.02
# transient simulation end time [sec]
# for PSR, the simulation time should be greater than the reactor residence time
fire_bomb.time = 0.2
# raise the initial reactor temperature to trigger ignition
fire_bomb.reset_initial_temperature(1500.0)

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

# connect the fuel stream to the reactor
fire_bomb.set_inlet(fuel)

# connect the air stream to the reactor
fire_bomb.set_inlet(air)

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances.

# set solver the tolerances
fire_bomb.tolerances = (1.0e-10, 1.0e-8)
# solver non-negative mode
fire_bomb.force_nonnegative = True
# transient solver max time step size [sec]
fire_bomb.set_solver_max_timestep_size(2.0e-3)
# set adaptive solution saving distance
fire_bomb.adaptive_solution_saving(mode=True, steps=60)
# use the legacy solver to get better coonvergence performance for
# small mechanism
fire_bomb.set_legacy_option(option=True)
# set ignition delay definition
fire_bomb.set_ignition_delay(method="T_inflection")

#############################
# Run the transient PSR model
# ============================
# Once all the necessary input parameters are provided, use the ``run()`` method
# to start the transient PSR simulation. After the simulation
# is completed successfully, use the ``get_ignition_delay`` method
# to extract the ignition delay time.

# set the start wall time
start_time = time.time()
# run the reactor model
runstatus = fire_bomb.run()
#
if runstatus == 0:
    # get ignition delay time
    delaytime = fire_bomb.get_ignition_delay()
    print(f"ignition delay time = {delaytime} [msec]")
else:
    # if get this, most likely the END time is too short
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)

runtime = time.time() - start_time
print(f"Total simulation duration: {runtime} [sec].")

# postprocess the solution profiles
fire_bomb.process_solution()

# get the number of data points in the solution profiles
n_points = fire_bomb.getnumbersolutionpoints()
print(f"number of solution time points = {n_points}")
# store solution profiles
# simulation time [sec]
sim_time = fire_bomb.get_solution_variable_profile("time")
# convert [sec] to [msec]
sim_time *= 1.0e3
# temperature solution profile [K]
temp_profile = fire_bomb.get_solution_variable_profile("temperature")
# total heat release rate from the gas-phase reactions [erg/sec]
gashrr_profile = fire_bomb.get_solution_variable_profile("gashrr")
# convert [erg] to [kJ]
gashrr_profile /= 1.0e3 * ck.ERGS_PER_JOULE
# gas phase CO solution profile [mass fraction]
co_profile = fire_bomb.get_solution_variable_profile("CO")
# gas phase CH4 solution profile [mass fraction]
ch4_profile = fire_bomb.get_solution_variable_profile("CH4")
# gas phase CO solution profile [mass fraction]
h2o_profile = fire_bomb.get_solution_variable_profile("H2O")
# gas phase NO solution profile [mass fraction]
no_profile = fire_bomb.get_solution_variable_profile("NO")

# clean up
ck.done()

##########################################
# Plot the transient PSR solution profiles
# ========================================
# Plot the solution profiles over the entire simulation time span.
# You could observe the ignition occurs around time = 0.002 [msec] by
# the spike in the heat release rate profile.
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("Transient PSR Solution", fontsize=16)
plt.subplot(221)
plt.plot(sim_time, temp_profile, "r-")
plt.ylabel("Temperature [K]")
plt.subplot(222)
plt.plot(sim_time, gashrr_profile, "b-")
plt.ylabel("Total Heat Release Rate [kJ/sec]")
plt.subplot(223)
plt.plot(sim_time, co_profile, "g-")
plt.plot(sim_time, ch4_profile, "r--")
plt.plot(sim_time, h2o_profile, "b:")
plt.legend(["CO", "CH4", "H2O"], loc="upper right")
plt.yscale("log")
plt.xlabel("Time [msec]")
plt.ylabel("Mass Fraction")
plt.subplot(224)
plt.plot(sim_time, no_profile, "b-")
plt.xlabel("Time [msec]")
plt.ylabel("NO Mass Fraction")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_trans_PSR_ignition.png", bbox_inches="tight")

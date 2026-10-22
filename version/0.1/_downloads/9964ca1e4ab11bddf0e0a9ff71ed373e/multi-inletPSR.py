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

r""".. _ref_multiple_inlets_PSR:

=============================================================
Determine the impact of residence time on combustion in a PSR
=============================================================

Ansys Chemkin offers some idealized reactor models commonly used for studying chemical
processes and for developing reaction mechanisms. The PSR (perfectly stirred reactor)
model is a steady-state 0-D model of the open perfectly mixed gas-phase reactor.
There is no limit on the number of inlets to the PSR. As soon as the inlet gases
enter the reactor, they are thoroughly mixed with the gas mixture inside. The PSR has
only one outlet, and the outlet gas is assumed to be exactly the same as
the gas mixture in the PSR. There are two basic types of PSR models:

- **constrained-pressure** (or set residence time)
- **constrained-volume**

By default, the PSR model is running under constant pressure. The PyChemkin PSR models
always require the connected inlets to be defined, that is, the total inlet flow rate
to the PSR is always known. Therefore, either the residence time or the reactor volume
is needed to satisfy the basic setup of the PSR model.

This example specifies the reactor volume of the PSR. The residence time is calculated
from the reactor volume and the total inlet volumetric flow rate.

For each type of PSR model, you can either specify the reactor temperature (as a fixed
value or by a piecewise-linear profile) or solve the energy conservation equation.
In total, you get four variations of the PSR model.

PSR models are mostly employed in chemical kinetics studies. By controlling the reactor
temperature, pressure, and/or residence time, you can gain knowledge about
the major intermediates of a complex chemical process and postulate possible
reaction pathways.

This example describes a parameter study of the influence of the PSR residence time
on the hydrogen combustion process. It uses two inlet streams, one for
the fuel mixture and the other for the air mixture. The fuel-to-air ratio inside
the PSR is determined by the mass or the volumetric flow rate ratio of
the two inlet streams.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_multi_inlet_PSR.png'

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
from ansys.chemkin.core.stirreactors.PSR import PSRSetVolumeEnergyConservation as Psr
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
# ``Hydrogen-Ammonia-NOx_chem_MFL2021.inp``.
# This mechanism is developed under Chemkin's
# **Model Fuel Library (MFL)** project.
# Like the rest of the MFL mechanisms, it is located in ``ModelFuelLibrary``
# in the ``/reaction/data`` directory of
# the standard Ansys Chemkin installation.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data" / "ModelFuelLibrary" / "Skeletal"
mechanism_dir = str(data_dir)
# create a chemistry set based on the gasoline 14 components mechanism
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

#################################
# Set up the fuel and air streams
# ===============================
# Instantiate a stream named ``fuel`` for the inlet stream containing
# the fuel mixture and a stream named ``air`` for the inlet stream
# containing the air mixture. The ``Inlet`` object is a mixture with
# the addition of the inlet flow rate.
#
# You specify inlet gas properties in the same way as you set up a mixture.
# Here, the ``fuel`` and ``air`` inlets are created separately. You can adjust
# their inlet volumetric flow rates to create the desired hydrogen-air mixture.
# You use the ``vol_flowrate()`` method to assign the inlet volumetric flow rate.
# In this project, the ``fuel`` and ``air`` inlets have fixed volumetric flow rates
# of 25 and 50 [cm3/sec], respectively.
#
# .. note ::
#   This equation is used to determine the hydrogen-to-oxygen molar ratio:
#
#   .. math ::
#       H_{2}:O_{2}=0.21[cm^3/cm^3]*25[cm^3/s]:0.21[cm^3/cm^3]*50[cm^3/s]=1:2
#
# .. note ::
#   The PSR residence time is specified indirectly through the reactor volume
#   in the parameter study.
#

# create the fuel inlet
fuel = Stream(MyGasMech, label="Fuel")
# set fuel composition
fuel.x = [("h2", 0.21), ("n2", 0.79)]
# setting pressure and temperature is not required in this case
fuel.pressure = ck.P_ATM
fuel.temperature = 450.0  # inlet temperature
# set inlet volumetric flow rate [cm3/sec]
fuel.vol_flowrate = 25.0

# create the oxidizer inlet: air
air = Stream(MyGasMech, label="Oxid")
air.x = [("o2", 0.21), ("n2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = fuel.pressure
air.temperature = fuel.temperature
# set inlet volumetric flow rate [cm3/sec]
air.vol_flowrate = 50.0

####################################################################
# Create the PSR to predict the gas composition of the outlet stream
# ==================================================================
# Use the ``PSRSetVolumeEnergyConservation()`` method to instantiate a PSR
# named ``combustor``. You must include the energy equation because the goal is to
# see how the residence time would affect the hydrogen combustion process.
# The ``combustor`` PSR is initiated with the parameter set to the ``fuel`` inlet,
# Hence, the gas property of the ``fuel`` inlet is applied
# as the estimated reactor condition of the ``combustor`` PSR. You can overwrite
# any estimated reactor conditions by using appropriate methods. For example,
# ``combustor.temperature = 2000.0`` changes the estimated reactor temperature
# from 450 (the temperature of the ``fuel`` inlet) to 2000 [K]. The ``volume()``
# method sets the reactor volume of the nominal case.

# create a PSR with fixed reactor volume and
# with the fuel inlet composition as the estimated reactor condition
combustor = Psr(fuel, label="tincan")

############################################
# Set up additional reactor model parameters
# ==========================================
# You must provide reactor parameters, solver controls, and output instructions
# before running the simulations. For the steady-state PSR, you must provide either
# the residence time or the reactor volume. You can also make changes to
# any estimated reactor conditions if desired.

# reset the estimated reactor temperature [K]
combustor.temperature = 2000.0
# set the reactor volume (cm3): required for PSRSetVolumeEnergyConservation model
combustor.volume = 200.0

###################################
# Connect the inlets to the reactor
# =================================
# You must connect at least one inlet to the open reactor. Use the ``set_inlet()``
# method to add a stream object to the PSR. Inversely, use the ``remove_inlet()``
# to disconnect an inlet from the PSR. Here two inlets, ``fuel`` and ``air``,
# are connected to the ``combustor`` PSR. The fuel-to-air ratio is controlled by
# the mass or the volumetric flow rate ratio of the ``fuel`` and the ``air`` inlets.

# add external inlets to the PSR
combustor.set_inlet(fuel)
combustor.set_inlet(air)

#####################
# Set solver controls
# ===================
# You can overwrite the default solver controls by using solver-related methods,
# such as those for tolerances. The following code changes the tolerances that
# the steady-state solver is to use for the steady-state search and changes
# the pseudo time stepping stages. Sometimes, during the iterations, some
# species mass fractions might become negative, causing the solver to report
# an error and stop. To overcome this issue, you can provide a small cushion to
# allow species mass fractions to go slightly negative by using the
# ``set_species_floor()`` method to reset the mass fraction floor value.

# reset the tolerances in the steady-state solver
combustor.steady_state_tolerances = (1.0e-9, 1.0e-6)
combustor.timestepping_tolerances = (1.0e-9, 1.0e-6)
# reset the gas species floor value in the steady-state solver
combustor.set_species_floor(-1.0e-10)

############################################
# Run the PSR residence time parameter study
# ==========================================
# The PSR residence time :math:`\tau` is calculated as follows:
#
# .. math ::
#
#   \tau =
#   \frac{reactor\text{ }volume}{total\text{ }volumetric\text{ }flow\text{ }rate}
#
# In this parameter study, the reactor volume is decreased from 200 to
# 160 [cm3]. Accordingly, :math:`\tau` is decreased from
# 2.67 to 2.13 [sec] because the total inlet volumetric flow rate is kept
# constant at 75 [cm3/sec]. Usually you want to run the burning cases
# (large residence times) first in the PSR parameter study.
#
# The ``process_solution`` method converts the result from each PSR run
# to a mixture. You can either overwrite the solution mixture or
# use a new one for each simulation result.

# reactor volume increment
delta_vol = -5
numbruns = 9
# solution arrays
residencetime = np.zeros(numbruns, dtype=np.double)
temp_ss_solution = np.zeros_like(residencetime, dtype=np.double)
# set the start wall time
start_time = time.time()
# loop over all inlet temperature values
for i in range(numbruns):
    # run the PSR model
    runstatus = combustor.run()
    # check run status
    if runstatus != 0:
        # Run failed.
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()
    # Run succeeded.
    print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
    # postprocess the solution profiles
    solnmixture = combustor.process_solution()
    # print the steady-state solution values
    print(f"steady-state temperature = {solnmixture.temperature} [K]")
    # solnmixture.list_composition(mode="mole")
    # store solution values
    # final reactor gas density [g/cm3]
    # density = solnmixture.RHO
    # final reactor mass [g]
    # mass = density * combustor.volume
    # PSR apparent residence time [sec]
    residencetime[i] = combustor.volume / combustor.net_vol_flowrate
    temp_ss_solution[i] = solnmixture.temperature
    # update reactor volume
    combustor.volume += delta_vol

# compute the total runtime
runtime = time.time() - start_time
print(f"Total simulation duration: {runtime} [sec] over {numbruns} runs")

##################################
# Plot the parameter study results
# ================================
# Plot the predicted PSR temperature against the residence time. You can observe
# that the hydrogen-air mixture ceases to burn when the residence time becomes
# too small. Gas turbine terminology refers to this as *blown out*
# because the fuel-air mixture gets blown out of the combustor before any
# significant chemical reaction can take place.
plt.plot(residencetime, temp_ss_solution, "bo-")
plt.xlabel("Apparent Residence Time [sec]")
plt.ylabel("Exit Gas Temperature [K]")
plt.title("PSR Solution")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_multi_inlet_PSR.png", bbox_inches="tight")

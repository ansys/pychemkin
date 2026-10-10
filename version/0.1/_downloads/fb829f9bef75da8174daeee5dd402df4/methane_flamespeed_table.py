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

r""".. _ref_flame_speed_table:

=============================================================================
Construct atmospheric  methane-air flame speed versus equivalence ratio table
=============================================================================

One of the prevailing use case of the *freely propagating premixed flame* model is
to build a *flame speed* table to be imported by another combustion simulation tools.
PyChemkin provides the flexibility to customize the data structure of
the flame speed table depending on the simulation goals and the tool. Furthermore,
over the years, the chemkin flame speed calculator has derived a set of default
solver settings that would greatly improve the convergence performance, especially
for those widely adopted hydrocarbon fuel combustion mechanisms. The required
input parameters the flame speed calculator are reduced to the composition of
the fuel-oxidizer mixture, the initial/inlet pressure and temperature,
and the calculation domain.

This tutorial shows the "minimal" effort to create a flame speed table of
CH\ :sub:`4`\ -air mixtures at the atmospheric pressure. The predicted flame speed
values are compared against the experimental data as a function of
the mixture equivalence ratio. Since the transport processes are critical for
flame calculations, the transport data must be included in the mechanism data
and preprocessed.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_flame_speed_table.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path
import time

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color
from ansys.chemkin.core.inlet import Stream  # external gaseous inlet
from ansys.chemkin.core.logger import logger

# Chemkin 1-D premixed freely propagating flame model (steady-state)
from ansys.chemkin.core.premixedflames.premixedflame import (
    FreelyPropagating as FlameSpeed,
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
# .. note::
#   The transport data *must* be included and preprocessed because
#   the transport processes, *convection and diffusion*, are important to
#   sustain the flame structure.
#

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# including the full file path is recommended
chemfile = str(mechanism_dir / "grimech30_chem.inp")
thermfile = str(mechanism_dir / "grimech30_thermo.dat")
tranfile = str(mechanism_dir / "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

##############################
# Preprocess the Chemistry Set
# ============================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()
if ierror != 0:
    print("Error: failed to preprocess the mechanism!")
    print(f"       error code = {ierror}")
    exit()

########################################################################
# Set up the CH\ :sub:`4`\ -air mixture for the flame speed calculation
# ======================================================================
# Instantiate a stream named ``premixed`` for the inlet gas mixture.
# This stream  is a mixture with the addition of the
# inlet flow rate. You can specify the inlet gas properties the same way you
# set up a ``Mixture``. Here the ``x_by_equivalence_ratio`` method is used.
# You create the ``fuel`` and the ``air`` mixtures first. Then define the
# *complete combustion product species* and provide the *additives* composition
# if applicable. And finally, during the parameter iteration runs, you can simply set
# different values to ``equivalenceratio`` to create different methane-air mixtures.
#

# create the fuel mixture
fuel = ck.Mixture(MyGasMech)
# set fuel composition: methane
fuel.x = [("CH4", 1.0)]
# setting pressure and temperature condition for the flame speed calculations
fuel.pressure = 1.0 * ck.P_ATM
fuel.temperature = 300.0  # inlet temperature

# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.x = ck.Air.x()
# setting pressure and temperature is not required in this case
air.pressure = fuel.pressure
air.temperature = fuel.temperature

# create the fuel-air Stream for the premixed flame speed calculation
premixed = Stream(MyGasMech, label="premixed")
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture.
# can also create an additives mixture here
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)  # no additives: all zeros

# setting pressure and temperature is not required in this case
premixed.pressure = fuel.pressure
premixed.temperature = fuel.temperature

# set estimated value of the flame mass flux [g/cm2-sec]
premixed.mass_flowrate = 0.4

# equivalence ratio for the first case
phi = 0.6
# create mixture by using the equivalence ratio
ierror = premixed.x_by_equivalence_ratio(
    MyGasMech, fuel.x, air.x, add_frac, products, equivalenceratio=phi
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
    print(
        "Error: failed to create the methane-air mixture "
        + "for equivalence ratio = "
        + str(phi)
    )
    exit()

##################################################
# Instantiate the laminar speed calculator
# ================================================
# Set up the *freely propagating premixed flame* model by using the stream
# representing the premixed fuel-oxidizer mixture (with the estimated
# mass flow rate value). When the flame speed is expected to be very small
# (< 10 [cm/sec]) or very large (> 300 [cm/sec]), it might be beneficial to
# set the ``mass_flowrate`` of this inlet to the mass flux [g/cm\ :sup:`2`\ -sec]
# based on the estimated flame speed value. There are many options and parameters
# related to the treatment of the species boundary condition, the transport properties.
# All the available options and parameters are described in the *Chemkin Input* manual.
#
# .. note::
#   The stream parameter used to instantiate a ``FlameSpeed`` object is
#   the properties of the unburned fuel-oxidizer mixture of which the *laminar
#   flame speed* will be determined.
#

flamespeedcalculator = FlameSpeed(premixed, label="premixed_methane")

###############################################
# Set up initial mesh and grid adaption options
# =============================================
# The ``end_poistion`` is a required input as it defines the length of
# the calculation domain. Typically, the length of the calculation domain is
# between 1 to 10 [cm]. For low pressure conditions, the flame thickness becomes
# wider and a larger calculation domain is required.
#

# set the maximum total number of grid points allowed in the calculation (optional)
# flamespeedcalculator.set_max_grid_points(150)
# define the calculation domain [cm]
flamespeedcalculator.end_position = 1.0

#####################################
# Run the flame speed parameter study
# ===================================
# Use the ``run()`` method to run the freely propagating premixed flame
# (flame speed) model. After the premixed flame calculation concludes successfully,
# use the ``process_solution()`` method to postprocess the solutions.
# The predicted laminar flame speed can be obtained by using the ``get_flame_speed()``
# method. You can create other property profiles by looping through
# the solution streams with proper ``Mixture`` methods. The parameter in this project
# is the equivalence ratio of the methane-air mixture. You can start
# the parameter run from the most fuel-lean or from the most fuel-rich case. Normally,
# the most "extreme" cases are difficult to converge. When running
# into these situations, start the parameter runs from the stoichiometric condition
# and go down the lean and/or the rich branch. Here the runs start from
# the most fuel-lean case (\ :math:`\phi = 0.6`\) and progress all the way to
# the most fuel-rich case (\ :math:`\phi = 1.6`\) in steps of 0.05.
#
# .. note::
#   - When the inlet stream condition is close to the flammability limit,
#     the flame speed calculation might fail. Remember that the reaction mechanism
#     (reaction rates, thermodynamic properties, and transport properties) and
#     the reactor model are *models* that contain assumptions and uncertainties.
#   - After complete the first run, you can use the ``continuation()`` method to
#     start the new runs from the solution of the previous run. However,
#     by doing this, the later runs will contain a lot of grid points accumulated
#     from all previous runs.
#   - Use the ``set_molefractions`` method to update the inlet gas composition
#     before each run. Similarly, use the ``pressure`` and the ``temperature``
#     methods to change the inlet condition.
#

# total number of parameter cases
points = 21
# equivalence ratio increment
delta_phi = 0.05
# create solution arrays
equival = np.zeros(points, dtype=np.double)
flamespeed = np.zeros_like(equival, dtype=np.double)
# set the start wall time
start_time = time.time()

# start the parameter study runs
for i in range(points):
    # run the flame speed calculation for this equivalence ratio
    status = flamespeedcalculator.run()
    if status != 0:
        print(
            Color.RED
            + "failed to calculate the laminar flame speed"
            + "for equivalence ratio = "
            + str(phi)
            + Color.END
        )
        exit()
    # get flame speed
    # postprocess the solutions
    flamespeedcalculator.process_solution()
    # save data
    equival[i] = phi
    # get flame speed
    flamespeed[i] = flamespeedcalculator.get_flame_speed()
    # print the predicted laminar flame speed
    print(
        f"methane-air equivalence ratio = {phi} :\n"
        + f"the predicted laminar flame speed = {flamespeed[i]} [cm/sec]"
    )
    #
    # update parameter
    phi += delta_phi
    # create mixture by using the equivalence ratio
    ierror = premixed.x_by_equivalence_ratio(
        MyGasMech, fuel.x, air.x, add_frac, products, equivalenceratio=phi
    )
    # check fuel-oxidizer mixture creation status
    if ierror != 0:
        print(
            "Error: failed to create the methane-air mixture ",
            "for equivalence ratio = ",
            str(phi),
        )
        exit()
    # update initial gas composition
    flamespeedcalculator.set_molefractions(premixed.x)

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"total simulation duration: {runtime} [sec]")
print()

# experimental data by Vagelopoulos
# equivalence ratios
data_equiv = [
    0.6126,
    0.6619,
    0.7109,
    0.7533,
    0.8268,
    0.9109,
    0.9826,
    1.0387,
    1.0901,
    1.1321,
    1.1858,
    1.2347,
    1.2695,
    1.325,
    1.3563,
    1.4279,
    1.4977,
]
# methane flame speeds at 1 atm
data_speed = [
    9.4434,
    12.7281,
    17.4088,
    21.0219,
    26.5237,
    33.2573,
    37.0347,
    38.677,
    38.8412,
    37.8558,
    34.9818,
    31.1223,
    25.9489,
    21.3504,
    17.2445,
    12.7281,
    9.7719,
]

###########################################
# Plot the premixed flame solution profiles
# =========================================
# Plot the predicted flame speeds against the experimental data
#

plt.plot(data_equiv, data_speed, label="data", linestyle="", marker="^", color="blue")
plt.plot(equival, flamespeed, label=MyGasMech.label, linestyle="-", color="blue")
plt.legend()
plt.ylabel("Flame Speed [cm/sec]")
plt.xlabel("Equivalence Ratio")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_flame_speed_table.png", bbox_inches="tight")

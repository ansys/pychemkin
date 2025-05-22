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
.. _ref_plug_flow_reactor:

===========================================
Simulate NO reduction in combustion exhaust
===========================================

The PFR (plug flow reactor) model is widely adopted in the chemical kinetic
studies of emission abatement processes because of its simplicity in flow treatment
as well as its resemblance of an exhaust duct or a tailpipe.

The PFR model is a steady-state open reactor because materials are allowed to
move in and out of the reactor. The solutions are obtained along the length of the reactor
as the inlet materials march toward the reactor exit. Typically, the pressure of the PFR
is fixed. However, Chemkin does permit the use of a pressure profile to enforce a certain
pressure gradient in the PFR.

Use the ``PlugFlowReactor_FixedTemperature()`` or ``PlugFlowReactor_EnergyConservation()``
method to create a PFR. Unlike the closed batch reactor model, which is instantiated by a mixture,
the open reactor model, such as the PFR model, is initiated by a stream, which is
simply a mixture with the addition of the inlet mass/volumetric
flow rate or velocity. You already know how to create a stream if you know how to create
a mixture. You can specify the inlet flow rate using one of these methods:
``velocity()``, ``mass_flowrate()``, ``vol_flowrate()``, or ``sccm()``
(standard cubic centimeters per minute).

This example shows how to use the Chemkin PFR model to study the
reduction of nitric oxide (NO) in the combustion exhaust by using the
CH\ :sub:`4` reburning process. This is achieved by injecting CH\ :sub:`4` into the
hot exhaust gas mixture. As the exhaust gas travels along the tubular reactor, the injected
CH\ :sub:`4` is oxidized by the NO to form N\ :sub:`2`\ , CO, and H\ :sub:`2`. This is why this NO reduction process is called *methane reburning*.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_plug_flow_reactor.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color

# chemkin plug flow reactor model
from ansys.chemkin.flowreactors.PFR import PlugFlowReactor_FixedTemperature
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True

########################
# Create a chemistry set
# ======================
# The mechanism to load is the C2 NOx mechanism for the combustion of C1-C2 hydrocarbons.
# The mechanism comes with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="C2 NOx")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")

######################################
# Preprocess the GRI 3.0 chemistry set
# ====================================

# preprocess the mechanism files
iError = MyGasMech.preprocess()

###################################
# Instantiate and set up the stream
# =================================
# Create a hot exhaust gas mixture by assigning the mole fractions of the
# species. A stream is simply a mixture with the addition of
# mass/volumetric flow rate or velocity. Here, the gas composition ``exhaust.X`` is given
# by a recipe consisting of the common species in the combustion exhaust, such as
# CO\ :sub:`2` and H\ :sub:`2`O and the CH\ :sub:`4` injected. The inlet velocity
# is given by ``exhaust.velocity = 26.815``.

# create the inlet (mixture + flow rate)
exhaust = Stream(MyGasMech)
# set inlet temperature [K]
exhaust.temperature = 1750.0
# set inlet/PFR pressure [atm]
exhaust.pressure = 0.83 * ck.Patm
# set inlet molar composition directly
exhaust.X = [
    ("CO", 0.0145),
    ("CO2", 0.0735),
    ("H2O", 0.1107),
    ("NO", 0.0012),
    ("N2", 0.7981),
    ("CH4", 0.002),
]
# set inlet velocity [cm/sec]
exhaust.velocity = 26.815
#

#####################################
# Create the plug flow reactor object
# ===================================
# Use the ``PlugFlowReactor_FixedTemperature()`` method to create a plug flow reactor.
# The required input parameter of the open reactor models in PyChemkin
# is the stream, while PyChemkin batch reactor models take a mixture
# as the input parameter. In this case, the exhaust stream is used to create a
# PFR named ``tubereactor``.
tubereactor = PlugFlowReactor_FixedTemperature(exhaust)

############################################
# Set up additional reactor model parameters
# ==========================================
# For the PFR, the required reactor parameters are the reactor diameter [cm] or the
# cross-sectional flow area [cm2] and the reactor length [cm].
# The mixture condition and inlet mass flow rate of the inlet are already defined
# by the stream when the ``tubereactor`` PFR is instantiated.

# set PFR diameter [cm]
tubereactor.diameter = 5.8431
# set PFR length [cm]
tubereactor.length = 5.0

#######################################
# Verify the inlet condition of the PFR
# =====================================
print(f"PFR inlet mass flow rate {tubereactor.mass_flowrate} [g/sec]")
print(f"PFR inlet velocity {tubereactor.velocity} [cm/sec]")
# show inlet gas composition of the PFR
print("PFR inlet gas compsition")
tubereactor.list_composition(mode="mass", bound=1.0e-8)

####################
# Set output options
# ==================
# You can turn on adaptive solution saving to resolve the steep variations in the solution
# profile. Here, additional solution data points are saved for every **100** internal solver steps.
#
# .. note::
#   By default, distance intervals for both print and save solution are 1/100 of the
#   reactor length. In this case, :math:`dt=time/100=0.001`\ . You can change them
#   to different values.

# set distance between saving solution
tubereactor.timestep_for_saving_solution = 0.0005
# turn on adaptive solution saving
tubereactor.adaptive_solution_saving(mode=True, steps=100)

#########################################
# Display the added parameters (keywords)
# =======================================
# Use the ``showkeywordinputlines()`` method to verify that the preceding
# parameters are correctly assigned to the reactor model.
tubereactor.showkeywordinputlines()

######################################################
# Run the parameter study to replicate the experiments
# ====================================================
# Use the ``run()`` method to start the plug flow simulation.
#
# .. note ::
#   You can use two ``time`` calls (one before the run and one after the run) to
#   get the simulation run time (wall time).
#

# set the start wall time
start_time = time.time()
# run the PFR model
runstatus = tubereactor.run()
# compute the total runtime
runtime = time.time() - start_time
# check run status
if runstatus != 0:
    # Run failed.
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    exit()
# run succeeded.
print(Color.GREEN + ">>> Run completed. <<<", end=Color.END)
print(f"Total simulation duration: {runtime * 1.0e3} [msec]")

##########################
# Postprocess the solution
# ========================
# The postprocessing step parses the solution and packages the solution values at each
# time point into a mixture. There are two ways to access the solution profiles:
#
# - The raw solution profiles (value as a function of distance) are available for distance,
#   temperature, pressure, volume, and species mass fractions.
#
# - The mixtures permit the use of all property and rate utilities to extract
#   information such as viscosity, density, and mole fractions.
#
# You can use the ``get_solution_variable_profile()`` method to get the raw solution profiles. You
# can get solution mixtures using either the ``get_solution_mixture_at_index()`` method for the
# solution mixture at the given saved location or the ``get_solution_mixture()`` method for the
# solution mixture at the given distance. (In this case, the mixture is constructed by interpolation.)
#
# You can get the outlet solution by simply grabbing the solution values at the very last point by
# using \ :math:`(outlet solution index) = (number of solutions) - 1`\.
#

# postprocess the solution profiles
tubereactor.process_solution()

# get the number of solution time points
solutionpoints = tubereactor.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")

# get the grid profile [cm]
xprofile = tubereactor.get_solution_variable_profile("time")
# get the temperature profile [K]
tempprofile = tubereactor.get_solution_variable_profile("temperature")
# get the CH4 mass fraction profile
YCH4profile = tubereactor.get_solution_variable_profile("CH4")
# get the CO mass fraction profile
YCOprofile = tubereactor.get_solution_variable_profile("CO")
# get the H2O mass fraction profile
YH2Oprofile = tubereactor.get_solution_variable_profile("H2O")
# get the H2 mass fraction profile
YH2profile = tubereactor.get_solution_variable_profile("H2")
# get the NO mass fraction profile
YNOprofile = tubereactor.get_solution_variable_profile("NO")
# get the N2 mass fraction profile
YN2profile = tubereactor.get_solution_variable_profile("N2")

# inlet NO mass fraction
YNO_inlet = exhaust.Y[MyGasMech.get_specindex("NO")]
# outlet grid index
xout_index = solutionpoints - 1
print("At the reactor inlet: x = 0 [xm]")
print(f"The NO mass fraction = {YNO_inlet}.")
print(f"At the reactor outlet: x = {xprofile[xout_index]} [cm]")
print(f"The NO mass fraction = {YNOprofile[xout_index]}.")
print(
    "The NO conversion rate = "
    + f"{(YNO_inlet - YNOprofile[xout_index]) / YNO_inlet * 100.0} %\n."
)

# more involved postprocessing using mixtures
#
# create arrays for the gas velocity profile
velocityprofile = np.zeros_like(xprofile, dtype=np.double)
# reactor mass flow rate (constant) [g/sec]
massflowrate = tubereactor.mass_flowrate
# reactor cross-section area [cm2]
areaflow = tubereactor.flowarea
print(f"Mass flow rate: {massflowrate} [g/sec]\nflow area: {areaflow} [cm2]")
# ratio
ratio = massflowrate / areaflow
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = tubereactor.get_solution_mixture_at_index(solution_index=i)
    # get gas density [g/cm3]
    den = solutionmixture.RHO
    # gas velocity [g]
    velocityprofile[i] = ratio / den

############################
# Plot the solution profiles
# ==========================
# Plot the species and the gas velocity profiles along the tubular reactor.
#
# You can see from the plots that CH\ :sub:`4` is mainly oxidized by NO to form
# N\ :sub:`2` and H\ :sub:`2`\ in the hot exhaust. The gas velocity accelerates because of
# the net increase in total number of moles of gas species due to the chemical reactions.
# You can conduct more in depth analyses of the reaction pathways of the CH\ :sub:`4` reburning
# mechanism by applying other PyChemkin utilities.
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle("Constant Temperature Plug-Flow Reactor", fontsize=16)
plt.subplot(221)
plt.plot(xprofile, YN2profile, "r-")
plt.ylabel("N2 Mass Fraction")
plt.subplot(222)
plt.plot(xprofile, YH2Oprofile, "b-", label="H2O")
plt.ylabel("H2O Mass Fraction")
plt.subplot(223)
plt.plot(xprofile, YNOprofile, "g-", label="NO")
plt.plot(xprofile, YCH4profile, "g:", label="CH4")
plt.plot(xprofile, YH2profile, "g--", label="H2")
plt.legend()
plt.xlabel("distance [cm]")
plt.ylabel("Mass Fraction")
plt.subplot(224)
plt.plot(xprofile, velocityprofile, "m-")
plt.xlabel("distance [cm]")
plt.ylabel("Gas Velocity [cm/sec]")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_plug_flow_reactor.png", bbox_inches="tight")

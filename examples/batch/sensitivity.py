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

r""".. _ref_brute_force_sensitivity:

========================================
Perform brute-force sensitivity analysis
========================================

One of the advantages of PyChemkin is customizability. You can easily build a
specialized workflow to facilitate your simulation goals.

This tutorial demonstrates how to create a purpose-built workflow with PyChemkin
by conducting a brute-force A-factor sensitivity analysis for ignition delay time
of a premixed natural gas-air mixture at given initial temperature and pressure.
Most Chemkin reactor models have the A-factor sensitivity analysis capability.
The catch is that the subject variable of the analysis must be a member of
the solution variables such as temperature, species mass fractions, and
mass flow rate. However, for derived variables such as the ignition delay time,
the built-in sensitivity analysis work as convenient. Thus, in this case,
you may want to resort to the brute-force method to obtain those
A-factor sensitivity coefficients with respect to the ignition delay time.

To conduct the brute-force A-factor sensitivity analysis, you will have to
repeat the three steps for every reaction in the mechanism one by one

    1.  perturb the A-factor (the Arrhenius pre-exponent parameters) of a reaction

    2.  obtain the ignition delay time with this perturbed A-factor by running
        a constant pressure batch reactor simulation

    3.  restore the A-factor to its original value

The normalized ignition delay time sensitivity coefficient of reaction :math:`j`
is the difference between the original and the perturbed ignition delay time values
divided by the size of the A-factor disturbance

.. math ::
    S_{j} = \\frac{(I_{j,ptb} - I_{j,org})}{(A_{j,ptb} - A_{j,org})/A_{j,org}}

"""

# sphinx_gallery_thumbnail_path = '_static/plot_sensitivity_analysis.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
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
# The mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation in the ``/reaction/data`` directory.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# create a chemistry set based on the diesel 14 components mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# including the full file path is recommended
MyGasMech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
MyGasMech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")

###################################
# Pre-process the ``Chemistry Set``
# =================================
ierror = MyGasMech.preprocess()
# check preprocess status
if ierror == 0:
    print("mechanism information:")
    print(f"number of gas species = {MyGasMech.kk:d}")
    print(f"number of gas reactions = {MyGasMech.ii_gas:d}")
else:
    # When a non-zero value is returned from the process, check the text output files
    # chem.out, tran.out, or summary.out for potential error messages about
    # the mechanism data.
    print(f"Preprocessing error encountered. Code = {ierror:d}.")
    print(f"see the summary file {MyGasMech.summaryfile} for details")
    exit()

####################################################################
# Set up gas mixtures based on the species in this ``Chemistry Set``
# ==================================================================
# Use the *equivalence ratio method* so that you can easily set up
# the premixed fuel-oxidizer mixture composition by assigning an
# *equivalence ratio* value. In this case, the fuel mixture consists
# of methane, ethane, and propane as the simulated "natural gas".
# The premixed air-fuel mixture has an equivalence ratio of 1.1.
oxid = ck.Mixture(MyGasMech)
# set mole fraction
oxid.x = [("O2", 1.0), ("N2", 3.76)]
oxid.temperature = 900
oxid.pressure = ck.P_ATM  # 1 atm

fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.x = [("C3H8", 0.1), ("CH4", 0.8), ("H2", 0.1)]
fuel.temperature = oxid.temperature
fuel.pressure = oxid.pressure

mixture = ck.Mixture(MyGasMech)
mixture.pressure = oxid.pressure
mixture.temperature = oxid.temperature
products = ["CO2", "H2O", "N2"]
add_frac = np.zeros(MyGasMech.kk, dtype=np.double)
# create the air-fuel mixture by using the equivalence ratio method
ierror = mixture.x_by_equivalence_ratio(
    MyGasMech, fuel.x, oxid.x, add_frac, products, equivalenceratio=1.1
)
# check fuel-oxidizer mixture creation status
if ierror != 0:
    print("Error: Failed to create the fuel-oxidizer mixture.")
    exit()

##############################
# List the mixture composition
# ============================
# list the composition of the premixed mixture for verification.
if ck.verbose():
    mixture.list_composition(mode="mole")

########################################
# Preparing for the sensitivity analysis
# ======================================
# You need to perform some preparation work before running
# the brute-force sensitivity analysis. These tasks include: making
# backup copies of the original Arrhenius rate parameters, set up a
# *constant pressure* batch reactor object to compute the ignition
# delay times, and establish the baseline ignition delay time value
# with the original mechanism.

##################################
# Get the original rate parameters
# ================================
# The first step is to save a copy of the Arrhenius rate parameters
# of all reactions. You can use the ``get_reaction_parameters`` method
# associated with the ``MyGasMech`` object. You can also verify the
# rate parameters by "screening" their values.
a_factor, beta, active_energy = MyGasMech.get_reaction_parameters()
if ck.verbose():
    for i in range(MyGasMech.ii_gas):
        print(f"reaction: {i + 1}")
        print(f"A  = {a_factor[i]}")
        print(f"B  = {beta[i]}")
        print(f"Ea = {active_energy[i]}\n")
        if np.isclose(0.0, a_factor[i], atol=1.0e-15):
            print("reaction pre-exponential factor = 0")
            exit()

################################################################
# Create the reactor object for ignition delay time calculations
# ==============================================================
# Use the ``GivenPressureBatchReactorEnergyConservation`` method to instantiate a
# *constant pressure batch reactor that also includes the energy equation*. You
# should use the ``mixture`` you just created.
MyCONP = GivenPressureBatchReactorEnergyConservation(mixture, label="CONP")
# show initial gas composition inside the reactor for verification
MyCONP.list_composition(mode="mole")

############################################
# Set up additional reactor model parameters
# ==========================================
# *Reactor parameters*, *solver controls*, and *output instructions* need to be
# provided before running the simulations. For a batch reactor,
# the *initial volume* and the *simulation end time* are required inputs.
# The ``set_ignition_delay`` method must be included for the reactor model to
# report the *ignition delay times* after the simulation is done.
# The *inflection points* definition is employed to detect
# the auto-ignition time because ``method="T_inflection"`` is specified.
# You can choose a different auto-ignition definition.
# Allow additional solution data point to be saved so that
# the predicted temperature profile can have enough resolution to provide
# more precise ignition delay time value. Here the adoptive solution saving
# is turned on by the ``adaptive_solution_saving`` method and the solution
# will be recorded for every **20** solver internal steps. Remember to set
# a simulation end time ``time`` that is long enough to catch
# the occurrence of auto-ignition.
#
# .. note::
#   By default, time intervals for both print and save solution are **1/100**
#   of the *simulation end time*. In this case :math:`dt=time/100=0.001`\ .
#   You can change them to different values.
#

# reactor volume [cm3]
MyCONP.volume = 10.0
# simulation end time [sec]
MyCONP.time = 2.0

# turn ON adaptive solution saving
MyCONP.adaptive_solution_saving(mode=True, steps=20)
# set ignition delay
MyCONP.set_ignition_delay(method="T_inflection")

# set tolerances in tuple: (absolute tolerance, relative tolerance)
MyCONP.tolerances = (1.0e-10, 1.0e-8)

# set the start wall time to get the total simulation run time
start_time = time.time()

###############################
# Establish the baseline result
# =============================
# Run the nominal case and get the baseline ignition delay time value
# by using the ``get_ignition_delay`` method.
runstatus = MyCONP.run()
#
if runstatus == 0:
    # get ignition delay time
    delaytime_org = MyCONP.get_ignition_delay()
    print(f"ignition delay time = {delaytime_org} [msec]")
else:
    # if get this, most likely the END time is too short
    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
    print("failed to find the ignition delay time of the nominal case")
    exit()

####################################
# Run the sensitivity analysis cases
# ==================================
# Now compute the "raw" A-factor sensitivity coefficients of ignition delay time.
# Firstly, you create an array ``IGsen`` to store the sensitivity coefficients,
# the size of ``ig_sen`` must be no less than the number of reactions in
# the mechanism ``MyGasMech``. Secondly, you introduce a small perturbation to
# the A-factor one reaction at a time by using the ``set_reaction_afactor`` method.
# The advantage of this method is that you do not need to preprocess
# the ``Chemistry Set`` every time you make a change to the rate parameter.
# Then you run the same batch reactor ``MyCONP`` to get the ignition delay time.
#
# Once the simulation is complete successfully, use the ``get_ignition_delay`` method
# to extract the ignition delay time. Compute the difference between this
# ignition delay time value (with altered A-factor) and the baseline value
# (from the original mechanism) and save the result to array ``ig_sen``. Remember to
# restore the A-factor to its original value before moving on to the next reaction.

# create sensitivity coefficient array
ig_sen = np.zeros(MyGasMech.ii_gas, dtype=np.double)
# set perturbation magnitude
perturb = 0.001  # increase by 0.1%
perturb_plus_1 = 1.0 + perturb
# loop over all reactions
for i in range(MyGasMech.ii_gas):
    a_new = a_factor[i] * perturb_plus_1
    # actual reaction index
    ireac = i + 1
    # update the A factor
    MyGasMech.set_reaction_afactor(ireac, a_new)
    # run the reactor model
    runstatus = MyCONP.run()
    #
    if runstatus == 0:
        # get ignition delay time
        delaytime = MyCONP.get_ignition_delay()
        print(f"ignition delay time = {delaytime} [msec]")
        # compute d(delaytime)
        ig_sen[i] = delaytime - delaytime_org
        # restore the A factor
        MyGasMech.set_reaction_afactor(ireac, a_factor[i])
    else:
        # if get this, most likely the END time is too short
        print(f"trouble finding ignition delay time for raection {ireac}")
        print(Color.RED + ">>> Run failed. <<<", end=Color.END)
        exit()

# compute and report the total runtime (wall time)
runtime = time.time() - start_time
print(f"\ntotal simulation time: {runtime} [sec] over {MyGasMech.ii_gas + 1} runs")

#################################################
# Compute the normalized sensitivity coefficients
# ===============================================
# Compute the normalized sensitivity coefficient = d(delaytime) * A[i] / d(A[i]).
ig_sen /= perturb

##################################
# Screen and rank the coefficients
# ================================
# Print top 5 positive and negative ignition delay time sensitivity coefficients
# to reveal the reactions of which the A-factor values have the strongest impact
# on the auto-ignition timing (positively or negatively). The ranking will change
# when the mixture composition or condition is changed.
top = 5
# rank the positive coefficients
posindex = np.argpartition(ig_sen, -top)[-top:]
poscoeffs = ig_sen[posindex]

# rank the negative coefficients
neg_ig_sen = np.negative(ig_sen)
negindex = np.argpartition(neg_ig_sen, -top)[-top:]
negcoeffs = ig_sen[negindex]

# print the top sensitivity coefficients
if ck.verbose():
    print("positive sensitivity coefficients")
    for i in range(top):
        print(f"reaction {posindex[i] + 1}: coefficient = {poscoeffs[i]}")
    print()
    print("negative sensitivity coefficients")
    for i in range(top):
        print(f"reaction {negindex[i] + 1}: coefficient = {negcoeffs[i]}")

##########################################
# Plot the ranked sensitivity coefficients
# ========================================
# Create plots to show the reactions whose A-factors have most positive
# and negative influence on the ignition delay time.
plt.rcParams.update({"figure.autolayout": True, "ytick.color": "blue"})
plt.subplots(2, 1, sharex="col", figsize=(10, 5))
# convert reaction # from integers to strings
rxnstring = []
for i in range(len(posindex)):
    # the array index starting from 0 so the actual reaction # = index + 1
    rxnstring.append(MyGasMech.get_gas_reaction_string(posindex[i] + 1))
# use horizontal bar chart
plt.subplot(211)
plt.barh(rxnstring, poscoeffs, color="orange", height=0.4)
plt.axvline(x=0, c="gray", lw=1)
# convert reaction # from integers to strings
rxnstring.clear()
fnegindex = np.flip(negindex)
fnegcoeffs = np.flip(negcoeffs)
for i in range(len(negindex)):
    # the array index starting from 0 so the actual reaction # = index + 1
    rxnstring.append(MyGasMech.get_gas_reaction_string(fnegindex[i] + 1))
plt.subplot(212)
plt.barh(rxnstring, fnegcoeffs, color="orange", height=0.4)
plt.axvline(x=0, c="gray", lw=1)
plt.xlabel("Sensitivity Coefficients")
plt.suptitle("Ignition Delay Time Sensitivity", fontsize=16)
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_sensitivity_analysis.png", bbox_inches="tight")

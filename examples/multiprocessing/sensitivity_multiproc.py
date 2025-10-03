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
.. _ref_brute_force_sensitivity_multi-processor:

===============================================================
Perform brute-force sensitivity analysis on multiple processors
===============================================================

The brutal force sensitivity analysis for ignition delay time shown in this example is based on the equation below:

    S[i] = ln(IDT_1/IDT2)/ln(perturb_factor^2)

IDT_1 is the ignition delay times when A factor of reaction i is multiplied by perturb_factor
IDT_2 is the ignition delay times when A factor of reaction i is multiplied by 1/perturb_factor
In literatures of combustion research, the perturb_factor is usually 2.0.

To conduct the brute-force A-factor sensitivity analysis, you will have to define the conditions in function
run_cases(), then define the mechanism files, perturb factor, and number of processors in function run_sens()
Running sensitivity on multiple processors is a more efficient approach, especially when the mechanism is large.
A known limit is that this script would skip reactions with negative A factors. Therefore, the sensitivity of
reactions in duplicate format with one or more negative A factors may not be calculated correctly.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_sensitivity_analysis_multiproc.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import os
import time
import multiprocessing as mp  #multi processors

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.batchreactors.batchreactor import (
    GivenVolumeBatchReactor_EnergyConservation,
)
from ansys.chemkin.logger import logger
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

#####################
# Run a list of cases
#====================
# This function runs a list of i cases. In each case, the A factor of reaction ids[i] will be multiplied by factors[i]
# Requires user input: mixture composition, pressure, temperature, reactor model parameters
def run_cases(work_path, chemfile, thermfile, ids, factors):
    # check working directory
    os.chdir(work_path)
    # set verbose mode
    ck.set_verbose(False)

######################################
# Create a chemistry set for all cases
# ====================================
    MyGasMech = ck.Chemistry(label="GRI 3.0")
    MyGasMech.chemfile = chemfile
    MyGasMech.thermfile = thermfile

###################################
# Pre-process the ``Chemistry Set``
# =================================
    iError = MyGasMech.preprocess()

####################################################################
# Set up gas mixtures based on the species in this ``Chemistry Set``
# ==================================================================
# Use the *equivalence ratio method* so that you can easily set up
# the premixed fuel-oxidizer mixture composition by assigning an
# *equivalence ratio* value. In this case, the fuel mixture consists
# of methane, ethane, and propane as the simulated "natural gas".
# The premixed air-fuel mixture has an equivalence ratio of 1.1.
    oxid = ck.Mixture(MyGasMech)
    oxid.X = [("O2", 1.0), ("N2", 3.76)]
    oxid.temperature = 900
    oxid.pressure = ck.Patm

    fuel = ck.Mixture(MyGasMech)
    fuel.X = [("C3H8", 0.1), ("CH4", 0.8), ("H2", 0.1)]
    fuel.temperature = oxid.temperature
    fuel.pressure = oxid.pressure

    equi = 1.1

    mixture = ck.Mixture(MyGasMech)
    mixture.pressure = oxid.pressure
    mixture.temperature = oxid.temperature
    products = ["CO2", "H2O", "N2"]
    add_frac = np.zeros(MyGasMech.KK, dtype=np.double)
    # create the air-fuel mixture by using the equivalence ratio method
    iError = mixture.X_by_Equivalence_Ratio(
        MyGasMech, fuel.X, oxid.X, add_frac, products, equivalenceratio=equi
    )
    # check fuel-oxidizer mixture creation status
    if iError != 0:
        print("Error: Failed to create the fuel-oxidizer mixture.")
        exit()

    # Create a file for sensitivity analysis results output
    f = open("sen_output.txt", "w")
    f.write(f"cases to run = {len(ids)}, from {ids[0] + 1} to {ids[-1] + 1}\n")
    f.flush()

################################################################
# Create the reactor object for ignition delay time calculations
# ==============================================================
# Use the ``GivenVolumeBatchReactor_EnergyConservation`` method to instantiate a
# *constant volume batch reactor that also includes the energy equation*. You
# should use the ``mixture`` you just created.
    MyCONV = GivenVolumeBatchReactor_EnergyConservation(mixture, label="CONV")
    # show initial gas composition inside the reactor for verification
    MyCONV.list_composition(mode="mole")

############################################
# Set up additional reactor model parameters
# ==========================================
    # reactor volume [cm3]
    MyCONV.volume = 10.0
    MyCONV.pressure = 20.0 * ck.Patm
    # simulation end time [sec]
    MyCONV.time = 0.1

    # turn ON adaptive solution saving
    MyCONV.adaptive_solution_saving(mode=True, steps=20)
    # set ignition delay
    MyCONV.set_ignition_delay(method="T_inflection")

    # set tolerances in tuple: (absolute tolerance, relative tolerance)
    MyCONV.tolerances = (1.0e-10, 1.0e-8)

###################################
# Change A factor and run the cases
#==================================
    # Get original A factors
    Afactor, Beta, ActiveEnergy = MyGasMech.get_reaction_parameters()
    # create sensitivity coefficient array
    idts = []
    # loop over all reactions in the array
    for i in range(len(ids)):
        # Skip negtive A factors
        if Afactor[ids[i]] < 0:
            delaytime = MyCONV.time
        else:
            Anew = Afactor[ids[i]] * factors[i]
            # actual reaction index
            ireac = ids[i] + 1
            # update the A factor
            MyGasMech.set_reaction_AFactor(ireac, Anew)
            # run the reactor model
            MyCONV.stop_output()
            runstatus = MyCONV.run()
            if runstatus == 0:
                # get ignition delay time
                delaytime = MyCONV.get_ignition_delay()
                # restore the A factor
                MyGasMech.set_reaction_AFactor(ireac, Afactor[ids[i]])
            else:
                # if get this, most likely the END time is too short
                print(f"trouble finding ignition delay time for reaction {ireac}")
                print(Color.RED + ">>> Run failed. <<<", end=Color.END)
                # Set ignition delay time to simulation time limit
                delaytime = MyCONV.time
        idts.append(delaytime)
        # write reaction id and ignitino delay time to output file
        f.write(f"{ids[i] + 1}, {delaytime}\n")
        f.flush()
    f.close()
    # Return the ignition delay times
    return idts

##############################################################
# Run brutal force sensitivity analysis on multiple processors
#=============================================================
# Requies user input: number of processors, perturb or perturb factor, mechanism files
def run_sens():
    #
    current_dir = os.getcwd()
    numProcessors = 6
    perturb = 1

    perturb_factor = 1.0 + perturb  # A factor is multiplied and divided by 2 in sensitivity analysis
    # set mechanism directory (the default Chemkin mechanism data directory)
    data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    mechanism_dir = data_dir
    # create a chemistry set based on the diesel 14 components mechanism
    MyGasMech = ck.Chemistry(label="GRI 3.0")
    # set mechanism input files
    # including the full file path is recommended
    MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
    MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")

    # set the start wall time to get the total simulation run time
    start_time = time.time()
    iError = MyGasMech.preprocess()
    if iError == 0:
        print("mechanism information:")
        print(f"number of gas species = {MyGasMech.KK:d}")
        print(f"number of gas reactions = {MyGasMech.IIGas:d}")
    else:
        # When a non-zero value is returned from the process, check the text output files
        # chem.out, tran.out, or summary.out for potential error messages about the mechanism data.
        print(f"Preprocessing error encountered. Code = {iError:d}.")
        print(f"see the summary file {MyGasMech.summaryfile} for details")
        exit()

    ids = np.arange(MyGasMech.IIGas)
    ids = np.repeat(ids, 2)
    factors = np.tile([perturb_factor, 1/perturb_factor], MyGasMech.IIGas)
    baseNumTask = len(ids) // numProcessors    # basic number of tasks per core
    numShift = len(ids) % numProcessors        # number of cores taking +1 task

    args = []
    offset = 0
    # make a new directory for each array of cases and run the cases on multiple processors
    for i in range(numProcessors):
        work_path = os.path.join(current_dir, "RUN_" + str(i + 1).zfill(4))
        if not os.path.exists(work_path):
            os.makedirs(work_path)
        if i < numShift:
            div = baseNumTask + 1
        else:
            div = baseNumTask
        args.append((work_path, MyGasMech.chemfile,MyGasMech.thermfile, ids[offset: offset + div], factors[offset: offset + div]))
        offset += div
    with mp.Pool(processes=numProcessors) as pool:
        results = pool.starmap(run_cases, args)
    flat_results = np.concatenate(results)

    sens = (np.log(flat_results[0::2])-np.log(flat_results[1::2]))/(np.log(perturb_factor * perturb_factor))

    top = 10
    # rank the positive coefficients
    posindex = np.argpartition(sens, -top)[-top:]
    poscoeffs = sens[posindex]
    pos_sorted_index = np.argsort(poscoeffs)
    pos_sorted = poscoeffs[pos_sorted_index]

    # rank the negative coefficients
    NegIGsen = np.negative(sens)
    negindex = np.argpartition(NegIGsen, -top)[-top:]
    negcoeffs = sens[negindex[::-1]]
    neg_sorted_index = np.argsort(negcoeffs)
    neg_sorted = negcoeffs[neg_sorted_index]

    # print the top sensitivity coefficients
    print("positive sensitivity coefficients")
    for i in range(top):
        print(f"reaction {posindex[np.flip(pos_sorted_index)[i]] + 1}: coefficient = {np.flip(pos_sorted)[i]}")
    print()
    print("negative sensitivity coefficients")
    for i in range(top):
        print(f"reaction {negindex[neg_sorted_index[i]] + 1}: coefficient = {neg_sorted[i]}")

    # compute and report the total runtime (wall time)
    runtime = time.time() - start_time
    print(f"\ntotal simulation time: {runtime} [sec] over {MyGasMech.IIGas + 1} runs")

##########################################
# Plot the ranked sensitivity coefficients
# ========================================
# Create plots to show the reactions whose A-factors have most positive
# and negative influence on the ignition delay time.
    plt.rcParams.update({"figure.autolayout": True, "ytick.color": "blue"})
    plt.subplots(2, 1, sharex="col", figsize=(10, 5))
    # convert reaction # from integers to strings
    rxnstring = []
    for i in range(top):
        # the array index starting from 0 so the actual reaction # = index + 1
        rxnid = posindex[pos_sorted_index[i]] + 1
        # add reaction index before reaction string
        rxnstring.append(str(rxnid) + ". " + MyGasMech.get_gas_reaction_string(rxnid))
    # use horizontal bar chart
    plt.subplot(211)
    plt.barh(rxnstring, pos_sorted, color="orange", height=0.4)
    plt.axvline(x=0, c="gray", lw=1)
    # convert reaction # from integers to strings
    rxnstring.clear()
    for i in range(top):
        # the array index starting from 0 so the actual reaction # = index + 1
        rxnid = negindex[neg_sorted_index[i]] + 1
        rxnstring.append(str(rxnid) + ". " + MyGasMech.get_gas_reaction_string(rxnid))
    plt.subplot(212)
    plt.barh(rxnstring, neg_sorted, color="orange", height=0.4)
    plt.axvline(x=0, c="gray", lw=1)
    plt.xlabel("Sensitivity Coefficients")
    plt.suptitle("Ignition Delay Time Sensitivity", fontsize=16)

    # plot results
    if interactive:
        plt.show()
    else:
        plt.savefig("plot_sensitivity_analysis_multiproc.png", bbox_inches="tight")
    return


if __name__ == "__main__":
    run_sens()


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

r""".. _ref_brute_force_sensitivity_multi-processor:

===============================================================
Perform brute-force sensitivity analysis on multiple processors
===============================================================

The brutal force sensitivity analysis for ignition delay time shown in this example is
based on the equation below:

.. math::

    S[i] = ln(IDT_1/IDT_2)/ln(perturb\_factor^2)

:math:`IDT_1` is the ignition delay times
when :math:`A` factor of reaction :math:`i` is multiplied by :math:`perturb\_factor`
:math:`IDT_2` is the ignition delay times
when :math:`A` factor of reaction :math:`i` is multiplied by :math:`1/perturb\_factor`

*In literatures of combustion research, the* :math:`perturb\_factor` *is usually 2.0.*

To conduct the brute-force A-factor sensitivity analysis, you will have to define
the conditions in function ``run_cases()``, then define the mechanism files,
perturb factor, and number of processors in function ``run_sens()``.
Running sensitivity on multiple processors is a more efficient approach, especially
when the mechanism is large. A known limit is that this script would skip reactions
with negative :math:`A` factors. Therefore, the sensitivity of reactions in
duplicate format with one or more negative :math:`A` factors may not be calculated
correctly.

This example is created by Kuiwen Zhang.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_sensitivity_analysis_multiproc.png'

###############################################
# Import PyChemkin package and start the logger
# =============================================

import multiprocessing as mp  # multi processors
import os
from pathlib import Path
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core import Color

# chemkin batch reactor models (transient)
from ansys.chemkin.core.batchreactors.batchreactor import (
    GivenVolumeBatchReactorEnergyConservation as BatchReactor,
)

# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True


#####################
# Run a list of cases
# ====================
# This function runs a list of :math:`i` cases. In each case, the :math:`A` factor of
# reaction :math:`ids[i]` will be multiplied by :math:`factors[i]`.


def run_cases(work_path: str, chemfile: str, thermfile: str, ids: list, factors: list):
    """
    Run a list of cases for sensitivity analysis.

    Parameters
    ----------
    work_path : str
        The working directory for the cases.
    chemfile : str
        The full path to the Chemkin gas chemistry file.
    thermfile : str
        The full path to the Chemkin thermodynamic data file.
    ids : list
        A list of reaction indices for which the A factor will be modified.
    factors : list
        A list of factors by which the A factor of the corresponding reaction
        will be multiplied.

    """
    # check working directory
    os.chdir(work_path)
    # set verbose mode
    ck.set_verbose(False)

    ######################################
    # Create a chemistry set for all cases
    # ====================================
    my_gas_mech = ck.Chemistry(label="GRI 3.0")
    my_gas_mech.chemfile = chemfile
    my_gas_mech.thermfile = thermfile

    ###################################
    # Pre-process the ``Chemistry Set``
    # =================================
    _ = my_gas_mech.preprocess()

    ####################################################################
    # Set up gas mixtures based on the species in this ``Chemistry Set``
    # ==================================================================
    # Use the *equivalence ratio method* so that you can easily set up
    # the premixed fuel-oxidizer mixture composition by assigning an
    # the *equivalence ratio* value. In this case, the fuel mixture consists
    # of methane, ethane, and propane as the simulated "natural gas".
    # The premixed air-fuel mixture has an equivalence ratio of 1.1.
    oxid = ck.Mixture(my_gas_mech)
    oxid.x = [("O2", 0.21), ("N2", 0.79)]
    oxid.temperature = 900
    oxid.pressure = ck.P_ATM

    fuel = ck.Mixture(my_gas_mech)
    fuel.x = [("C3H8", 0.1), ("CH4", 0.8), ("H2", 0.1)]
    fuel.temperature = oxid.temperature
    fuel.pressure = oxid.pressure

    equi = 1.1

    mixture = ck.Mixture(my_gas_mech)
    mixture.pressure = oxid.pressure
    mixture.temperature = oxid.temperature
    products = ["CO2", "H2O", "N2"]
    add_frac = np.zeros(my_gas_mech.kk, dtype=np.double)
    # create the air-fuel mixture by using the equivalence ratio method
    ierror = mixture.x_by_equivalence_ratio(
        my_gas_mech, fuel.x, oxid.x, add_frac, products, equivalenceratio=equi
    )
    # check fuel-oxidizer mixture creation status
    if ierror != 0:
        print("Error: Failed to create the fuel-oxidizer mixture.")
        exit()

    ################################################################
    # Create the reactor object for ignition delay time calculations
    # ==============================================================
    # Use ``GivenVolumeBatchReactorEnergyConservation`` to instantiate a
    # *constant volume batch reactor that also includes the energy equation*. You
    # should use the ``mixture`` you just created.
    conv_bomb = BatchReactor(mixture, label="CONV")
    # show initial gas composition inside the reactor for verification
    conv_bomb.list_composition(mode="mole")

    # ##########################################
    # Set up additional reactor model parameters
    # ==========================================
    conv_bomb.volume = 10.0
    conv_bomb.pressure = 20.0 * ck.P_ATM
    # simulation end time [sec]
    conv_bomb.time = 0.1

    # turn ON adaptive solution saving
    conv_bomb.adaptive_solution_saving(mode=True, steps=20)
    # set ignition delay
    conv_bomb.set_ignition_delay(method="T_inflection")

    # set tolerances in tuple: (absolute tolerance, relative tolerance)
    conv_bomb.tolerances = (1.0e-10, 1.0e-8)

    ###################################
    # Change A factor and run the cases
    # ==================================
    # Get original A factors
    afactor, beta, active_energy = my_gas_mech.get_reaction_parameters()
    # create sensitivity coefficient array
    idts = []
    # Create a file for sensitivity analysis results output
    with Path("sen_output.txt").open("w", encoding="utf-8") as f:
        f.write(f"cases to run = {len(ids)}, from {ids[0] + 1} to {ids[-1] + 1}\n")
        f.flush()
        # loop over all reactions in the array
        for i in range(len(ids)):
            # Skip negative A factors
            if afactor[ids[i]] < 0:
                delaytime = conv_bomb.time
            else:
                afactor_new = afactor[ids[i]] * factors[i]
                # actual reaction index
                ireac = ids[i] + 1
                # update the A factor
                my_gas_mech.set_reaction_afactor(ireac, afactor_new)
                # run the reactor model
                conv_bomb.stop_output()
                runstatus = conv_bomb.run()
                if runstatus == 0:
                    # get ignition delay time
                    delaytime = conv_bomb.get_ignition_delay()
                    # restore the A factor
                    my_gas_mech.set_reaction_afactor(ireac, afactor[ids[i]])
                else:
                    # if get this, most likely the END time is too short
                    print(f"trouble finding ignition delay time for reaction {ireac}")
                    print(Color.RED + ">>> Run failed. <<<", end=Color.END)
                    # Set ignition delay time to simulation time limit
                    delaytime = conv_bomb.time
            idts.append(delaytime)
            # write reaction id and ignitino delay time to output file
            f.write(f"{ids[i] + 1}, {delaytime}\n")
            f.flush()
    # f.close()  # No need to close explicitly when using 'with' context
    # Return the ignition delay times
    return idts


##############################################################
# Run brutal force sensitivity analysis on multiple processors
# =============================================================
def run_sens():
    """
    Run brutal force sensitivity analysis on multiple processors.

    This function sets up the mechanism, divides the tasks among multiple processors,
    and computes the sensitivity coefficients for the ignition delay times.
    """
    #
    current_dir = str(Path.cwd())
    num_processors = 6
    perturb = 1

    perturb_factor = (
        1.0 + perturb
    )  # A factor is multiplied and divided by 2 in sensitivity analysis
    # set mechanism directory (the default Chemkin mechanism data directory)
    data_dir = Path(ck.ansys_dir) / "reaction" / "data"
    mechanism_dir = data_dir
    # create a chemistry set based on the diesel 14 components mechanism
    my_gas_mech = ck.Chemistry(label="GRI 3.0")
    # set mechanism input files
    # including the full file path is recommended
    my_gas_mech.chemfile = str(mechanism_dir / "grimech30_chem.inp")
    my_gas_mech.thermfile = str(mechanism_dir / "grimech30_thermo.dat")

    # set the start wall time to get the total simulation run time
    start_time = time.time()
    ierror = my_gas_mech.preprocess()
    if ierror == 0:
        print("mechanism information:")
        print(f"number of gas species = {my_gas_mech.kk:d}")
        print(f"number of gas reactions = {my_gas_mech.ii_gas:d}")
    else:
        # When a non-zero value is returned from the process,
        # check the text output files
        # chem.out, tran.out, or summary.out for potential error messages
        # about the mechanism data.
        print(f"Preprocessing error encountered. Code = {ierror:d}.")
        print(f"see the summary file {my_gas_mech.summaryfile} for details")
        exit()

    ids = np.arange(my_gas_mech.ii_gas)
    ids = np.repeat(ids, 2)
    factors = np.tile([perturb_factor, 1 / perturb_factor], my_gas_mech.ii_gas)
    base_num_task = len(ids) // num_processors  # basic number of tasks per core
    num_shift = len(ids) % num_processors  # number of cores taking +1 task

    args = []
    offset = 0
    # make a new directory for each array of cases and run the cases on
    # multiple processors
    for i in range(num_processors):
        work_path = Path(current_dir) / "RUN_" / str(i + 1).zfill(4)
        if not work_path.exists():
            work_path.mkdir(parents=True, exist_ok=True)
        if i < num_shift:
            div = base_num_task + 1
        else:
            div = base_num_task
        args.append(
            (
                work_path,
                my_gas_mech.chemfile,
                my_gas_mech.thermfile,
                ids[offset : offset + div],
                factors[offset : offset + div],
            )
        )
        offset += div
    with mp.Pool(processes=num_processors) as pool:
        results = pool.starmap(run_cases, args)
    flat_results = np.concatenate(results)

    sens = (np.log(flat_results[0::2]) - np.log(flat_results[1::2])) / (
        np.log(perturb_factor * perturb_factor)
    )

    top = 10
    # rank the positive coefficients
    posindex = np.argpartition(sens, -top)[-top:]
    poscoeffs = sens[posindex]
    pos_sorted_index = np.argsort(poscoeffs)
    pos_sorted = poscoeffs[pos_sorted_index]

    # rank the negative coefficients
    neg_ig_sen = np.negative(sens)
    negindex = np.argpartition(neg_ig_sen, -top)[-top:]
    negcoeffs = sens[negindex[::-1]]
    neg_sorted_index = np.argsort(negcoeffs)
    neg_sorted = negcoeffs[neg_sorted_index]

    # print the top sensitivity coefficients
    print("positive sensitivity coefficients")
    for i in range(top):
        print(
            f"reaction {posindex[np.flip(pos_sorted_index)[i]] + 1}: "
            f"coefficient = {np.flip(pos_sorted)[i]}"
        )
    print()
    print("negative sensitivity coefficients")
    for i in range(top):
        print(
            f"reaction {negindex[neg_sorted_index[i]] + 1}: "
            f"coefficient = {neg_sorted[i]}"
        )

    # compute and report the total runtime (wall time)
    runtime = time.time() - start_time
    print(
        f"\ntotal simulation time: {runtime} [sec] over {my_gas_mech.ii_gas + 1} runs"
    )

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
        rxnstring.append(str(rxnid) + ". " + my_gas_mech.get_gas_reaction_string(rxnid))
    # use horizontal bar chart
    plt.subplot(211)
    plt.barh(rxnstring, pos_sorted, color="orange", height=0.4)
    plt.axvline(x=0, c="gray", lw=1)
    # convert reaction # from integers to strings
    rxnstring.clear()
    for i in range(top):
        # the array index starting from 0 so the actual reaction # = index + 1
        rxnid = negindex[neg_sorted_index[i]] + 1
        rxnstring.append(str(rxnid) + ". " + my_gas_mech.get_gas_reaction_string(rxnid))
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

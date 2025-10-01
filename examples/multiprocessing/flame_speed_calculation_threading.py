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
.. _ref_flame_speed_threading:

==============================================================
Setting up a multi-thread laminar flame speed parameter study
==============================================================

One of the prevailing use case of the *freely propagating premixed flame* model is
to build a *flame speed* table to be imported by another combustion simulation tools.
PyChemkin provides the flexibility to customize the data structure of the flame speed table
depending on the simulation goals and the tool. Furthermore, over the years, the chemkin
flame speed calculator has derived a set of default solver settings that would greatly improve
the convergence performance, especially for those widely adopted hydrocarbon fuel
combustion mechanisms. The required input parameters the flame speed calculator are reduced
to the composition of the fuel-oxidizer mixture, the initial/inlet pressure and temperature,
and the calculation domain.

This tutorial shows the steps of setting up a flame speed parameter study for
CH\ :sub:`4`\ -air mixtures at the 5 atmosphere pressure. The predicted flame speed values are
compared against the experimental data as a function of the mixture equivalence ratio.
The parameter study is performed in the multi-thread mode by using the ``threading`` package.

Since the transport processes are critical for flame calculations, the transport data must be
included in the mechanism data and preprocessed.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_flame_speed_threading.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os
import threading
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger

# Chemkin 1-D premixed freely propagating flame model (steady-state)
from ansys.chemkin.premixedflames.premixedflame import FreelyPropagating as FlameSpeed
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

#######################################
# Create a flame speed calculator class
# =====================================
# Create a local class that wraps around the actual ``FlameSpeed`` class
# to make the setup of the multi-thread flame speed calculation parameter study
# more convenient.


class FlameSpeedCalculator:
    """
    Laminar flame speed calculator with fixed set up parameters
    """

    def __init__(self, fresh_mixture: Stream, index: int):
        """
        Laminar flame speed calculator that instantiates a FlameSpeed object with
        the given fresh (unburnt) mixture condition.

        Parameters
        ----------
            fresh_mixture: Mixture object
                the initial/fresh/unburnt condition
            index: integer
                run index of this flame speed calculator
        """
        # instantiate a flame speed object
        # set up the run and working directory name
        name = "Flame_Speed_" + str(index)
        # instantiate the FlameSpeed object for this run
        self.FScalculator = FlameSpeed(fresh_mixture, label=name)
        # set the required premixed flame model parameters
        #
        # set the maximum total number of grid points allowed in the calculation (optional)
        # self.FScalculator.set_max_grid_points(150)
        # define the calculation domain [cm]
        self.FScalculator.end_position = 1.0
        # set the root directory
        self.root_dir = os.getcwd()
        # set the working directory
        self.work_dir = os.path.join(self.root_dir, name)
        # run status
        self.runstatus = -100
        # calculated laminar flame speed [cm/sec]
        self.flame_speed = 0.0

    def run(self):
        """
        Run the flame speed calculation in a separate working directory
        """
        # create or clean up the working directory for this run
        if os.path.isdir(self.work_dir):
            # directory exists
            for f in os.listdir(self.work_dir):
                file_path = os.path.join(self.work_dir, f)
                if os.path.isfile(file_path):
                    # remove any existing file
                    try:
                        os.remove(file_path)
                    except OSError as e:
                        print(f"Error removing {file_path}: {e}")
        else:
            # create a new directory
            os.mkdir(self.work_dir)
        # change to the working directory for this run
        os.chdir(self.work_dir)
        # run the flame speed calculation
        self.runstatus = self.FScalculator.run()
        # extract the laminar flame speed from the solution
        if self.runstatus == 0:
            # postprocess the solutions
            self.FScalculator.process_solution()
            # get the flame speed value [cm/sec]
            # because the memory is shared, it must be done as soon as the run is finished
            self.flame_speed = self.FScalculator.get_flame_speed()
        # go back the root directory
        os.chdir(self.root_dir)

    def get_flame_speed(self) -> float:
        """
        Get the predicted laminar flame speed

        Returns
        -------
            flame_speed: double
                predicted laminar flame speed [cm/sec]
        """
        return self.flame_speed


############################################################
# Set up the flame speed parameter study for multi-threading
# ==========================================================
# Create a list of ``FlameSpeedCalculator`` objects with different
# initial methane-air equivalence ratios from 0.6 to 1.6. Each object
# represents one parameter study case and will be assigned to its own thread
# when the parameter study is executed.


def prepare_multi_thread_runs() -> dict[float, FlameSpeedCalculator]:
    """
    Set up the parameter study runs for multi-threading.

    Returns
    -------
        flame_speed_runs: list of FlameSpeedCalculator objects
            flame speed calculation cases
    """
    ##########################################
    # Create an instance of the Chemistry Set
    # ========================================
    # The mechanism loaded is the GRI 3.0 mechanism for methane combustion.
    # The mechanism and its associated data files come with the standard Ansys Chemkin
    # installation under the subdirectory *"/reaction/data"*.
    #

    # set mechanism directory (the default Chemkin mechanism data directory)
    data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    mechanism_dir = data_dir
    # including the full file path is recommended
    chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
    thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
    tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
    # create a chemistry set based on GRI 3.0
    MyGasMech = ck.Chemistry(
        chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0"
    )

    ##############################
    # Preprocess the Chemistry Set
    # ============================

    # preprocess the mechanism files
    iError = MyGasMech.preprocess()
    if iError != 0:
        print("Error: failed to preprocess the mechanism!")
        print(f"       error code = {iError}")
        exit()

    ########################################################################
    # Set up the CH\ :sub:`4`\ -air mixture for the flame speed calculation
    # ======================================================================
    # Instantiate a stream named ``premixed`` for the inlet gas mixture.
    # This stream  is a mixture with the addition of the
    # inlet flow rate. You can specify the inlet gas properties the same way you
    # set up a ``Mixture``. Here the ``X_by_Equivalence_Ratio`` method is used.
    # You create the ``fuel`` and the ``air`` mixtures first. Then define the
    # *complete combustion product species* and provide the *additives* composition
    # if applicable. And finally, during the parameter iteration runs, you can simply set
    # different values to ``equivalenceratio`` to create different methane-air mixtures.
    #

    # create the fuel mixture
    fuel = ck.Mixture(MyGasMech)
    # set fuel composition: methane
    fuel.X = [("CH4", 1.0)]
    # setting pressure and temperature condition for the flame speed calculations
    fuel.pressure = 5.0 * ck.Patm
    fuel.temperature = 300.0  # inlet temperature

    # create the oxidizer mixture: air
    air = ck.Mixture(MyGasMech)
    air.X = ck.Air.X()
    # setting pressure and temperature is not required in this case
    air.pressure = fuel.pressure
    air.temperature = fuel.temperature

    # create the fuel-air Stream for the premixed flame speed calculation
    premixed = Stream(MyGasMech, label="premixed")
    # products from the complete combustion of the fuel mixture and air
    products = ["CO2", "H2O", "N2"]
    # species mole fractions of added/inert mixture. can also create an additives mixture here
    add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros

    # setting pressure and temperature is not required in this case
    premixed.pressure = fuel.pressure
    premixed.temperature = fuel.temperature

    # set estimated value of the flame mass flux [g/cm2-sec]
    premixed.mass_flowrate = 0.4

    # Set up the flame speed parameter study for multi-threading
    # equivalence ratio for the first case
    phi = 0.6
    # total number of parameter cases
    points = 21
    # equivalence ratio increment
    delta_phi = 0.05
    # set up flame speed calculation runs
    FlameSpeedRuns: dict[float, FlameSpeedCalculator] = {}
    for i in range(points):
        # create mixture by using the equivalence ratio
        iError = premixed.X_by_Equivalence_Ratio(
            MyGasMech, fuel.X, air.X, add_frac, products, equivalenceratio=phi
        )
        # check fuel-oxidizer mixture creation status
        if iError != 0:
            print(
                "Error: failed to create the methane-air mixture "
                + "for equivalence ratio = "
                + str(phi)
            )
            exit()
        # create a flame speed calculation instance
        FlameSpeedRuns[phi] = FlameSpeedCalculator(premixed, index=i)
        # update parameter
        phi += delta_phi

    # return the job setup parameters
    return FlameSpeedRuns


########################################
# Set up and start the multi-thread runs
# ======================================
# Use the ``Thread()`` method to assign the flame speed runs.
# In this project, each flame speed run/case has its own thread.
# The `target` parameter of ``Thread()`` method should be the
# ``run()`` method of the ``FlameSpeedCalculator``.
# Use the ``start()`` method to initiate the threads,
# and the ``join()`` method to sync the threads after they
# are done.
#

FlameSpeedRuns = prepare_multi_thread_runs()
# set the start wall time
start_time = time.time()
threads = []
for phi, fsc in FlameSpeedRuns.items():
    t = threading.Thread(target=fsc.run())
    threads.append(t)
    # start each thread
    t.start()
# wait for all threads to finish
for t in threads:
    t.join()

# compute the total runtime
runtime = time.time() - start_time
print()
print(f"total simulation duration: {runtime} [sec]")
print()

#################################
# Get the parameter study results
# ===============================
# Use the ``get_flame_speed()`` method to get the calculated
# laminar flame speed values from each run. Then set up the experimental
# flame speed data for comparison.
#

points = len(FlameSpeedRuns.keys())
equival = np.zeros(points, dtype=np.double)
flamespeed = np.zeros_like(equival, dtype=np.double)
for i, case in enumerate(FlameSpeedRuns.items()):
    # equivalence ratio
    equival[i] = case[0]
    # flame speed calculator case
    fsc = case[1]
    flamespeed[i] = fsc.get_flame_speed()

# experimental data by Kochar
# equivalence ratios
data_equiv = [
    0.7005,
    0.8007,
    0.9009,
    1.001,
    1.1032,
    1.2014,
    1.3014,
]
# methane flame speeds at 5 atm
data_speed = [
    6.906,
    12.0094,
    15.9072,
    19.2376,
    19.6601,
    15.8274,
    10.2925,
]

###########################################
# Plot the premixed flame solution profiles
# =========================================
# Plot the predicted flame speeds against the experimental data.

plt.plot(data_equiv, data_speed, label="data", linestyle="", marker="^", color="blue")
plt.plot(equival, flamespeed, label="GRI 3.0", linestyle="-", color="blue")
plt.legend()
plt.ylabel("Flame Speed [cm/sec]")
plt.xlabel("Equivalence Ratio")
# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_flame_speed_threading.png", bbox_inches="tight")

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
.. _ref_flame_speed_pooling:

=================================================================
Setting up a multi-process laminar flame speed parameter study
=================================================================

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
The parameter study is performed in the multi-process (or multi-core) mode by using ``ProcessPoolExecutor``
from the ``concurrent.futures`` package. If the computer has enough number of CPU cores, the parameter
study can be run in parallel with each case running on a designated CPU core.

Since the transport processes are critical for flame calculations, the transport data must be
included in the mechanism data and preprocessed.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_flame_speed_pooling.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import os
import time

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger
from ansys.chemkin.utilities import workingFolders

# Chemkin 1-D premixed freely propagating flame model (steady-state)
from ansys.chemkin.premixedflames.premixedflame import FreelyPropagating as FlameSpeed
from concurrent.futures import ProcessPoolExecutor
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
# to make the setup of the multi-process flame speed calculation parameter study
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
        # set up the run and working direcotry name
        self.name = "Flame_Speed_" + str(index)
        # instantiate the FlameSpeed object for this run
        self.FScalculator = FlameSpeed(fresh_mixture, label=self.name)
        # set the required premixed flame model parameters
        #
        # set the maximum total number of grid points allowed in the calculation (optional)
        # self.FScalculator.set_max_grid_points(150)
        # define the calculation domain [cm]
        self.FScalculator.end_position = 1.0
        # set the root directory
        self.root_dir = os.getcwd()
        # set the working directory
        self.work_dir = os.path.join(self.root_dir, self.name)
        # run status
        self.runstatus = -100
        # calculated laminar flame speed [cm/sec]
        self.flame_speed = 0.0

    def run(self):
        """
        Run the flame speed calculation in a separate working directory
        """
        # run the flame speed calculation
        self.runstatus = self.FScalculator.run()
        # extract the laminar flame speed from the solution
        if self.runstatus == 0:
            # postprocess the solutions
            self.FScalculator.process_solution()
            # get the flame speed value [cm/sec]
            # because the memory is shared, it must be done as soon as the run is finished
            self.flame_speed = self.FScalculator.get_flame_speed()

    def get_flame_speed(self) -> float:
        """
        Get the predicted laminar flame speed

        Returns
        -------
            flame_speed: double
                predicted laminar flame speed [cm/sec]
        """
        return self.flame_speed


#############################################################
# Set up the flame speed parameter study for multi-processing
# ===========================================================
# Create a list of ``FlameSpeedCalculator`` objects with different
# initial methane-air equivalence ratios. Each object
# represents one parameter study case and will be run on an designated
# cpu core when the parameter study is executed.

def flame_speed_run(case: tuple[int, float]) -> tuple[float, float]:
    """
    Set up the parameter study runs for multi-processing.

    Parameters
    ----------
        case: tuple (integer, double)
            flame speed calculation case condition: case index and the equivalence ratio

    Returns
    -------
        phi: double
            the parameter: fresh mixture equivalence ratio
        flame_speed: double
            predicted laminar flame speed of the fresh mixture [cm/sec]
    """
##########################################
# Create an instance of the Chemistry Set
# ========================================
# The mechanism loaded is the GRI 3.0 mechanism for methane combustion.
# The mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"/reaction/data"*.
#
    # case index
    index = case[0]
    # fresh mixture equivalence ratio
    phi = case[1]
    # create and change to the working directory for this run
    name = "Flame_Speed_" + str(index)
    work_folder = workingFolders(name, current_dir)
    # set mechanism directory (the default Chemkin mechanism data directory)
    data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    mechanism_dir = data_dir
    # including the full file path is recommended
    chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
    thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
    tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
    # create a chemistry set based on GRI 3.0
    MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

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
    this_run = FlameSpeedCalculator(premixed, index=index)
    # run the case
    this_run.run()
    # change back to the original top folder
    work_folder.done()
    #
    return phi, this_run.flame_speed

#########################################
# Set up and start the multi-process runs
# =======================================
# Use the ``ProcessPoolExecutor()`` method to assign the flame speed runs.
# In this project, each flame speed run/case runs on an exclusive cpu core.
# The first parameter of ``map()`` method should be ``flame_speed_run()``.
# Make the ``flame_speed_run()`` method return the required parameter and results
# to make post-processing the results from the parameter study easier.
#
# .. note::
#   ``if __name__ == '__main__'`` is required when ``multiprocessing``
#   package is used.
#

if __name__ == '__main__':
    # number of available cpu cores
    numb_cores = max(os.cpu_count(), 1)
    # set the number of worker cpu cores to be used by the parameter study
    numb_workers = 14
    numb_workers = max(numb_workers, 1)
    numb_workers = min(numb_workers, numb_cores - 2)
    # Set up the flame speed parameter study for multi-processing
    # equivalence ratio for the first case
    phi = 0.6
    # total number of parameter cases
    numb_cases = 21
    # equivalence ratio increment
    delta_phi = 0.05
    # set up flame speed calculation runs with different equivalence ratios
    FlameSpeed_cases: tuple[int, float] = []
    for i in range(numb_cases):
        # create mixture by using the equivalence ratio
        FlameSpeed_cases.append((i, phi))
        # update parameter
        phi += delta_phi
    # start the multi-process
    case_id = 0
    # set the start wall time
    start_time = time.time()
    #
    numb_workers = min(numb_workers, numb_cases)
    equiv: float = []
    flame_speed: float = []
    with ProcessPoolExecutor(max_workers=numb_workers) as e:
        for ret_value in e.map(flame_speed_run, FlameSpeed_cases):
            # results returned by the flame speed calculator
            # equivalence ratio
            equiv.append(ret_value[0])
            # predicted laminar flame speed [cm/sec]
            flame_speed.append(ret_value[1])

    # compute the total runtime
    runtime = time.time() - start_time
    print()
    print(f"total simulation duration: {runtime} [sec]")
    print()

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
    plt.plot(equiv, flame_speed, label="GRI 3.0", linestyle="-", color="blue")
    plt.legend()
    plt.ylabel("Flame Speed [cm/sec]")
    plt.xlabel("Equivalence Ratio")
# plot results
    if interactive:
        plt.show()
    else:
        plt.savefig("plot_flame_speed_pooling.png", bbox_inches="tight")
